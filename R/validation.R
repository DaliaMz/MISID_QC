.softmax_scores <- function(scores, temperature = 1) {
  scaled <- scores / temperature
  scaled <- scaled - apply(scaled, 1L, max)
  probabilities <- exp(scaled)
  probabilities / rowSums(probabilities)
}

.fit_temperature <- function(scores, truth) {
  class_index <- match(truth, colnames(scores))
  objective <- function(log_temperature) {
    probabilities <- .softmax_scores(scores, exp(log_temperature))
    -mean(log(pmax(probabilities[cbind(seq_along(class_index), class_index)],
                   .Machine$double.xmin)))
  }
  exp(stats::optimize(objective, interval = c(-5, 5))$minimum)
}

.safe_divide <- function(numerator, denominator) {
  if (denominator == 0) NA_real_ else numerator / denominator
}

.binary_metrics <- function(positive, flagged, score = NULL) {
  tp <- sum(positive & flagged)
  tn <- sum(!positive & !flagged)
  fp <- sum(!positive & flagged)
  fn <- sum(positive & !flagged)
  result <- data.frame(
    n = length(positive),
    sensitivity = .safe_divide(tp, tp + fn),
    specificity = .safe_divide(tn, tn + fp),
    false_positive_rate = .safe_divide(fp, fp + tn),
    false_negative_rate = .safe_divide(fn, fn + tp),
    precision = .safe_divide(tp, tp + fp),
    negative_predictive_value = .safe_divide(tn, tn + fn),
    accuracy = .safe_divide(tp + tn, length(positive))
  )
  if (!is.null(score) && any(positive) && any(!positive)) {
    ranks <- rank(score, ties.method = "average")
    result$auroc <- (sum(ranks[positive]) -
                       sum(positive) * (sum(positive) + 1) / 2) /
      (sum(positive) * sum(!positive))
    ordering <- order(score, decreasing = TRUE)
    ranked_positive <- positive[ordering]
    precision_at_rank <- cumsum(ranked_positive) / seq_along(ranked_positive)
    result$auprc <- mean(precision_at_rank[ranked_positive])
  } else {
    result$auroc <- NA_real_
    result$auprc <- NA_real_
  }
  result
}

.select_lower_threshold <- function(values, target_specificity) {
  unname(stats::quantile(values, probs = 1 - target_specificity, type = 1,
                         na.rm = TRUE))
}

.unknown_validation <- function(x, y, temperature, known_threshold, shrinkage,
                                known_confidence, progress = FALSE) {
  classes <- sort(unique(y))
  if (length(classes) < 3L) {
    if (isTRUE(progress)) {
      message("Skipping held-out-individual validation: at least three individuals are needed.")
    }
    return(list(
      metrics = NULL,
      threshold = NA_real_,
      n = 0L,
      status = "insufficient_evidence"
    ))
  }
  if (isTRUE(progress)) {
    message(sprintf("Validating held-out individuals across %d individuals.",
                    length(classes)))
  }
  unknown_confidence <- numeric()
  individual_progress_interval <- max(1L, ceiling(length(classes) / 10L))
  for (i in seq_along(classes)) {
    held_out <- classes[i]
    train <- y != held_out
    test <- !train
    model <- fit_shrinkage_lda(x[train, , drop = FALSE], y[train],
                               shrinkage = shrinkage)
    scores <- predict_shrinkage_lda(model, x[test, , drop = FALSE])
    confidence <- apply(.softmax_scores(scores, temperature), 1L, max)
    unknown_confidence <- c(unknown_confidence, confidence)
    if (isTRUE(progress) &&
        (i %% individual_progress_interval == 0L || i == length(classes))) {
      message(sprintf("Held-out-individual validation: %d/%d individuals.",
                      i, length(classes)))
    }
  }
  confidence <- c(known_confidence, unknown_confidence)
  positive <- c(rep(FALSE, length(known_confidence)),
                rep(TRUE, length(unknown_confidence)))
  flags <- confidence < known_threshold
  list(
    metrics = .binary_metrics(positive, flags, score = -confidence),
    threshold = known_threshold,
    n = length(unknown_confidence),
    status = "estimated"
  )
}
