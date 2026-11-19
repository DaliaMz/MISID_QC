# Raw-score margins computed without sorting every row. `top_margin` is the
# top minus second-best score; `annotation_margin` is the top minus the
# observed label's score (0 when the observed label is the top prediction).
.score_margins <- function(scores, y) {
  rows <- seq_len(nrow(scores))
  top_cell <- cbind(rows, max.col(scores, ties.method = "first"))
  top_score <- scores[top_cell]
  masked <- scores
  masked[top_cell] <- -Inf
  second_score <- masked[cbind(rows, max.col(masked, ties.method = "first"))]
  observed_score <- scores[cbind(rows, match(y, colnames(scores)))]
  predicted <- colnames(scores)[top_cell[, 2L]]
  list(
    predicted = predicted,
    correct = predicted == y,
    top_score = top_score,
    observed_score = observed_score,
    top_margin = top_score - second_score,
    annotation_margin = top_score - observed_score
  )
}

# Hold-out cutoff based on margins from correctly predicted samples.
.holdout_threshold <- function(scores, correct, coverage,
                               reference = c("all_alternatives", "top")) {
  reference <- match.arg(reference)
  if (!any(correct)) {
    return(Inf)
  }
  scores <- scores[correct, , drop = FALSE]
  top_cell <- cbind(seq_len(nrow(scores)), max.col(scores, ties.method = "first"))
  margins <- scores[top_cell] - scores
  margins[top_cell] <- NA_real_
  if (reference == "top") {
    margins <- apply(margins, 1L, min, na.rm = TRUE)
  }
  unname(stats::quantile(margins, probs = coverage, type = 1, na.rm = TRUE))
}

# Exchange labels between random pairs of samples from different individuals,
# which preserves per-individual sample counts.
.swap_labels <- function(y, n_swaps) {
  y_sim <- y
  available <- seq_along(y)
  for (s in seq_len(n_swaps)) {
    candidates <- available
    if (length(candidates) < 2L) break
    i <- candidates[sample.int(length(candidates), 1L)]
    others <- candidates[y[candidates] != y[i]]
    if (!length(others)) break
    j <- others[sample.int(length(others), 1L)]
    y_sim[c(i, j)] <- y[c(j, i)]
    available <- setdiff(available, c(i, j))
  }
  y_sim
}

# Re-runs out-of-fold scoring on labels with known, simulated swaps, applying
# the same confidently-misassigned hold-out as the real data, so agreement
# with the pre-swap labels can be measured.
.simulate_mismatches <- function(x, y, rate, n_folds, shrinkage,
                                 margin_coverage, holdout_reference) {
  y_sim <- .swap_labels(y, max(1L, round(rate * length(y) / 2)))
  if (all(y_sim == y)) {
    return(NULL)
  }
  y_true <- y
  scores <- .simulation_oof_scores(x, y_sim, n_folds, shrinkage)
  margins <- .score_margins(scores, y_sim)
  threshold <- .holdout_threshold(scores, margins$correct, margin_coverage,
                                  holdout_reference)
  held_out <- !margins$correct & margins$annotation_margin > threshold
  data.frame(
    sample_id = rownames(x),
    true_id = y_true,
    simulated_id = y_sim,
    predicted_id = margins$predicted,
    top_margin = margins$top_margin,
    top_score = margins$top_score,
    annotation_margin = margins$annotation_margin,
    actually_misannotated = y_sim != y_true,
    confidently_misassigned = held_out,
    agreement_simulated = margins$predicted == y_sim,
    agreement_real = margins$predicted == y_true,
    stringsAsFactors = FALSE
  )
}

# Risk-coverage table. A sample is retained when its top margin is at or above
# the threshold. Cumulative sums evaluate every distinct cutoff in one pass; the
# table is then thinned to at most `n_thresholds` evenly spaced rows.
.risk_coverage <- function(top_margin, agreement, agreement_real = NULL,
                           n_thresholds = 100L) {
  ordering <- order(top_margin, decreasing = TRUE)
  sorted_margin <- top_margin[ordering]
  last_tie <- !duplicated(sorted_margin, fromLast = TRUE)
  retained <- seq_along(ordering)[last_tie]
  table <- data.frame(
    threshold = sorted_margin[last_tie],
    coverage = retained / length(ordering),
    n_retained = retained,
    agreement = cumsum(agreement[ordering])[last_tie] / retained,
    agreement_real = if (is.null(agreement_real)) {
      NA_real_
    } else {
      cumsum(agreement_real[ordering])[last_tie] / retained
    }
  )
  if (nrow(table) > n_thresholds) {
    table <- table[unique(round(seq(1, nrow(table),
                                    length.out = n_thresholds))), ,
                   drop = FALSE]
  }
  rownames(table) <- NULL
  table
}

.safe_divide <- function(numerator, denominator) {
  if (denominator == 0) NA_real_ else numerator / denominator
}

.estimate_detection_errors <- function(flagged, sensitivity, specificity, source) {
  result <- list(
    status = "unavailable", reason = NULL, n_samples = length(flagged),
    n_flagged = sum(flagged), sensitivity = sensitivity,
    specificity = specificity, source = source,
    prevalence = NA_real_, scenarios = NULL
  )
  if (!source %in% c("simulated_swaps_cross_fitted", "simulated_swaps_in_sample") ||
      length(sensitivity) != 1L || length(specificity) != 1L ||
      !is.finite(sensitivity) || !is.finite(specificity)) {
    result$reason <- "Simulated swap sensitivity and specificity are unavailable."
    return(result)
  }
  separation <- sensitivity + specificity - 1
  if (separation <= 0) {
    result$reason <- "Simulated sensitivity and specificity cannot distinguish errors."
    return(result)
  }
  flag_rate <- mean(flagged)
  false_positive_rate <- 1 - specificity
  tolerance <- 16 * .Machine$double.eps
  if (!is.finite(flag_rate) ||
      flag_rate < false_positive_rate - tolerance ||
      flag_rate > sensitivity + tolerance) {
    result$reason <- "The observed flag rate is incompatible with the simulated rates."
    return(result)
  }
  prevalence <- max(0, min(1, (flag_rate - false_positive_rate) / separation))
  rates <- c(max(0, prevalence - 0.025), prevalence,
             min(1, prevalence + 0.025))
  result$status <- "available"
  result$prevalence <- prevalence
  result$scenarios <- data.frame(
    scenario = c("lower", "estimate", "upper"),
    prevalence = rates,
    expected_false_positives = length(flagged) * (1 - rates) * (1 - specificity),
    expected_false_negatives = length(flagged) * rates * (1 - sensitivity)
  )
  result
}

.format_detection_error_estimate <- function(estimate) {
  if (is.null(estimate)) {
    return("Estimated detection errors unavailable (not present in this result).")
  }
  if (!identical(estimate$status, "available")) {
    return(paste("Estimated detection errors unavailable:", estimate$reason))
  }
  scenarios <- estimate$scenarios
  center <- scenarios[scenarios$scenario == "estimate", ]
  c(
    sprintf("Estimated misannotation rate: %.1f%% (assumed range %.1f%%-%.1f%%).",
            100 * center$prevalence, 100 * scenarios$prevalence[1L],
            100 * scenarios$prevalence[3L]),
    sprintf("Projected false positives: %.1f (range %.1f-%.1f).",
            center$expected_false_positives,
            min(scenarios$expected_false_positives), max(scenarios$expected_false_positives)),
    sprintf("Projected missed misannotations: %.1f (range %.1f-%.1f).",
            center$expected_false_negatives,
            min(scenarios$expected_false_negatives), max(scenarios$expected_false_negatives)),
    "Exploratory projections, not confidence intervals. Swaps may not represent",
    "real errors or unknown individuals."
  )
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

# Cutoff separating `negatives` from `positives`. With direction "upper" a
# value is flagged when it exceeds the cutoff, with "lower" when it falls below
# it. Each supplied negative population must independently meet the target.
.select_threshold <- function(negatives, positives, target_specificity,
                              direction = c("upper", "lower"),
                              additional_negatives = NULL) {
  direction <- match.arg(direction)
  negatives <- negatives[is.finite(negatives)]
  positives <- positives[is.finite(positives)]
  additional_negatives <- additional_negatives[is.finite(additional_negatives)]
  if (!length(positives)) {
    cutoff <- if (direction == "upper") {
      unname(stats::quantile(negatives, probs = target_specificity, type = 1))
    } else {
      .select_lower_threshold(negatives, target_specificity)
    }
    if (length(additional_negatives)) {
      extra_cutoff <- if (direction == "upper") {
        unname(stats::quantile(additional_negatives, probs = target_specificity, type = 1))
      } else {
        .select_lower_threshold(additional_negatives, target_specificity)
      }
      cutoff <- if (direction == "upper") max(cutoff, extra_cutoff) else
        min(cutoff, extra_cutoff)
    }
    return(cutoff)
  }
  values <- sort(unique(c(negatives, positives, additional_negatives)))
  candidates <- c(values[1L] - 1, (values[-1L] + values[-length(values)]) / 2,
                  values[length(values)] + 1)
  if (direction == "upper") {
    specificity <- vapply(candidates, function(c) mean(negatives <= c), 1)
    extra_specificity <- if (length(additional_negatives)) {
      vapply(candidates, function(c) mean(additional_negatives <= c), 1)
    } else rep(1, length(candidates))
    sensitivity <- vapply(candidates, function(c) mean(positives > c), 1)
  } else {
    specificity <- vapply(candidates, function(c) mean(negatives >= c), 1)
    extra_specificity <- if (length(additional_negatives)) {
      vapply(candidates, function(c) mean(additional_negatives >= c), 1)
    } else rep(1, length(candidates))
    sensitivity <- vapply(candidates, function(c) mean(positives < c), 1)
  }
  eligible <- which(specificity >= target_specificity &
                      extra_specificity >= target_specificity)
  best <- eligible[sensitivity[eligible] == max(sensitivity[eligible])]
  best <- best[specificity[best] == max(specificity[best])]
  tied <- candidates[best]
  # Margins span orders of magnitude and mostly sit near zero, so the middle
  # tied candidate by rank sits low in the gap; use the geometric midpoint.
  tied <- tied[tied > 0]
  if (length(tied) > 1L) {
    return(sqrt(min(tied) * max(tied)))
  }
  candidates[best][ceiling(length(best) / 2)]
}

# Evaluate swaps and unchanged samples using cutoffs chosen from other repeats'
# simulated samples and the real reference margins.
.crossfit_swap_detection <- function(negatives, simulation,
                                     target_specificity) {
  swaps <- simulation[simulation$actually_misannotated &
                        is.finite(simulation$annotation_margin), ,
                      drop = FALSE]
  by_repeat <- split(swaps$annotation_margin, swaps$repeat_id)
  by_repeat <- by_repeat[lengths(by_repeat) > 0L]
  if (length(by_repeat) < 2L) {
    return(NULL)
  }
  unchanged <- simulation[!simulation$actually_misannotated &
                            is.finite(simulation$annotation_margin), , drop = FALSE]
  if (!nrow(unchanged) ||
      any(vapply(names(by_repeat), function(r) {
        !any(unchanged$repeat_id != r)
      }, logical(1)))) {
    return(NULL)
  }
  per_repeat <- lapply(names(by_repeat), function(r) {
    threshold <- .select_threshold(
      negatives, unlist(by_repeat[names(by_repeat) != r], use.names = FALSE),
      target_specificity, "upper",
      additional_negatives = unchanged$annotation_margin[unchanged$repeat_id != r]
    )
    held_unchanged <- unchanged$annotation_margin[unchanged$repeat_id == r]
    c(sensitivity = mean(by_repeat[[r]] > threshold),
      n_unchanged = length(held_unchanged),
      n_correct = sum(held_unchanged <= threshold))
  })
  per_repeat <- do.call(rbind, per_repeat)
  sizes <- lengths(by_repeat)
  list(
    sensitivity = sum(per_repeat[, "sensitivity"] * sizes) / sum(sizes),
    sensitivity_sd = stats::sd(per_repeat[, "sensitivity"]),
    simulated_specificity = .safe_divide(sum(per_repeat[, "n_correct"]),
                                         sum(per_repeat[, "n_unchanged"])),
    n_simulated_negatives = sum(per_repeat[, "n_unchanged"]),
    n_repeats = nrow(per_repeat)
  )
}

# Top margins and nearest-class distances of simulated strangers. For each held-out individual only up to
# `n_strays` of their samples (never all of them) stay in the data, relabelled
# as random other enrolled individuals, as a stranger's samples would be in a
# real cohort; their remaining samples are dropped. Each stray is scored by a
# model that never saw it but was trained with the other strays under their
# wrong labels. With `n_folds`, individuals are processed in groups.
.unknown_validation <- function(x, y, shrinkage, exclude = rep(FALSE, length(y)),
                                progress = FALSE, n_folds = NULL,
                                n_strays = 2L) {
  classes <- sort(unique(y))
  if (length(classes) < 3L) {
    if (isTRUE(progress)) {
      message("Skipping held-out-individual validation: at least three individuals are needed.")
    }
    return(list(n = 0L, top_margin = numeric(), distance = numeric(),
                status = "insufficient_evidence"))
  }
  n_groups <- if (is.null(n_folds)) {
    length(classes)
  } else {
    min(as.integer(n_folds), length(classes))
  }
  # Every group must leave at least two enrolled individuals to train on.
  while (ceiling(length(classes) / n_groups) > length(classes) - 2L) {
    n_groups <- n_groups + 1L
  }
  group <- sample(rep(seq_len(n_groups), length.out = length(classes)))
  names(group) <- classes
  if (isTRUE(progress)) {
    message(sprintf("Validating held-out individuals in %d groups across %d individuals.",
                    n_groups, length(classes)))
  }
  unknown_margin <- numeric()
  unknown_distance <- numeric()
  progress_interval <- max(1L, ceiling(n_groups / 10L))
  for (g in seq_len(n_groups)) {
    members <- classes[group == g]
    in_group <- group[y] == g
    strays <- integer()
    for (individual in members) {
      idx <- which(y == individual & !exclude)
      if (length(idx) >= 2L) {
        strays <- c(strays, idx[sample.int(length(idx),
                                           min(n_strays, length(idx) - 1L))])
      }
    }
    if (length(strays)) {
      others <- classes[!classes %in% members]
      wrong <- vapply(strays, function(i) sample(others, 1L), character(1))
      half <- if (length(strays) == 1L) 1L else
        sample(rep(1:2, length.out = length(strays)))
      base <- !in_group & !exclude
      for (h in unique(half)) {
        scored <- strays[half == h]
        model <- fit_shrinkage_lda(
          rbind(x[base, , drop = FALSE], x[strays[half != h], , drop = FALSE]),
          c(y[base], wrong[half != h]), shrinkage = shrinkage
        )
        scores <- predict_shrinkage_lda(model, x[scored, , drop = FALSE])
        unknown_margin <- c(
          unknown_margin, .score_margins(scores, y[scored])$top_margin
        )
        unknown_distance <- c(
          unknown_distance,
          .nearest_class_distance(model, x[scored, , drop = FALSE])
        )
      }
    }
    if (isTRUE(progress) &&
        (g %% progress_interval == 0L || g == n_groups)) {
      message(sprintf("Held-out-individual validation: %d/%d groups.",
                      g, n_groups))
    }
  }
  list(n = length(unknown_margin), top_margin = unknown_margin,
       distance = unknown_distance,
       status = if (length(unknown_margin)) "estimated" else
         "insufficient_evidence")
}

# Robust lower cutoff on the log top margin of enrolled samples. It uses no
# unknown samples and no single extreme sample, and a clean cohort flags
# almost nothing.
.unknown_margin_threshold <- function(known_margin, mad_multiplier) {
  log_margin <- log(pmax(known_margin, .Machine$double.eps))
  exp(stats::median(log_margin) - mad_multiplier * stats::mad(log_margin))
}

# Out-of-fold scores for a simulated label set; NULL means leave-one-out.
.simulation_oof_scores <- function(x, y, n_folds, shrinkage) {
  classes <- sort(unique(y))
  scores <- matrix(NA_real_, nrow(x), length(classes),
                   dimnames = list(rownames(x), classes))
  if (is.null(n_folds)) {
    fold_id <- seq_len(nrow(x))
  } else {
    fold_id <- integer(nrow(x))
    for (individual in classes) {
      samples <- which(y == individual)
      fold_id[samples] <- sample(rep(seq_len(n_folds), length.out = length(samples)))
    }
  }
  folds <- unique(fold_id)
  for (i in seq_along(folds)) {
    fold <- folds[i]
    test <- fold_id == fold
    train <- !test
    model <- fit_shrinkage_lda(x[train, , drop = FALSE], y[train],
                               shrinkage = shrinkage)
    scores[test, ] <- predict_shrinkage_lda(model, x[test, , drop = FALSE])
  }
  scores
}

# Lowest top-margin cutoff at which agreement among retained samples
# (margin at or above the cutoff) reaches `min_accuracy`.
.abstention_threshold <- function(known_margin, known_correct, min_accuracy,
                                  n_candidates = 200L) {
  candidates <- sort(unique(known_margin))
  if (length(candidates) > n_candidates) {
    candidates <- candidates[unique(round(seq(1, length(candidates),
                                              length.out = n_candidates)))]
  }
  accuracy <- vapply(candidates, function(cutoff) {
    keep <- known_margin >= cutoff
    if (!any(keep)) NA_real_ else mean(known_correct[keep])
  }, numeric(1))
  ok <- !is.na(accuracy) & accuracy >= min_accuracy
  list(threshold = if (any(ok)) min(candidates[ok]) else Inf)
}

.correction_threshold <- function(simulation, min_accuracy) {
  empty <- list(threshold = Inf, n_eligible = 0L, n_retained = 0L,
                coverage = NA_real_, accuracy = NA_real_)
  if (is.null(simulation)) return(empty)
  eligible <- simulation$actually_misannotated |
    !simulation$agreement_simulated
  selected <- simulation[eligible, , drop = FALSE]
  empty$n_eligible <- nrow(selected)
  if (!nrow(selected)) return(empty)
  empty$coverage <- 0
  curve <- .risk_coverage(selected$top_margin, selected$agreement_real,
                          n_thresholds = nrow(selected))
  passing <- which(curve$agreement >= min_accuracy)
  if (!length(passing)) return(empty)
  chosen <- passing[which.min(curve$threshold[passing])]
  list(threshold = curve$threshold[chosen], n_eligible = nrow(selected),
       n_retained = curve$n_retained[chosen], coverage = curve$coverage[chosen],
       accuracy = curve$agreement[chosen])
}
