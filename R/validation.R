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
.holdout_threshold <- function(scores, correct, coverage) {
  if (!any(correct)) {
    return(Inf)
  }
  scores <- scores[correct, , drop = FALSE]
  top_cell <- cbind(seq_len(nrow(scores)), max.col(scores, ties.method = "first"))
  margins <- scores[top_cell] - scores
  margins[top_cell] <- NA_real_
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
                                 margin_coverage) {
  y_sim <- .swap_labels(y, max(1L, round(rate * length(y) / 2)))
  if (all(y_sim == y)) {
    return(NULL)
  }
  y_true <- y
  scores <- .simulation_oof_scores(x, y_sim, n_folds, shrinkage)
  margins <- .score_margins(scores, y_sim)
  threshold <- .holdout_threshold(scores, margins$correct, margin_coverage)
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

# Check whether covariates contain enough independent variation to be useful
# fixed effects or linear covariates. Sparse categorical levels and constant
# numeric covariates add no reliable adjustment and can make residualization
# unstable.
.validate_covariates <- function(covariates) {
  covariates <- as.data.frame(covariates)
  if (!ncol(covariates) || nrow(covariates) < 2L) {
    stop("Covariates must provide at least two samples.")
  }
  for (j in seq_len(ncol(covariates))) {
    name <- names(covariates)[j]
    value <- covariates[[j]]
    if (is.numeric(value)) {
      if (any(!is.finite(value))) {
        stop("Numeric covariates must contain only finite values.")
      }
      if (length(unique(value)) < 2L) {
        stop(sprintf("Numeric covariate '%s' is constant and cannot be residualized.", name))
      }
      if (length(unique(value)) < 5L) {
        warning(sprintf(
          "Numeric covariate '%s' has few distinct values; it may behave like a categorical covariate.",
          name
        ), call. = FALSE)
      }
    } else {
      if (anyNA(value) || any(!nzchar(as.character(value)))) {
        stop("Covariates must contain no missing or empty values.")
      }
      counts <- table(as.character(value))
      if (length(counts) < 2L) {
        stop(sprintf("Categorical covariate '%s' has only one level.", name))
      }
      if (length(counts) >= nrow(covariates) - 1L) {
        stop(sprintf(
          "Categorical covariate '%s' has too many levels for the number of samples.",
          name
        ))
      }
      sparse <- counts[counts <= 2L]
      if (length(sparse)) {
        warning(sprintf(
          "Categorical covariate '%s' has sparse level(s) represented by two or fewer samples: %s.",
          name,
          paste(sprintf("%s (n=%d)", names(sparse), sparse), collapse = ", ")
        ), call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

# Residualize sample-by-feature expression against fixed-effect covariates.
# Covariates are removed without using individual identities, swaps, or
# unknown-sample labels. Character, logical, and factor columns become fixed
# categorical effects; numeric columns remain linear effects.
.covariate_residualize <- function(x, covariates) {
  covariates <- as.data.frame(covariates)
  if (!ncol(covariates) || nrow(covariates) != nrow(x)) {
    stop("Covariates must provide one row per sample.")
  }
  .validate_covariates(covariates)
  covariate_frame <- covariates
  for (j in seq_len(ncol(covariates))) {
    value <- covariates[[j]]
    if (is.numeric(value)) {
      if (any(!is.finite(value))) {
        stop("Numeric covariates must contain only finite values.")
      }
      covariate_frame[[j]] <- as.numeric(value)
    } else {
      if (anyNA(value) || any(!nzchar(as.character(value)))) {
        stop("Covariates must contain no missing or empty values.")
      }
      covariate_frame[[j]] <- factor(as.character(value))
    }
  }
  design <- stats::model.matrix(~ ., data = covariate_frame)
  fit <- stats::lm.fit(design, as.matrix(x))
  residual <- fit$residuals
  residual <- as.matrix(residual)
  dimnames(residual) <- dimnames(x)
  storage.mode(residual) <- "double"
  residual
}

.detection_verdict <- function(sensitivity, specificity, target_specificity,
                               enough_known, swap_simulation) {
  if (!enough_known || !swap_simulation ||
      !is.finite(sensitivity) || !is.finite(specificity)) {
    return("insufficient_evidence")
  }
  shortfalls <- c(0.8 - sensitivity, target_specificity - specificity)
  if (all(shortfalls <= 0)) return("appears_feasible")
  if (any(shortfalls > 0.1 + 16 * .Machine$double.eps)) return("unreliable")
  "below_target"
}

.estimate_detection_errors <- function(flagged, sensitivity, specificity, source) {
  result <- list(
    status = "unavailable", reason = NULL, n_samples = length(flagged),
    n_flagged = sum(flagged), sensitivity = sensitivity,
    specificity = specificity, source = source,
    prevalence = NA_real_, scenarios = NULL
  )
  if (!identical(source, "simulated_swaps_in_sample") ||
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
    return(list(n = 0L, distance = numeric(),
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
  unknown_distance <- numeric()
  records <- list()
  score_rows <- list()
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
        margins <- .score_margins(scores, wrong[half == h])
        distances <- .nearest_class_distance(model, x[scored, , drop = FALSE],
                                              scores)
        unknown_distance <- c(unknown_distance, distances)
        records[[length(records) + 1L]] <- data.frame(
          sample_id = rownames(x)[scored],
          original_id = y[scored],
          simulated_id = wrong[half == h],
          predicted_id = margins$predicted,
          top_candidate_score = margins$top_score,
          simulated_id_score = margins$observed_score,
          annotation_margin = margins$annotation_margin,
          top_margin = margins$top_margin,
          nearest_class_distance = distances,
          holdout_group = g,
          scoring_fold = h,
          stringsAsFactors = FALSE
        )
        aligned_scores <- matrix(
          NA_real_, nrow(scores), length(classes),
          dimnames = list(rownames(scores), classes)
        )
        aligned_scores[, colnames(scores)] <- scores
        score_rows[[length(score_rows) + 1L]] <- aligned_scores
      }
    }
    if (isTRUE(progress) &&
        (g %% progress_interval == 0L || g == n_groups)) {
      message(sprintf("Held-out-individual validation: %d/%d groups.",
                      g, n_groups))
    }
  }
  list(n = length(unknown_distance), distance = unknown_distance,
       records = if (length(records)) do.call(rbind, records) else NULL,
       identity_scores = if (length(score_rows)) do.call(rbind, score_rows) else NULL,
       status = if (length(unknown_distance)) "estimated" else
         "insufficient_evidence")
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

.correction_threshold <- function(simulation, min_accuracy) {
  empty <- list(threshold = Inf, n_eligible = 0L, n_retained = 0L,
                coverage = NA_real_, accuracy = NA_real_,
                repeat_accuracy = numeric(),
                repeat_accuracy_range = c(NA_real_, NA_real_))
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
  threshold <- curve$threshold[chosen]
  retained <- selected$top_margin >= threshold
  repeat_accuracy <- if ("repeat_id" %in% names(selected)) {
    repeat_ids <- selected$repeat_id[retained]
    retained_rows <- selected[retained, , drop = FALSE]
    vapply(
      split(retained_rows, repeat_ids),
      function(rows) mean(rows$agreement_real),
      numeric(1)
    )
  } else numeric()
  repeat_accuracy_range <- if (length(repeat_accuracy)) {
    c(min(repeat_accuracy), max(repeat_accuracy))
  } else c(NA_real_, NA_real_)
  list(threshold = threshold, n_eligible = nrow(selected),
       n_retained = curve$n_retained[chosen], coverage = curve$coverage[chosen],
       accuracy = curve$agreement[chosen],
       repeat_accuracy = repeat_accuracy,
       repeat_accuracy_range = repeat_accuracy_range)
}

.correction_gate <- function(simulation, correction, n_input,
                             min_accuracy, min_correction_samples,
                             repeat_tolerance) {
  result <- list(
    status = "insufficient_evidence",
    reason = "No usable correction calibration.",
    threshold = correction$threshold,
    n_input = n_input,
    n_input_required = 20L,
    n_eligible = correction$n_eligible,
    n_retained = correction$n_retained,
    n_retained_required = min_correction_samples,
    coverage = correction$coverage,
    accuracy = correction$accuracy,
    target_accuracy = min_accuracy,
    repeat_accuracy = correction$repeat_accuracy,
    repeat_accuracy_range = correction$repeat_accuracy_range,
    repeat_tolerance = repeat_tolerance
  )
  if (is.null(simulation) || !is.finite(correction$threshold) ||
      correction$n_eligible < 1L) {
    return(result)
  }
  if (n_input < 20L) {
    result$reason <- "Fewer than 20 input samples are available."
    return(result)
  }
  if (correction$n_retained < min_correction_samples) {
    result$reason <- sprintf(
      "Fewer than %d simulated rows remain at the correction threshold.",
      min_correction_samples
    )
    return(result)
  }
  if (!is.finite(correction$accuracy) ||
      correction$accuracy < min_accuracy) {
    result$status <- "unreliable"
    result$reason <- "Retained correction accuracy misses the target."
    return(result)
  }
  repeat_lower <- if (length(correction$repeat_accuracy)) {
    min(correction$repeat_accuracy)
  } else NA_real_
  if (!is.finite(repeat_lower)) {
    result$reason <- "Correction accuracy is unavailable by simulation repeat."
    return(result)
  }
  if (repeat_lower < min_accuracy - repeat_tolerance) {
    result$status <- "unreliable"
    result$reason <- sprintf(
      "The lowest repeat accuracy is more than %.1f percentage points below target.",
      100 * repeat_tolerance
    )
    return(result)
  }
  result$status <- "supported"
  result$reason <- "Correction calibration meets the sample, row, accuracy, and repeat gates."
  result
}

.confident_prediction <- function(predicted, observed, margin, cutoff) {
  predicted == observed | margin >= cutoff
}
