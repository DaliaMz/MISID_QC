#' Assess repeated-measures RNA-seq sample annotations
#'
#' The expression matrix must contain transformed, finite values with features
#' in rows and samples in columns. Each enrolled individual needs at least two
#' samples. Samples are scored out of fold: an evaluated sample is excluded
#' from its own fitted reference. This initial implementation uses all supplied
#' features and does not normalize or select features.
#'
#' Class probabilities are temperature-scaled probabilities across enrolled
#' identities. They are not probabilities that an observed annotation is
#' correct. The function therefore leaves
#' `annotation_correctness_probability` missing; simulated error prevalence
#' does not identify a real-data annotation-error prevalence.
#'
#' @param expr Numeric expression matrix, features by samples, with unique
#'   sample names in `colnames(expr)`.
#' @param metadata Data frame containing sample and observed identity columns.
#' @param sample_col Name of the sample identifier column.
#' @param id_col Name of the observed individual identifier column.
#' @param shrinkage Shrinkage fraction for the diagonal covariance target.
#' @param target_specificity Target specificity used to select the lower-tail
#'   thresholds for annotation support and unknown detection.
#' @param assess_feasibility Whether to return task-specific feasibility
#'   assessments based on internal validation.
#' @return A list with sample-level `predictions`, raw `identity_scores`,
#'   temperature-scaled `identity_probabilities`, `performance`,
#'   `calibration`, `feasibility`, and `thresholds` components.
#' @export
check_annotations <- function(expr, metadata, sample_col, id_col,
                              shrinkage = 0.1, target_specificity = 0.95,
                              assess_feasibility = TRUE) {
  expr <- as.matrix(expr)
  storage.mode(expr) <- "double"
  if (length(dim(expr)) != 2L || nrow(expr) < 1L || ncol(expr) < 1L ||
      any(!is.finite(expr))) {
    stop("expr must be a non-empty matrix containing finite numeric values.")
  }
  if (is.null(colnames(expr)) || anyNA(colnames(expr)) ||
      any(!nzchar(colnames(expr))) || anyDuplicated(colnames(expr))) {
    stop("expr must have unique, non-empty sample names in its columns.")
  }
  if (!is.data.frame(metadata) ||
      !all(c(sample_col, id_col) %in% names(metadata))) {
    stop("metadata must contain the columns named by sample_col and id_col.")
  }
  sample_ids <- as.character(metadata[[sample_col]])
  observed_ids <- as.character(metadata[[id_col]])
  if (anyNA(sample_ids) || any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids) || anyNA(observed_ids) ||
      any(!nzchar(observed_ids))) {
    stop("Sample and individual identifiers must be unique/non-missing as applicable.")
  }
  if (!setequal(sample_ids, colnames(expr))) {
    stop("metadata sample identifiers must match expr column names exactly.")
  }
  metadata <- metadata[match(colnames(expr), sample_ids), , drop = FALSE]
  y <- as.character(metadata[[id_col]])
  class_counts <- table(y)
  if (length(class_counts) < 2L || any(class_counts < 2L)) {
    stop("At least two individuals with at least two samples each are required.")
  }
  if (length(target_specificity) != 1L ||
      !is.finite(target_specificity) ||
      target_specificity <= 0 || target_specificity >= 1) {
    stop("target_specificity must be a single number between 0 and 1.")
  }

  x <- t(expr)
  rownames(x) <- colnames(expr)
  oof_scores <- matrix(
    NA_real_, nrow = nrow(x), ncol = length(class_counts),
    dimnames = list(rownames(x), names(class_counts))
  )
  for (i in seq_len(nrow(x))) {
    train <- seq_len(nrow(x)) != i
    model <- fit_shrinkage_lda(x[train, , drop = FALSE], y[train],
                               shrinkage = shrinkage)
    oof_scores[i, ] <- predict_shrinkage_lda(
      model, x[i, , drop = FALSE]
    )[1L, ]
  }

  temperature <- .fit_temperature(oof_scores, y)
  class_probabilities <- .softmax_scores(oof_scores, temperature)
  class_index <- match(y, colnames(oof_scores))
  observed_support <- class_probabilities[
    cbind(seq_along(y), class_index)
  ]
  predicted_index <- max.col(oof_scores, ties.method = "first")
  predicted_id <- colnames(oof_scores)[predicted_index]
  ordered_scores <- t(apply(oof_scores, 1L, sort, decreasing = TRUE))
  top_margin <- if (ncol(oof_scores) > 1L) {
    ordered_scores[, 1L] - ordered_scores[, 2L]
  } else {
    rep(NA_real_, nrow(oof_scores))
  }
  support_threshold <- .select_lower_threshold(
    observed_support, target_specificity
  )
  flagged <- observed_support < support_threshold

  correct_identity <- predicted_id == y
  known_max_probability <- apply(class_probabilities, 1L, max)
  unknown_threshold <- .select_lower_threshold(
    known_max_probability, target_specificity
  )
  wrong_indices <- lapply(seq_along(y), function(i) {
    setdiff(seq_along(class_counts), class_index[i])
  })
  false_labels <- unlist(lapply(seq_along(y), function(i) {
    colnames(oof_scores)[wrong_indices[[i]]]
  }), use.names = FALSE)
  false_sample_index <- rep(seq_along(y), lengths(wrong_indices))
  false_support <- class_probabilities[
    cbind(false_sample_index, match(false_labels, colnames(oof_scores)))
  ]
  misannotation_positive <- c(rep(FALSE, length(y)),
                              rep(TRUE, length(false_support)))
  misannotation_flagged <- c(flagged, false_support < support_threshold)
  misannotation_score <- c(-observed_support, -false_support)
  misannotation_metrics <- .binary_metrics(
    misannotation_positive, misannotation_flagged, misannotation_score
  )
  identity_recovery <- if (length(false_sample_index)) {
    mean(predicted_id[false_sample_index] == y[false_sample_index])
  } else {
    NA_real_
  }
  swap_samples <- integer()
  swap_labels <- character()
  for (i in seq_len(length(class_counts) - 1L)) {
    for (j in seq.int(i + 1L, length(class_counts))) {
      first_samples <- which(y == names(class_counts)[i])
      second_samples <- which(y == names(class_counts)[j])
      n_pairs <- min(length(first_samples), length(second_samples))
      swap_samples <- c(
        swap_samples, first_samples[seq_len(n_pairs)],
        second_samples[seq_len(n_pairs)]
      )
      swap_labels <- c(
        swap_labels, rep(names(class_counts)[j], n_pairs),
        rep(names(class_counts)[i], n_pairs)
      )
    }
  }
  swap_support <- class_probabilities[
    cbind(swap_samples, match(swap_labels, colnames(oof_scores)))
  ]
  pairwise_swap_metrics <- .binary_metrics(
    c(rep(FALSE, length(y)), rep(TRUE, length(swap_support))),
    c(flagged, swap_support < support_threshold),
    c(-observed_support, -swap_support)
  )
  top_three <- vapply(seq_along(y), function(i) {
    ranked <- names(sort(oof_scores[i, ], decreasing = TRUE))
    y[i] %in% ranked[seq_len(min(3L, length(ranked)))]
  }, logical(1))

  unknown <- .unknown_validation(
    x, y, temperature, unknown_threshold, shrinkage, known_max_probability
  )
  calibration <- list(
    method = "temperature_scaling",
    temperature = temperature,
    multiclass_brier_score = mean(rowSums((class_probabilities -
      outer(y, colnames(oof_scores), `==`) * 1)^2)),
    multiclass_log_loss = -mean(log(pmax(
      observed_support, .Machine$double.xmin
    ))),
    note = paste(
      "Calibrates candidate identity probabilities, not annotation correctness;",
      "no real-world error prevalence is assumed."
    )
  )
  predictions <- data.frame(
    sample_id = rownames(x),
    observed_id = y,
    predicted_id = predicted_id,
    observed_id_score = oof_scores[cbind(seq_along(y), class_index)],
    top_candidate_score = apply(oof_scores, 1L, max),
    top_vs_observed_score_margin = apply(oof_scores, 1L, max) -
      oof_scores[cbind(seq_along(y), class_index)],
    observed_id_probability = observed_support,
    top_candidate_probability = known_max_probability,
    top_score_margin = top_margin,
    annotation_correctness_probability = NA_real_,
    suspicious = flagged,
    potential_unknown_individual = known_max_probability < unknown_threshold,
    most_likely_alternative = ifelse(
      flagged & known_max_probability >= unknown_threshold & predicted_id != y,
      predicted_id, NA_character_
    ),
    candidate_correction = NA_character_,
    stringsAsFactors = FALSE
  )
  performance <- list(
    correct_identity_recovery_accuracy = mean(correct_identity),
    top_three_identity_recovery_accuracy = mean(top_three),
    misannotation_detection = misannotation_metrics,
    identity_recovery_accuracy = identity_recovery,
    pairwise_swap_detection = pairwise_swap_metrics,
    pairwise_swap_identity_recovery_accuracy = mean(
      predicted_id[swap_samples] == y[swap_samples]
    ),
    unknown_individual_detection = unknown$metrics,
    unknown_simulation_n = unknown$n
  )
  feasibility <- NULL
  if (isTRUE(assess_feasibility)) {
    enough_known <- length(y) >= 20L
    detection_status <- if (!enough_known) {
      "insufficient_evidence"
    } else if (is.finite(misannotation_metrics$sensitivity) &&
               misannotation_metrics$sensitivity >= 0.8 &&
               misannotation_metrics$specificity >= target_specificity) {
      "appears_feasible"
    } else {
      "unreliable"
    }
    recovery_status <- if (!enough_known) {
      "insufficient_evidence"
    } else if (identity_recovery >= 0.8) {
      "appears_feasible"
    } else {
      "unreliable"
    }
    unknown_status <- if (unknown$n < 10L) {
      "insufficient_evidence"
    } else if (is.finite(unknown$metrics$sensitivity) &&
               unknown$metrics$sensitivity >= 0.8) {
      "appears_feasible"
    } else {
      "unreliable"
    }
    feasibility <- list(
      misannotation_detection = detection_status,
      identity_recovery = recovery_status,
      unknown_individual_detection = unknown_status,
      warning = paste(
        "Internal estimates are exploratory: threshold selection and",
        "performance use the same out-of-fold predictions."
      )
    )
    reliable_correction <- detection_status == "appears_feasible" &&
      recovery_status == "appears_feasible"
  } else {
    reliable_correction <- FALSE
  }
  predictions$candidate_correction <- ifelse(
      flagged & reliable_correction &
        !predictions$potential_unknown_individual & predicted_id != y,
      predicted_id,
      NA_character_
  )
  reliable_assessment <- !is.null(feasibility) &&
    feasibility$misannotation_detection == "appears_feasible"
  reliable_unknown <- !is.null(feasibility) &&
    feasibility$unknown_individual_detection == "appears_feasible"
  predictions$annotation_status <- ifelse(
    !reliable_assessment,
    "insufficient_evidence",
    ifelse(predictions$potential_unknown_individual & !reliable_unknown,
           "insufficient_evidence",
           ifelse(predictions$potential_unknown_individual,
                  "potential_unknown_individual",
           ifelse(!flagged, "no_evidence_of_incorrect_annotation",
                  ifelse(!is.na(predictions$candidate_correction),
                         "suspicious_with_candidate_correction",
                         "suspicious_without_reliable_alternative"))))
  )

  list(
    predictions = predictions,
    identity_scores = oof_scores,
    identity_probabilities = class_probabilities,
    performance = performance,
    calibration = calibration,
    feasibility = feasibility,
    thresholds = list(
      observed_id_probability = support_threshold,
      unknown_max_probability = unknown_threshold,
      target_specificity = target_specificity
    )
  )
}
