#' Assess repeated-measures RNA-seq sample annotations
#'
#' The expression matrix must contain transformed, finite values with features
#' in rows and samples in columns. Each enrolled individual needs at least two
#' samples. Samples are scored out of fold: an evaluated sample is excluded
#' from its own fitted reference. This initial implementation uses all supplied
#' features and does not normalize or select features.
#'
#' Evidence is expressed as margins between raw shrinkage-LDA scores: the
#' top-versus-second margin and the annotation margin (top versus observed
#' label). No probabilities are produced. Annotation margins are not
#' probabilities that an observed annotation is correct. The function therefore
#' leaves
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
#'   annotation margin threshold.
#' @param assess_feasibility Whether to return task-specific feasibility
#'   assessments based on internal validation.
#' @param progress Whether to report progress updates. Defaults to `TRUE` in
#'   interactive sessions and `FALSE` otherwise.
#' @param cv_folds `NULL` (default) uses leave-one-sample-out validation;
#'   otherwise, an integer from 2 to the smallest individual sample count
#'   selects the number of stratified cross-validation folds. Fold assignment
#'   is randomized; call `set.seed()` beforehand for reproducible results.
#' @param min_agreement Minimum fraction of samples whose out-of-fold top
#'   predicted identity equals the observed identity. Below this, the function
#'   stops before calibration.
#' @param margin_coverage Fraction of correctly predicted samples that must
#'   fall at or below the top-versus-second raw score margin threshold. Samples
#'   predicted incorrectly with an annotation margin (top versus observed
#'   score) above this threshold are treated as confidently misassigned and
#'   excluded from the threshold reference set.
#' @param refit_clean Whether to refit the out-of-fold models after removing
#'   confidently misassigned samples from training, and use those second-pass
#'   scores for calibration and all downstream metrics. Individuals with fewer
#'   than two remaining clean samples are dropped (with all their samples) from
#'   the second pass and listed in `calibration$dropped_individuals`.
#'   `calibration$first_pass_agreement` and `calibration$oof_agreement` report
#'   agreement before and after the refit; the minimum-agreement gate applies
#'   to the first pass.
#' @param mismatch_rate Fraction of clean samples whose labels are swapped
#'   pairwise in a simulation that gives known-truth agreement for the
#'   risk-coverage curve. `NULL` or 0 disables the simulation.
#' @param min_accuracy Target agreement among samples that receive a
#'   prediction. The abstention cutoff is the lowest top-versus-second margin
#'   at which simulated true-label agreement (observed agreement when the
#'   simulation is off) among retained samples reaches this value. Abstaining
#'   samples get `predicted_id = NA`, are never flagged as suspicious and have
#'   status `"no_confident_prediction"`; if the target is never reached, all
#'   samples abstain.
#' @param unknown_mad Samples whose top margin is below
#'   `median - unknown_mad * MAD` of the log top margin of enrolled samples
#'   abstain and are marked as potential unknown individuals. No share of samples is forced to be unknown.
#' @return A list with sample-level `predictions`, raw `identity_scores`,
#'   `performance`, `calibration`, `risk_coverage` (see
#'   [plot_risk_coverage()]), `mismatch_simulation`, `held_out` (count, fraction and IDs of samples held out as
#'   confidently misassigned, from the first pass), `feasibility`, and
#'   `thresholds` components.
#' @export
check_annotations <- function(expr, metadata, sample_col, id_col,
                              shrinkage = 0.1, target_specificity = 0.95,
                              assess_feasibility = TRUE,
                              progress = interactive(), cv_folds = NULL,
                              min_agreement = 0.7, margin_coverage = 0.95,
                              refit_clean = TRUE, mismatch_rate = 0.1,
                              min_accuracy = 0.9, unknown_mad = 3) {
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
  if (!is.null(cv_folds) &&
      (!is.numeric(cv_folds) || is.complex(cv_folds) ||
       length(cv_folds) != 1L || !is.finite(cv_folds) ||
       cv_folds != floor(cv_folds) || cv_folds < 2L ||
       cv_folds > min(class_counts))) {
    stop("cv_folds must be NULL or an integer from 2 to the smallest individual sample count.")
  }
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("progress must be TRUE or FALSE.")
  }

  x <- t(expr)
  rownames(x) <- colnames(expr)
  oof_scores <- matrix(
    NA_real_, nrow = nrow(x), ncol = length(class_counts),
    dimnames = list(rownames(x), names(class_counts))
  )
  n_folds <- if (is.null(cv_folds)) nrow(x) else as.integer(cv_folds)
  fold_id <- if (is.null(cv_folds)) {
    seq_len(nrow(x))
  } else {
    sample_folds <- integer(nrow(x))
    for (individual in names(class_counts)) {
      individual_samples <- which(y == individual)
      sample_folds[individual_samples] <- sample(rep(
        seq_len(n_folds), length.out = length(individual_samples)
      ))
    }
    sample_folds
  }
  if (isTRUE(progress)) {
    if (is.null(cv_folds)) {
      message(sprintf("Fitting leave-one-sample-out models for %d samples.",
                      nrow(x)))
    } else {
      message(sprintf("Fitting stratified %d-fold CV for %d samples.",
                      n_folds, nrow(x)))
    }
  }
  fold_progress_interval <- max(1L, ceiling(n_folds / 10L))
  for (fold in seq_len(n_folds)) {
    test <- which(fold_id == fold)
    train <- fold_id != fold
    model <- fit_shrinkage_lda(x[train, , drop = FALSE], y[train],
                               shrinkage = shrinkage)
    oof_scores[test, ] <- predict_shrinkage_lda(
      model, x[test, , drop = FALSE]
    )
    if (isTRUE(progress) &&
        (fold %% fold_progress_interval == 0L || fold == n_folds)) {
      if (is.null(cv_folds)) {
        message(sprintf("Leave-one-sample-out models: %d/%d samples.",
                        fold, n_folds))
      } else {
        message(sprintf("Cross-validation folds: %d/%d.", fold, n_folds))
      }
    }
  }

  class_index <- match(y, colnames(oof_scores))
  margins <- .score_margins(oof_scores, y)
  predicted_id <- margins$predicted
  correct_identity <- margins$correct
  oof_agreement <- mean(correct_identity)
  if (isTRUE(progress)) {
    message(sprintf(
      "Out-of-fold agreement between predicted and observed identities: %.1f%% (minimum %.1f%%).",
      100 * oof_agreement, 100 * min_agreement
    ))
  }
  if (oof_agreement < min_agreement) {
    stop(sprintf(
      paste("Out-of-fold identity agreement (%.1f%%) is below the required",
            "minimum (%.1f%%); calibration was not attempted."),
      100 * oof_agreement, 100 * min_agreement
    ))
  }

  top_margin <- margins$top_margin
  annotation_margin <- margins$annotation_margin
  top_margin_threshold <- unname(stats::quantile(
    top_margin[correct_identity], probs = margin_coverage, type = 1
  ))
  confidently_misassigned <- !correct_identity &
    annotation_margin > top_margin_threshold
  if (isTRUE(progress)) {
    message(sprintf(
      "Holding out %d of %d samples (%.1f%%) as confidently misassigned.",
      sum(confidently_misassigned), nrow(x),
      100 * mean(confidently_misassigned)
    ))
  }
  first_pass_agreement <- oof_agreement
  held_out_ids <- rownames(x)[confidently_misassigned]
  n_samples_first_pass <- nrow(x)
  dropped_individuals <- character()
  if (isTRUE(refit_clean) && any(confidently_misassigned)) {
    clean_counts <- table(factor(y[!confidently_misassigned],
                                 levels = names(class_counts)))
    dropped_individuals <- names(clean_counts)[clean_counts < 2L]
    use <- !(y %in% dropped_individuals)
    if (length(setdiff(names(class_counts), dropped_individuals)) < 2L) {
      stop("Fewer than two individuals have at least two clean samples; cannot refit.")
    }
    if (isTRUE(progress)) {
      message(sprintf(
        "Refitting on clean samples (%d individuals dropped: %d samples).",
        length(dropped_individuals), sum(!use)
      ))
    }
    x <- x[use, , drop = FALSE]
    y <- y[use]
    confidently_misassigned <- confidently_misassigned[use]
    class_counts <- table(y)
    oof_scores <- .refit_oof_scores(
      x, y, confidently_misassigned, if (is.null(cv_folds)) NULL else n_folds,
      shrinkage
    )
    class_index <- match(y, colnames(oof_scores))
    margins <- .score_margins(oof_scores, y)
    predicted_id <- margins$predicted
    correct_identity <- margins$correct
    oof_agreement <- mean(correct_identity)
    top_margin <- margins$top_margin
    annotation_margin <- margins$annotation_margin
  }
  if (isTRUE(progress)) {
    message("Setting the annotation margin threshold and summarizing validation metrics.")
  }
  # Confidently misassigned samples are probable misannotations; keeping them
  # in the reference set would inflate the thresholds and mask the errors.
  reference <- !confidently_misassigned
  # Every wrong identity as a false label: margin = top score - that score.
  is_observed <- col(oof_scores) == class_index
  false_margin <- (margins$top_score - oof_scores)[!is_observed]
  margin_threshold <- .select_threshold(
    annotation_margin[reference], false_margin, target_specificity, "upper"
  )
  unknown <- .unknown_validation(
    x, y, shrinkage, exclude = confidently_misassigned, progress = progress
  )
  simulation <- NULL
  clean_counts <- table(y[reference])
  if (!is.null(mismatch_rate) && mismatch_rate > 0 &&
      sum(clean_counts >= 2L) >= 2L) {
    clean <- reference & y %in% names(clean_counts)[clean_counts >= 2L]
    if (isTRUE(progress)) {
      message(sprintf(
        "Simulating label swaps (rate %.0f%%) on %d clean samples.",
        100 * mismatch_rate, sum(clean)
      ))
    }
    simulation <- .simulate_mismatches(
      x[clean, , drop = FALSE], y[clean], mismatch_rate,
      if (is.null(cv_folds)) NULL else n_folds, shrinkage, margin_coverage,
      refit_clean
    )
  }
  risk_coverage <- cbind(
    source = "observed",
    .risk_coverage(top_margin, correct_identity),
    stringsAsFactors = FALSE
  )
  if (!is.null(simulation)) {
    risk_coverage <- rbind(risk_coverage, cbind(
      source = "simulated",
      .risk_coverage(simulation$top_margin, simulation$agreement_simulated,
                     simulation$agreement_real),
      stringsAsFactors = FALSE
    ))
  }

  # Abstention: no prediction when the top margin is below a cutoff (ambiguous
  # between enrolled individuals). The cutoff is the lowest margin at which
  # agreement among retained samples reaches min_accuracy, with simulated
  # true-label agreement preferred because observed labels may be wrong. The
  # unknown cutoff is a robust lower bound on enrolled top margins; unknowns
  # always abstain too.
  known <- if (is.null(simulation)) {
    list(margin = top_margin[reference], correct = correct_identity[reference])
  } else {
    list(margin = simulation$top_margin, correct = simulation$agreement_real)
  }
  abstain_threshold <- .abstention_threshold(
    known$margin, known$correct, min_accuracy
  )$threshold
  unknown_threshold <- .unknown_margin_threshold(top_margin[reference],
                                                 unknown_mad)
  low_margin <- top_margin < max(abstain_threshold, unknown_threshold)
  unknown_flag <- top_margin < unknown_threshold
  confident_prediction <- !low_margin
  abstain_reason_counts <- c(
    ambiguous_margin = sum(low_margin & !unknown_flag),
    potential_unknown_individual = sum(unknown_flag)
  )
  unknown_metrics <- if (unknown$n) {
    .binary_metrics(
      c(rep(FALSE, sum(reference)), rep(TRUE, unknown$n)),
      c(unknown_flag[reference], unknown$top_margin < unknown_threshold),
      -c(top_margin[reference], unknown$top_margin)
    )
  }
  if (isTRUE(progress)) {
    message(sprintf(
      "Abstaining on %d of %d samples (margin < %.3g; %d potential unknown, margin < %.3g) to keep agreement >= %.0f%%.",
      sum(!confident_prediction), length(y), abstain_threshold,
      sum(unknown_flag), unknown_threshold, 100 * min_accuracy
    ))
  }
  # Samples without a confident prediction are never flagged: a small top
  # margin makes disagreement with the observed label uninformative.
  flagged <- confident_prediction & annotation_margin > margin_threshold

  misannotation_metrics <- .binary_metrics(
    c(rep(FALSE, sum(reference)), rep(TRUE, length(false_margin))),
    c(flagged[reference],
      confident_prediction[row(oof_scores)[!is_observed]] &
        false_margin > margin_threshold),
    c(annotation_margin[reference], false_margin)
  )
  # Each sample is scored against every wrong label equally often, so recovery
  # equals overall agreement.
  identity_recovery <- mean(correct_identity)
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
  swap_margin <- margins$top_score[swap_samples] -
    oof_scores[cbind(swap_samples, match(swap_labels, colnames(oof_scores)))]
  pairwise_swap_metrics <- .binary_metrics(
    c(rep(FALSE, sum(reference)), rep(TRUE, length(swap_margin))),
    c(flagged[reference],
      confident_prediction[swap_samples] & swap_margin > margin_threshold),
    c(annotation_margin[reference], swap_margin)
  )
  top_three <- 1L + rowSums(oof_scores > margins$observed_score) <= 3L

  calibration <- list(
    method = "score_margin_threshold",
    margin_threshold = margin_threshold,
    exclusion_margin_threshold = top_margin_threshold,
    abstain_margin_threshold = abstain_threshold,
    prediction_coverage = mean(confident_prediction),
    abstain_reason_counts = abstain_reason_counts,
    min_accuracy = min_accuracy,
    first_pass_agreement = first_pass_agreement,
    oof_agreement = oof_agreement,
    dropped_individuals = dropped_individuals,
    n_confidently_misassigned = sum(confidently_misassigned),
    mismatch_simulation = if (is.null(simulation)) NULL else list(
      rate = mismatch_rate,
      n_samples = nrow(simulation),
      n_misannotated = sum(simulation$actually_misannotated),
      n_confidently_misassigned = sum(simulation$confidently_misassigned),
      agreement_simulated = mean(simulation$agreement_simulated),
      agreement_real = mean(simulation$agreement_real)
    ),
    note = paste(
      "Flags samples whose top-versus-observed raw score margin exceeds a",
      "cutoff that keeps at least the target specificity among reference",
      "samples while maximizing sensitivity to simulated false labels;",
      "not an annotation-correctness probability."
    )
  )
  predictions <- data.frame(
    sample_id = rownames(x),
    observed_id = y,
    predicted_id = ifelse(confident_prediction, predicted_id, NA_character_),
    confident_prediction = confident_prediction,
    observed_id_score = margins$observed_score,
    top_candidate_score = margins$top_score,
    top_vs_observed_score_margin = annotation_margin,
    top_score_margin = top_margin,
    annotation_score_margin = annotation_margin,
    confidently_misassigned = confidently_misassigned,
    annotation_correctness_probability = NA_real_,
    suspicious = flagged,
    potential_unknown_individual = unknown_flag,
    most_likely_alternative = ifelse(
      flagged & predicted_id != y,
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
    unknown_individual_detection = unknown_metrics,
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
    } else if (is.finite(unknown_metrics$sensitivity) &&
               unknown_metrics$sensitivity >= 0.8) {
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
        "performance use the same out-of-fold predictions. Thresholds and",
        "specificity exclude confidently misassigned samples. Margin-based",
        "flagging has no probability interpretation. Margin and unknown",
        "thresholds are tuned on the same simulated positives used to estimate",
        "sensitivity, which makes sensitivity optimistic."
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
    rep(!reliable_assessment, nrow(predictions)),
    "insufficient_evidence",
    ifelse(predictions$potential_unknown_individual & !reliable_unknown,
           "insufficient_evidence",
           ifelse(predictions$potential_unknown_individual,
                  "potential_unknown_individual",
           ifelse(!confident_prediction, "no_confident_prediction",
           ifelse(!flagged, "no_evidence_of_incorrect_annotation",
                  ifelse(!is.na(predictions$candidate_correction),
                         "suspicious_with_candidate_correction",
                         "suspicious_without_reliable_alternative")))))
  )

  list(
    predictions = predictions,
    identity_scores = oof_scores,
    performance = performance,
    calibration = calibration,
    risk_coverage = risk_coverage,
    mismatch_simulation = simulation,
    held_out = list(
      n_confidently_misassigned = length(held_out_ids),
      n_samples = n_samples_first_pass,
      fraction = length(held_out_ids) / n_samples_first_pass,
      sample_ids = held_out_ids
    ),
    feasibility = feasibility,
    thresholds = list(
      annotation_margin = margin_threshold,
      unknown_margin = unknown_threshold,
      abstain_margin = abstain_threshold,
      target_specificity = target_specificity
    )
  )
}
