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
#' @param target_specificity Minimum specificity for each of the reference
#'   and unchanged simulated samples when selecting the annotation threshold.
#'   Also used for the simulated unknown-individual cutoff.
#' @param assess_feasibility Whether to return task-specific feasibility
#'   assessments based on internal validation.
#' @param progress Whether to report progress updates. Defaults to `TRUE` in
#'   interactive sessions and `FALSE` otherwise.
#' @param cv_folds `NULL` (default) uses leave-one-sample-out validation;
#'   otherwise, an integer from 2 to the smallest individual sample count
#'   selects the number of stratified cross-validation folds. Fold assignment
#'   is randomized; call `set.seed()` beforehand for reproducible results.
#' @param min_agreement Minimum fraction of samples whose out-of-fold top
#'   predicted identity equals the observed identity. If agreement is below
#'   this value, the function stops before calibration. Defaults to 0.6;
#'   a warning is issued when agreement is below 0.75.
#' @param margin_coverage Quantile used to set the first-pass holdout cutoff.
#'   The reference values are chosen by `holdout_reference`.
#' @param holdout_reference Whether to use the top-versus-second margin
#'   (`"top"`) or pool the top-versus-each-alternative margins
#'   (`"all_alternatives"`) from correctly predicted samples. Samples predicted
#'   incorrectly with an annotation margin (top versus observed score) above
#'   this threshold are treated as confidently misassigned and excluded from
#'   the threshold reference set. The same choice is applied to simulated
#'   mismatches. The original out-of-fold models are not refit after holdout.
#' @param mismatch_rate Fraction of reference samples whose labels are swapped
#'   pairwise in the simulation. `"held_out"` uses the fraction of all input
#'   samples held out in the first pass, applied to the simulation reference
#'   samples. Swaps are rounded to whole pairs, with at least one pair for
#'   `"held_out"`, even when no samples were held out. The actual simulated
#'   fraction may differ. `NULL` or numeric 0 disables the simulation.
#' @param mismatch_repeats Number of independent label-swap simulations whose
#'   results are pooled to tune the annotation cutoff. With at least two
#'   repeats, detection sensitivity and simulated specificity are cross-fitted:
#'   each repeat is scored with a cutoff tuned on the other repeats' swapped
#'   and unchanged samples and the reference negatives. More repeats give a less extreme
#'   cutoff and a steadier sensitivity estimate at the cost of roughly
#'   proportional simulation time. With 1, sensitivity is in-sample.
#' @param diagnostics Whether to include intermediate gate data in the result.
#'   This is intended for understanding and plotting the validation steps; it
#'   can add substantially to the result size.
#' @param min_accuracy Target agreement for selecting the correction cutoff.
#'   The reference includes simulated swaps and unchanged samples whose
#'   predicted identity differs from their simulated label. Agreeing predictions
#'   are returned in `predicted_id_clean`; disagreements are returned there
#'   only at or above the correction cutoff. Without a usable cutoff,
#'   disagreements have `predicted_id_clean = NA`. `predicted_id` always
#'   retains the raw prediction.
#' @param detect_unknown If `FALSE`, unknown-individual detection is switched
#'   off: no stranger validation is run, no sample is marked as a potential
#'   unknown individual.
#'   `feasibility$unknown_individual_detection` is `"not_assessed"`.
#' @param unknown_threshold_method How the unknown-individual cutoff on the top
#'   statistic is chosen. `"simulated"` (default) uses the Mahalanobis
#'   distance (under the shrinkage covariance) from a sample to its nearest
#'   enrolled individual's mean, with the cutoff chosen so at most
#'   `1 - target_specificity` of enrolled reference samples exceed it while
#'   catching the most simulated strangers; with no simulated strangers it is
#'   the `target_specificity` quantile of the enrolled distances. `"mad"`
#'   restores the earlier rule, a cutoff on the top-versus-second margin of
#'   `median - unknown_mad * MAD` of the log top margin of enrolled samples,
#'   which uses no simulated strangers. Samples beyond the cutoff are marked
#'   as potential unknown individuals and cannot receive candidate corrections;
#'   the flag does not hide `predicted_id_clean`.
#'   No share of samples is forced to be unknown.
#' @param unknown_mad Multiplier used only when
#'   `unknown_threshold_method = "mad"`.
#' @param unknown_stray_samples Unknown-individual validation pretends each
#'   enrolled individual is a stranger. Up to this many of their samples (never
#'   all of them) are kept in the data relabelled as random other individuals,
#'   as a stranger's samples would be in a real cohort, and the rest are
#'   dropped. Each kept sample is scored by a model that never saw it, giving a
#'   margin for an unknown person.
#' @param grouped_unknown_validation By default the stranger simulation is
#'   repeated for each individual separately. If `TRUE`
#'   and `cv_folds` is set, individuals are removed in `cv_folds` groups
#'   instead, so far fewer models are fit, at the cost of slightly optimistic
#'   (too large) unknown margins. With the default
#'   `unknown_threshold_method = "simulated"` it also affects the unknown
#'   cutoff.
#' @return A list with sample-level `predictions`, raw `identity_scores`,
#'   `performance`, `calibration`, `risk_coverage` (see
#'   [plot_risk_coverage()]), `mismatch_simulation`, `held_out` (count, fraction and IDs of samples held out as
#'   confidently misassigned, from the first pass), `feasibility`, and
#'   `thresholds` components. `predictions$potential_unknown_individual` is a
#'   separate flag, not an `annotation_status` value. Predictions also include
#'   the nearest-class distance and the unknown-detection score and threshold.
#'   `predictions$annotation_status` is descriptive (flag, confidence,
#'   agreement) and does not depend on `feasibility`.
#'   All samples, including held-out samples, are flagged as suspicious only
#'   when their annotation margin exceeds the calibrated threshold. Every input
#'   sample has a row.
#'   `performance$misannotation_detection` includes simulated-swap sensitivity,
#'   `simulated_specificity` on unchanged simulation samples and `specificity`
#'   on the original reference samples. The detection verdict uses simulated
#'   sensitivity and simulated specificity; without simulated specificity it
#'   is `"insufficient_evidence"`. With one repeat both rates are in-sample.
#'   Precision, accuracy and ranking metrics in that table use the reference
#'   negatives and final cutoff.
#'   `performance$estimated_detection_errors` estimates misannotation
#'   prevalence from the observed suspicious-flag rate and simulated swap
#'   sensitivity and specificity. It projects false positives and missed
#'   misannotations at that estimate and at +/- 2.5 percentage points,
#'   bounded to 0--100%. These are exploratory projections, not confidence
#'   intervals or sample-level probabilities; simulated swaps may not
#'   represent real errors or unknown individuals. It is unavailable when
#'   simulated rates are missing, uninformative, or incompatible with the
#'   observed flag rate. Set `progress = TRUE` to display the summary.
#'   `performance$all_alternatives_detection`
#'   keeps the every-wrong-label diagnostic, and
#'   `performance$recovery_accuracy_misannotated` and
#'   `performance$correction_validation` report recovery on simulated
#'   misannotated samples and the precision and recall of proposed corrections.
#'   The correction cutoff is the lowest top margin reaching `min_accuracy`
#'   on that filtered simulation reference. It includes swapped samples even
#'   when the prediction equals their simulated label, so this accuracy is not
#'   the precision of proposed corrections. The result records its threshold,
#'   retained count, coverage and in-sample accuracy in `calibration$correction_margin`.
#'   Identity recovery appears feasible when that cutoff exists and the input
#'   has at least 20 samples; without swaps it is insufficient evidence.
#'   Candidates must also pass the unknown and suspicious gates. Simulated
#'   correction validation applies the margin-based unknown gate for `"mad"`,
#'   but not the default distance-based unknown gate because simulated swap
#'   distances are not recorded; it also does not apply feasibility. There is
#'   no minimum retained simulated count; even one correct row
#'   can meet the target. Repeated rows are dependent, and selecting and
#'   evaluating the cutoff on the same simulations can be optimistic.
#' @export
check_annotations <- function(expr, metadata, sample_col, id_col,
                              shrinkage = 0.1, target_specificity = 0.95,
                              assess_feasibility = TRUE,
                              progress = interactive(), cv_folds = NULL,
                              min_agreement = 0.6, margin_coverage = 0.025,
                              holdout_reference = c("all_alternatives", "top"),
                              mismatch_rate = 0.1,
                              min_accuracy = 0.9, detect_unknown = TRUE,
                              unknown_threshold_method = c("simulated", "mad"),
                              unknown_mad = 3,
                              grouped_unknown_validation = FALSE,
                              unknown_stray_samples = 2L,
                              mismatch_repeats = 3L,
                              diagnostics = FALSE) {
  unknown_threshold_method <- match.arg(unknown_threshold_method)
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
  if (!is.numeric(min_agreement) || is.complex(min_agreement) ||
      length(min_agreement) != 1L || !is.finite(min_agreement) ||
      min_agreement < 0 || min_agreement > 1) {
    stop("min_agreement must be a single number between 0 and 1.")
  }
  if (!is.numeric(margin_coverage) || length(margin_coverage) != 1L ||
      !is.finite(margin_coverage) || margin_coverage < 0 ||
      margin_coverage > 1) {
    stop("margin_coverage must be a single number between 0 and 1.")
  }
  holdout_reference <- match.arg(holdout_reference)
  if (!is.null(mismatch_rate) &&
      !(is.character(mismatch_rate) && length(mismatch_rate) == 1L &&
        !is.na(mismatch_rate) && mismatch_rate == "held_out") &&
      !(is.numeric(mismatch_rate) && !is.complex(mismatch_rate) &&
        length(mismatch_rate) == 1L && is.finite(mismatch_rate) &&
        mismatch_rate >= 0 && mismatch_rate <= 1)) {
    stop('mismatch_rate must be NULL, "held_out", or a single number between 0 and 1.')
  }
  if (!is.numeric(mismatch_repeats) || length(mismatch_repeats) != 1L ||
      !is.finite(mismatch_repeats) || mismatch_repeats < 1L ||
      mismatch_repeats != floor(mismatch_repeats)) {
    stop("mismatch_repeats must be a positive integer.")
  }
  mismatch_repeats <- as.integer(mismatch_repeats)
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
  if (!is.logical(diagnostics) || length(diagnostics) != 1L ||
      is.na(diagnostics)) {
    stop("diagnostics must be TRUE or FALSE.")
  }

  x <- t(expr)
  rownames(x) <- colnames(expr)
  oof_scores <- matrix(
    NA_real_, nrow = nrow(x), ncol = length(class_counts),
    dimnames = list(rownames(x), names(class_counts))
  )
  nearest_distance <- rep(NA_real_, nrow(x))
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
    nearest_distance[test] <- .nearest_class_distance(
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
      "Run stopped because first-pass agreement (%.1f%%) is below min_agreement (%.1f%%).",
      100 * oof_agreement, 100 * min_agreement
    ), call. = FALSE)
  }
  if (oof_agreement < 0.75) {
    warning(sprintf(
      paste("First-pass agreement (%.1f%%) is below 75%%; calibration may be",
            "off and results may be unreliable."),
      100 * oof_agreement
    ), call. = FALSE)
  }
  low_agreement <- FALSE

  top_margin <- margins$top_margin
  annotation_margin <- margins$annotation_margin
  top_margin_threshold <- .holdout_threshold(
    oof_scores, correct_identity, margin_coverage, holdout_reference
  )
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
  mismatch_rate_held_out <- identical(mismatch_rate, "held_out")
  simulation_rate <- if (mismatch_rate_held_out) {
    length(held_out_ids) / n_samples_first_pass
  } else mismatch_rate
  first_pass_diagnostics <- if (diagnostics) {
    data.frame(
      sample_id = rownames(x),
      observed_id = y,
      predicted_id = predicted_id,
      correct = correct_identity,
      top_margin = top_margin,
      annotation_margin = annotation_margin,
      exclusion_margin_threshold = top_margin_threshold,
      held_out = confidently_misassigned,
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }
  first_pass_scores <- if (diagnostics) oof_scores else NULL
  if (isTRUE(progress)) {
    message("Setting the annotation margin threshold and summarizing validation metrics.")
  }
  # Confidently misassigned samples are probable misannotations; keeping them
  # in the reference set would inflate the thresholds and mask the errors.
  reference <- !confidently_misassigned
  # Every wrong identity as a false label: margin = top score - that score.
  is_observed <- col(oof_scores) == class_index
  false_margin <- (margins$top_score - oof_scores)[!is_observed]
  validate_unknown <- function() .unknown_validation(
    x, y, shrinkage, exclude = confidently_misassigned, progress = progress,
    n_strays = as.integer(unknown_stray_samples),
    n_folds = if (is.null(cv_folds) || !isTRUE(grouped_unknown_validation)) {
      NULL
    } else {
      n_folds
    }
  )
  unknown <- if (!detect_unknown) {
    list(n = 0L, top_margin = numeric())
  } else {
    validate_unknown()
  }
  simulation <- NULL
  clean_counts <- table(y[reference])
  if (!is.null(simulation_rate) &&
      (simulation_rate > 0 || mismatch_rate_held_out) &&
      sum(clean_counts >= 2L) >= 2L) {
    clean <- reference & y %in% names(clean_counts)[clean_counts >= 2L]
    if (isTRUE(progress)) {
      message(sprintf(
        "Simulating label swaps (rate %.0f%%) on %d clean samples.",
        100 * simulation_rate, sum(clean)
      ))
    }
    draws <- lapply(seq_len(mismatch_repeats), function(r) {
      if (isTRUE(progress) && mismatch_repeats > 1L) {
        message(sprintf("Label-swap simulation %d/%d.", r, mismatch_repeats))
      }
      draw <- .simulate_mismatches(
        x[clean, , drop = FALSE], y[clean], simulation_rate,
        if (is.null(cv_folds)) NULL else n_folds, shrinkage, margin_coverage,
        holdout_reference
      )
      if (!is.null(draw)) draw$repeat_id <- r
      draw
    })
    draws <- draws[!vapply(draws, is.null, logical(1))]
    if (length(draws)) {
      simulation <- do.call(rbind, draws)
    }
  }
  swap_margin_positives <- if (!is.null(simulation)) {
    simulation$annotation_margin[simulation$actually_misannotated]
  } else {
    numeric()
  }
  swap_margin_negatives <- if (!is.null(simulation)) {
    simulation$annotation_margin[!simulation$actually_misannotated]
  } else {
    numeric()
  }
  margin_threshold <- .select_threshold(
    annotation_margin[reference],
    if (length(swap_margin_positives)) swap_margin_positives else false_margin,
    target_specificity, "upper", additional_negatives = swap_margin_negatives
  )
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

  correction_threshold <- .correction_threshold(simulation, min_accuracy)
  confident_prediction <- .confident_prediction(
    predicted_id, y, top_margin, correction_threshold$threshold
  )
  unknown_threshold <- if (!detect_unknown) {
    -Inf
  } else if (unknown_threshold_method == "mad") {
    .unknown_margin_threshold(top_margin[reference], unknown_mad)
  } else {
    .select_threshold(nearest_distance[reference], unknown$distance,
                      target_specificity, "upper")
  }
  # "mad" flags a low top margin; "simulated" flags a large distance to the
  # nearest enrolled individual.
  use_distance <- detect_unknown && unknown_threshold_method == "simulated"
  unknown_statistic <- if (use_distance) nearest_distance else top_margin
  stranger_statistic <- if (use_distance) unknown$distance else
    unknown$top_margin
  is_unknown <- function(statistic) {
    if (!detect_unknown) rep(FALSE, length(statistic))
    else if (use_distance) statistic > unknown_threshold
    else statistic < unknown_threshold
  }
  unknown_flag <- is_unknown(unknown_statistic)
  unknown_metrics <- if (unknown$n) {
    .binary_metrics(
      c(rep(FALSE, sum(reference)), rep(TRUE, unknown$n)),
      c(unknown_flag[reference], is_unknown(stranger_statistic)),
      (if (use_distance) 1 else -1) *
        c(unknown_statistic[reference], stranger_statistic)
    )
  }
  if (isTRUE(progress)) {
    message(sprintf(
      "%d of %d samples lack a clean prediction (disagreement below correction margin %.3g); %d potential unknown (%s %.3g).",
      sum(!confident_prediction), length(y), correction_threshold$threshold,
      sum(unknown_flag), if (use_distance) "distance >" else "margin <",
      unknown_threshold
    ))
  }
  # Flagging depends on the annotation margin alone: a clearly beaten observed
  # label is suspicious even when the top margin is too small to name the
  # replacement, which needs a confident prediction.
  flagged <- annotation_margin > margin_threshold

  # Every wrong label as a false label: a broad robustness diagnostic that
  # also includes easily confused identities real swaps rarely involve.
  all_alternatives_metrics <- .binary_metrics(
    c(rep(FALSE, sum(reference)), rep(TRUE, length(false_margin))),
    c(flagged[reference], false_margin > margin_threshold),
    c(annotation_margin[reference], false_margin)
  )
  # Estimate simulated sensitivity and specificity on held-out repeats where possible.
  swap_simulation <- !is.null(simulation) && any(simulation$actually_misannotated)
  if (swap_simulation) {
    swap_margins <- simulation$annotation_margin[simulation$actually_misannotated]
    misannotation_metrics <- .binary_metrics(
      c(rep(FALSE, sum(reference)), rep(TRUE, length(swap_margins))),
      c(flagged[reference], swap_margins > margin_threshold),
      c(annotation_margin[reference], swap_margins)
    )
    crossfit <- .crossfit_swap_detection(
      annotation_margin[reference], simulation, target_specificity
    )
    detection_source <- if (is.null(crossfit)) {
      "simulated_swaps_in_sample"
    } else {
      "simulated_swaps_cross_fitted"
    }
  } else {
    misannotation_metrics <- all_alternatives_metrics
    crossfit <- NULL
    detection_source <- "all_wrong_identities"
  }
  misannotation_metrics$sensitivity_in_sample <- misannotation_metrics$sensitivity
  sim_negatives <- swap_margin_negatives
  sim_negatives <- sim_negatives[is.finite(sim_negatives)]
  misannotation_metrics$simulated_specificity_in_sample <- if (length(sim_negatives)) {
    mean(sim_negatives <= margin_threshold)
  } else NA_real_
  misannotation_metrics$simulated_specificity <- if (is.null(crossfit)) {
    misannotation_metrics$simulated_specificity_in_sample
  } else crossfit$simulated_specificity
  misannotation_metrics$n_simulated_negatives <- if (is.null(crossfit)) {
    length(sim_negatives)
  } else crossfit$n_simulated_negatives
  if (!is.null(crossfit)) {
    misannotation_metrics$sensitivity <- crossfit$sensitivity
    misannotation_metrics$false_negative_rate <- 1 - crossfit$sensitivity
  }
  misannotation_metrics$sensitivity_sd <- if (is.null(crossfit)) {
    NA_real_
  } else {
    crossfit$sensitivity_sd
  }
  misannotation_metrics$n_repeats_validated <- if (is.null(crossfit)) {
    0L
  } else {
    crossfit$n_repeats
  }
  misannotation_metrics$source <- detection_source
  estimated_detection_errors <- .estimate_detection_errors(
    flagged, misannotation_metrics$sensitivity,
    misannotation_metrics$simulated_specificity, detection_source
  )
  if (isTRUE(progress)) {
    message(paste(.format_detection_error_estimate(estimated_detection_errors),
                  collapse = "\n"))
  }
  # Overall agreement: each sample is scored against every wrong label equally
  # often, so this equals the all-alternatives recovery.
  identity_recovery <- mean(correct_identity)
  # Recovery on the samples that actually need a correction: simulated swaps
  # whose top prediction equals the original label, and the quality of the
  # corrections that would be proposed.
  recovery_misannotated <- if (swap_simulation) {
    mean(simulation$agreement_real[simulation$actually_misannotated])
  } else {
    NA_real_
  }
  correction_validation <- if (!is.null(simulation)) {
    proposed <- simulation$annotation_margin > margin_threshold &
      simulation$top_margin >= max(
        correction_threshold$threshold,
        if (use_distance || !detect_unknown) -Inf else unknown_threshold
      ) &
      !simulation$agreement_simulated
    swapped <- simulation$actually_misannotated
    list(
      n_proposed = sum(proposed),
      precision = .safe_divide(sum(simulation$agreement_real[proposed]),
                               sum(proposed)),
      recall = .safe_divide(sum(simulation$agreement_real[proposed & swapped]),
                            sum(swapped))
    )
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
  swap_margin <- margins$top_score[swap_samples] -
    oof_scores[cbind(swap_samples, match(swap_labels, colnames(oof_scores)))]
  pairwise_swap_metrics <- .binary_metrics(
    c(rep(FALSE, sum(reference)), rep(TRUE, length(swap_margin))),
    c(flagged[reference], swap_margin > margin_threshold),
    c(annotation_margin[reference], swap_margin)
  )
  top_three <- 1L + rowSums(oof_scores > margins$observed_score) <= 3L

  calibration <- list(
    method = "score_margin_threshold",
    margin_threshold = margin_threshold,
    exclusion_margin_threshold = top_margin_threshold,
    holdout_reference = holdout_reference,
    correction_margin = correction_threshold,
    prediction_coverage = mean(confident_prediction),
    min_accuracy = min_accuracy,
    first_pass_agreement = first_pass_agreement,
    low_agreement = low_agreement,
    oof_agreement = oof_agreement,
    dropped_individuals = character(),
    n_confidently_misassigned = sum(confidently_misassigned),
    mismatch_simulation = if (is.null(simulation)) NULL else list(
      rate = simulation_rate,
      rate_mode = if (mismatch_rate_held_out) "held_out" else "fixed",
      actual_rate = mean(simulation$actually_misannotated),
      repeats = length(unique(simulation$repeat_id)),
      n_samples = nrow(simulation),
      n_misannotated = sum(simulation$actually_misannotated),
      n_confidently_misassigned = sum(simulation$confidently_misassigned),
      agreement_simulated = mean(simulation$agreement_simulated),
      agreement_real = mean(simulation$agreement_real)
    ),
    note = paste(
      "Flags samples whose top-versus-observed raw score margin exceeds a",
      "cutoff that keeps at least the target specificity separately among",
      "reference and unchanged simulated samples, when available, while",
      "maximizing sensitivity to simulated false labels;",
      "not an annotation-correctness probability."
    )
  )
  predictions <- data.frame(
    sample_id = rownames(x),
    observed_id = y,
    predicted_id = predicted_id,
    predicted_id_clean = ifelse(confident_prediction, predicted_id,
                                NA_character_),
    confident_prediction = confident_prediction,
    observed_id_score = margins$observed_score,
    top_candidate_score = margins$top_score,
    top_vs_observed_score_margin = annotation_margin,
    top_score_margin = top_margin,
    annotation_score_margin = annotation_margin,
    nearest_class_distance = nearest_distance,
    unknown_detection_score = if (detect_unknown) unknown_statistic else
      rep(NA_real_, length(y)),
    unknown_detection_threshold = if (detect_unknown) {
      rep(unknown_threshold, length(y))
    } else {
      rep(NA_real_, length(y))
    },
    confidently_misassigned = confidently_misassigned,
    annotation_correctness_probability = NA_real_,
    suspicious = flagged,
    potential_unknown_individual = unknown_flag,
    most_likely_alternative = ifelse(
      flagged & predicted_id != y,
      predicted_id, NA_character_
    ),
    candidate_correction = NA_character_,
    excluded_from_refit = FALSE,
    stringsAsFactors = FALSE
  )
  performance <- list(
    correct_identity_recovery_accuracy = mean(correct_identity),
    top_three_identity_recovery_accuracy = mean(top_three),
    misannotation_detection = misannotation_metrics,
    estimated_detection_errors = estimated_detection_errors,
    all_alternatives_detection = all_alternatives_metrics,
    identity_recovery_accuracy = identity_recovery,
    recovery_accuracy_misannotated = recovery_misannotated,
    correction_validation = correction_validation,
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
    detection_status <- if (!enough_known || !swap_simulation ||
                            !is.finite(misannotation_metrics$simulated_specificity) ||
                            !is.finite(misannotation_metrics$sensitivity)) {
      "insufficient_evidence"
    } else if (is.finite(misannotation_metrics$sensitivity) &&
               misannotation_metrics$sensitivity >= 0.8 &&
               misannotation_metrics$simulated_specificity >= target_specificity) {
      "appears_feasible"
    } else {
      "unreliable"
    }
    recovery_status <- if (!enough_known || !swap_simulation ||
                           correction_threshold$n_eligible == 0L) {
      "insufficient_evidence"
    } else if (is.finite(correction_threshold$threshold)) {
      "appears_feasible"
    } else {
      "unreliable"
    }
    unknown_status <- if (!detect_unknown) {
      "not_assessed"
    } else if (unknown$n < 10L) {
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
        "performance use the same out-of-fold predictions. Reference thresholds",
        "and specificity exclude confidently misassigned samples. Margin-based",
        "flagging has no probability interpretation. Detection sensitivity and",
        "simulated specificity are cross-fitted by repeat when at least two are",
        "available (see performance$misannotation_detection$source); the",
        "unknown cutoff (unless unknown_threshold_method = \"mad\") is",
        "tuned on the same simulated strangers used to estimate unknown",
        "sensitivity, which makes that sensitivity optimistic. Correction",
        "cutoff selection and accuracy use the same simulated swaps, which",
        "can make the accuracy optimistic, especially with few retained rows."
      )
    )
    # Corrections need validated identity recovery; the suspicious flag itself
    # does not depend on feasibility.
    reliable_correction <- recovery_status == "appears_feasible"
  } else {
    reliable_correction <- FALSE
  }
  predictions$candidate_correction <- ifelse(
      flagged & reliable_correction &
        confident_prediction & !unknown_flag &
        predicted_id != y,
      predicted_id,
      NA_character_
  )
  # Descriptive only: feasibility is reported in `feasibility`, not folded in.
  predictions$annotation_status <- ifelse(
    !is.na(predictions$candidate_correction),
           "suspicious_with_candidate_correction",
    ifelse(flagged, "suspicious_without_reliable_alternative",
    ifelse(!confident_prediction & predicted_id != y,
           "no_confident_prediction",
    ifelse(confident_prediction & predicted_id != y,
           "weak_disagreement",
           "no_evidence_of_incorrect_annotation"))))
  result <- list(
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
      unknown_margin = if (use_distance || !detect_unknown) -Inf else
        unknown_threshold,
      unknown_distance = if (use_distance) unknown_threshold else Inf,
      correction_margin = correction_threshold$threshold,
      target_specificity = target_specificity
    )
  )
  if (diagnostics) {
    false_cells <- which(!is_observed, arr.ind = TRUE)
    false_identity_margins <- data.frame(
      sample_id = rownames(oof_scores)[false_cells[, 1L]],
      alternative_id = colnames(oof_scores)[false_cells[, 2L]],
      margin = margins$top_score[false_cells[, 1L]] -
        oof_scores[false_cells],
      stringsAsFactors = FALSE
    )
    simulated_positive_margins <- if (!is.null(simulation)) {
      simulated <- simulation[simulation$actually_misannotated, ,
                              drop = FALSE]
      data.frame(
        sample_id = simulated$sample_id,
        identity = simulated$simulated_id,
        population = "simulated_swap_positive",
        margin = simulated$annotation_margin,
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
    annotation_threshold_data <- data.frame(
      sample_id = rownames(oof_scores)[reference],
      identity = y[reference],
      population = "reference_negative",
      margin = annotation_margin[reference],
      stringsAsFactors = FALSE
    )
    if (!is.null(simulated_positive_margins) &&
        nrow(simulated_positive_margins)) {
      simulated_negatives <- simulation[!simulation$actually_misannotated, ,
                                        drop = FALSE]
      simulated_negative_margins <- data.frame(
        sample_id = simulated_negatives$sample_id,
        identity = simulated_negatives$simulated_id,
        population = rep("simulated_unchanged_negative", nrow(simulated_negatives)),
        margin = simulated_negatives$annotation_margin,
        stringsAsFactors = FALSE
      )
      annotation_threshold_data <- rbind(
        annotation_threshold_data, simulated_positive_margins,
        simulated_negative_margins
      )
    } else {
      false_identity_margins$population <- "all_wrong_identity_alternatives"
      names(false_identity_margins)[names(false_identity_margins) ==
                                     "alternative_id"] <- "identity"
      annotation_threshold_data <- rbind(
        annotation_threshold_data,
        false_identity_margins[, names(annotation_threshold_data),
                               drop = FALSE]
      )
    }
    annotation_threshold_data$threshold <- margin_threshold
    annotation_threshold_data$threshold_source <- if (
      !is.null(simulated_positive_margins) &&
      nrow(simulated_positive_margins)
    ) "simulated_swaps" else "all_wrong_identities"

    final_diagnostics <- data.frame(
      sample_id = rownames(x),
      observed_id = y,
      predicted_id = predicted_id,
      correct = correct_identity,
      top_margin = top_margin,
      annotation_margin = annotation_margin,
      stringsAsFactors = FALSE
    )
    first_holdout_by_id <- setNames(first_pass_diagnostics$held_out,
                                    first_pass_diagnostics$sample_id)
    final_diagnostics$confidently_misassigned_first_pass <-
      unname(first_holdout_by_id[final_diagnostics$sample_id])
    final_diagnostics$confident_prediction <- predictions$confident_prediction
    final_diagnostics$suspicious <- predictions$suspicious
    final_diagnostics$potential_unknown_individual <-
      predictions$potential_unknown_individual
    final_diagnostics$candidate_correction <- predictions$candidate_correction
    final_diagnostics$annotation_status <- predictions$annotation_status

    correction_eligible <- if (is.null(simulation)) logical() else
      simulation$actually_misannotated | !simulation$agreement_simulated
    correction_curve <- if (any(correction_eligible)) {
      .risk_coverage(simulation$top_margin[correction_eligible],
                     simulation$agreement_real[correction_eligible])
    } else NULL
    result$diagnostics <- list(
      first_pass = first_pass_diagnostics,
      first_pass_scores = first_pass_scores,
      final_pass = final_diagnostics,
      annotation_threshold = annotation_threshold_data,
      annotation_threshold_source = if (
        !is.null(simulated_positive_margins) &&
        nrow(simulated_positive_margins)
      ) "simulated_swaps" else "all_wrong_identities",
      correction = list(threshold = correction_threshold$threshold,
                        min_accuracy = min_accuracy,
                        curve = correction_curve),
      unknown = list(
        enabled = detect_unknown,
        threshold = unknown_threshold,
        statistic = if (use_distance) "nearest_class_distance" else
          "top_margin",
        known_values = unknown_statistic[reference],
        simulated_unknown_values = stranger_statistic,
        metrics = unknown_metrics
      ),
      gates = list(
        min_agreement = min_agreement,
        first_pass_agreement = first_pass_agreement,
        low_agreement = low_agreement,
        target_specificity = target_specificity,
        feasibility = feasibility,
        performance = performance
      )
    )
  }
  result
}
