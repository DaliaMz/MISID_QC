# Manual, step-by-step version of check_annotations() (R/check_annotations.R).
#
# Run it interactively from the repository root, section by section, and
# inspect each intermediate object. Every variable keeps the name it has inside
# check_annotations(), so the decisions can be compared line by line.
#
# Before running, define `expr` (features x samples, finite, transformed) and
# `metadata` (data frame) in your session, or edit the INPUTS section.

# ---- 0. Setup ---------------------------------------------------------------
source("R/lda.R")
source("R/validation.R")

# ---- 1. Inputs and parameters ----------------------------------------------
# expr <- readRDS("path/to/expr.rds")
# metadata <- read.csv("path/to/metadata.csv")
sample_col <- "sample_id"
id_col <- "individual_id"

shrinkage <- 0.1
target_specificity <- 0.95
cv_folds <- NULL            # NULL = leave-one-sample-out
min_agreement <- 0.7
margin_coverage <- 0.025
holdout_reference <- "all_alternatives" # or "top"
refit_clean <- TRUE
mismatch_rate <- 0.1
mismatch_repeats <- 3L
min_accuracy <- 0.9
detect_unknown <- TRUE
unknown_threshold_method <- "simulated"  # or "mad"
unknown_mad <- 3          # only used with "mad"
grouped_unknown_validation <- FALSE
unknown_stray_samples <- 2L
stop_on_low_agreement <- TRUE
progress <- TRUE
set.seed(123)

# ---- 2. Validate and align inputs ------------------------------------------
expr <- as.matrix(expr)
storage.mode(expr) <- "double"
stopifnot(length(dim(expr)) == 2L, all(is.finite(expr)))
stopifnot(!is.null(colnames(expr)), !anyDuplicated(colnames(expr)))
stopifnot(is.data.frame(metadata), all(c(sample_col, id_col) %in% names(metadata)))
sample_ids <- as.character(metadata[[sample_col]])
observed_ids <- as.character(metadata[[id_col]])
stopifnot(!anyNA(sample_ids), !anyDuplicated(sample_ids), !anyNA(observed_ids))
stopifnot(setequal(sample_ids, colnames(expr)))
metadata <- metadata[match(colnames(expr), sample_ids), , drop = FALSE]
y <- as.character(metadata[[id_col]])
class_counts <- table(y)
stopifnot(length(class_counts) >= 2L, all(class_counts >= 2L))
if (!is.null(cv_folds)) stopifnot(cv_folds >= 2L, cv_folds <= min(class_counts))
print(class_counts)

# ---- 3. First-pass out-of-fold identity scores -----------------------------
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
for (fold in seq_len(n_folds)) {
  test <- which(fold_id == fold)
  train <- fold_id != fold
  model <- fit_shrinkage_lda(x[train, , drop = FALSE], y[train],
                             shrinkage = shrinkage)
  oof_scores[test, ] <- predict_shrinkage_lda(model, x[test, , drop = FALSE])
  nearest_distance[test] <- .nearest_class_distance(
    model, x[test, , drop = FALSE]
  )
  if (progress && (fold %% max(1L, ceiling(n_folds / 10L)) == 0L ||
                   fold == n_folds)) {
    message(sprintf("Fold %d/%d", fold, n_folds))
  }
}

# ---- 4. Margins and the minimum-agreement gate -----------------------------
class_index <- match(y, colnames(oof_scores))
margins <- .score_margins(oof_scores, y)
predicted_id <- margins$predicted
correct_identity <- margins$correct
oof_agreement <- mean(correct_identity)
if (oof_agreement < 0.6 && stop_on_low_agreement) {
  stop(sprintf(paste("Run stopped because first-pass agreement (%.1f%%) is",
                     "below 60%%. To proceed anyway, set",
                     "stop_on_low_agreement = FALSE."), 100 * oof_agreement))
}
if (oof_agreement < 0.75) {
  warning(sprintf(paste("First-pass agreement (%.1f%%) is below 75%%;",
                        "calibration may be off and results may be",
                        "unreliable."), 100 * oof_agreement))
}
low_agreement <- oof_agreement < min_agreement
top_margin <- margins$top_margin
annotation_margin <- margins$annotation_margin
cat(sprintf("First-pass agreement %.3f (min %.3f) -> low_agreement = %s\n",
            oof_agreement, min_agreement, low_agreement))

# ---- 5. First-pass hold-out of confidently misassigned samples -------------
top_margin_threshold <- SampleAnnotationLikelihood:::.holdout_threshold(
  oof_scores, correct_identity, margin_coverage, holdout_reference
)
confidently_misassigned <- !correct_identity &
  annotation_margin > top_margin_threshold
first_pass_agreement <- oof_agreement
held_out_ids <- rownames(x)[confidently_misassigned]
first_pass_all <- data.frame(
  sample_id = rownames(x), observed_id = y, predicted_id = predicted_id,
  observed_id_score = margins$observed_score,
  top_candidate_score = margins$top_score,
  top_margin = top_margin, annotation_margin = annotation_margin,
  confidently_misassigned = confidently_misassigned,
  stringsAsFactors = FALSE
)
first_pass_oof_scores <- oof_scores
cat(sprintf("Exclusion threshold %.4g; holding out %d of %d samples\n",
            top_margin_threshold, sum(confidently_misassigned), nrow(x)))
print(first_pass_all[confidently_misassigned, ])

# ---- 6. Refit on clean samples ---------------------------------------------
dropped_individuals <- character()
if (refit_clean && any(confidently_misassigned) && !low_agreement) {
  clean_counts <- table(factor(y[!confidently_misassigned],
                               levels = names(class_counts)))
  dropped_individuals <- names(clean_counts)[clean_counts < 2L]
  use <- !(y %in% dropped_individuals)
  stopifnot(length(setdiff(names(class_counts), dropped_individuals)) >= 2L)
  x <- x[use, , drop = FALSE]
  y <- y[use]
  confidently_misassigned <- confidently_misassigned[use]
  class_counts <- table(y)
  oof_scores <- .refit_oof_scores(
    x, y, confidently_misassigned, if (is.null(cv_folds)) NULL else n_folds,
    shrinkage, progress = progress
  )
  nearest_distance <- attr(oof_scores, "nearest_distance")
  class_index <- match(y, colnames(oof_scores))
  margins <- .score_margins(oof_scores, y)
  predicted_id <- margins$predicted
  correct_identity <- margins$correct
  oof_agreement <- mean(correct_identity)
  top_margin <- margins$top_margin
  annotation_margin <- margins$annotation_margin
}
cat("Dropped individuals:", dropped_individuals, "\n")
cat(sprintf("Post-refit agreement %.3f\n", oof_agreement))

# ---- 7. Reference set and false-label margins ------------------------------
reference <- !confidently_misassigned
is_observed <- col(oof_scores) == class_index
false_margin <- (margins$top_score - oof_scores)[!is_observed]

# ---- 8. Simulated unknown individuals --------------------------------------
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
} else if (low_agreement) {
  tryCatch(validate_unknown(), error = function(e) {
    list(n = 0L, top_margin = numeric(), status = "insufficient_evidence")
  })
} else {
  validate_unknown()
}
str(unknown)

# ---- 9. Simulated label swaps ----------------------------------------------
simulation <- NULL
clean_counts <- table(y[reference])
if (!low_agreement && !is.null(mismatch_rate) && mismatch_rate > 0 &&
    sum(clean_counts >= 2L) >= 2L) {
  clean <- reference & y %in% names(clean_counts)[clean_counts >= 2L]
  draws <- lapply(seq_len(mismatch_repeats), function(r) {
    draw <- .simulate_mismatches(
      x[clean, , drop = FALSE], y[clean], mismatch_rate,
      if (is.null(cv_folds)) NULL else n_folds, shrinkage, margin_coverage,
      refit_clean, holdout_reference
    )
    if (!is.null(draw)) draw$repeat_id <- r
    draw
  })
  draws <- draws[!vapply(draws, is.null, logical(1))]
  if (length(draws)) simulation <- do.call(rbind, draws)
}
if (!is.null(simulation)) print(table(simulation$repeat_id,
                                      simulation$actually_misannotated))

# ---- 10. Annotation (flagging) margin threshold ----------------------------
swap_margin_positives <- if (!is.null(simulation)) {
  simulation$annotation_margin[simulation$actually_misannotated]
} else {
  numeric()
}
margin_threshold <- .select_threshold(
  annotation_margin[reference],
  if (length(swap_margin_positives)) swap_margin_positives else false_margin,
  target_specificity, "upper"
)
cat(sprintf("Annotation margin threshold %.4g (positives: %s)\n",
            margin_threshold,
            if (length(swap_margin_positives)) "simulated swaps"
            else "all wrong identities"))

# ---- 11. Risk-coverage ------------------------------------------------------
risk_coverage <- cbind(
  source = "observed", .risk_coverage(top_margin, correct_identity),
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

# ---- 12. Abstention and unknown-individual thresholds ----------------------
known <- if (is.null(simulation)) {
  list(margin = top_margin[reference], correct = correct_identity[reference])
} else {
  list(margin = simulation$top_margin, correct = simulation$agreement_real)
}
abstain_threshold <- .abstention_threshold(
  known$margin, known$correct, min_accuracy
)$threshold
abstention_candidates <- sort(unique(known$margin))
if (length(abstention_candidates) > 200L) {
  abstention_candidates <- abstention_candidates[
    unique(round(seq(1, length(abstention_candidates), length.out = 200L)))
  ]
}
abstention_curve <- data.frame(
  threshold = abstention_candidates,
  n_retained = vapply(abstention_candidates, function(cutoff) {
    sum(known$margin >= cutoff)
  }, integer(1)),
  agreement = vapply(abstention_candidates, function(cutoff) {
    mean(known$correct[known$margin >= cutoff])
  }, numeric(1))
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
low_margin <- top_margin < abstain_threshold | unknown_flag
confident_prediction <- !low_margin
abstain_reason_counts <- c(
  ambiguous_margin = sum(low_margin & !unknown_flag),
  potential_unknown_individual = sum(unknown_flag)
)
unknown_metrics <- if (unknown$n) {
  .binary_metrics(
    c(rep(FALSE, sum(reference)), rep(TRUE, unknown$n)),
    c(unknown_flag[reference], is_unknown(stranger_statistic)),
    (if (use_distance) 1 else -1) *
      c(unknown_statistic[reference], stranger_statistic)
  )
}
cat(sprintf("abstain < %.4g, unknown < %.4g\n", abstain_threshold,
            unknown_threshold))
print(abstain_reason_counts)

# ---- 13. Flagging and detection metrics ------------------------------------
flagged <- annotation_margin > margin_threshold | confidently_misassigned
all_alternatives_metrics <- .binary_metrics(
  c(rep(FALSE, sum(reference)), rep(TRUE, length(false_margin))),
  c(flagged[reference], false_margin > margin_threshold),
  c(annotation_margin[reference], false_margin)
)
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
if (!is.null(crossfit)) {
  misannotation_metrics$sensitivity <- crossfit$sensitivity
  misannotation_metrics$false_negative_rate <- 1 - crossfit$sensitivity
}
misannotation_metrics$sensitivity_sd <-
  if (is.null(crossfit)) NA_real_ else crossfit$sensitivity_sd
misannotation_metrics$n_repeats_validated <-
  if (is.null(crossfit)) 0L else crossfit$n_repeats
misannotation_metrics$source <- detection_source
str(misannotation_metrics)

# ---- 14. Identity recovery and correction validation -----------------------
identity_recovery <- mean(correct_identity)
recovery_misannotated <- if (swap_simulation) {
  mean(simulation$agreement_real[simulation$actually_misannotated])
} else {
  NA_real_
}
correction_validation <- if (!is.null(simulation)) {
  proposed <- simulation$annotation_margin > margin_threshold &
    simulation$top_margin >= max(
      abstain_threshold,
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
str(correction_validation)

# ---- 15. Pairwise swap detection -------------------------------------------
swap_samples <- integer()
swap_labels <- character()
for (i in seq_len(length(class_counts) - 1L)) {
  for (j in seq.int(i + 1L, length(class_counts))) {
    first_samples <- which(y == names(class_counts)[i])
    second_samples <- which(y == names(class_counts)[j])
    n_pairs <- min(length(first_samples), length(second_samples))
    swap_samples <- c(swap_samples, first_samples[seq_len(n_pairs)],
                      second_samples[seq_len(n_pairs)])
    swap_labels <- c(swap_labels, rep(names(class_counts)[j], n_pairs),
                     rep(names(class_counts)[i], n_pairs))
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

# ---- 16. Calibration and performance summaries -----------------------------
calibration <- list(
  method = "score_margin_threshold",
  margin_threshold = margin_threshold,
  exclusion_margin_threshold = top_margin_threshold,
  abstain_margin_threshold = abstain_threshold,
  prediction_coverage = mean(confident_prediction),
  abstain_reason_counts = abstain_reason_counts,
  min_accuracy = min_accuracy,
  first_pass_agreement = first_pass_agreement,
  low_agreement = low_agreement,
  oof_agreement = oof_agreement,
  dropped_individuals = dropped_individuals,
  n_confidently_misassigned = sum(confidently_misassigned),
  mismatch_simulation = if (is.null(simulation)) NULL else list(
    rate = mismatch_rate,
    repeats = length(unique(simulation$repeat_id)),
    n_samples = nrow(simulation),
    n_misannotated = sum(simulation$actually_misannotated),
    n_confidently_misassigned = sum(simulation$confidently_misassigned),
    agreement_simulated = mean(simulation$agreement_simulated),
    agreement_real = mean(simulation$agreement_real)
  )
)
performance <- list(
  correct_identity_recovery_accuracy = mean(correct_identity),
  top_three_identity_recovery_accuracy = mean(top_three),
  misannotation_detection = misannotation_metrics,
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

# ---- 17. Per-sample predictions (before feasibility) -----------------------
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
  confidently_misassigned = confidently_misassigned,
  annotation_correctness_probability = NA_real_,
  suspicious = flagged,
  potential_unknown_individual = unknown_flag,
  most_likely_alternative = ifelse(flagged & predicted_id != y,
                                   predicted_id, NA_character_),
  candidate_correction = NA_character_,
  excluded_from_refit = FALSE,
  stringsAsFactors = FALSE
)

# ---- 18. Feasibility assessment --------------------------------------------
assess_feasibility <- TRUE
feasibility <- NULL
if (assess_feasibility) {
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
  recovery_estimate <- if (is.finite(recovery_misannotated)) {
    recovery_misannotated
  } else {
    identity_recovery
  }
  recovery_status <- if (!enough_known) {
    "insufficient_evidence"
  } else if (recovery_estimate >= 0.8) {
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
    unknown_individual_detection = unknown_status
  )
  reliable_correction <- recovery_status == "appears_feasible"
} else {
  reliable_correction <- FALSE
}
str(feasibility)

# ---- 19. Candidate corrections and annotation status -----------------------
predictions$candidate_correction <- ifelse(
  flagged & reliable_correction & confident_prediction & predicted_id != y,
  predicted_id, NA_character_
)
predictions$annotation_status <- ifelse(
  !flagged & predictions$potential_unknown_individual,
  "potential_unknown_individual",
  ifelse(!is.na(predictions$candidate_correction),
         "suspicious_with_candidate_correction",
  ifelse(flagged, "suspicious_without_reliable_alternative",
  ifelse(!confident_prediction & predicted_id != y, "no_confident_prediction",
  ifelse(confident_prediction & predicted_id != y, "weak_disagreement",
         "no_evidence_of_incorrect_annotation"))))
)

# ---- 20. Low-agreement override --------------------------------------------
if (low_agreement) {
  predictions$predicted_id_clean <- NA_character_
  predictions$confident_prediction <- FALSE
  predictions$suspicious <- NA
  predictions$potential_unknown_individual <- if (detect_unknown) NA else FALSE
  predictions$most_likely_alternative <- NA_character_
  predictions$candidate_correction <- NA_character_
  predictions$annotation_status <- "insufficient_evidence"
  if (!is.null(feasibility)) {
    feasibility$misannotation_detection <- "insufficient_evidence"
    feasibility$identity_recovery <- "insufficient_evidence"
    if (detect_unknown) feasibility$unknown_individual_detection <-
      "insufficient_evidence"
  }
}

# ---- 21. Re-add samples of individuals dropped at the refit ----------------
if (length(dropped_individuals)) {
  gone <- first_pass_all[first_pass_all$observed_id %in% dropped_individuals, ,
                         drop = FALSE]
  held_gone <- gone$confidently_misassigned
  predictions <- rbind(predictions, data.frame(
    sample_id = gone$sample_id,
    observed_id = gone$observed_id,
    predicted_id = gone$predicted_id,
    predicted_id_clean = NA_character_,
    confident_prediction = FALSE,
    observed_id_score = gone$observed_id_score,
    top_candidate_score = gone$top_candidate_score,
    top_vs_observed_score_margin = gone$annotation_margin,
    top_score_margin = gone$top_margin,
    annotation_score_margin = gone$annotation_margin,
    confidently_misassigned = held_gone,
    annotation_correctness_probability = NA_real_,
    suspicious = if (low_agreement) NA else held_gone,
    potential_unknown_individual = if (detect_unknown) NA else FALSE,
    most_likely_alternative = NA_character_,
    candidate_correction = NA_character_,
    annotation_status = ifelse(held_gone,
                               "suspicious_without_reliable_alternative",
                               "insufficient_evidence"),
    excluded_from_refit = TRUE,
    stringsAsFactors = FALSE
  ))
  predictions <- predictions[match(colnames(expr), predictions$sample_id), ,
                             drop = FALSE]
  rownames(predictions) <- NULL
}

# ---- 22. Final result -------------------------------------------------------
result <- list(
  predictions = predictions,
  identity_scores = oof_scores,
  performance = performance,
  calibration = calibration,
  risk_coverage = risk_coverage,
  mismatch_simulation = simulation,
  held_out = list(
    n_confidently_misassigned = length(held_out_ids),
    n_samples = nrow(first_pass_all),
    fraction = length(held_out_ids) / nrow(first_pass_all),
    sample_ids = held_out_ids
  ),
  feasibility = feasibility,
  thresholds = list(
    annotation_margin = margin_threshold,
    unknown_margin = if (use_distance || !detect_unknown) -Inf else
      unknown_threshold,
    unknown_distance = if (use_distance) unknown_threshold else Inf,
    abstain_margin = abstain_threshold,
    target_specificity = target_specificity
  )
)
print(table(result$predictions$annotation_status, useNA = "ifany"))
