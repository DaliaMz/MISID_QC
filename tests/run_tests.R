library(SampleAnnotationLikelihood)

set.seed(12)
ids <- rep(c("donor_a", "donor_b", "donor_c"), each = 4)
centers <- rbind(
  donor_a = c(-3, 0, 0, 0, 0, 0),
  donor_b = c(0, 3, 0, 0, 0, 0),
  donor_c = c(0, 0, 3, 0, 0, 0)
)
samples <- paste0("sample_", seq_along(ids))
expression <- vapply(seq_along(ids), function(i) {
  centers[ids[i], ] + rnorm(ncol(centers), sd = 0.25)
}, numeric(ncol(centers)))
rownames(expression) <- paste0("gene_", seq_len(nrow(expression)))
colnames(expression) <- samples
metadata <- data.frame(
  sample = rev(samples),
  individual = rev(ids),
  stringsAsFactors = FALSE
)

model <- fit_shrinkage_lda(t(expression), ids, shrinkage = 0.2)
scores <- predict_shrinkage_lda(model, t(expression))
stopifnot(is.matrix(scores), nrow(scores) == length(ids), ncol(scores) == 3L)
stopifnot(all(is.finite(scores)))
stopifnot(isTRUE(all.equal(
  unname(scores),
  unname(t(expression) %*% t(model$beta) +
           matrix(model$intercepts, nrow = length(ids), ncol = 3L,
                  byrow = TRUE))
)))

holdout_scores <- matrix(c(10, 9, 0), nrow = 1L,
                         dimnames = list("sample", c("a", "b", "c")))
holdout_threshold <- get(".holdout_threshold",
                         asNamespace("SampleAnnotationLikelihood"))
stopifnot(holdout_threshold(holdout_scores, TRUE, 0.75, "top") == 1)
stopifnot(holdout_threshold(holdout_scores, TRUE, 0.75,
                            "all_alternatives") == 10)

correction_threshold <- get(".correction_threshold",
                            asNamespace("SampleAnnotationLikelihood"))
simulated <- data.frame(
  top_margin = c(5, 4, 3, 100, 100),
  actually_misannotated = c(TRUE, TRUE, FALSE, FALSE, FALSE),
  agreement_simulated = c(FALSE, TRUE, FALSE, TRUE, TRUE),
  agreement_real = c(TRUE, FALSE, FALSE, TRUE, TRUE)
)
cutoff <- correction_threshold(simulated, 0.9)
stopifnot(cutoff$threshold == 5, cutoff$n_eligible == 3L,
          cutoff$n_retained == 1L, cutoff$coverage == 1/3, cutoff$accuracy == 1)
stopifnot(correction_threshold(simulated, 0.5)$threshold == 4)
simulated$top_margin[2] <- 5
stopifnot(is.infinite(correction_threshold(simulated, 0.9)$threshold),
          is.infinite(correction_threshold(NULL, 0.9)$threshold))

progress_messages <- character()
withCallingHandlers(
  result <- check_annotations(expression, metadata, "sample", "individual",
                              progress = TRUE),
  message = function(condition) {
    progress_messages <<- c(progress_messages, conditionMessage(condition))
    invokeRestart("muffleMessage")
  }
)
stopifnot(any(grepl("Leave-one-sample-out models", progress_messages)))
stopifnot(any(grepl("held-out individuals", progress_messages,
                    ignore.case = TRUE)))
stopifnot(any(grepl("Estimated detection errors unavailable|Estimated misannotation rate",
                    progress_messages)))
stopifnot(nrow(result$predictions) == length(ids))
stopifnot(all(c("potential_unknown_individual", "nearest_class_distance",
                "unknown_detection_score", "unknown_detection_threshold") %in%
              names(result$predictions)))
stopifnot(!any(result$predictions$annotation_status ==
                 "potential_unknown_individual"))
stopifnot(isTRUE(all.equal(result$predictions$unknown_detection_score,
                           result$predictions$nearest_class_distance)))
stopifnot(all(result$predictions$unknown_detection_threshold ==
              result$thresholds$unknown_distance))
stopifnot(identical(dim(result$identity_scores), c(length(ids), 3L)))
stopifnot(identical(result$predictions$sample_id, samples))
stopifnot(all(is.na(
  result$predictions$annotation_correctness_probability
)))
stopifnot(all(c("misannotation_detection", "identity_recovery",
                "unknown_individual_detection") %in%
              names(result$feasibility)))
stopifnot("pairwise_swap_detection" %in% names(result$performance))
stopifnot(length(result$calibration$margin_threshold) == 1L)
stopifnot(all(is.na(result$predictions$predicted_id_clean[
  !result$predictions$confident_prediction])))
stopifnot(!anyNA(result$predictions$predicted_id))
stopifnot(!any(!is.na(result$predictions$candidate_correction) &
                 !result$predictions$confident_prediction))
stopifnot(length(result$calibration$abstain_reason_counts) == 2L)
stopifnot(identical(result$calibration$holdout_reference,
                    "all_alternatives"))
top_reference <- check_annotations(
  expression, metadata, "sample", "individual", progress = FALSE,
  holdout_reference = "top"
)
stopifnot(identical(top_reference$calibration$holdout_reference, "top"))
stopifnot(inherits(
  try(check_annotations(expression, metadata, "sample", "individual",
                        progress = FALSE, holdout_reference = "invalid"),
      silent = TRUE),
  "try-error"
))
stopifnot(all(c("observed", "simulated") %in% result$risk_coverage$source))
stopifnot(all(is.finite(result$risk_coverage$coverage)))
pdf(NULL)
stopifnot(is.data.frame(plot_risk_coverage(result)))
invisible(dev.off())
stopifnot(all(is.na(result$predictions$candidate_correction)))

# Detection is validated on simulated swaps; sensitivity is cross-fitted when
# at least two repeats are available and in-sample otherwise.
stopifnot(all(c("all_alternatives_detection",
                "recovery_accuracy_misannotated",
                "correction_validation") %in% names(result$performance)))
stopifnot(result$performance$misannotation_detection$source %in%
            c("simulated_swaps_cross_fitted", "simulated_swaps_in_sample"))
stopifnot(result$performance$misannotation_detection$n_repeats_validated <= 3L)

estimate_errors <- SampleAnnotationLikelihood:::.estimate_detection_errors
estimate <- estimate_errors(c(rep(TRUE, 12), rep(FALSE, 88)), 0.8, 0.95,
                            "simulated_swaps_cross_fitted")
stopifnot(estimate$status == "available",
          isTRUE(all.equal(estimate$prevalence, 7 / 75)),
          isTRUE(all.equal(estimate$scenarios$prevalence,
                           7 / 75 + c(-0.025, 0, 0.025))),
          isTRUE(all.equal(estimate$scenarios$expected_false_positives[2], 68 / 15)),
          isTRUE(all.equal(estimate$scenarios$expected_false_negatives[2], 28 / 15)))
boundary <- estimate_errors(rep(FALSE, 100), 0.8, 1, "simulated_swaps_in_sample")
stopifnot(boundary$status == "available",
          identical(boundary$scenarios$prevalence[1], 0))
boundary <- estimate_errors(c(rep(TRUE, 5), rep(FALSE, 95)), 0.8, 0.95,
                            "simulated_swaps_in_sample")
stopifnot(boundary$status == "available", boundary$prevalence == 0)
boundary <- estimate_errors(rep(TRUE, 100), 1, 0.95, "simulated_swaps_in_sample")
stopifnot(boundary$status == "available",
          identical(boundary$scenarios$prevalence[3], 1))
stopifnot(estimate_errors(rep(FALSE, 100), 0.8, 0.95,
                          "simulated_swaps_cross_fitted")$status == "unavailable",
          estimate_errors(rep(FALSE, 100), 0.05, 0.95,
                          "simulated_swaps_in_sample")$status == "unavailable",
          estimate_errors(rep(FALSE, 100), 0.8, NA_real_,
                          "simulated_swaps_in_sample")$status == "unavailable",
          estimate_errors(rep(FALSE, 100), 0.8, 0.95,
                          "all_wrong_identities")$status == "unavailable")
estimate <- result$performance$estimated_detection_errors
stopifnot(estimate$n_flagged == sum(result$predictions$suspicious),
          estimate$source == result$performance$misannotation_detection$source)
if (estimate$status == "available") {
  stopifnot(isTRUE(all.equal(
    estimate$prevalence,
    (mean(result$predictions$suspicious) - (1 - estimate$specificity)) /
      (estimate$sensitivity + estimate$specificity - 1)
  )))
}

# Feasibility never overwrites sample statuses.
no_feasibility <- check_annotations(
  expression, metadata, "sample", "individual", progress = FALSE,
  assess_feasibility = FALSE
)
stopifnot(!any(no_feasibility$predictions$annotation_status ==
                 "insufficient_evidence"))
stopifnot(!any(result$predictions$annotation_status == "insufficient_evidence"))

swap_data <- data.frame(
  annotation_margin = c(10, 11, 12, 13, 14, 15),
  actually_misannotated = TRUE,
  repeat_id = rep(1:3, each = 2)
)
stopifnot(is.null(SampleAnnotationLikelihood:::.crossfit_swap_detection(
  negatives = c(0, 0.1, 0.2, 0.3), swap_data, target_specificity = 0.95
)))
stopifnot(is.null(SampleAnnotationLikelihood:::.crossfit_swap_detection(
  c(0, 0.1), swap_data[swap_data$repeat_id == 1L, ], 0.95
)))
swap_data <- data.frame(
  annotation_margin = c(10, 11, 0, 9, 12, 13, 0, 20, 14, 15, 0, 3),
  actually_misannotated = rep(c(TRUE, TRUE, FALSE, FALSE), 3),
  repeat_id = rep(1:3, each = 4)
)
crossfit <- SampleAnnotationLikelihood:::.crossfit_swap_detection(
  c(0, 0.1, 0.2, 0.3), swap_data, 0.95
)
stopifnot(crossfit$sensitivity == 1 / 3,
          crossfit$simulated_specificity == 5 / 6,
          crossfit$n_simulated_negatives == 6L)
without_negatives <- swap_data[swap_data$actually_misannotated, ]
stopifnot(is.null(SampleAnnotationLikelihood:::.crossfit_swap_detection(
  c(0, 0.1, 0.2, 0.3), without_negatives, 0.95
)))
select_threshold <- SampleAnnotationLikelihood:::.select_threshold
ref_negatives <- c(0, 0.1, 0.2, 0.3)
sim_negatives <- c(0, 9, 0, 20, 0, 3)
positives <- c(10, 11, 12, 13, 14, 15)
reference_only <- select_threshold(ref_negatives, positives, 0.95, "upper")
dual <- select_threshold(ref_negatives, positives, 0.95, "upper",
                         additional_negatives = sim_negatives)
stopifnot(mean(sim_negatives <= reference_only) < 0.95,
          mean(ref_negatives <= dual) >= 0.95,
          mean(sim_negatives <= dual) >= 0.95,
          dual == 21,
          select_threshold(ref_negatives, numeric(), 0.95, "upper",
                           additional_negatives = sim_negatives) == 20,
          select_threshold(c(3, 4), numeric(), 0.95, "lower",
                           additional_negatives = c(0, 1)) == 0)

set.seed(12)
diagnostic_result <- check_annotations(
  expression, metadata, "sample", "individual", progress = FALSE,
  diagnostics = TRUE
)
stopifnot(all(c("first_pass", "first_pass_scores", "final_pass",
                "annotation_threshold", "abstention", "unknown", "gates") %in%
              names(diagnostic_result$diagnostics)))
stopifnot(identical(
  diagnostic_result$diagnostics$first_pass$sample_id, samples
))
stopifnot(setequal(
  diagnostic_result$held_out$sample_ids,
  diagnostic_result$diagnostics$first_pass$sample_id[
    diagnostic_result$diagnostics$first_pass$held_out
  ]
))
stopifnot(nrow(diagnostic_result$diagnostics$annotation_threshold) > 0L)
stopifnot("simulated_unchanged_negative" %in%
            diagnostic_result$diagnostics$annotation_threshold$population,
          diagnostic_result$performance$misannotation_detection$specificity >= 0.95,
          diagnostic_result$performance$misannotation_detection$simulated_specificity_in_sample >= 0.95)
stopifnot(nrow(diagnostic_result$diagnostics$abstention$curve) > 0L)
holdout_reference <- SampleAnnotationLikelihood:::.report_holdout_reference(diagnostic_result)
stopifnot(isTRUE(all.equal(
  unname(quantile(holdout_reference, 0.025, type = 1)),
  diagnostic_result$calibration$exclusion_margin_threshold
)))
top_result <- check_annotations(
  expression, metadata, "sample", "individual", progress = FALSE,
  holdout_reference = "top", diagnostics = TRUE
)
top_reference <- SampleAnnotationLikelihood:::.report_holdout_reference(top_result)
stopifnot(length(top_reference) == sum(top_result$diagnostics$first_pass$correct),
          isTRUE(all.equal(unname(quantile(top_reference, 0.025, type = 1)),
                           top_result$calibration$exclusion_margin_threshold)))

recovery_curves <- SampleAnnotationLikelihood:::.report_recovery_curves(
  data.frame(top_margin = c(1, 2, 3, 4),
             actually_misannotated = c(TRUE, FALSE, FALSE, FALSE),
             agreement_simulated = c(FALSE, FALSE, TRUE, TRUE),
             agreement_real = c(TRUE, FALSE, TRUE, TRUE))
)
full_coverage <- recovery_curves[recovery_curves$x == 1, ]
stopifnot(full_coverage$y[full_coverage$group == "All simulated samples"] == 0.75,
          full_coverage$y[full_coverage$group ==
                            "Swaps and unchanged disagreements"] == 0.5,
          nrow(SampleAnnotationLikelihood:::.report_recovery_curves(NULL)) == 0L)
unchanged <- data.frame(top_margin = c(1, 2),
                        actually_misannotated = FALSE,
                        agreement_simulated = TRUE, agreement_real = TRUE)
stopifnot(all(SampleAnnotationLikelihood:::.report_recovery_curves(unchanged)$group ==
                "All simulated samples"))

report_file <- tempfile(fileext = ".pdf")
stopifnot(identical(annotation_report(diagnostic_result, report_file), report_file))
stopifnot(file.exists(report_file), file.info(report_file)$size > 0)
unlink(report_file)
report_file <- tempfile(fileext = ".pdf")
annotation_report(no_feasibility, report_file)
stopifnot(file.exists(report_file), file.info(report_file)$size > 0)
unlink(report_file)
older_result <- no_feasibility
older_result$performance$estimated_detection_errors <- NULL
report_file <- tempfile(fileext = ".pdf")
annotation_report(older_result, report_file)
stopifnot(file.exists(report_file), file.info(report_file)$size > 0)
unlink(report_file)
unknown_off <- check_annotations(
  expression, metadata, "sample", "individual", progress = FALSE,
  detect_unknown = FALSE, diagnostics = TRUE
)
report_file <- tempfile(fileext = ".pdf")
annotation_report(unknown_off, report_file)
stopifnot(file.exists(report_file), file.info(report_file)$size > 0)
unlink(report_file)
curve <- SampleAnnotationLikelihood:::.report_rate_curve(
  c(1, 2, 2, 4), c(0, 1, 2, 3), 2
)
at_cutoff <- curve[curve$x == 2, ]
stopifnot(at_cutoff$rate[at_cutoff$metric == "Sensitivity"] == 0.25,
          at_cutoff$rate[at_cutoff$metric == "Specificity"] == 0.75)
curve <- SampleAnnotationLikelihood:::.report_rate_curve(
  c(1, 2, 2, 4), c(0, 1, 2, 3), 2, "lower"
)
at_cutoff <- curve[curve$x == 2, ]
stopifnot(at_cutoff$rate[at_cutoff$metric == "Sensitivity"] == 0.25,
          at_cutoff$rate[at_cutoff$metric == "Specificity"] == 0.5)
report_file <- tempfile(fileext = ".pdf")
set.seed(123)
seed_before <- .Random.seed
annotation_report(diagnostic_result, report_file)
stopifnot(identical(seed_before, .Random.seed))
pdf_bytes <- readBin(report_file, "raw", n = file.info(report_file)$size)
stopifnot(length(grepRaw(charToRaw("/Type /Page"), pdf_bytes, all = TRUE)) == 7L)
unlink(report_file)

cv_progress_messages <- character()
withCallingHandlers(
  cv_result <- check_annotations(expression, metadata, "sample", "individual",
                                 cv_folds = 3L, progress = TRUE),
  message = function(condition) {
    cv_progress_messages <<- c(cv_progress_messages,
                               conditionMessage(condition))
    invokeRestart("muffleMessage")
  }
)
stopifnot(any(grepl("stratified 3-fold CV", cv_progress_messages)))
stopifnot(any(grepl("Cross-validation folds: 3/3", cv_progress_messages)))
stopifnot(identical(dim(cv_result$identity_scores), c(length(ids), 3L)))
stopifnot(all(is.finite(cv_result$identity_scores)))
stopifnot(inherits(
  try(check_annotations(expression, metadata, "sample", "individual",
                         cv_folds = 5L), silent = TRUE),
  "try-error"
))

bad_metadata <- metadata
bad_metadata$sample[1L] <- "not_a_sample"
stopifnot(inherits(
  try(check_annotations(expression, bad_metadata, "sample", "individual"),
       silent = TRUE),
  "try-error"
))

stopifnot(result$calibration$mismatch_simulation$repeats == 3L)
stopifnot(result$calibration$mismatch_simulation$rate == 0.1)
for (disabled_rate in list(0, NULL)) {
  disabled <- check_annotations(
    expression, metadata, "sample", "individual", progress = FALSE,
    mismatch_rate = disabled_rate, detect_unknown = FALSE
  )
  stopifnot(is.null(disabled$mismatch_simulation),
            is.null(disabled$calibration$mismatch_simulation),
            disabled$performance$estimated_detection_errors$status == "unavailable")
}
set.seed(12)
held_zero <- check_annotations(
  expression, metadata, "sample", "individual", progress = FALSE,
  mismatch_rate = "held_out", mismatch_repeats = 1, detect_unknown = FALSE
)
stopifnot(held_zero$held_out$n_confidently_misassigned == 0L,
          held_zero$calibration$mismatch_simulation$rate == 0,
          held_zero$calibration$mismatch_simulation$rate_mode == "held_out",
          sum(held_zero$mismatch_simulation$actually_misannotated) == 2L,
          held_zero$calibration$mismatch_simulation$actual_rate ==
            mean(held_zero$mismatch_simulation$actually_misannotated))
for (invalid in list("invalid", NA_character_, -0.1, 1.1, NA_real_,
                     c(0.1, 0.2), numeric())) {
  stopifnot(inherits(try(
    check_annotations(expression, metadata, "sample", "individual",
                      progress = FALSE, mismatch_rate = invalid),
    silent = TRUE), "try-error"))
}
single <- check_annotations(expression, metadata, "sample", "individual",
                            progress = FALSE, mismatch_repeats = 1L)
stopifnot(single$calibration$mismatch_simulation$repeats == 1L)
stopifnot(inherits(
  try(check_annotations(expression, metadata, "sample", "individual",
                        progress = FALSE, mismatch_repeats = 0),
      silent = TRUE),
  "try-error"
))

scrambled <- metadata
scrambled$individual <- rep(c("donor_a", "donor_b", "donor_c"), 4)
low <- suppressWarnings(check_annotations(
  expression, scrambled, "sample", "individual", progress = FALSE,
  min_agreement = 0, mismatch_rate = 0, detect_unknown = FALSE
))
agreement <- low$calibration$first_pass_agreement
stopifnot(agreement < 0.6, !low$calibration$low_agreement)
stopped <- tryCatch(
  check_annotations(expression, scrambled, "sample", "individual",
                    progress = FALSE, mismatch_rate = 0, detect_unknown = FALSE),
  error = function(e) e
)
stopifnot(inherits(stopped, "error"),
          grepl("below min_agreement \\(60.0%\\)", conditionMessage(stopped)))
at_boundary <- suppressWarnings(check_annotations(
  expression, scrambled, "sample", "individual", progress = FALSE,
  min_agreement = agreement, mismatch_rate = 0, detect_unknown = FALSE
))
stopifnot(at_boundary$calibration$first_pass_agreement == agreement)
stopped <- tryCatch(
  check_annotations(expression, scrambled, "sample", "individual",
                    progress = FALSE, min_agreement = 0.9,
                    mismatch_rate = 0, detect_unknown = FALSE),
  error = function(e) e
)
stopifnot(inherits(stopped, "error"),
          grepl("below min_agreement \\(90.0%\\)", conditionMessage(stopped)))
set.seed(2)
boundary_expr <- matrix(rnorm(40), nrow = 4)
colnames(boundary_expr) <- paste0("boundary_", seq_len(10))
boundary_meta <- data.frame(sample = colnames(boundary_expr),
                            individual = rep(c("a", "b"), each = 5))
boundary <- suppressWarnings(check_annotations(
  boundary_expr, boundary_meta, "sample", "individual", progress = FALSE,
  mismatch_rate = 0, detect_unknown = FALSE
))
stopifnot(boundary$calibration$first_pass_agreement == 0.6)
warning_at_boundary <- tryCatch(
  check_annotations(boundary_expr, boundary_meta, "sample", "individual",
                    progress = FALSE, mismatch_rate = 0, detect_unknown = FALSE),
  warning = function(w) w
)
stopifnot(inherits(warning_at_boundary, "warning"),
          grepl("below 75%", conditionMessage(warning_at_boundary)))
stopifnot(inherits(try(
  check_annotations(boundary_expr, boundary_meta, "sample", "individual",
                    progress = FALSE, min_agreement = 0.61,
                    mismatch_rate = 0, detect_unknown = FALSE),
  silent = TRUE), "try-error"))
for (invalid in list(-0.1, 1.01, NA_real_, Inf, numeric(), c(0, 1))) {
  stopifnot(inherits(try(
    check_annotations(expression, metadata, "sample", "individual",
                      progress = FALSE, min_agreement = invalid),
    silent = TRUE), "try-error"))
}

# Held-out samples remain scored and use the same cutoff as every other sample.
set.seed(3)
K <- 6; n <- 5
true_ids <- rep(paste0("d", 1:K), each = n)
cent <- diag(K) * 4
X <- vapply(seq_along(true_ids), function(i) {
  cent[match(true_ids[i], paste0("d", 1:K)), ] + rnorm(K, sd = 0.3)
}, numeric(K))
colnames(X) <- paste0("s", seq_along(true_ids))
rownames(X) <- paste0("g", 1:K)
obs_ids <- true_ids
obs_ids[2:5] <- "d2"; obs_ids[c(6, 11)] <- "d1"
set.seed(4)
held_result <- check_annotations(
  X, data.frame(sample = colnames(X), individual = obs_ids), "sample",
  "individual", cv_folds = NULL, min_agreement = 0, progress = FALSE,
  mismatch_repeats = 2, margin_coverage = 0, detect_unknown = FALSE
)
selected <- with(held_result$mismatch_simulation,
                 actually_misannotated | simulated_id != predicted_id)
retained <- selected &
  held_result$mismatch_simulation$top_margin >= held_result$thresholds$correction_margin
correction_calibration <- held_result$calibration$correction_margin
stopifnot(correction_calibration$n_eligible == sum(selected),
          correction_calibration$n_retained == sum(retained))
if (any(retained)) {
  stopifnot(correction_calibration$accuracy ==
              mean(held_result$mismatch_simulation$agreement_real[retained]),
            correction_calibration$accuracy >= held_result$calibration$min_accuracy)
}
stopifnot(all(held_result$predictions$top_score_margin[
  !is.na(held_result$predictions$candidate_correction)] >=
    held_result$thresholds$correction_margin))
sim <- held_result$mismatch_simulation
proposed <- sim$annotation_margin > held_result$thresholds$annotation_margin &
  sim$top_margin >= max(held_result$thresholds$abstain_margin,
                        held_result$thresholds$correction_margin) &
  sim$simulated_id != sim$predicted_id
stopifnot(held_result$performance$correction_validation$n_proposed == sum(proposed))
if (any(proposed)) {
  stopifnot(held_result$performance$correction_validation$precision ==
              mean(sim$agreement_real[proposed]))
}
stopifnot(held_result$feasibility$identity_recovery ==
            if (is.finite(correction_calibration$threshold))
              "appears_feasible" else "unreliable")
without_swaps <- check_annotations(
  X, data.frame(sample = colnames(X), individual = obs_ids), "sample",
  "individual", min_agreement = 0, progress = FALSE,
  mismatch_rate = 0, detect_unknown = FALSE
)
stopifnot(without_swaps$feasibility$identity_recovery == "insufficient_evidence",
          is.infinite(without_swaps$thresholds$correction_margin),
          all(is.na(without_swaps$predictions$candidate_correction)))
dp <- held_result$predictions
set.seed(4)
held_rate_result <- check_annotations(
  X, data.frame(sample = colnames(X), individual = obs_ids), "sample",
  "individual", cv_folds = NULL, min_agreement = 0, progress = FALSE,
  mismatch_repeats = 2, margin_coverage = 0, detect_unknown = FALSE,
  mismatch_rate = "held_out"
)
held_rate <- held_rate_result$held_out$fraction
stopifnot(held_rate > 0,
          held_rate_result$calibration$mismatch_simulation$rate == held_rate,
          held_rate_result$calibration$mismatch_simulation$actual_rate ==
            mean(held_rate_result$mismatch_simulation$actually_misannotated))
set.seed(4)
fixed_rate_result <- check_annotations(
  X, data.frame(sample = colnames(X), individual = obs_ids), "sample",
  "individual", cv_folds = NULL, min_agreement = 0, progress = FALSE,
  mismatch_repeats = 2, margin_coverage = 0, detect_unknown = FALSE,
  mismatch_rate = held_rate
)
stopifnot(identical(held_rate_result$mismatch_simulation,
                    fixed_rate_result$mismatch_simulation))
stopifnot(identical(dp$sample_id, colnames(X)))
stopifnot(!("refit_clean" %in% names(formals(check_annotations))))
stopifnot(all(c("nearest_class_distance", "unknown_detection_score",
                "unknown_detection_threshold") %in% names(dp)))
stopifnot(!any(dp$excluded_from_refit),
          !length(held_result$calibration$dropped_individuals),
          length(held_result$held_out$sample_ids) > 0L,
          any(dp$confidently_misassigned & !dp$suspicious),
          any(dp$confidently_misassigned & dp$suspicious),
          identical(dp$suspicious,
                    dp$annotation_score_margin > held_result$thresholds$annotation_margin))
clean_counts <- table(dp$observed_id[!dp$confidently_misassigned])
simulation_ids <- dp$sample_id[!dp$confidently_misassigned &
                                dp$observed_id %in% names(clean_counts)[clean_counts >= 2L]]
stopifnot(all(vapply(split(held_result$mismatch_simulation$sample_id,
                           held_result$mismatch_simulation$repeat_id),
                     function(ids) setequal(ids, simulation_ids), logical(1))))
metrics <- held_result$performance$misannotation_detection
stopifnot(identical(held_result$feasibility$misannotation_detection,
                    if (metrics$sensitivity >= 0.8 &&
                        metrics$simulated_specificity >= 0.95) {
                      "appears_feasible"
                    } else "unreliable"))
set.seed(2)
single_held <- check_annotations(
  X, data.frame(sample = colnames(X), individual = obs_ids), "sample",
  "individual", min_agreement = 0, progress = FALSE, mismatch_repeats = 1,
  margin_coverage = 0, detect_unknown = FALSE
)
single_metrics <- single_held$performance$misannotation_detection
stopifnot(identical(single_metrics$source, "simulated_swaps_in_sample"),
          identical(single_metrics$simulated_specificity,
                    single_metrics$simulated_specificity_in_sample))
set.seed(2)
no_swaps <- check_annotations(
  X, data.frame(sample = colnames(X), individual = obs_ids), "sample",
  "individual", min_agreement = 0, progress = FALSE, mismatch_rate = 0,
  margin_coverage = 0, detect_unknown = FALSE
)
stopifnot(is.na(no_swaps$performance$misannotation_detection$simulated_specificity),
          identical(no_swaps$feasibility$misannotation_detection,
                    "insufficient_evidence"))
