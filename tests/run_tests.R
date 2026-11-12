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
crossfit <- SampleAnnotationLikelihood:::.crossfit_swap_detection(
  negatives = c(0, 0.1, 0.2, 0.3), swap_data, target_specificity = 0.95
)
stopifnot(crossfit$n_repeats == 3L, crossfit$sensitivity == 1)
stopifnot(is.null(SampleAnnotationLikelihood:::.crossfit_swap_detection(
  c(0, 0.1), swap_data[swap_data$repeat_id == 1L, ], 0.95
)))

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
stopifnot(nrow(diagnostic_result$diagnostics$abstention$curve) > 0L)

report_file <- tempfile(fileext = ".pdf")
stopifnot(identical(annotation_report(diagnostic_result, report_file), report_file))
stopifnot(file.exists(report_file), file.info(report_file)$size > 0)
unlink(report_file)
report_file <- tempfile(fileext = ".pdf")
annotation_report(no_feasibility, report_file)
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
many_calls <- diagnostic_result
many_calls$predictions <- many_calls$predictions[rep(seq_len(nrow(result$predictions)),
                                                   length.out = 23L), ]
many_calls$predictions$suspicious <- TRUE
report_file <- tempfile(fileext = ".pdf")
annotation_report(many_calls, report_file)
pdf_bytes <- readBin(report_file, "raw", n = file.info(report_file)$size)
# Eight pages plus the PDF page-tree object: overview, five plots, and two tables.
stopifnot(length(grepRaw(charToRaw("/Type /Page"), pdf_bytes, all = TRUE)) == 9L)
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
single <- check_annotations(expression, metadata, "sample", "individual",
                            progress = FALSE, mismatch_repeats = 1L)
stopifnot(single$calibration$mismatch_simulation$repeats == 1L)
stopifnot(inherits(
  try(check_annotations(expression, metadata, "sample", "individual",
                        progress = FALSE, mismatch_repeats = 0),
      silent = TRUE),
  "try-error"
))

low <- suppressWarnings(check_annotations(
  expression, metadata, "sample", "individual", progress = FALSE,
  min_agreement = 1.01, stop_on_low_agreement = FALSE
))
stopifnot(low$calibration$low_agreement)
stopifnot(all(low$predictions$annotation_status == "insufficient_evidence"))
stopifnot(all(is.na(low$predictions$predicted_id_clean)))
stopifnot(!anyNA(low$predictions$predicted_id))
scrambled <- metadata
scrambled$individual <- sample(scrambled$individual)
stopifnot(inherits(try(
  check_annotations(expression, scrambled, "sample", "individual",
                    progress = FALSE),
  silent = TRUE), "try-error"))
stopifnot(inherits(tryCatch(
  check_annotations(expression, metadata, "sample", "individual",
                    progress = FALSE, min_agreement = 1.01,
                    stop_on_low_agreement = FALSE),
  warning = function(w) w), "warning"))

# Individuals dropped at the refit and held-out samples still get a row and a call.
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
dropped <- check_annotations(
  X, data.frame(sample = colnames(X), individual = obs_ids), "sample",
  "individual", cv_folds = NULL, min_agreement = 0, progress = FALSE,
  mismatch_repeats = 2
)
dp <- dropped$predictions
stopifnot(identical(dp$sample_id, colnames(X)))
stopifnot(all(c("nearest_class_distance", "unknown_detection_score",
                "unknown_detection_threshold") %in% names(dp)))
stopifnot(all(is.na(dp$unknown_detection_score[dp$excluded_from_refit])))
stopifnot("d1" %in% dropped$calibration$dropped_individuals)
stopifnot(all(dp$excluded_from_refit == (dp$observed_id %in% dropped$calibration$dropped_individuals)))
stopifnot(all(dp$suspicious[dp$sample_id %in% dropped$held_out$sample_ids]))
stopifnot(all(dp$annotation_status[dp$excluded_from_refit & dp$suspicious] ==
                "suspicious_without_reliable_alternative"))
