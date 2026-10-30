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
stopifnot(all(c("observed", "simulated") %in% result$risk_coverage$source))
stopifnot(all(is.finite(result$risk_coverage$coverage)))
pdf(NULL)
stopifnot(is.data.frame(plot_risk_coverage(result)))
invisible(dev.off())
stopifnot(all(is.na(result$predictions$candidate_correction)))

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
  min_agreement = 1.01
))
stopifnot(low$calibration$low_agreement)
stopifnot(all(low$predictions$annotation_status == "insufficient_evidence"))
stopifnot(all(is.na(low$predictions$predicted_id_clean)))
stopifnot(!anyNA(low$predictions$predicted_id))
stopifnot(inherits(tryCatch(
  check_annotations(expression, metadata, "sample", "individual",
                    progress = FALSE, min_agreement = 1.01),
  warning = function(w) w), "warning"))
