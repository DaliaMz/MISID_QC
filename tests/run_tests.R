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

result <- check_annotations(expression, metadata, "sample", "individual")
stopifnot(nrow(result$predictions) == length(ids))
stopifnot(identical(result$predictions$sample_id, samples))
stopifnot(all(is.na(
  result$predictions$annotation_correctness_probability
)))
stopifnot(all(c("misannotation_detection", "identity_recovery",
                "unknown_individual_detection") %in%
              names(result$feasibility)))
stopifnot(length(result$calibration$temperature) == 1L)

bad_metadata <- metadata
bad_metadata$sample[1L] <- "not_a_sample"
stopifnot(inherits(
  try(check_annotations(expression, bad_metadata, "sample", "individual"),
       silent = TRUE),
  "try-error"
))
