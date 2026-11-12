score_unknown <- function(predictions, truth, enabled) {
  if (!is.logical(enabled) || length(enabled) != 1L || is.na(enabled)) {
    stop("detect_unknown must be TRUE or FALSE to score unknown detection.")
  }
  if (!all(c("sample", "is_unknown") %in% names(truth)) ||
      !is.logical(truth$is_unknown) || anyNA(truth$is_unknown) ||
      anyNA(truth$sample) || anyDuplicated(truth$sample)) {
    stop("Truth must contain unique sample IDs and non-missing logical is_unknown.")
  }
  if (!all(c("sample_id", "potential_unknown_individual") %in% names(predictions)) ||
      anyNA(predictions$sample_id) || anyDuplicated(predictions$sample_id) ||
      !is.logical(predictions$potential_unknown_individual)) {
    stop("Predictions must contain unique sample IDs and logical unknown flags.")
  }
  idx <- match(predictions$sample_id, truth$sample)
  if (anyNA(idx)) stop("Prediction sample IDs do not match truth.")
  actual <- truth$is_unknown[idx]
  flag <- predictions$potential_unknown_individual
  assessed <- !is.na(flag) & enabled
  n_unknown <- sum(truth$is_unknown)
  n_known <- sum(!truth$is_unknown)
  assessed_unknown <- sum(assessed & actual)
  assessed_known <- sum(assessed & !actual)
  list(
    unknown_n_true = n_unknown,
    unknown_n_known = n_known,
    unknown_n_assessed_true = assessed_unknown,
    unknown_n_assessed_known = assessed_known,
    unknown_sens_true = if (assessed_unknown) {
      sum(flag[assessed & actual]) / assessed_unknown
    } else NA_real_,
    unknown_spec_true = if (assessed_known) {
      sum(!flag[assessed & !actual]) / assessed_known
    } else NA_real_
  )
}

rescore_unknown_rows <- function(files, dirty_dir) {
  lapply(files, function(file) {
    row <- data.table::as.data.table(readRDS(file))
    if (nrow(row) != 1L || !"status" %in% names(row)) {
      stop("Expected one result row with status in ", file)
    }
    if (is.na(row$status)) stop("Missing result status in ", file)
    if (row$status == "error") return(row)
    dirty_file <- file.path(dirty_dir, basename(file))
    if (!file.exists(dirty_file)) stop("Missing saved result: ", dirty_file)
    saved <- readRDS(dirty_file)
    if (is.null(saved$result$predictions) || is.null(saved$truth) ||
        is.null(saved$check_args$detect_unknown)) {
      stop("Saved result lacks predictions, truth, or detect_unknown: ", dirty_file)
    }
    metrics <- score_unknown(saved$result$predictions, saved$truth,
                             saved$check_args$detect_unknown)
    for (name in names(metrics)) data.table::set(row, j = name, value = metrics[[name]])
    row
  })
}
