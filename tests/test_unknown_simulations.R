# Run from the repository root with Rscript tests/test_unknown_simulations.R.
library(data.table)

synthetic <- new.env(parent = globalenv())
sys.source("simulations/Simulate_full_dataset.R", envir = synthetic)
synthetic_grid_counts <- synthetic$grid[, .N, by = .(
  type, var_label, nuisance, n_per_label, misid_rate, rep
)]
stopifnot(nrow(synthetic_grid_counts) == 360L,
          all(synthetic_grid_counts$N == 2L),
          setequal(synthetic$grid$n_unknown_individuals, c(0L, 1L)),
          setequal(synthetic$grid[ n_unknown_individuals == 1L, scenario], 1:360),
          setequal(synthetic$grid[n_unknown_individuals == 0L, scenario], 361:720),
          !anyDuplicated(synthetic$grid$scenario))
covariate_grid_env <- new.env(parent = globalenv())
covariate_grid_env$covariate_nuisance_levels <-
  c("none", "categorical", "numeric", "both")
covariate_grid_env$n_unknown_individuals <- 0L
sys.source("simulations/Simulate_full_dataset.R", envir = covariate_grid_env)
stopifnot(nrow(covariate_grid_env$grid) == 1440L,
          setequal(covariate_grid_env$grid$covariate_nuisance,
                   c("none", "categorical", "numeric", "both")),
          all(covariate_grid_env$grid$n_tissues == 4L),
          all(covariate_grid_env$grid$var_label +
              covariate_grid_env$grid$var_discrete +
              covariate_grid_env$grid$var_cont +
              covariate_grid_env$grid$var_tissue +
              covariate_grid_env$grid$var_numeric_covariate <= 1))
for (type in c("gaussian", "rnaseq", "olink")) {
  set.seed(42)
  sim <- synthetic$simulate_dataset(
    type, 0.2, 0.1, 0.1, n_per_label = 6, n_labels = 5,
    n_features = 20, misid_rate = 0, n_unknown_individuals = 2,
    unknown_samples_per_individual = 2
  )
  meta <- sim$meta
  known <- meta[is_unknown == FALSE]
  unknown <- meta[is_unknown == TRUE]
  stopifnot(nrow(known) == 30L, nrow(unknown) == 4L,
            uniqueN(unknown$label_true) == 2L,
            all(table(unknown$label_true) == 2L),
            all(unknown$is_misid), !any(known$is_misid),
            all(unknown$label_observed %in% known$label_true),
            !any(unknown$label_true %in% meta$label_observed),
            all(!is.na(meta$tissue)), all(nzchar(meta$tissue)),
            all(is.finite(meta$time)),
            identical(rownames(sim$X), meta$sample))
}
set.seed(42)
clean <- synthetic$simulate_dataset(
  "gaussian", 0.2, 0, 0, n_per_label = 6, n_labels = 5,
  n_features = 20, misid_rate = 0, n_unknown_individuals = 0
)
stopifnot(nrow(clean$meta) == 30L, !any(clean$meta$is_unknown),
          !any(clean$meta$is_misid),
          all(!is.na(clean$meta$tissue)), all(nzchar(clean$meta$tissue)),
          all(is.finite(clean$meta$time)))
small_grid <- data.table(type = "gaussian", var_label = 0.5,
                         var_discrete = 0, var_cont = 0, n_per_label = 6,
                         misid_rate = 0, seed = 42L, scenario = 1L,
                         n_unknown_individuals = 1L,
                         unknown_samples_per_individual = 2L)
generated <- synthetic$simulate_grid(small_grid, n_labels = 5, n_features = 20)[[1L]]
stopifnot(sum(generated$meta$is_unknown) == 2L,
          generated$params$n_unknown_individuals == 1L)
set.seed(43)
categorical_covariate <- synthetic$simulate_dataset(
  "gaussian", 0.2, 0.1, 0.1, n_per_label = 6, n_labels = 5,
  n_features = 20, misid_rate = 0, n_unknown_individuals = 1,
  unknown_samples_per_individual = 2, n_tissues = 3L, var_tissue = 0.2
)
stopifnot(uniqueN(categorical_covariate$meta$tissue) == 3L,
          all(!is.na(categorical_covariate$meta$tissue)),
          all(is.finite(categorical_covariate$meta$time)))
set.seed(44)
numeric_covariate <- synthetic$simulate_dataset(
  "gaussian", 0.2, 0.1, 0.1, n_per_label = 6, n_labels = 5,
  n_features = 20, misid_rate = 0, n_unknown_individuals = 1,
  unknown_samples_per_individual = 2, var_numeric_covariate = 0.2
)
stopifnot(all(is.finite(numeric_covariate$meta$time)),
          diff(range(numeric_covariate$meta$time)) > 0)
set.seed(45)
both_covariates <- synthetic$simulate_dataset(
  "gaussian", 0.2, 0.1, 0.1, n_per_label = 6, n_labels = 5,
  n_features = 20, misid_rate = 0, n_unknown_individuals = 1,
  unknown_samples_per_individual = 2, n_tissues = 3L,
  var_tissue = 0.2, var_numeric_covariate = 0.2
)
stopifnot(uniqueN(both_covariates$meta$tissue) == 3L,
          all(is.finite(both_covariates$meta$time)))
stopifnot(inherits(try(synthetic$simulate_dataset(
  "gaussian", 0.2, 0.1, 0.1, n_per_label = 6,
  n_tissues = 1L, var_tissue = 0.2
), silent = TRUE), "try-error"))
stopifnot(inherits(try(synthetic$simulate_dataset(
  "gaussian", 0.2, 0.1, 0.1, n_per_label = 6,
  var_tissue = 0.8
), silent = TRUE), "try-error"))
for (file in c("R/lda.R", "R/validation.R", "R/check_annotations.R")) {
  sys.source(file, envir = synthetic)
}
set.seed(42)
result <- suppressWarnings(synthetic$check_annotations(
  t(generated$X), as.data.frame(generated$meta), "sample", "label_observed",
  progress = FALSE, min_agreement = 0,
  mismatch_rate = 0.1, detect_unknown = TRUE
))
stopifnot(nrow(result$predictions) == nrow(generated$meta),
          !is.null(result$performance$unknown_individual_detection))
set.seed(42)
candidate <- synthetic$check_annotations(
  t(generated$X), as.data.frame(generated$meta), "sample", "label_observed",
  progress = FALSE, min_agreement = 0, mismatch_rate = 0.1,
  detect_unknown = TRUE, diagnostics = TRUE
)
stopifnot(nrow(candidate$predictions) == nrow(generated$meta),
          identical(candidate$diagnostics$unknown$statistic,
                    "nearest_class_distance"),
          all(is.finite(candidate$predictions$unknown_detection_score)),
          length(candidate$diagnostics$unknown$known_values) > 0L,
          length(candidate$diagnostics$unknown$simulated_unknown_values) > 0L)
covariate_meta <- as.data.frame(generated$meta)
covariate_meta$tissue <- rep(c("tissue_a", "tissue_b"),
                             length.out = nrow(covariate_meta))
set.seed(42)
covariate_result <- synthetic$check_annotations(
  t(generated$X), covariate_meta, "sample", "label_observed",
  progress = FALSE, min_agreement = 0, mismatch_rate = 0.1,
  detect_unknown = TRUE, covariate_col = "tissue"
)
stopifnot(
  identical(covariate_result$covariate_adjustment$method, "residualize"),
  identical(covariate_result$covariate_adjustment$covariates$tissue$type,
            "categorical"),
  setequal(covariate_result$covariate_adjustment$covariates$tissue$levels,
           c("tissue_a", "tissue_b")),
  nrow(covariate_result$predictions) == nrow(covariate_meta)
)
set.seed(46)
both_covariate_result <- synthetic$check_annotations(
  t(both_covariates$X), as.data.frame(both_covariates$meta),
  "sample", "label_observed", progress = FALSE, min_agreement = 0,
  mismatch_rate = 0.1, detect_unknown = TRUE,
  covariate_col = c("tissue", "time")
)
stopifnot(
  identical(names(both_covariate_result$covariate_adjustment$covariates),
            c("tissue", "time")),
  identical(both_covariate_result$covariate_adjustment$covariates$tissue$type,
            "categorical"),
  identical(both_covariate_result$covariate_adjustment$covariates$time$type,
            "numeric")
)
covariate_x <- matrix(c(1, 3, 5, 7, 2, 4, 6, 8), nrow = 4L)
covariate_residual <- synthetic$.covariate_residualize(
  covariate_x,
  data.frame(tissue = c("a", "a", "b", "b"))
)
stopifnot(
  max(abs(colMeans(covariate_residual[c(1, 2), , drop = FALSE]))) < 1e-10,
  max(abs(colMeans(covariate_residual[c(3, 4), , drop = FALSE]))) < 1e-10
)
missing_covariate_meta <- covariate_meta
missing_covariate_meta$tissue[1] <- NA_character_
stopifnot(inherits(try(synthetic$check_annotations(
  t(generated$X), missing_covariate_meta, "sample", "label_observed",
  progress = FALSE, covariate_col = "tissue"
), silent = TRUE), "try-error"))
stopifnot(inherits(try(synthetic$check_annotations(
  t(generated$X), covariate_meta, "sample", "label_observed",
  progress = FALSE, covariate_col = "missing_column"
), silent = TRUE), "try-error"))
stopifnot(inherits(try(synthetic$check_annotations(
  t(generated$X), as.data.frame(generated$meta), "sample", "label_observed",
  progress = FALSE, mismatch_rate = 0
), silent = TRUE), "try-error"))
stopifnot(inherits(try(synthetic$simulate_dataset(
  "gaussian", 0.2, 0, 0, n_per_label = 6,
  unknown_samples_per_individual = 7
), silent = TRUE), "try-error"))

gtex <- new.env(parent = globalenv())
samples <- CJ(donor = seq_len(30), tissue = seq_len(8))
samples[, sample := sprintf("GTEX-%04d-%02d", donor, tissue)]
gtex$sample_info <- samples[, .(
  SAMPID = sample, SMTS = paste0("tissue", tissue),
  SMTSD = paste0("tissue", tissue)
)]
set.seed(1)
gtex$expr_vst <- matrix(rnorm(20L * nrow(samples)), nrow = 20L,
                        dimnames = list(NULL, samples$sample))
gtex$min_tissue_samples <- 20L
gtex$n_unknown_individuals <- 2L
gtex$n_unknown_tissues <- 1L
gtex$unknown_samples_per_individual <- 2L
sys.source("simulations/Simulate_from_GTEX_grid.R", envir = gtex)
scenario <- data.table(task = "individual", n_labels = 10L, n_per_label = 4,
                       set = "random", misid_rate = 0, seed = 42L,
                       scenario = 1L)
meta <- gtex$build_scenario(scenario)
known <- meta[is_unknown == FALSE]
unknown <- meta[is_unknown == TRUE]
stopifnot(nrow(known) == 40L, nrow(unknown) == 4L,
          uniqueN(unknown$label_true) == 2L,
          all(table(unknown$label_true) == 2L),
          all(unknown$is_misid), !any(known$is_misid),
          all(unknown$label_observed %in% known$label_true),
          !any(unknown$label_true %in% meta$label_observed),
          !anyDuplicated(meta$sample),
          isTRUE(all.equal(as.data.frame(meta),
                           as.data.frame(gtex$build_scenario(scenario)),
                           check.attributes = FALSE)))
saved <- gtex$simulate_gtex_grid(scenario)[[1L]]
stopifnot(saved$params$n_unknown_individuals == 2L,
          saved$params$unknown_samples_per_individual == 2L)
gtex$n_unknown_individuals <- 0L
stopifnot(!any(gtex$build_scenario(scenario)$is_unknown))
matched_individual <- copy(scenario)
matched_individual[, `:=`(
  matching = "paired",
  n_labels = 10L,
  n_per_label = 4L
)]
gtex$n_unknown_individuals <- 0L
matched_individual_meta <- gtex$build_scenario(matched_individual)
stopifnot(nrow(matched_individual_meta) == 40L,
          all(matched_individual_meta[, uniqueN(SMTSD), by = ID][
            , V1] == 4L),
          uniqueN(matched_individual_meta$SMTSD) <= 4L)
gtex$n_unknown_individuals <- 100L
stopifnot(inherits(try(gtex$build_scenario(scenario), silent = TRUE),
                   "try-error"))
gtex$n_unknown_individuals <- 2L
tissue <- copy(scenario)
tissue[, `:=`(task = "tissue", n_labels = 3L, n_per_label = 4,
        n_unknown_tissues = 1L)]
tissue_meta <- gtex$build_scenario(tissue)
tissue_unknown <- tissue_meta[is_unknown == TRUE]
tissue_known <- tissue_meta[is_unknown == FALSE]
stopifnot(nrow(tissue_unknown) == 2L,
      uniqueN(tissue_unknown$label_true) == 1L,
      !any(tissue_unknown$label_true %in% tissue_known$label_true),
      all(tissue_unknown$is_misid),
      all(tissue_unknown$label_observed %in% tissue_known$label_true),
      !any(tissue_unknown$label_true %in% tissue_meta$label_observed),
      uniqueN(tissue_unknown$ID) == 2L)
matched_tissue <- copy(tissue)
matched_tissue[, `:=`(matching = "paired", n_unknown_tissues = 0L)]
gtex$n_unknown_tissues <- 0L
matched_tissue_meta <- gtex$build_scenario(matched_tissue)
matched_tissue_ids <- matched_tissue_meta[, .(ID = unique(ID))]
common_across_tissues <- matched_tissue_meta[
  , .(n_tissues = uniqueN(SMTSD)), by = ID]
stopifnot(nrow(matched_tissue_meta) ==
            uniqueN(matched_tissue_meta$SMTS) *
            matched_tissue$n_per_label,
          all(common_across_tissues$n_tissues > 0L),
          all(matched_tissue_meta[, .(n = uniqueN(ID)), by = SMTS][, n] ==
              matched_tissue$n_per_label),
          length(Reduce(intersect,
                        split(matched_tissue_meta$ID,
                              matched_tissue_meta$SMTS))) ==
            matched_tissue$n_per_label)
tissue_saved <- gtex$simulate_gtex_grid(tissue)[[1L]]
stopifnot(tissue_saved$params$unknown_samples_per_tissue == 2L,
          sum(tissue_saved$meta$is_unknown) == 2L)
tissue[, n_unknown_tissues := 0L]
stopifnot(!any(gtex$build_scenario(tissue)$is_unknown))
gtex$n_unknown_tissues <- 100L
stopifnot(inherits(try(gtex$build_scenario(transform(tissue, n_unknown_tissues = 100L)),
             silent = TRUE), "try-error"))
gtex_grid_counts <- gtex$grid[, .N, by = .(task, n_unknown_individuals,
                       n_unknown_tissues)]
stopifnot(setequal(gtex_grid_counts$task, c("individual", "tissue")),
      nrow(gtex_grid_counts) == 4L,
      all(gtex_grid_counts[task == "individual",
                 n_unknown_tissues == 0L]),
      setequal(gtex_grid_counts[task == "individual", n_unknown_individuals],
           c(0L, 2L)),
      all(gtex_grid_counts[task == "tissue",
                 n_unknown_individuals == 0L]),
      setequal(gtex_grid_counts[task == "tissue", n_unknown_tissues],
           c(0L, 1L)),
          !anyDuplicated(gtex$grid$scenario))
gtex$matching_modes <- c("random", "paired")
sys.source("simulations/Simulate_from_GTEX_grid.R", envir = gtex)
stopifnot(setequal(gtex$grid$matching, c("random", "paired")),
          length(unique(gtex$grid$scenario)) == nrow(gtex$grid))

source("simulations/score_unknown.R")
truth <- data.frame(sample = letters[1:6],
                    is_unknown = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))
predictions <- data.frame(sample_id = letters[c(6, 1, 4, 2, 5)],
                          potential_unknown_individual = c(FALSE, TRUE, TRUE, FALSE, NA))
scored <- score_unknown(predictions, truth, TRUE)
stopifnot(scored$unknown_n_true == 3, scored$unknown_n_known == 3,
          scored$unknown_n_assessed_true == 2, scored$unknown_n_assessed_known == 2,
          scored$unknown_sens_true == 0.5, scored$unknown_spec_true == 0.5)
disabled <- score_unknown(predictions, truth, FALSE)
stopifnot(is.na(disabled$unknown_sens_true), is.na(disabled$unknown_spec_true),
          disabled$unknown_n_assessed_true == 0, disabled$unknown_n_assessed_known == 0)
no_unknown <- score_unknown(predictions, transform(truth, is_unknown = FALSE), TRUE)
stopifnot(is.na(no_unknown$unknown_sens_true), no_unknown$unknown_spec_true == 0.5)
stopifnot(inherits(try(score_unknown(predictions, truth[, "sample", drop = FALSE], TRUE),
                       silent = TRUE), "try-error"))

tmp <- tempfile("unknown-rescore-")
dir.create(tmp)
dir.create(file.path(tmp, "dirty"))
row_file <- file.path(tmp, "scenario_1__cfg.rds")
saveRDS(data.table(status = "ok"), row_file)
saveRDS(list(result = list(predictions = predictions), truth = truth,
             check_args = list(detect_unknown = TRUE)),
        file.path(tmp, "dirty", basename(row_file)))
rescored <- rescore_unknown_rows(row_file, file.path(tmp, "dirty"))[[1L]]
stopifnot(rescored$unknown_sens_true == 0.5, rescored$unknown_spec_true == 0.5)
unlink(file.path(tmp, "dirty", basename(row_file)))
unlink(file.path(tmp, "dirty"))
unlink(row_file)
unlink(tmp)
cat("Unknown simulation tests passed.\n")
