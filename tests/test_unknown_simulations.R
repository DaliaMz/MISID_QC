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
            identical(rownames(sim$X), meta$sample))
}
set.seed(42)
clean <- synthetic$simulate_dataset(
  "gaussian", 0.2, 0, 0, n_per_label = 6, n_labels = 5,
  n_features = 20, misid_rate = 0, n_unknown_individuals = 0
)
stopifnot(nrow(clean$meta) == 30L, !any(clean$meta$is_unknown),
          !any(clean$meta$is_misid))
small_grid <- data.table(type = "gaussian", var_label = 0.5,
                         var_discrete = 0, var_cont = 0, n_per_label = 6,
                         misid_rate = 0, seed = 42L, scenario = 1L,
                         n_unknown_individuals = 1L,
                         unknown_samples_per_individual = 2L)
generated <- synthetic$simulate_grid(small_grid, n_labels = 5, n_features = 20)[[1L]]
stopifnot(sum(generated$meta$is_unknown) == 2L,
          generated$params$n_unknown_individuals == 1L)
for (file in c("R/lda.R", "R/validation.R", "R/check_annotations.R")) {
  sys.source(file, envir = synthetic)
}
set.seed(42)
result <- suppressWarnings(synthetic$check_annotations(
  t(generated$X), as.data.frame(generated$meta), "sample", "label_observed",
  progress = FALSE, min_agreement = 0, stop_on_low_agreement = FALSE,
  mismatch_rate = 0, detect_unknown = TRUE
))
stopifnot(nrow(result$predictions) == nrow(generated$meta),
          !is.null(result$performance$unknown_individual_detection))
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
stopifnot(nrow(tissue_unknown) == 4L,
      uniqueN(tissue_unknown$label_true) == 1L,
      !any(tissue_unknown$label_true %in% tissue_known$label_true),
      all(tissue_unknown$is_misid),
      all(tissue_unknown$label_observed %in% tissue_known$label_true),
      !any(tissue_unknown$label_true %in% tissue_meta$label_observed),
      uniqueN(tissue_unknown$ID) == 4L)
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
cat("Unknown simulation tests passed.\n")
