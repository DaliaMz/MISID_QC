custom_lib <- "/scratch/project/palpant_scratch/DMizikovsky/R_packages"
.libPaths(c(custom_lib, .libPaths()))

# Run check_annotations() on every scenario in results/gtex_scenarios.rds
# (made by Simulate_from_GTEX_grid.R) and score it against the known truth.
# Same outputs as Run_across_simulated_data.R, for the GTEx scenarios.
#
# Usage (repo root, package installed with R CMD INSTALL .): submit as a Slurm
# job array with simulations/run_gtex.sh, or run a slice by hand:
#   Rscript simulations/Run_across_GTEX_simulations.R task_id=1 n_tasks=20
# Task k of n handles scenarios k, k+n, k+2n, ... (interleaved so slow
# scenarios spread evenly). Without arguments it runs all scenarios in one
# process. Inside a Slurm array, task_id defaults to SLURM_ARRAY_TASK_ID.
# Combine the results afterwards with simulations/Aggregate_gtex_results.R.
#
# Arguments (key=value, all optional):
#   task_id, n_tasks  which slice to run (defaults: 1 of 1)
#   scenarios_file    GTEx scenarios      (default simulations/results/gtex_scenarios.rds)
#   expr              VST matrix, genes x samples (default: the cluster path below)
#   out               output directory    (default simulations/results)
#   n_features        most variable genes kept per dataset (default 2000)
#   subset            data.table expression on the scenario parameters
#   scenarios         R expression for scenario numbers, e.g. 1:20
#   cv_folds, mismatch_repeats, target_specificity
#                     comma-separated values to compare; every combination is
#                     run, e.g.  cv_folds=5,loo mismatch_repeats=1,3
#                     target_specificity=0.9,0.95,0.99   ("loo" = leave-one-out)
#                     defaults: 5, 1, 0.95
#   force             TRUE to recompute scenarios that are already done
#
# Outputs (one file per scenario, so tasks never write to the same file):
#   <out>/rows_gtex/scenario_<n>__<config>.rds   summary row (parameters + metrics)
#   <out>/dirty_gtex/scenario_<n>__<config>.rds  full check_annotations() result
#                                      plus the truth (disposable)
# Finished scenarios are skipped on a rerun, so a requeued or timed-out array
# task resumes where it stopped. Each scenario is seeded, so a result does not
# depend on how the scenarios are split across tasks.
suppressMessages({library(data.table); library(SampleAnnotationLikelihood)})

# ---- Settings -----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
opt <- list(scenarios_file = "simulations/results/gtex_scenarios.rds",
            expr = "/scratch/project_mnt/S0138/DMizikovsky/misid/GTEX/GTEx_v11_VST.rds",
            out = "simulations/results", n_features = "2000",
            task_id = Sys.getenv("SLURM_ARRAY_TASK_ID", "1"), n_tasks = "1",
            subset = "", scenarios = "", cv_folds = "5", mismatch_repeats = "1", target_specificity = "0.95",
            force = "FALSE")
for (a in args[grepl("=", args)]) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  opt[[kv[1]]] <- paste(kv[-1], collapse = "=")
}
check_args <- list(detect_unknown = FALSE)   # fixed arguments for check_annotations()
# Settings to compare: every combination is run on every scenario, on the same
# seed, so differences between configs are paired. cv_folds: a number, or
# "loo" for leave-one-sample-out (cv_folds = NULL in the package). A number
# is capped at the smallest observed label count of the scenario.
as_vals <- function(x) strsplit(x, ",", fixed = TRUE)[[1]]
configs <- CJ(cv_folds = as_vals(opt$cv_folds),
              mismatch_repeats = as.integer(as_vals(opt$mismatch_repeats)),
              target_specificity = as.numeric(as_vals(opt$target_specificity)),
              sorted = FALSE)
configs[, cfg := sprintf("cv%s_mr%d_ts%s", cv_folds, mismatch_repeats, target_specificity)]
task_id <- as.integer(opt$task_id); n_tasks <- as.integer(opt$n_tasks)
stopifnot(task_id >= 1L, task_id <= n_tasks)
dirty_dir <- file.path(opt$out, "dirty_gtex")
rows_dir <- file.path(opt$out, "rows_gtex")
for (d in c(dirty_dir, rows_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Scoring helpers (same as Run_across_simulated_data.R) ----------------------

auroc <- function(score, positive) {
  n1 <- sum(positive); n0 <- sum(!positive)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  (sum(rank(score)[positive]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
ratio <- function(a, b) if (b > 0) a / b else NA_real_
all_statuses <- c("insufficient_evidence", "potential_unknown_individual",
                  "suspicious_with_candidate_correction",
                  "suspicious_without_reliable_alternative", "no_confident_prediction",
                  "weak_disagreement", "no_evidence_of_incorrect_annotation")

score_scenario <- function(r, meta, seconds) {
  p <- r$predictions
  m <- match(p$sample_id, meta$sample)
  bad <- meta$is_misid[m]
  truth <- as.character(meta$label_true[m])
  low <- isTRUE(r$calibration$low_agreement)
  flagged <- if (low) rep(FALSE, nrow(p)) else p$suspicious
  tp <- sum(flagged & bad); fp <- sum(flagged & !bad)
  fn <- sum(!flagged & bad); tn <- sum(!flagged & !bad)
  held <- r$held_out$sample_ids
  held_bad <- meta$is_misid[match(held, meta$sample)]
  corr <- !is.na(p$candidate_correction)
  conf <- !is.na(p$predicted_id_clean)
  sens <- ratio(tp, tp + fn); spec <- ratio(tn, tn + fp)
  rep_m <- r$performance$misannotation_detection
  verdict <- r$feasibility$misannotation_detection
  st <- table(factor(p$annotation_status, levels = all_statuses))
  out <- data.table(
    status = if (low) "low_agreement" else "ok",
    n_samples = nrow(p), n_labels = length(unique(p$observed_id)),
    n_misid = sum(bad),
    first_pass_agreement = r$calibration$first_pass_agreement,
    auroc_annotation_margin = auroc(p$annotation_score_margin, bad),
    auroc_top_margin = auroc(-p$top_score_margin, bad),
    tp = tp, fp = fp, fn = fn, tn = tn,
    sens_true = sens, spec_true = spec, precision_true = ratio(tp, tp + fp),
    sens_reported = rep_m$sensitivity, spec_reported = rep_m$specificity,
    sens_gap = rep_m$sensitivity - sens, spec_gap = rep_m$specificity - spec,
    verdict = verdict,
    verdict_correct = if (low) NA else (verdict == "appears_feasible") ==
      (isTRUE(sens >= 0.8) && isTRUE(spec >= 0.95)),
    recovery_verdict = r$feasibility$identity_recovery,
    n_held = length(held), n_held_misid = sum(held_bad),
    held_precision = ratio(sum(held_bad), length(held)),
    held_recall = ratio(sum(held_bad), sum(bad)),
    n_corrections = sum(corr),
    correction_accuracy = ratio(sum(p$candidate_correction[corr] == truth[corr]), sum(corr)),
    correction_recall = ratio(sum(p$candidate_correction[corr & bad] == truth[corr & bad]), sum(bad)),
    prediction_coverage = mean(conf),
    clean_accuracy_true = ratio(sum(p$predicted_id_clean[conf] == truth[conf]), sum(conf)),
    raw_accuracy_true = mean(p$predicted_id == truth),
    margin_threshold = r$thresholds$annotation_margin,
    abstain_threshold = r$thresholds$abstain_margin,
    exclusion_threshold = r$calibration$exclusion_margin_threshold,
    seconds = seconds
  )
  cbind(out, as.data.table(as.list(setNames(as.integer(st), paste0("n_", names(st))))))
}

# ---- Running one scenario -----------------------------------------------------

# Features x samples matrix with the most variable genes within the dataset.
slice_expression <- function(expr, samples, n_features) {
  x <- expr[, samples, drop = FALSE]
  v <- apply(x, 1, var)
  x[order(v, decreasing = TRUE)[seq_len(min(n_features, nrow(x)))], , drop = FALSE]
}

run_config <- function(sim, id, cfg, meta, X, md) {
  f <- file.path(dirty_dir, sprintf("scenario_%s__%s.rds", id, cfg$cfg))
  folds <- if (cfg$cv_folds == "loo") NULL else
    min(as.integer(cfg$cv_folds), min(table(md$label)))
  warn <- character()
  set.seed(sim$params$seed)   # same draws for every config of a scenario
  t0 <- Sys.time()
  args <- c(list(expr = X, metadata = md, sample_col = "sample", id_col = "label",
                 cv_folds = folds, progress = FALSE,
                 mismatch_repeats = cfg$mismatch_repeats,
                 target_specificity = cfg$target_specificity), check_args)
  r <- tryCatch(
    withCallingHandlers(do.call(check_annotations, args),
      warning = function(w) {
        warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning")
      }),
    error = function(e) e)
  seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(r, "error")) {
    saveRDS(list(error = conditionMessage(r), truth = sim$meta, params = sim$params), f)
    return(data.table(status = "error", error = conditionMessage(r), seconds = seconds))
  }
  res <- score_scenario(r, meta, seconds)
  res[, warnings := paste(unique(warn), collapse = " | ")]
  saveRDS(list(result = r, truth = sim$meta, params = sim$params, config = cfg,
               check_args = check_args, cv_folds = folds), f)
  res
}

row_file <- function(id, cfg) file.path(rows_dir, sprintf("scenario_%s__%s.rds", id, cfg))

# Runs the configs that are not finished yet for one scenario.
run_scenario <- function(sim, id, expr) {
  todo <- if (as.logical(opt$force)) configs$cfg else
    configs$cfg[!file.exists(row_file(id, configs$cfg))]
  if (!length(todo)) return(invisible())
  meta <- as.data.frame(sim$meta)
  meta$label_true <- as.character(meta$label_true)
  X <- slice_expression(expr, meta$sample, as.integer(opt$n_features))
  md <- data.frame(sample = meta$sample, label = as.character(meta$label_observed))
  for (k in todo) {
    cfg <- configs[cfg == k]
    res <- tryCatch(run_config(sim, id, cfg, meta, X, md),
                    error = function(e) data.table(status = "error", error = conditionMessage(e)))
    saveRDS(cbind(params[scenario == id], cfg, res), row_file(id, k))
    cat(sprintf("scenario %d [%s]: %s\n", id, k, res$status))
  }
}

# ---- Run ----------------------------------------------------------------------

sims <- readRDS(opt$scenarios_file)
params <- rbindlist(lapply(sims, `[[`, "params"), fill = TRUE)
keep <- rep(TRUE, nrow(params))
if (nzchar(opt$subset)) keep <- keep & params[, eval(parse(text = opt$subset))]
if (nzchar(opt$scenarios)) keep <- keep & params$scenario %in% eval(parse(text = opt$scenarios))
# Slice by position among the selected scenarios, so every task gets a share.
mine <- which(keep)
mine <- mine[(seq_along(mine) - 1L) %% n_tasks == task_id - 1L]
if (!as.logical(opt$force)) {
  # drop scenarios whose configs are all finished
  all_done <- vapply(params$scenario[mine], function(id) {
    all(file.exists(row_file(id, configs$cfg)))
  }, logical(1))
  mine <- mine[!all_done]
}
cat(sprintf("task %d/%d: %d scenarios x %d configs (%d scenarios in the grid)\n",
            task_id, n_tasks, length(mine), nrow(configs), nrow(params)))

if (length(mine)) {
  cat("Loading expression:", opt$expr, "\n")
  expr <- readRDS(opt$expr)
  needed <- unique(unlist(lapply(sims[mine], function(s) s$meta$sample)))
  if (!all(needed %in% colnames(expr))) {
    stop("Some scenario samples are missing from the expression matrix.")
  }
  expr <- expr[, needed, drop = FALSE]  # keep only what this task uses

  for (i in mine) {
    run_scenario(sims[[i]], params$scenario[i], expr)
  }
}
cat("task", task_id, "done\n")
