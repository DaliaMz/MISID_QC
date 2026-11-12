# Run check_annotations() on every scenario in results/simulated_datasets_v2.rds
# (made by Simulate_full_dataset.R) and score it against the known truth.
#
# Usage (repo root, package installed with R CMD INSTALL .):
#   Rscript simulations/Run_across_simulated_data.R
#   Rscript simulations/Run_across_simulated_data.R cores=4 subset="rep == 1"
#   Rscript simulations/Run_across_simulated_data.R scenarios=1:20 force=TRUE
# Or as a Slurm job array with simulations/run_sims.sh. Task k of n handles
# scenarios k, k+n, k+2n, ... (interleaved so slow scenarios spread evenly) and
# writes one row file per scenario and config under <out>/rows, so tasks never
# write to the same file. After all tasks finish, build the summary with
#   Rscript simulations/Run_across_simulated_data.R combine=TRUE
# A single process (n_tasks=1) combines its own results at the end.
#
# Arguments (key=value, all optional):
#   input      simulated datasets        (default simulations/results/simulated_datasets_v2.rds)
#   out        output directory          (default simulations/results)
#   cores      parallel workers inside one process (default 1)
#   task_id, n_tasks  which slice to run (defaults: SLURM_ARRAY_TASK_ID, 1 of 1)
#   combine    TRUE to only combine the finished row files into the summary
#   rescore_unknown TRUE with combine=TRUE to score saved dirty results again
#   subset     data.table expression on the scenario parameters, e.g. "rep == 1"
#   scenarios  R expression for scenario numbers, e.g. 1:20
#   cv_folds, mismatch_repeats, target_specificity, margin_coverage,
#   holdout_reference
#              comma-separated values to compare; every combination is run, e.g.
#              cv_folds=5,loo mismatch_repeats=1,3 target_specificity=0.9,0.95,0.99
#              margin_coverage=0.01,0.025,0.05
#              holdout_reference=all_alternatives,top
#              ("loo" = leave-one-out). Defaults: 5, 3, 0.95, 0.025
#   detect_unknown TRUE/FALSE; comma-separated values to compare (default FALSE)
#   force      TRUE to recompute scenarios that already have a row file
#   min_agreement  agreement gate passed to check_annotations() (default 0 =
#              off; first_pass_agreement is recorded to apply cutoffs later)
#
# Outputs:
#   <out>/summary_v2.rds, .csv   one row per scenario and config (parameters + metrics)
#   <out>/rows/scenario_<n>__<config>.rds  one scored row; the resume point
#   <out>/dirty/scenario_<n>__<config>.rds the full check_annotations() result plus the
#                                truth for that scenario. These are large and
#                                disposable: delete the folder at any time.
# Finished scenarios and configs are skipped on a rerun (use force=TRUE to redo
# them, or a new out= directory after changing the package), so an interrupted
# run can be resumed.
suppressMessages({library(data.table); library(SampleAnnotationLikelihood)})
source("simulations/score_unknown.R")

# ---- Settings -----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
opt <- list(input = "simulations/results/simulated_datasets_v2.rds",
            out = "simulations/results", cores = "1", subset = "", scenarios = "",
            task_id = Sys.getenv("SLURM_ARRAY_TASK_ID", "1"), n_tasks = "1",
            combine = "FALSE", rescore_unknown = "FALSE",
            cv_folds = "5", mismatch_repeats = "3", target_specificity = "0.95",
            margin_coverage = "0.025",
            holdout_reference = "all_alternatives",
            detect_unknown = "FALSE",
            min_agreement = "0", force = "FALSE")
for (a in args[grepl("=", args)]) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (kv[1] == "refit_clean") stop("refit_clean has been removed; models are not refit.")
  opt[[kv[1]]] <- paste(kv[-1], collapse = "=")
}
# min_agreement = 0 disables the agreement gate so every scenario is fully
# processed; first_pass_agreement is recorded to apply any cutoff afterwards.
# Settings to compare: every combination is run on every scenario, on the same
# seed, so differences between configs are paired. cv_folds: a number, or
# "loo" for leave-one-sample-out (cv_folds = NULL in the package). A number
# is capped at the smallest observed label count of the scenario.
as_vals <- function(x) strsplit(x, ",", fixed = TRUE)[[1]]
configs <- CJ(cv_folds = as_vals(opt$cv_folds),
              mismatch_repeats = as.integer(as_vals(opt$mismatch_repeats)),
              target_specificity = as.numeric(as_vals(opt$target_specificity)),
              margin_coverage = as.numeric(as_vals(opt$margin_coverage)),
              holdout_reference = as_vals(opt$holdout_reference),
              detect_unknown = as.logical(as_vals(opt$detect_unknown)),
              sorted = FALSE)
# The suffix distinguishes cutoff-only results from historical runs.
          configs[, cfg := sprintf("cv%s_mr%d_ts%s%s%s%s_cutoffonly", cv_folds, mismatch_repeats,
                         target_specificity,
                         ifelse(margin_coverage == 0.025, "", paste0("_mc", margin_coverage)),
                         ifelse(holdout_reference == "all_alternatives", "",
                                paste0("_hr", holdout_reference)),
                    ifelse(detect_unknown, "_du", ""))]
dirty_dir <- file.path(opt$out, "dirty")
rows_dir <- file.path(opt$out, "rows")
for (d in c(dirty_dir, rows_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
task_id <- as.integer(opt$task_id); n_tasks <- as.integer(opt$n_tasks)
stopifnot(task_id >= 1L, task_id <= n_tasks)

# ---- Scoring helpers ----------------------------------------------------------

auroc <- function(score, positive) {
  n1 <- sum(positive); n0 <- sum(!positive)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  (sum(rank(score)[positive]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
ratio <- function(a, b) if (b > 0) a / b else NA_real_
or_na <- function(x) if (is.null(x) || !length(x)) NA_real_ else x
all_statuses <- c("insufficient_evidence", "potential_unknown_individual",
                  "suspicious_with_candidate_correction",
                  "suspicious_without_reliable_alternative", "no_confident_prediction",
                  "weak_disagreement", "no_evidence_of_incorrect_annotation")

score_scenario <- function(r, meta, seconds, detect_unknown) {
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
  perf <- r$performance
  det <- rep_m
  all_alt <- perf$all_alternatives_detection
  corr_v <- perf$correction_validation
  rec_true <- ratio(sum(p$predicted_id[bad] == truth[bad]), sum(bad))
  verdict <- r$feasibility$misannotation_detection
  true_detection_feasible <- if (!any(bad)) {
    isTRUE(spec >= 0.95)
  } else {
    isTRUE(sens >= 0.8) && isTRUE(spec >= 0.95)
  }
  st <- table(factor(p$annotation_status, levels = all_statuses))
  out <- data.table(
    status = if (low) "low_agreement" else "ok",
    n_samples = nrow(p), n_labels = length(unique(p$observed_id)),
    n_misid = sum(bad),
    first_pass_agreement = r$calibration$first_pass_agreement,
    # separation: how well the raw annotation margin ranks real mislabels
    auroc_annotation_margin = auroc(p$annotation_score_margin, bad),
    # correction: among truly mislabelled samples, does a larger top margin rank
    # the predictions that recover the true identity higher?
    auroc_top_margin = auroc(p$top_score_margin[bad],
                             p$predicted_id[bad] == truth[bad]),
    tp = tp, fp = fp, fn = fn, tn = tn,
    sens_true = sens, spec_true = spec, precision_true = ratio(tp, tp + fp),
    sens_reported = rep_m$sensitivity, spec_reported = rep_m$specificity,
    sens_gap = rep_m$sensitivity - sens, spec_gap = rep_m$specificity - spec,
    verdict = verdict,
    # Sensitivity is undefined for clean scenarios, so score those on specificity.
    verdict_correct = if (low) NA else (verdict == "appears_feasible") ==
      true_detection_feasible,
    recovery_verdict = r$feasibility$identity_recovery,
    recovery_verdict_correct = if (low || is.na(rec_true)) NA else
      (r$feasibility$identity_recovery == "appears_feasible") == (rec_true >= 0.8),
    # detection benchmark behind the verdict: where sensitivity came from, its
    # in-sample value and spread, and the every-wrong-label diagnostic
    detection_source = or_na(det$source),
    n_repeats_validated = or_na(det$n_repeats_validated),
    sens_in_sample = or_na(det$sensitivity_in_sample),
    sens_sd = or_na(det$sensitivity_sd),
    sens_all_alt = or_na(all_alt$sensitivity),
    spec_all_alt = or_na(all_alt$specificity),
    sens_all_alt_gap = or_na(all_alt$sensitivity) - sens,
    # recovery on the truly misannotated samples: reported (simulated swaps)
    # versus true
    recovery_true_misid = rec_true,
    recovery_reported_misid = or_na(perf$recovery_accuracy_misannotated),
    recovery_gap = or_na(perf$recovery_accuracy_misannotated) - rec_true,
    # proposed corrections: reported (simulated swaps) versus true
    n_corrections_reported = or_na(corr_v$n_proposed),
    correction_precision_reported = or_na(corr_v$precision),
    correction_recall_reported = or_na(corr_v$recall),
    correction_precision_gap = or_na(corr_v$precision) -
      ratio(sum(p$candidate_correction[corr] == truth[corr]), sum(corr)),
    correction_recall_gap = or_na(corr_v$recall) -
      ratio(sum(p$candidate_correction[corr & bad] == truth[corr & bad]), sum(bad)),
    n_held = length(held), n_held_misid = sum(held_bad),
    held_precision = ratio(sum(held_bad), length(held)),
    # Held-out samples and their denominators cover all input samples.
    n_samples_input = nrow(meta), n_misid_input = sum(meta$is_misid),
    held_recall = ratio(sum(held_bad), sum(meta$is_misid)),
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
  cbind(out, as.data.table(score_unknown(p, meta, detect_unknown)),
        as.data.table(as.list(setNames(as.integer(st), paste0("n_", names(st))))))
}

run_config <- function(sim, id, cfg, meta, X, md) {
  f <- file.path(dirty_dir, sprintf("scenario_%s__%s.rds", id, cfg$cfg))
  folds <- if (cfg$cv_folds == "loo") NULL else
    min(as.integer(cfg$cv_folds), min(table(md$label)))
  warn <- character()
  set.seed(sim$params$seed)   # same draws for every config of a scenario
  t0 <- Sys.time()
  check_args <- list(detect_unknown = cfg$detect_unknown,
                     min_agreement = as.numeric(opt$min_agreement))
  args <- c(list(expr = X, metadata = md, sample_col = "sample", id_col = "label",
                 cv_folds = folds, progress = FALSE,
                 mismatch_repeats = cfg$mismatch_repeats,
                 margin_coverage = cfg$margin_coverage,
                 holdout_reference = cfg$holdout_reference,
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
  res <- score_scenario(r, meta, seconds, cfg$detect_unknown)
  res[, warnings := paste(unique(warn), collapse = " | ")]
  saveRDS(list(result = r, truth = sim$meta, params = sim$params, config = cfg,
               check_args = check_args, cv_folds = folds), f)
  res
}

row_file <- function(id, cfg) file.path(rows_dir, sprintf("scenario_%s__%s.rds", id, cfg))

# Runs the configs that are not finished yet for one scenario and saves one row
# file per config, so array tasks never write to the same file.
run_scenario <- function(sim, id, g) {
  todo <- if (as.logical(opt$force)) configs$cfg else
    configs$cfg[!file.exists(row_file(id, configs$cfg))]
  if (!length(todo)) return(invisible())
  meta <- as.data.frame(sim$meta)
  meta$label_true <- as.character(meta$label_true)
  X <- t(sim$X)   # the simulator returns samples x features
  md <- data.frame(sample = meta$sample, label = as.character(meta$label_observed))
  for (k in todo) {
    cfg <- configs[cfg == k]
    res <- tryCatch(run_config(sim, id, cfg, meta, X, md),
                    error = function(e) data.table(status = "error", error = conditionMessage(e)))
    saveRDS(cbind(g, cfg, res), row_file(id, k))
    cat(sprintf("scenario %d [%s]: %s\n", id, k, res$status))
  }
  invisible()
}

# Combines every finished row file into the summary table.
combine_rows <- function() {
  files <- list.files(rows_dir, pattern = "^scenario_.*\\.rds$", full.names = TRUE)
  if (!length(files)) {
    cat("No finished rows in", rows_dir, "\n")
    return(invisible())
  }
  rows <- if (as.logical(opt$rescore_unknown)) {
    rescore_unknown_rows(files, dirty_dir)
  } else {
    lapply(files, readRDS)
  }
  new <- rbindlist(rows, fill = TRUE)
  setorder(new, scenario, cfg)
  summary_file <- file.path(opt$out, "summary_v2.rds")
  saveRDS(new, summary_file)
  fwrite(new, sub("\\.rds$", ".csv", summary_file))
  cat(sprintf("\nCombined %d rows into %s\n", nrow(new), summary_file))

  cat("\nStatus counts:\n"); print(new[, .N, by = status])
  ok <- new[status == "ok"]
  if (nrow(ok)) {
    cat("\nBy config and mislabel rate (ok runs):\n")
    print(ok[, .(runs = .N, sens = mean(sens_true, na.rm = TRUE), fp = mean(fp),
                 spec = mean(spec_true), auroc = mean(auroc_annotation_margin, na.rm = TRUE),
                 sens_gap = mean(sens_gap, na.rm = TRUE),
                 verdict_ok = mean(verdict_correct, na.rm = TRUE)),
             by = .(cfg, misid_rate)][order(cfg, misid_rate)])
  }
  invisible(new)
}

# ---- Run ----------------------------------------------------------------------

if (!as.logical(opt$combine)) {
  sims <- readRDS(opt$input)
  params <- rbindlist(lapply(sims, `[[`, "params"), fill = TRUE)
  params[, scenario := as.integer(names(sims))]
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
  cat(sprintf("task %d/%d: %d scenarios x %d configs (%d scenarios in the grid, %s cores)\n",
              task_id, n_tasks, length(mine), nrow(configs), nrow(params), opt$cores))

  run_one <- function(j) run_scenario(sims[[j]], params$scenario[j], params[j])
  # Keep only this task's scenarios in memory; the full file is large.
  sims <- sims[mine]; params <- params[mine]; invisible(gc())
  res <- parallel::mclapply(seq_along(mine), run_one, mc.cores = as.integer(opt$cores),
                            mc.preschedule = FALSE)
  failed <- vapply(res, inherits, logical(1), "try-error")
  if (any(failed)) cat(sum(failed), "worker(s) crashed; rerun to retry them\n")
  cat("task", task_id, "done\n")
}

# A single process combines its own results; in a job array, run
# `combine=TRUE` once after all tasks finish.
if (as.logical(opt$combine) || n_tasks == 1L) combine_rows()
