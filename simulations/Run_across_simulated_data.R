# Run check_annotations() on every scenario in results/simulated_datasets_v2.rds
# (made by Simulate_full_dataset.R) and score it against the known truth.
#
# Usage (repo root, package installed with R CMD INSTALL .):
#   Rscript simulations/Run_across_simulated_data.R
#   Rscript simulations/Run_across_simulated_data.R cores=4 subset="rep == 1"
#   Rscript simulations/Run_across_simulated_data.R scenarios=1:20 force=TRUE
#
# Arguments (key=value, all optional):
#   input      simulated datasets        (default simulations/results/simulated_datasets_v2.rds)
#   out        output directory          (default simulations/results)
#   cores      parallel workers          (default 1)
#   subset     data.table expression on the scenario parameters, e.g. "rep == 1"
#   scenarios  R expression for scenario numbers, e.g. 1:20
#   cv_folds, mismatch_repeats, target_specificity
#              comma-separated values to compare; every combination is run, e.g.
#              cv_folds=5,loo mismatch_repeats=1,3 target_specificity=0.9,0.95,0.99
#              ("loo" = leave-one-out). Defaults: 5, 1, 0.95
#   force      TRUE to recompute scenarios that already have a dirty file
#
# Outputs:
#   <out>/summary_v2.rds, .csv   one row per scenario (parameters + metrics)
#   <out>/dirty/scenario_<n>__<config>.rds the full check_annotations() result plus the
#                                truth for that scenario. These are large and
#                                disposable: delete the folder at any time.
# Finished scenarios are skipped on a rerun (use force=TRUE to redo them), so
# an interrupted run can be resumed.
suppressMessages({library(data.table); library(SampleAnnotationLikelihood)})

# ---- Settings -----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
opt <- list(input = "simulations/results/simulated_datasets_v2.rds",
            out = "simulations/results", cores = "1", subset = "", scenarios = "",
            cv_folds = "5", mismatch_repeats = "1", target_specificity = "0.95",
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
dirty_dir <- file.path(opt$out, "dirty")
dir.create(dirty_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Scoring helpers ----------------------------------------------------------

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
    # separation: how well the raw annotation margin ranks real mislabels
    auroc_annotation_margin = auroc(p$annotation_score_margin, bad),
    auroc_top_margin = auroc(-p$top_score_margin, bad),
    tp = tp, fp = fp, fn = fn, tn = tn,
    sens_true = sens, spec_true = spec, precision_true = ratio(tp, tp + fp),
    sens_reported = rep_m$sensitivity, spec_reported = rep_m$specificity,
    sens_gap = rep_m$sensitivity - sens, spec_gap = rep_m$specificity - spec,
    verdict = verdict,
    # a verdict is right if it matches whether the true performance is usable
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

# Runs the configs that are not finished yet for one scenario; returns rows.
run_scenario <- function(sim, id) {
  done <- file.exists(file.path(dirty_dir, sprintf("scenario_%s__%s.rds", id, configs$cfg)))
  todo <- if (as.logical(opt$force)) configs$cfg else configs$cfg[!done]
  meta <- as.data.frame(sim$meta)
  meta$label_true <- as.character(meta$label_true)
  X <- t(sim$X)   # the simulator returns samples x features
  md <- data.frame(sample = meta$sample, label = as.character(meta$label_observed))
  rbindlist(lapply(todo, function(k) {
    cfg <- configs[cfg == k]
    res <- tryCatch(run_config(sim, id, cfg, meta, X, md),
                    error = function(e) data.table(status = "error", error = conditionMessage(e)))
    cat(sprintf("scenario %d [%s]: %s\n", id, k, res$status))
    cbind(cfg, res)
  }), fill = TRUE)
}

# ---- Run ----------------------------------------------------------------------

sims <- readRDS(opt$input)
params <- rbindlist(lapply(sims, `[[`, "params"), fill = TRUE)
params[, scenario := as.integer(names(sims))]
keep <- rep(TRUE, nrow(params))
if (nzchar(opt$subset)) keep <- keep & params[, eval(parse(text = opt$subset))]
if (nzchar(opt$scenarios)) keep <- keep & params$scenario %in% eval(parse(text = opt$scenarios))
todo <- which(keep)
if (!as.logical(opt$force)) {
  # drop scenarios whose configs are all finished
  all_done <- vapply(params$scenario[todo], function(id) {
    all(file.exists(file.path(dirty_dir, sprintf("scenario_%s__%s.rds", id, configs$cfg))))
  }, logical(1))
  todo <- todo[!all_done]
}
cat(sprintf("%d scenarios x %d configs selected (%d scenarios in the grid, %s cores)\n",
            length(todo), nrow(configs), nrow(params), opt$cores))

run_one <- function(i) {
  g <- params[i]
  cbind(g, run_scenario(sims[[i]], g$scenario))
}
rows <- parallel::mclapply(todo, run_one, mc.cores = as.integer(opt$cores),
                           mc.preschedule = FALSE)
failed <- vapply(rows, inherits, logical(1), "try-error")
if (any(failed)) cat(sum(failed), "worker(s) crashed; rerun to retry them\n")
new <- rbindlist(rows[!failed], fill = TRUE)

# Merge with earlier runs so a resumed or partial run keeps its summary.
summary_file <- file.path(opt$out, "summary_v2.rds")
if (file.exists(summary_file)) {
  old <- readRDS(summary_file)
  new <- rbindlist(list(old[!paste(scenario, cfg) %in% paste(new$scenario, new$cfg)], new), fill = TRUE)
}
setorder(new, scenario, cfg)
saveRDS(new, summary_file)
fwrite(new, sub("\\.rds$", ".csv", summary_file))

# ---- Quick look ---------------------------------------------------------------

cat("\nStatus counts:\n"); print(new[, .N, by = status])
ok <- new[status == "ok"]
if (nrow(ok)) {
  cat("\nBy config and mislabel rate (ok runs):\n")
  print(ok[, .(runs = .N, sens = mean(sens_true, na.rm = TRUE), fp = mean(fp),
               spec = mean(spec_true), auroc = mean(auroc_annotation_margin, na.rm = TRUE),
               sens_gap = mean(sens_gap, na.rm = TRUE), verdict_ok = mean(verdict_correct, na.rm = TRUE)),
           by = .(cfg, misid_rate)][order(cfg, misid_rate)])
}
