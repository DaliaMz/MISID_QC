custom_lib <- "/scratch/project/palpant_scratch/DMizikovsky/R_packages"
.libPaths(c(custom_lib, .libPaths()))

# Combine the per-scenario rows written by Run_across_GTEX_simulations.R.
#   Rscript simulations/Aggregate_gtex_results.R [out=simulations/results]
suppressMessages(library(data.table))
out <- "simulations/results"
for (a in commandArgs(trailingOnly = TRUE)) {
  if (startsWith(a, "out=")) out <- sub("^out=", "", a)
}
files <- list.files(file.path(out, "rows_gtex"), pattern = "\\.rds$", full.names = TRUE)
if (!length(files)) stop("No rows found in ", file.path(out, "rows_gtex"))
all <- rbindlist(lapply(files, readRDS), fill = TRUE)
setorder(all, scenario, cfg)
saveRDS(all, file.path(out, "summary_gtex.rds"))
fwrite(all, file.path(out, "summary_gtex.csv"))
cat(nrow(all), "scenarios aggregated\n")

cat("\nStatus counts:\n"); print(all[, .N, by = status])
ok <- all[status == "ok"]
if (nrow(ok)) {
  cat("\nBy config, task, tissue set and mislabel rate (ok runs):\n")
  print(ok[, .(runs = .N, sens = mean(sens_true, na.rm = TRUE), fp = mean(fp),
               spec = mean(spec_true), auroc = mean(auroc_annotation_margin, na.rm = TRUE),
               sens_gap = mean(sens_gap, na.rm = TRUE),
               verdict_ok = mean(verdict_correct, na.rm = TRUE)),
           by = .(cfg, task, set, misid_rate)][order(cfg, task, set, misid_rate)])
}
