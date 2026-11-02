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

# ---- Diagnostic plots --------------------------------------------------------

plot_data <- copy(ok[is.finite(auroc_annotation_margin) &
                     is.finite(sens_true) & is.finite(spec_true)])
if (nrow(plot_data)) {
  cfg_levels <- sort(unique(plot_data$cfg))
  cfg_colors <- grDevices::hcl.colors(length(cfg_levels), "Dark 3")
  misid_levels <- sort(unique(plot_data$misid_rate))
  misid_colors <- grDevices::hcl.colors(length(misid_levels), "Set 2")
  plot_data[, cfg_color := cfg_colors[match(cfg, cfg_levels)]]
  plot_data[, misid_color := misid_colors[match(misid_rate, misid_levels)]]

  high_auc_low_sens <- plot_data[auroc_annotation_margin >= 0.9 & sens_true < 0.8]
  draw_diagnostics <- function() {
    old_par <- graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1))
    on.exit(graphics::par(old_par))

    graphics::plot(plot_data$auroc_annotation_margin, plot_data$sens_true,
                   xlim = c(0, 1), ylim = c(0, 1), pch = 16,
                   col = plot_data$misid_color,
                   xlab = "AUROC (continuous annotation-margin ranking)",
                   ylab = "Sensitivity at the selected cutoff",
                   main = "Ranking quality vs detected mislabels",
                   sub = "Each point is one scenario; colors show mislabel rate")
    graphics::abline(v = 0.9, h = 0.8, col = "grey45", lty = 2)
    graphics::legend("bottomright", legend = paste("Mislabel rate", misid_levels),
                     col = misid_colors, pch = 16, bty = "n", cex = 0.8)
    graphics::mtext(sprintf("%d scenarios have AUROC >= 0.90 and sensitivity < 0.80",
                            nrow(high_auc_low_sens)), side = 3, line = 0.3, cex = 0.8)

    graphics::plot(plot_data$sens_true, plot_data$spec_true,
                   xlim = c(0, 1), ylim = c(0, 1), pch = 16,
                   col = plot_data$cfg_color,
                   xlab = "Sensitivity (true mislabels detected)",
                   ylab = "Observed specificity (correct labels not flagged)",
                   main = "Operating points by configuration",
                   sub = "Each point is one scenario; task and tissue set are pooled")
    target_specificities <- suppressWarnings(as.numeric(sub(".*_ts", "", cfg_levels)))
    target_specificities <- unique(target_specificities[is.finite(target_specificities)])
    if (length(target_specificities)) {
      graphics::abline(h = target_specificities, col = "grey85", lty = 3)
    }
    graphics::abline(v = 0.8, h = 0.95, col = "grey45", lty = 2)
    graphics::legend("bottomright", legend = cfg_levels, col = cfg_colors,
                     pch = 16, bty = "n", cex = 0.7, title = "Config")
  }

  if (interactive()) draw_diagnostics()
  plot_file <- file.path(out, "summary_gtex_diagnostics.pdf")
  grDevices::pdf(plot_file, width = 11, height = 7.5)
  draw_diagnostics()
  grDevices::dev.off()
  cat("\nDiagnostic plots written to", plot_file, "\n")
} else {
  cat("\nNo successful rows with finite AUROC, sensitivity, and specificity; skipping diagnostic plots.\n")
}

# ---- Stratified performance report -------------------------------------------

report_data <- copy(all)
report_data[, setting := sprintf("%s labels / %s samples per label",
                                 n_labels, n_per_label)]
group_cols <- c("cfg", "task", "set", "n_labels", "n_per_label",
                "misid_rate", "setting")

status_rates <- list(
  list(value = "ok", title = "Actionable result fraction",
       ylab = "Fraction with status = ok"),
  list(value = "low_agreement", title = "Low-agreement gate rate",
       ylab = "Fraction with low agreement"),
  list(value = "error", title = "Execution error rate",
       ylab = "Fraction with status = error")
)
status_summaries <- lapply(status_rates, function(x) {
  report_data[, .(mean = mean(status == x$value)), by = group_cols]
})

summarize_metric <- function(data, metric) {
  data[status == "ok", {
    x <- get(metric)
    valid <- is.finite(x)
    .(mean = if (any(valid)) mean(x[valid]) else NA_real_,
      lower = if (any(valid)) min(x[valid]) else NA_real_,
      upper = if (any(valid)) max(x[valid]) else NA_real_)
  }, by = group_cols]
}

# One page per config / task / tissue set, with one panel per metric.
# Settings are nudged sideways and use different symbols and line types so
# lines with similar values stay distinguishable.
draw_report_page <- function(cfg_value, task_value, set_value, panels) {
  base_data <- report_data[cfg == cfg_value & task == task_value & set == set_value]
  settings <- sort(unique(base_data$setting))
  n_set <- length(settings)
  setting_colors <- grDevices::hcl.colors(max(n_set, 3L), "Dark 3")[seq_len(n_set)]
  setting_pch <- rep_len(c(16, 17, 15, 18, 1, 2, 0, 5), n_set)
  setting_lty <- rep_len(1:4, n_set)
  x_range <- range(base_data$misid_rate, finite = TRUE)
  if (diff(x_range) == 0) x_range <- x_range + c(-0.02, 0.02)
  step <- diff(x_range) * 0.015
  offsets <- (seq_len(n_set) - (n_set + 1) / 2) * step
  x_lim <- x_range + c(-1, 1) * (max(abs(offsets), 0) + step)

  graphics::par(mfrow = c(2, 4), mar = c(4, 4, 3, 1), oma = c(2, 0, 3, 0))
  for (p in panels) {
    data <- p$data[cfg == cfg_value & task == task_value & set == set_value]
    target <- p$target
    if (p$column == "spec_true") {
      target <- unique(report_data[cfg == cfg_value]$target_specificity)
      if (length(target) != 1L || !is.finite(target)) target <- NA_real_
    }
    graphics::plot(NA, xlim = x_lim, ylim = c(0, 1), xlab = "Mislabel rate",
                   ylab = p$ylab, main = p$title, cex.main = 0.95)
    if (is.finite(target)) graphics::abline(h = target, col = "grey45", lty = 2)
    for (j in seq_len(n_set)) {
      curve <- data[setting == settings[j]][order(misid_rate)]
      curve <- curve[is.finite(mean)]
      if (!nrow(curve)) next
      x <- curve$misid_rate + offsets[j]
      if (p$show_range) {
        has_range <- is.finite(curve$lower) & is.finite(curve$upper) &
          curve$lower < curve$upper
        if (any(has_range)) {
          graphics::arrows(x[has_range], curve$lower[has_range],
                           x[has_range], curve$upper[has_range],
                           angle = 90, code = 3, length = 0.02,
                           col = grDevices::adjustcolor(setting_colors[j], 0.6))
        }
      }
      graphics::lines(x, curve$mean, type = "b", pch = setting_pch[j],
                      lty = setting_lty[j], col = setting_colors[j], cex = 0.9)
    }
    if (!any(is.finite(data$mean))) {
      graphics::text(mean(x_range), 0.5, "No finite values", cex = 0.8)
    }
  }
  graphics::plot.new()
  graphics::legend("center", legend = settings, col = setting_colors,
                   pch = setting_pch, lty = setting_lty, bty = "n", cex = 0.85,
                   title = "Labels / samples per label")
  graphics::mtext(sprintf("Task: %s  |  Tissue set: %s  |  Config: %s",
                          task_value, set_value, cfg_value),
                  outer = TRUE, cex = 1.1)
  graphics::mtext(paste("Points are means over replicates (positions nudged",
                        "sideways to avoid overlap); whiskers are min-max.",
                        "Metric panels use status = ok rows only."),
                  side = 1, outer = TRUE, line = 0.5, cex = 0.75)
}

stratified_file <- file.path(out, "summary_gtex_performance_by_setting.pdf")
draw_stratified_report <- function() {
  panels <- c(
    Map(function(x, s) list(data = s, column = x$value, title = x$title,
                            ylab = x$ylab, show_range = FALSE, target = NA_real_),
        status_rates, status_summaries),
    lapply(list(
      list(column = "sens_true", title = "Mislabel detection", ylab = "Sensitivity",
           target = 0.8),
      list(column = "spec_true", title = "Correct-label protection", ylab = "Specificity",
           target = NA_real_),
      list(column = "auroc_annotation_margin", title = "Ranking quality",
           ylab = "AUROC", target = NA_real_),
      list(column = "correction_accuracy", title = "Candidate correction quality",
           ylab = "Correction accuracy", target = NA_real_)
    ), function(m) c(list(data = summarize_metric(report_data, m$column),
                          show_range = TRUE), m))
  )
  pages <- unique(report_data[, .(cfg, task, set)])[order(cfg, task, set)]
  for (i in seq_len(nrow(pages))) {
    draw_report_page(pages$cfg[i], pages$task[i], pages$set[i], panels)
  }
}

if (interactive()) {
  # Pause between pages so each can be inspected; click or press Enter to advance.
  old_ask <- grDevices::devAskNewPage(TRUE)
  draw_stratified_report()
  grDevices::devAskNewPage(old_ask)
}
grDevices::pdf(stratified_file, width = 11, height = 8.5)
draw_stratified_report()
grDevices::dev.off()
cat("\nStratified performance report written to", stratified_file, "\n")
