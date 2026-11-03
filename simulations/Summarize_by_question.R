custom_lib <- "/scratch/project/palpant_scratch/DMizikovsky/R_packages"
.libPaths(c(custom_lib, .libPaths()))

# One PDF organised by the questions the simulations are meant to answer
# (see simulations/README.md), plus the key numbers printed to the console.
#
#   Rscript simulations/Summarize_by_question.R
#   Rscript simulations/Summarize_by_question.R \
#     input=simulations/results/summary_v2.rds,simulations/results/summary_gtex.rds \
#     out=simulations/results cfg=cv5_mr3_mc0.025_ts0.95
#
# Arguments (key=value, all optional):
#   input  comma-separated summary files; each gets its own section
#          (default: summary_v2.rds and summary_gtex.rds under out, if present).
#          Make them with Run_across_simulated_data.R combine=TRUE and
#          Aggregate_gtex_results.R.
#   out    output directory (default simulations/results)
#   cfg    config used for Q1-Q4, Q6 (default: the config with the most rows)
#
# Q1  How well does it work (true sensitivity, specificity, AUROC, corrections)?
# Q2  How far are its own estimates from the truth (reported minus true)?
# Q3  Is the verdict right, and at what first-pass agreement does it break down?
# Q4  Does the first-pass holdout catch the mislabels?
# Q5  Do the tuning parameters matter (needs more than one config)?
# Q7  False positives on clean data are in the rate-0 points of Q1.
suppressMessages(library(data.table))

opt <- list(out = "simulations/results", input = "", cfg = "")
for (a in commandArgs(trailingOnly = TRUE)) {
  if (grepl("=", a, fixed = TRUE)) {
    kv <- strsplit(a, "=", fixed = TRUE)[[1]]
    opt[[kv[1]]] <- paste(kv[-1], collapse = "=")
  }
}
inputs <- if (nzchar(opt$input)) strsplit(opt$input, ",", fixed = TRUE)[[1]] else
  file.path(opt$out, c("summary_v2.rds", "summary_gtex.rds"))
inputs <- inputs[file.exists(inputs)]
if (!length(inputs)) stop("No summary files found; run the combine/aggregate steps first.")

agreement_breaks <- c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0001)
agreement_labels <- c("<0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", ">=0.9")

# ---- Plot helpers ---------------------------------------------------------------

as_unit <- function(v) pmin(1, pmax(0, v))
gap_limits <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return(c(-1, 1))
  m <- max(abs(v), 0.1)
  c(-m, m) * 1.05
}

# Per-scenario points (jittered) and group means across the mislabel rate.
rate_panel <- function(d, y, ylab, main, ylim = c(0, 1), hline = NULL) {
  d <- d[is.finite(get(y))]
  if (!nrow(d)) {
    graphics::plot.new(); graphics::title(main); graphics::text(0.5, 0.5, "No data")
    return(invisible())
  }
  rates <- sort(unique(d$misid_rate))
  groups <- sort(unique(d$group))
  cols <- grDevices::hcl.colors(length(groups), "Dark 3")
  d[, x := match(misid_rate, rates)]
  graphics::plot(range(d$x) + c(-0.4, 0.4), ylim, type = "n", xaxt = "n",
                 xlab = "Simulated mislabel rate", ylab = ylab, main = main)
  graphics::axis(1, at = seq_along(rates), labels = rates)
  if (!is.null(hline)) graphics::abline(h = hline, col = "grey55", lty = 2)
  for (k in seq_along(groups)) {
    g <- d[d$group == groups[k], ]
    graphics::points(jitter(g$x, 0.5), g[[y]], pch = 16, cex = 0.6,
                     col = grDevices::adjustcolor(cols[k], 0.4))
    m <- g[, .(v = mean(get(y))), by = x][order(x)]
    graphics::lines(m$x, m$v, col = cols[k], lwd = 2)
  }
  graphics::legend("bottomleft", legend = groups, col = cols, lwd = 2, bty = "n",
                   cex = 0.6)
}

# Gap (reported minus true) against first-pass agreement, with binned means.
gap_panel <- function(d, y, main) {
  d <- d[is.finite(get(y)) & is.finite(first_pass_agreement)]
  graphics::plot(d$first_pass_agreement, d[[y]], pch = 16, cex = 0.7,
                 col = grDevices::adjustcolor("steelblue", 0.5), ylim = gap_limits(d[[y]]),
                 xlim = c(min(0.3, d$first_pass_agreement), 1),
                 xlab = "First-pass agreement", ylab = "Reported minus true", main = main)
  graphics::abline(h = 0, col = "grey40")
  graphics::abline(v = 0.7, col = "grey70", lty = 3)
  if (nrow(d)) {
    d[, bin := cut(first_pass_agreement, agreement_breaks, labels = agreement_labels)]
    m <- d[, .(x = mean(first_pass_agreement), v = mean(get(y))), by = bin][order(x)]
    graphics::lines(m$x, m$v, col = "firebrick", lwd = 2, type = "b", pch = 15)
  }
}

# ---- Sections -------------------------------------------------------------------

summarise_source <- function(path) {
  all <- readRDS(path)
  label <- basename(path)
  all <- as.data.table(all)
  cat("\n==========", label, "==========\n")
  cat(nrow(all), "rows;", sum(all$status == "error", na.rm = TRUE), "errors\n")
  ok <- all[status %in% c("ok", "low_agreement")]
  if (!nrow(ok)) { cat("No usable rows.\n"); return(invisible()) }
  ok[, group := if ("var_label" %in% names(ok)) sprintf("label var %s, n/label %s", var_label, n_per_label)
                else sprintf("%s, %s labels", task, n_labels)]
  cfgs <- ok[, .N, by = cfg][order(-N)]
  cfg_main <- if (nzchar(opt$cfg)) opt$cfg else cfgs$cfg[1]
  cat("Config for Q1-Q4:", cfg_main, "\n")
  d <- ok[cfg == cfg_main]

  # Q1
  graphics::par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1), oma = c(0, 0, 2, 0))
  rate_panel(d, "sens_true", "True sensitivity", "Q1: mislabels detected", hline = 0.8)
  rate_panel(d, "spec_true", "True specificity", "Q1/Q7: correct labels not flagged (rate 0 = false positives)", hline = 0.95)
  rate_panel(d, "auroc_annotation_margin", "AUROC", "Q1: annotation-margin ranking")
  rate_panel(d, "correction_accuracy", "Correction accuracy", "Q1: proposed corrections that are right")
  graphics::mtext(paste(label, "- Q1: how well does it work?"), outer = TRUE, font = 2)

  # Q2
  graphics::par(mfrow = c(2, 2))
  gap_panel(d, "sens_gap", "Q2: sensitivity, reported - true")
  gap_panel(d, "spec_gap", "Q2: specificity, reported - true")
  gap_panel(d, "recovery_gap", "Q2: recovery on swapped samples, reported - true")
  gap_panel(d, "correction_precision_gap", "Q2: correction precision, reported - true")
  graphics::mtext(paste(label, "- Q2: how far are its estimates from the truth? (red = binned mean)"),
                  outer = TRUE, font = 2)

  # Q3
  d[, bin := cut(first_pass_agreement, agreement_breaks, labels = agreement_labels)]
  graphics::par(mfrow = c(1, 2))
  for (v in c("verdict_correct", "recovery_verdict_correct")) {
    if (!v %in% names(d)) next
    m <- d[!is.na(get(v)), .(rate = mean(get(v)), n = .N), by = bin][order(bin)]
    mids <- graphics::barplot(m$rate, names.arg = as.character(m$bin), ylim = c(0, 1),
                              xlab = "First-pass agreement", ylab = "Verdict matches the truth",
                              main = paste("Q3:", sub("_", " ", v)), col = "steelblue")
    graphics::text(mids, m$rate, paste0("n=", m$n), pos = 3, cex = 0.7)
  }
  graphics::mtext(paste(label, "- Q3: where does the verdict break down?"), outer = TRUE, font = 2)

  # Q4
  graphics::par(mfrow = c(1, 2))
  rate_panel(d, "held_precision", "Held out that are truly mislabelled", "Q4: holdout precision")
  rate_panel(d, "held_recall", "Mislabelled that are held out", "Q4: holdout recall")
  graphics::mtext(paste(label, "- Q4: does the first-pass holdout work?"), outer = TRUE, font = 2)

  # Q5
  if (nrow(cfgs) > 1L) {
    graphics::par(mfrow = c(2, 2))
    for (v in c("sens_true", "sens_gap", "spec_gap", "correction_precision_gap")) {
      vals <- ok[is.finite(get(v))]
      if (!nrow(vals)) { graphics::plot.new(); next }
      graphics::boxplot(vals[[v]] ~ vals$cfg, las = 2, cex.axis = 0.5, xlab = "",
                        ylab = v, main = paste("Q5:", v), col = "grey85")
      if (grepl("gap", v)) graphics::abline(h = 0, col = "grey40")
    }
    graphics::mtext(paste(label, "- Q5: do the tuning parameters matter?"), outer = TRUE, font = 2)
  }

  # Console summary of the key numbers.
  cat("\nQ2/Q3 by first-pass agreement (config", cfg_main, "):\n")
  print(d[, .(n = .N, sens_gap = round(mean(sens_gap, na.rm = TRUE), 3),
              spec_gap = round(mean(spec_gap, na.rm = TRUE), 3),
              recovery_gap = round(mean(recovery_gap, na.rm = TRUE), 3),
              verdict_correct = round(mean(verdict_correct, na.rm = TRUE), 3)),
          by = bin][order(bin)])
  cat("\nQ1/Q7 by mislabel rate:\n")
  print(d[, .(n = .N, sens = round(mean(sens_true, na.rm = TRUE), 3),
              spec = round(mean(spec_true, na.rm = TRUE), 3),
              auroc = round(mean(auroc_annotation_margin, na.rm = TRUE), 3),
              held_precision = round(mean(held_precision, na.rm = TRUE), 3),
              held_recall = round(mean(held_recall, na.rm = TRUE), 3)),
          by = misid_rate][order(misid_rate)])
  invisible()
}

pdf_file <- file.path(opt$out, "summary_by_question.pdf")
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
grDevices::pdf(pdf_file, width = 11, height = 8)
for (p in inputs) summarise_source(p)
grDevices::dev.off()
cat("\nPlots written to", pdf_file, "\n")
