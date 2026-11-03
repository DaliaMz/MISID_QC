# Source this file, then call:
#   plot_first_pass_truth(result, truth, file = "first-pass-truth.pdf")
# `result` must come from check_annotations(..., diagnostics = TRUE); `truth`
# needs `sample` and logical `is_misid` (simulation truth). Returns the
# merged first-pass table invisibly.

plot_first_pass_truth <- function(result, truth, file = NULL,
                                  sample_col = "sample",
                                  misid_col = "is_misid") {
  first <- result$diagnostics$first_pass
  if (is.null(first)) {
    stop("Run check_annotations(..., diagnostics = TRUE) first.")
  }
  m <- match(first$sample_id, truth[[sample_col]])
  if (anyNA(m)) stop("Some first-pass samples are missing from `truth`.")
  first$truly_misannotated <- as.logical(truth[[misid_col]][m])
  if ("label_true" %in% names(truth)) {
    first$true_id <- as.character(truth$label_true[m])
    first$predicted_is_true <- first$predicted_id == first$true_id
  }
  thr <- first$exclusion_margin_threshold[1L]

  if (!is.null(file)) {
    grDevices::pdf(file, width = 11, height = 8.5, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  col_clean <- grDevices::adjustcolor("steelblue", 0.6)
  col_bad <- grDevices::adjustcolor("firebrick", 0.7)
  pt_col <- ifelse(first$truly_misannotated, col_bad, col_clean)
  logx <- function(v) log10(pmax(v, .Machine$double.eps))
  thr_line <- function(axis_log = TRUE) {
    if (is.finite(thr) && thr > 0) {
      v <- if (axis_log) log10(thr) else thr
      graphics::abline(v = v, col = "black", lwd = 2, lty = 2)
    }
  }

  # ---- Page 0: margin distributions, predicted correctly vs not --------------
  graphics::par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
  col_a <- grDevices::adjustcolor("grey30", 0.5)
  col_b <- grDevices::adjustcolor("darkorange", 0.6)
  dens_overlay <- function(a, b, xlab, main, lab_a, lab_b) {
    da <- if (length(a) > 1L) stats::density(a)
    db <- if (length(b) > 1L) stats::density(b)
    if (is.null(da) && is.null(db)) {
      graphics::plot.new(); graphics::text(0.5, 0.5, "Too few samples"); return()
    }
    graphics::plot(NA, xlim = range(c(a, b)), ylim = range(c(da$y, db$y)),
                   xlab = xlab, ylab = "Density", main = main)
    if (!is.null(da)) graphics::polygon(da, col = col_a, border = NA)
    if (!is.null(db)) graphics::polygon(db, col = col_b, border = NA)
    graphics::legend("topright", bty = "n", fill = c(col_a, col_b),
                     legend = c(sprintf("%s (n = %d)", lab_a, length(a)),
                                sprintf("%s (n = %d)", lab_b, length(b))))
  }
  ok <- first$correct
  tm_ok <- logx(first$top_margin[ok]); tm_bad <- logx(first$top_margin[!ok])
  dens_overlay(tm_ok, tm_bad, "log10(top - second score)", "Top margin",
               "predicted = observed", "predicted != observed")
  thr_line()

  graphics::plot(NA, xlim = range(logx(first$top_margin)), ylim = c(0, 1),
                 xlab = "log10(top - second score)",
                 ylab = "Cumulative fraction", main = "Top margin ECDF")
  if (length(tm_ok) > 1L) graphics::lines(stats::ecdf(tm_ok), col = "grey30",
                                          lwd = 2, do.points = FALSE)
  if (length(tm_bad) > 1L) graphics::lines(stats::ecdf(tm_bad), col = "darkorange",
                                           lwd = 2, do.points = FALSE)
  thr_line()
  graphics::legend("bottomright", bty = "n", lwd = 2,
                   col = c("grey30", "darkorange"),
                   legend = c("predicted = observed", "predicted != observed"))

  dis0 <- !ok
  wrong <- first$truly_misannotated
  dens_overlay(logx(first$annotation_margin[dis0 & !wrong]),
               logx(first$annotation_margin[dis0 & wrong]),
               "log10(top - observed-label score)",
               "Annotation margin (disagreeing only)",
               "label actually correct", "label actually wrong")
  thr_line()
  dens_overlay(logx(first$top_margin[dis0 & !wrong]),
               logx(first$top_margin[dis0 & wrong]),
               "log10(top - second score)",
               "Top margin (disagreeing only)",
               "label actually correct", "label actually wrong")

  # ---- Page 1: how the cutoff is derived and what it holds out ---------------
  graphics::par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))

  # A: top margins of agreeing samples define the cutoff.
  agree <- first$correct
  tm <- logx(first$top_margin[agree])
  graphics::hist(
    tm, breaks = 40, col = "grey80", border = "white",
    main = sprintf("A. Top margin, predicted = observed (n = %d)", sum(agree)),
    xlab = "log10(top - second score)"
  )
  thr_line()
  graphics::legend("topright", bty = "n", lty = 2, lwd = 2,
                   legend = sprintf("cutoff = %.3g (low quantile)", thr))

  # B: annotation margins of disagreeing samples by true status.
  dis <- !first$correct
  if (any(dis)) {
    am <- logx(first$annotation_margin[dis])
    brk <- seq(min(am), max(am) + 1e-9, length.out = 41)
    h_bad <- graphics::hist(am[first$truly_misannotated[dis]], breaks = brk,
                            plot = FALSE)
    h_ok <- graphics::hist(am[!first$truly_misannotated[dis]], breaks = brk,
                           plot = FALSE)
    graphics::barplot(
      rbind(h_ok$counts, h_bad$counts), beside = FALSE, space = 0,
      col = c(col_clean, col_bad), border = "white",
      names.arg = rep("", length(h_ok$counts)),
      main = sprintf("B. Disagreeing samples (n = %d)", sum(dis)),
      xlab = "log10(top - observed-label score)", ylab = "Samples"
    )
    ticks <- pretty(range(brk), 6)
    at <- (ticks - min(brk)) / diff(brk)[1L]
    graphics::axis(1, at = at, labels = ticks)
    if (is.finite(thr) && thr > 0) {
      graphics::abline(v = (log10(thr) - min(brk)) / diff(brk)[1L],
                       lwd = 2, lty = 2)
    }
    graphics::legend("topright", bty = "n", fill = c(col_clean, col_bad),
                     legend = c("label actually correct", "label actually wrong"))
  } else {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No disagreeing samples")
  }

  # C: top vs annotation margin for disagreeing samples.
  if (any(dis)) {
    d <- first[dis, , drop = FALSE]
    graphics::plot(
      logx(d$top_margin), logx(d$annotation_margin),
      col = ifelse(d$truly_misannotated, col_bad, col_clean),
      pch = ifelse(d$held_out, 17, 1), cex = 1.2,
      xlab = "log10(top - second score)",
      ylab = "log10(top - observed-label score)",
      main = "C. Disagreeing samples: triangles = held out"
    )
    graphics::abline(h = log10(thr), lwd = 2, lty = 2)
    graphics::legend("topleft", bty = "n", pch = c(16, 16, 17, 1),
                     col = c(col_clean, col_bad, "black", "black"),
                     legend = c("label actually correct", "label actually wrong",
                                "held out", "kept"))
  } else {
    graphics::plot.new()
  }

  # D: outcome table.
  held <- first$held_out
  bad <- first$truly_misannotated
  tp <- sum(held & bad); fp <- sum(held & !bad)
  fn <- sum(!held & bad); tn <- sum(!held & !bad)
  rt <- function(a, b) if (b > 0) a / b else NA_real_
  graphics::plot.new()
  graphics::title("D. Holdout vs truth")
  lines <- c(
    sprintf("Samples: %d   truly misannotated: %d", nrow(first), sum(bad)),
    sprintf("Exclusion cutoff (top margin quantile): %.4g", thr),
    "",
    sprintf("Held out & wrong label (TP):      %d", tp),
    sprintf("Held out & correct label (FP):    %d", fp),
    sprintf("Kept & wrong label (FN):          %d", fn),
    sprintf("Kept & correct label (TN):        %d", tn),
    "",
    sprintf("Holdout precision: %.3f", rt(tp, tp + fp)),
    sprintf("Holdout recall:    %.3f", rt(tp, tp + fn)),
    sprintf("Wrong labels not disagreeing at all: %d",
            sum(bad & first$correct)),
    sprintf("Wrong labels disagreeing but below cutoff: %d",
            sum(bad & !first$correct & !held))
  )
  graphics::text(0, seq(0.95, by = -0.075, length.out = length(lines)),
                 lines, adj = 0, family = "mono", cex = 0.9)

  # ---- Page 2: every sample, both margins, coloured by truth ----------------
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  ord <- order(first$truly_misannotated)  # draw wrong labels last
  graphics::plot(
    logx(first$top_margin[ord]), logx(first$annotation_margin[ord] + 1e-3),
    col = pt_col[ord], pch = ifelse(first$held_out[ord], 17, 16), cex = 1.1,
    xlab = "log10(top - second score)",
    ylab = "log10(annotation margin + 0.001); 0 margin = agreeing",
    main = "All samples (first pass)"
  )
  graphics::abline(h = log10(thr + 1e-3), lwd = 2, lty = 2)
  graphics::legend("topleft", bty = "n", pch = c(16, 16, 17),
                   col = c(col_clean, col_bad, "black"),
                   legend = c("label actually correct", "label actually wrong",
                              "held out (triangle)"))

  groups <- list(
    "correct label" = first$top_margin[!bad],
    "wrong label" = first$top_margin[bad]
  )
  groups <- groups[lengths(groups) > 0L]
  graphics::boxplot(lapply(groups, logx), col = c(col_clean, col_bad)[
    seq_along(groups)], ylab = "log10(top - second score)",
    main = "Top margin by true label status")

  # ---- Page 3: the samples behind the dots ------------------------------------
  keep <- which(!first$correct | first$held_out | bad)
  if (length(keep)) {
    tab <- first[keep, intersect(
      c("sample_id", "observed_id", "true_id", "predicted_id",
        "truly_misannotated", "held_out", "top_margin", "annotation_margin"),
      names(first)), drop = FALSE]
    tab <- tab[order(-tab$annotation_margin), , drop = FALSE]
    graphics::par(mfrow = c(1, 1))
    graphics::plot.new()
    graphics::title("Disagreeing, held-out or wrong-label samples (top 35)")
    top <- utils::head(tab, 35L)
    num <- vapply(top, is.numeric, logical(1))
    top[num] <- lapply(top[num], function(v) signif(v, 4))
    txt <- utils::capture.output(print(top, row.names = FALSE))
    graphics::text(0, seq(0.98, by = -0.026, length.out = length(txt)), txt,
                   adj = 0, family = "mono", cex = 0.45)
  }

  invisible(first)
}
