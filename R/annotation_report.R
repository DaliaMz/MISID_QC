#' Write a printable annotation QC report
#'
#' Summarizes the first-pass holdout, detection and unknown thresholds, recovery
#' coverage, and final sample calls. Run [check_annotations()] with
#' `diagnostics = TRUE` for the calibration plots. Without diagnostics, the
#' report explicitly marks those plots as unavailable.
#'
#' Internal validation uses simulated errors and is not ground truth or a
#' sample-level probability. First-pass agreement is with observed labels.
#'
#' @param result A result from [check_annotations()].
#' @param file Path to the PDF to create.
#' @return The PDF path, invisibly.
#' @export
annotation_report <- function(result, file = "annotation_report.pdf") {
  if (!is.list(result) || !is.data.frame(result$predictions) ||
      is.null(result$calibration) || is.null(result$held_out) ||
      is.null(result$performance) || is.null(result$thresholds)) {
    stop("result must be an output from check_annotations().")
  }
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("file must be a non-empty PDF path.")
  }

  p <- result$predictions
  diagnostics <- result$diagnostics
  metric <- function(x) {
    if (length(x) != 1L || is.na(x) || !is.finite(x)) "not available" else
      sprintf("%.1f%%", 100 * x)
  }
  count <- function(x) sum(x, na.rm = TRUE)
  verdict <- function(x) if (is.null(x)) "not assessed" else x
  f <- result$feasibility
  status <- function(name) verdict(if (is.null(f)) NULL else f[[name]])
  empty_panel <- function(title, text) {
    graphics::plot.new()
    graphics::title(main = title)
    graphics::text(0.5, 0.5, text)
  }
  title_page <- function(title, subtitle = NULL) {
    graphics::plot.new()
    graphics::title(main = title, sub = subtitle)
  }
  # A signed logarithm preserves direction, zero, and the threshold comparison.
  transform_score <- function(x) sign(x) * log10(1 + abs(x))
  distribution <- function(groups, cutoff, title, direction, xlab) {
    groups <- lapply(groups, function(x) transform_score(x[is.finite(x)]))
    groups <- groups[lengths(groups) > 0L]
    if (!length(groups)) {
      empty_panel(title, "No finite calibration values")
      return(invisible())
    }
    values <- unlist(groups, use.names = FALSE)
    threshold <- if (length(cutoff) == 1L && is.finite(cutoff)) {
      transform_score(cutoff)
    } else numeric()
    limits <- range(c(values, threshold))
    if (diff(limits) == 0) limits <- limits + c(-0.5, 0.5)
    graphics::plot(NA, xlim = limits, ylim = c(0.5, length(groups) + 0.5),
                   xlab = xlab, ylab = "", yaxt = "n", main = title)
    graphics::axis(2, at = seq_along(groups),
                   labels = paste0(names(groups), " (n=", lengths(groups), ")"),
                   las = 1, cex.axis = 0.8)
    for (i in seq_along(groups)) {
      graphics::points(groups[[i]], rep(i, length(groups[[i]])),
                       pch = 16, cex = 0.55,
                       col = grDevices::adjustcolor(c("#286090", "#D2691E")[1L + (i - 1L) %% 2L], 0.45))
    }
    if (length(threshold)) graphics::abline(v = threshold, lty = 2, col = "firebrick", lwd = 2)
    graphics::mtext(paste("Dashed = cutoff;", direction), side = 3, line = 0.2, cex = 0.75)
  }

  grDevices::pdf(file, width = 11, height = 8.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(5, 5, 4, 2))
  title_page("Annotation QC overview", "Internal estimates are exploratory, not probabilities or ground truth")
  det <- result$performance$misannotation_detection
  rows <- c(
    sprintf("Samples: %d    First-pass agreement with observed labels: %s",
            nrow(p), metric(result$calibration$first_pass_agreement)),
    sprintf("First-pass held out: %d of %d",
            result$held_out$n_confidently_misassigned, result$held_out$n_samples),
    sprintf("Suspicious: %d    Candidate corrections: %d    Unknown flags: %d",
            count(p$suspicious), count(!is.na(p$candidate_correction)),
            count(p$potential_unknown_individual)),
    sprintf("Unavailable suspicious flags: %d    Unavailable unknown flags: %d",
            sum(is.na(p$suspicious)), sum(is.na(p$potential_unknown_individual))),
    sprintf("Detection feasibility: %s    Recovery feasibility: %s",
            status("misannotation_detection"), status("identity_recovery")),
    sprintf("Unknown feasibility: %s", status("unknown_individual_detection")),
    sprintf("Detection sensitivity: %s    Specificity: %s",
            metric(det$sensitivity), metric(det$specificity)),
    sprintf("Sensitivity source: %s; validated repeats: %s",
            if (is.null(det$source)) "not available" else det$source,
            if (is.null(det$n_repeats_validated)) "not available" else det$n_repeats_validated),
    sprintf("Recovery on simulated mislabels: %s",
            metric(result$performance$recovery_accuracy_misannotated)),
    "Counts exclude unavailable flags; unknown and suspicious flags can overlap."
  )
  graphics::text(0.04, seq(0.88, 0.16, length.out = length(rows)),
                 rows, adj = 0, cex = 0.98)

  first <- diagnostics$first_pass
  if (is.null(first)) {
    empty_panel("First-pass holdout", "Unavailable: run with diagnostics = TRUE")
  } else {
    graphics::par(mar = c(5, 5, 4, 2))
    values <- transform_score(first$annotation_margin)
    cutoff <- first$exclusion_margin_threshold[1L]
    lim <- range(c(values[is.finite(values)],
                   if (is.finite(cutoff)) transform_score(cutoff)), finite = TRUE)
    if (!all(is.finite(lim))) {
      empty_panel("First-pass holdout", "No finite first-pass margins")
    } else {
      if (diff(lim) == 0) lim <- lim + c(-0.5, 0.5)
      group <- ifelse(first$correct, 1, 2)
      graphics::plot(values, group, xlim = lim, ylim = c(0.5, 2.5), yaxt = "n",
                     xlab = "Signed log10(1 + annotation margin)",
                     ylab = "", main = "First-pass holdout",
                     pch = ifelse(first$held_out, 17, 16),
                     col = ifelse(first$held_out, "firebrick", "#286090"))
      graphics::axis(2, 1:2, c("Agrees with observed", "Disagrees"), las = 1, cex.axis = 0.8)
      if (is.finite(cutoff)) graphics::abline(v = transform_score(cutoff), lty = 2, col = "firebrick")
      graphics::mtext("Triangles: held out; dashed: exclusion cutoff (applied to disagreements)",
                      side = 3, line = 0.2, cex = 0.75)
    }
  }

  annotation <- diagnostics$annotation_threshold
  if (is.null(annotation)) {
    empty_panel("Detection calibration", "Unavailable: run with diagnostics = TRUE")
  } else {
    groups <- split(annotation$margin, annotation$population)
    graphics::par(mar = c(5, 11, 4, 2))
    distribution(groups, result$thresholds$annotation_margin,
                 "Detection calibration", "larger margins are suspicious",
                 "Signed log10(1 + annotation margin)")
    graphics::mtext(paste("Source:", diagnostics$annotation_threshold_source,
                          "| simulated errors are not ground truth"),
                    side = 1, line = 3, cex = 0.75)
  }

  graphics::par(mar = c(5, 5, 4, 2))
  if (is.data.frame(result$risk_coverage) && nrow(result$risk_coverage)) {
    plot_risk_coverage(result)
    graphics::mtext(sprintf("Selected top-margin cutoff: %s; prediction coverage: %s",
                            format(result$thresholds$abstain_margin, digits = 3),
                            metric(result$calibration$prediction_coverage)),
                    side = 3, line = 0.3, cex = 0.75)
  } else {
    empty_panel("Recovery and coverage", "Risk-coverage data unavailable")
  }

  unknown <- diagnostics$unknown
  graphics::par(mar = c(5, 11, 4, 2))
  if (is.null(unknown)) {
    empty_panel("Unknown detection", "Unavailable: run with diagnostics = TRUE")
  } else if (!isTRUE(unknown$enabled)) {
    empty_panel("Unknown detection", "Unknown detection was disabled")
  } else {
    distance <- identical(unknown$statistic, "nearest_class_distance")
    distribution(list("Reference" = unknown$known_values,
                      "Simulated unknowns" = unknown$simulated_unknown_values),
                 unknown$threshold, "Unknown detection calibration",
                 if (distance) "larger distances are flagged" else "smaller margins are flagged",
                 if (distance) "Signed log10(1 + nearest-class distance)" else
                   "Signed log10(1 + top margin)")
    graphics::mtext("Simulated unknowns; scores are not probabilities",
                    side = 1, line = 3, cex = 0.75)
  }

  graphics::par(mar = c(5, 12, 4, 2))
  counts <- table(p$annotation_status, useNA = "ifany")
  graphics::barplot(counts, horiz = TRUE, las = 1, col = "#286090",
                    xlab = "Number of samples", main = "Final annotation statuses",
                    cex.names = 0.75, xlim = c(0, max(counts, 1) * 1.15))
  graphics::mtext("Unknown flags are separate and may overlap with these statuses",
                  side = 3, line = 0.2, cex = 0.75)

  show <- (!is.na(p$suspicious) & p$suspicious) |
    (!is.na(p$potential_unknown_individual) & p$potential_unknown_individual)
  calls <- p[show, , drop = FALSE]
  if (!nrow(calls)) {
    empty_panel("Flagged samples", "No suspicious or unknown-flagged samples")
  } else {
    pages <- split(seq_len(nrow(calls)), ceiling(seq_len(nrow(calls)) / 22L))
    for (i in seq_along(pages)) {
      title_page(sprintf("Flagged samples (%d/%d)", i, length(pages)),
                 "Unknown flags and suspicious flags are separate")
      d <- calls[pages[[i]], , drop = FALSE]
      columns <- list(d$sample_id, d$observed_id,
                      ifelse(is.na(d$candidate_correction), "-", d$candidate_correction),
                      ifelse(is.na(d$suspicious), "NA", ifelse(d$suspicious, "yes", "no")),
                      ifelse(is.na(d$potential_unknown_individual), "NA",
                             ifelse(d$potential_unknown_individual, "yes", "no")))
      x <- c(0.02, 0.28, 0.49, 0.72, 0.87)
      headings <- c("Sample", "Observed", "Correction", "Suspicious", "Unknown")
      graphics::text(x, 0.93, headings, adj = 0, font = 2, cex = 0.8)
      for (j in seq_along(columns)) {
        graphics::text(x[j], seq(0.88, by = -0.035, length.out = nrow(d)),
                       as.character(columns[[j]]), adj = 0, cex = 0.72)
      }
    }
  }
  invisible(file)
}
