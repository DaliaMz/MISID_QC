# Source this file, then call:
#   plot_check_annotation_gates(result, file = "annotation-gates.pdf")
# `result` must come from check_annotations(..., diagnostics = TRUE).

plot_check_annotation_gates <- function(result, file = NULL) {
  diagnostics <- result$diagnostics
  if (is.null(diagnostics)) {
    stop("Run check_annotations(..., diagnostics = TRUE) first.")
  }
  if (!is.null(file)) {
    grDevices::pdf(file, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  first <- diagnostics$first_pass
  final <- diagnostics$final_pass
  finite_hist <- function(x, breaks = 30L, ...) {
    x <- x[is.finite(x)]
    if (!length(x)) {
      graphics::plot.new()
      graphics::text(0.5, 0.5, "No finite values")
    } else {
      graphics::hist(x, breaks = breaks, ...)
    }
  }
  add_threshold <- function(value, label) {
    if (length(value) == 1L && is.finite(value)) {
      graphics::abline(v = value, col = "firebrick", lwd = 2, lty = 2)
      graphics::legend("topright", legend = label, col = "firebrick",
                       lty = 2, lwd = 2, bty = "n")
    } else {
      graphics::mtext(label, side = 3, col = "firebrick", cex = 0.8)
    }
  }

  gate <- diagnostics$gates
  thr <- first$exclusion_margin_threshold[1L]
  status <- sprintf(
    "First-pass agreement %.1f%% (required %.1f%%): %s",
    100 * gate$first_pass_agreement, 100 * gate$min_agreement,
    if (gate$low_agreement) "below gate, downstream calls suppressed" else "passed"
  )
  grp <- factor(ifelse(first$correct, "Predicted = observed",
                       "Predicted != observed"),
                levels = c("Predicted = observed", "Predicted != observed"))
  cols <- c("steelblue", "firebrick")
  margin_box <- function(values, ylab, main, cutoff = NULL) {
    keep <- is.finite(values)
    y <- log10(1 + values[keep])
    g <- grp[keep]
    graphics::boxplot(
      y ~ g, col = grDevices::adjustcolor(cols, 0.4), outline = FALSE,
      ylim = range(c(y, 0), finite = TRUE), xlab = "", ylab = ylab,
      main = main, drop = FALSE
    )
    graphics::stripchart(
      y ~ g, vertical = TRUE, method = "jitter", add = TRUE, pch = 16,
      cex = 0.6, col = grDevices::adjustcolor(cols, 0.5)
    )
    graphics::mtext(sprintf("n = %d", as.integer(table(g))),
                    side = 1, at = 1:2, line = 2.5, cex = 0.8)
    if (!is.null(cutoff) && is.finite(cutoff)) {
      graphics::abline(h = log10(1 + cutoff), lty = 2, lwd = 2)
    }
  }
  graphics::par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
  margin_box(first$top_margin, "log10(1 + top - second score)",
             "First-pass top margin", thr)
  margin_box(first$annotation_margin,
             "log10(1 + top - observed-label score)",
             "First-pass annotation margin", thr)
  graphics::mtext(paste(status, "| dashed = exclusion cutoff"),
                  outer = TRUE, cex = 0.9,
                  col = if (gate$low_agreement) "firebrick" else "darkgreen")
  graphics::par(oma = c(0, 0, 0, 0))

  graphics::par(mfrow = c(1, 2))
  finite_hist(
    first$top_margin[first$correct],
    main = "Correct first-pass predictions",
    xlab = "Top score - second score",
    col = "grey80", border = "white"
  )
  add_threshold(first$exclusion_margin_threshold[1L],
                "Exclusion cutoff: low quantile of correct top margins")
  finite_hist(
    first$annotation_margin[!first$correct],
    main = "Disagreeing first-pass predictions",
    xlab = "Top score - observed-label score",
    col = "grey80", border = "white"
  )
  add_threshold(first$exclusion_margin_threshold[1L],
                "Compared with annotation margin")

  annotation <- diagnostics$annotation_threshold
  graphics::par(mfrow = c(1, 1))
  groups <- split(annotation$margin, annotation$population)
  finite_groups <- lapply(groups, function(x) x[is.finite(x)])
  finite_groups <- finite_groups[lengths(finite_groups) > 0L]
  if (!length(finite_groups)) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No annotation-margin values")
  } else {
    pooled <- unlist(finite_groups, use.names = FALSE)
    breaks <- graphics::hist(pooled, plot = FALSE, breaks = 35L)$breaks
    colors <- grDevices::adjustcolor(
      c("grey50", "darkorange", "steelblue", "purple"),
      alpha.f = 0.45
    )
    graphics::hist(
      finite_groups[[1L]], breaks = breaks, col = colors[1L],
      border = "white", xlim = range(breaks),
      main = "Gate 2: annotation-margin threshold",
      xlab = "Top score - tested label score"
    )
    if (length(finite_groups) > 1L) {
      for (i in 2:length(finite_groups)) {
        graphics::hist(finite_groups[[i]], breaks = breaks, add = TRUE,
                       col = colors[(i - 1L) %% length(colors) + 1L],
                       border = "white")
      }
    }
    add_threshold(result$thresholds$annotation_margin,
                  paste("Selected from", diagnostics$annotation_threshold_source))
    graphics::legend(
      "topright", legend = names(finite_groups),
      fill = colors[seq_along(finite_groups)], bty = "n", cex = 0.8
    )
  }

  abstention <- diagnostics$abstention
  curve <- abstention$curve
  graphics::par(mfrow = c(1, 2))
  if (nrow(curve)) {
    graphics::plot(
      curve$threshold, curve$agreement, type = "b", pch = 16,
      xlab = "Top score - second score cutoff", ylab = "Agreement",
      ylim = c(0, 1), main = "Gate 3: abstention accuracy"
    )
    graphics::abline(h = abstention$min_accuracy, col = "grey40", lty = 3)
    add_threshold(abstention$threshold, "Selected abstention cutoff")
    graphics::plot(
      curve$threshold, curve$n_retained / max(curve$n_retained),
      type = "b", pch = 16, xlab = "Top score - second score cutoff",
      ylab = "Fraction retained", ylim = c(0, 1),
      main = paste("Coverage; source:", abstention$source)
    )
    add_threshold(abstention$threshold, "Selected abstention cutoff")
  } else {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No abstention data")
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No coverage data")
  }

  unknown <- diagnostics$unknown
  graphics::par(mfrow = c(1, 1))
  known <- unknown$known_margins[is.finite(unknown$known_margins)]
  strangers <- unknown$simulated_unknown_margins[
    is.finite(unknown$simulated_unknown_margins)
  ]
  if (length(known) || length(strangers)) {
    all_margins <- c(known, strangers)
    breaks <- graphics::hist(all_margins, plot = FALSE, breaks = 35L)$breaks
    if (length(known)) {
      graphics::hist(
        known, breaks = breaks, col = grDevices::adjustcolor("steelblue", 0.5),
        border = "white", xlab = "Top score - second score",
        main = "Gate 4: unknown-individual margin"
      )
    } else {
      graphics::hist(
        strangers, breaks = breaks,
        col = grDevices::adjustcolor("darkorange", 0.5),
        border = "white", xlab = "Top score - second score",
        main = "Gate 4: unknown-individual margin"
      )
      strangers <- numeric()
    }
    if (length(strangers)) {
      graphics::hist(strangers, breaks = breaks, add = TRUE,
                     col = grDevices::adjustcolor("darkorange", 0.5),
                     border = "white")
    }
    add_threshold(unknown$threshold,
                  if (unknown$enabled) "Unknown cutoff" else
                    "Unknown detection disabled")
    graphics::legend(
      "topright", legend = c("Enrolled reference", "Simulated strangers"),
      fill = grDevices::adjustcolor(c("steelblue", "darkorange"), 0.5),
      bty = "n"
    )
  } else {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "Unknown-individual gate not assessed")
  }

  performance <- gate$performance
  misannotation <- performance$misannotation_detection
  feasibility <- gate$feasibility
  detection_status <- if (is.null(feasibility)) "not assessed" else
    feasibility$misannotation_detection
  recovery_status <- if (is.null(feasibility)) "not assessed" else
    feasibility$identity_recovery
  unknown_status <- if (is.null(feasibility)) "not assessed" else
    feasibility$unknown_individual_detection
  graphics::par(mfrow = c(2, 2))
  graphics::plot(
    c(0, 1), c(0, 1), type = "n", xlim = c(0.5, 2.5), ylim = c(0, 1),
    xaxt = "n", xlab = "", ylab = "Rate",
    main = paste("Gate 5: misannotation", detection_status)
  )
  graphics::axis(1, at = 1:2, labels = c("Sensitivity", "Specificity"))
  graphics::abline(h = c(0.8, gate$target_specificity),
                   col = c("firebrick", "darkorange"), lty = 2)
  graphics::segments(1, 0, 1, misannotation$sensitivity,
                     lwd = 8, col = "steelblue")
  graphics::segments(2, 0, 2, misannotation$specificity,
                     lwd = 8, col = "steelblue")
  graphics::mtext(
    sprintf("%s; %d repeat(s)", sub("_", " ", misannotation$source),
            misannotation$n_repeats_validated),
    side = 3, cex = 0.7
  )

  recovery <- c(
    Overall = performance$identity_recovery_accuracy,
    `Swapped samples` = performance$recovery_accuracy_misannotated
  )
  graphics::plot(
    c(0, 1), c(0, 1), type = "n", xlim = c(0.5, 2.5), ylim = c(0, 1),
    xaxt = "n", xlab = "", ylab = "Accuracy",
    main = paste("Identity recovery", recovery_status)
  )
  graphics::axis(1, at = 1:2, labels = names(recovery))
  graphics::abline(h = 0.8, col = "firebrick", lty = 2)
  for (i in which(is.finite(recovery))) {
    graphics::segments(i, 0, i, recovery[[i]], lwd = 8, col = "steelblue")
  }

  unknown_sensitivity <- if (is.null(unknown$metrics)) NA_real_ else
    unknown$metrics$sensitivity
  graphics::plot(
    c(0, 1), c(0, 1), type = "n", xlim = c(0.5, 1.5), ylim = c(0, 1),
    xaxt = "n", xlab = "", ylab = "Sensitivity",
    main = paste("Unknown detection", unknown_status,
                 sprintf("(n = %d)", performance$unknown_simulation_n))
  )
  graphics::axis(1, at = 1, labels = "Simulated strangers")
  graphics::abline(h = 0.8, col = "firebrick", lty = 2)
  if (is.finite(unknown_sensitivity)) {
    graphics::segments(1, 0, 1, unknown_sensitivity,
                       lwd = 8, col = "steelblue")
  }
  pairwise <- performance$pairwise_swap_detection
  graphics::plot(
    c(0, 1), c(0, 1), type = "n", xlim = c(0.5, 2.5), ylim = c(0, 1),
    xaxt = "n", xlab = "", ylab = "Rate",
    main = "Pairwise-swap diagnostic (not a feasibility gate)"
  )
  graphics::axis(1, at = 1:2, labels = c("Sensitivity", "Specificity"))
  graphics::abline(h = c(0.8, gate$target_specificity),
                   col = c("firebrick", "darkorange"), lty = 2)
  graphics::segments(1, 0, 1, pairwise$sensitivity,
                     lwd = 8, col = "steelblue")
  graphics::segments(2, 0, 2, pairwise$specificity,
                     lwd = 8, col = "steelblue")

  invisible(list(
    first_pass = first,
    final_pass = final,
    annotation_threshold = annotation,
    abstention = abstention,
    unknown = unknown,
    gates = gate
  ))
}
