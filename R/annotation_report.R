utils::globalVariables(c("x", "y", "rate", "metric", "group", "held"))

.report_rate_curve <- function(positives, negatives, cutoff, direction = "upper") {
  positives <- sort(positives[is.finite(positives)])
  negatives <- sort(negatives[is.finite(negatives)])
  values <- sort(unique(c(positives, negatives, cutoff[is.finite(cutoff)])))
  if (!length(values)) return(data.frame(x = numeric(), rate = numeric(), metric = character()))
  span <- max(diff(range(values)), abs(values), 1)
  thresholds <- sort(unique(c(values, min(values) - span * 0.01,
                              max(values) + span * 0.01)))
  if (length(thresholds) > 400L) {
    thresholds <- sort(unique(c(thresholds[unique(round(seq(1, length(thresholds),
                                                            length.out = 400L)))],
                                cutoff[is.finite(cutoff)])))
  }
  fraction <- function(v, lower) {
    if (!length(v)) return(rep(NA_real_, length(thresholds)))
    if (lower) {
      findInterval(thresholds, v, left.open = TRUE) / length(v)
    } else {
      1 - findInterval(thresholds, v) / length(v)
    }
  }
  lower <- identical(direction, "lower")
  sensitivity <- fraction(positives, lower)
  specificity <- 1 - fraction(negatives, lower)
  data.frame(x = rep(thresholds, 2L),
             rate = c(sensitivity, specificity),
             metric = rep(c("Sensitivity", "Specificity"), each = length(thresholds)))
}

.report_holdout_reference <- function(result) {
  first <- result$diagnostics$first_pass
  scores <- result$diagnostics$first_pass_scores
  if (is.null(first) || is.null(scores)) return(numeric())
  scores <- scores[first$correct, , drop = FALSE]
  if (!nrow(scores)) return(numeric())
  top <- cbind(seq_len(nrow(scores)), max.col(scores, ties.method = "first"))
  margins <- scores[top] - scores
  margins[top] <- NA_real_
  as.vector(margins[!is.na(margins)])
}

.report_recovery_curves <- function(simulation) {
  empty <- data.frame(x = numeric(), y = numeric(), group = character())
  if (is.null(simulation) || !nrow(simulation)) return(empty)
  eligible <- simulation$actually_misannotated |
    !simulation$agreement_simulated
  curve <- function(keep, label) {
    if (!any(keep)) return(empty)
    risk <- .risk_coverage(simulation$top_margin[keep],
                           simulation$agreement_simulated[keep],
                           simulation$agreement_real[keep])
    data.frame(x = risk$coverage, y = risk$agreement_real, group = label)
  }
  rbind(curve(rep(TRUE, nrow(simulation)), "All simulated samples"),
        curve(eligible, "Swaps and unchanged disagreements"))
}

.report_detection_summary <- function(result) {
  sim <- result$mismatch_simulation
  det <- result$performance$misannotation_detection
  cutoff <- result$thresholds$annotation_margin
  positives <- if (is.null(sim)) numeric() else
    sim$annotation_margin[sim$actually_misannotated]
  negatives <- if (is.null(sim)) numeric() else
    sim$annotation_margin[!sim$actually_misannotated]
  annotation <- result$diagnostics$annotation_threshold
  reference <- if (is.null(annotation)) numeric() else
    annotation$margin[annotation$population == "reference_negative"]
  pct <- function(x) {
    if (length(x) != 1L || !is.finite(x)) "not available" else
      sprintf("%.2f%%", 100 * x)
  }
  flagged_rate <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x) || length(cutoff) != 1L || is.na(cutoff)) {
      return(NA_real_)
    }
    mean(x > cutoff)
  }
  source <- "In-sample simulated-swap verdict rates"
  c(
    sprintf("Detection targets: sensitivity 80%%; specificity %s.",
            pct(result$thresholds$target_specificity)),
    sprintf("In-sample at cutoff %s: sensitivity %s (swapped n=%d);",
            format(cutoff, digits = 4), pct(flagged_rate(positives)),
            sum(is.finite(positives))),
    sprintf("simulated specificity %s (unchanged n=%d); reference specificity %s.",
            pct(1 - flagged_rate(negatives)), sum(is.finite(negatives)),
            pct(1 - flagged_rate(reference))),
    sprintf("%s: sensitivity %s; simulated specificity %s.",
            source, pct(det$sensitivity), pct(det$simulated_specificity)),
    "Rows may reuse samples across repeats; pre-swap labels are assumed correct.",
    "Estimated errors using the simulated verdict rates above:",
    .format_detection_error_estimate(result$performance$estimated_detection_errors)
  )
}

#' Write a printable annotation QC report
#'
#' Writes a six-page report covering the overview, first-pass holdout,
#' misannotation calibration, recovery coverage, unknown detection, and actual
#' sample distributions. Recovery coverage shows both all simulated samples
#' and the filtered correction-calibration samples, each relative to its own
#' population. The marked correction cutoff is selected from the filtered
#' samples. Set `diagnostics = TRUE` in [check_annotations()] for
#' the detailed plots; unavailable sections are labelled.
#' The misannotation page shows detection targets, in-sample and verdict
#' simulation rates, and exploratory estimated detection errors together.
#'
#' Internal validation uses simulated errors and is not ground truth or a
#' sample-level probability. First-pass agreement is with observed labels.
#' Estimated detection errors are exploratory projections, not confidence
#' intervals; simulated swaps may not represent real errors or unknown people.
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
  diag <- result$diagnostics
  sim <- result$mismatch_simulation
  det <- result$performance$misannotation_detection
  unknown <- if (is.null(diag)) NULL else diag$unknown
  first <- if (is.null(diag)) NULL else diag$first_pass
  pct <- function(x) {
    if (length(x) != 1L || is.na(x) || !is.finite(x)) "not available" else
      sprintf("%.1f%%", 100 * x)
  }
  num <- function(x) {
    if (length(x) != 1L || is.na(x)) "not available" else format(x, digits = 4)
  }
  verdict <- function(name) {
    x <- result$feasibility[[name]]
    if (is.null(x)) "not assessed" else gsub("_", " ", x, fixed = TRUE)
  }
  rate <- function(x, cutoff, direction = "upper") {
    x <- x[is.finite(x)]
    if (!length(x) || length(cutoff) != 1L || is.na(cutoff)) return(NA_real_)
    if (direction == "upper") mean(x > cutoff) else mean(x < cutoff)
  }
  theme <- ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14),
                   plot.subtitle = ggplot2::element_text(size = 9),
                   plot.caption = ggplot2::element_text(hjust = 0, size = 9),
                   legend.position = "bottom")
  empty <- function(title, message) {
    ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = message, size = 4) +
      ggplot2::labs(title = title) + ggplot2::theme_void() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14))
  }
  cutoff_line <- function(plot, cutoff, horizontal = FALSE) {
    if (length(cutoff) != 1L || !is.finite(cutoff)) return(plot)
    if (horizontal) plot + ggplot2::geom_hline(yintercept = cutoff,
                                                colour = "firebrick", linetype = "dashed") else
      plot + ggplot2::geom_vline(xintercept = cutoff,
                                 colour = "firebrick", linetype = "dashed")
  }
  draw_page <- function(plots, rows, cols) {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(rows, cols)))
    for (i in seq_along(plots)) {
      print(plots[[i]], vp = grid::viewport(layout.pos.row = (i - 1L) %/% cols + 1L,
                                            layout.pos.col = (i - 1L) %% cols + 1L))
    }
    grid::popViewport()
  }
  curve_plot <- function(positives, negatives, cutoff, direction, title, caption, xlab) {
    curve <- .report_rate_curve(positives, negatives, cutoff, direction)
    curve <- curve[is.finite(curve$rate), , drop = FALSE]
    if (!nrow(curve)) return(empty(title, "No calibration data available"))
    plot <- ggplot2::ggplot(curve, ggplot2::aes(x = x, y = rate, colour = metric)) +
      ggplot2::geom_step(linewidth = 0.7,
                         direction = if (direction == "lower") "vh" else "hv") +
      ggplot2::scale_colour_manual(values = c(Sensitivity = "#c35420", Specificity = "#27658b")) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(title = title, x = xlab, y = "Rate",
                    colour = NULL, caption = caption) + theme
    cutoff_line(plot, cutoff)
  }
  histogram <- function(values, cutoff, title, caption, xlab) {
    values <- values[is.finite(values)]
    if (!length(values)) return(empty(title, paste("No assessed values.", caption)))
    data <- data.frame(x = values)
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = x)) +
      ggplot2::geom_histogram(bins = 30, fill = "#4286a5", colour = "white") +
      ggplot2::labs(title = title, x = xlab, y = "Samples",
                    caption = caption) + theme
    cutoff_line(plot, cutoff)
  }

  grDevices::pdf(file, width = 11.7, height = 8.3, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  # 1. Overview: verdicts and sample-level outcomes.
  statuses <- as.data.frame(table(p$annotation_status, useNA = "ifany"))
  names(statuses) <- c("group", "y")
  statuses$group <- gsub("_", " ", statuses$group, fixed = TRUE)
  status_plot <- ggplot2::ggplot(statuses, ggplot2::aes(x = group, y = y)) +
    ggplot2::geom_col(fill = "#4286a5") + ggplot2::coord_flip() +
    ggplot2::labs(title = "Final sample statuses", x = NULL, y = "Samples",
                  caption = "Unknown flags are separate and may overlap.") + theme
  flagged <- which((!is.na(p$suspicious) & p$suspicious) |
                     (!is.na(p$potential_unknown_individual) & p$potential_unknown_individual))
  calls <- p[head(flagged, 6L), , drop = FALSE]
  call_lines <- if (!nrow(calls)) "No samples require inspection." else {
    paste(calls$sample_id, "| observed:", calls$observed_id,
          "| correction:", ifelse(is.na(calls$candidate_correction), "-", calls$candidate_correction),
          "| suspicious:", ifelse(is.na(calls$suspicious), "NA", calls$suspicious),
          "| unknown:", ifelse(is.na(calls$potential_unknown_individual), "NA",
                               calls$potential_unknown_individual))
  }
  overview <- c(
    sprintf("Samples: %d | First-pass agreement: %s | Held out: %d",
            nrow(p), pct(result$calibration$first_pass_agreement),
            result$held_out$n_confidently_misassigned),
    sprintf("Detection: %s | sensitivity %s | simulated specificity %s",
            verdict("misannotation_detection"), pct(det$sensitivity), pct(det$simulated_specificity)),
    .format_detection_error_estimate(result$performance$estimated_detection_errors),
    sprintf("Detection cutoff: %s | target specificity: %s",
            num(result$thresholds$annotation_margin), pct(result$thresholds$target_specificity)),
    sprintf("Recovery: %s | filtered simulated accuracy %s (n=%s) | correction cutoff %s",
            verdict("identity_recovery"), pct(result$calibration$correction_margin$accuracy),
            num(result$calibration$correction_margin$n_retained),
            num(result$thresholds$correction_margin)),
    sprintf("Correction gate: %s",
            if (!is.null(result$correction_gate)) result$correction_gate$reason
            else "Correction gate result unavailable."),
    sprintf("Unknown: %s | sensitivity %s | specificity %s",
            verdict("unknown_individual_detection"),
            pct(result$performance$unknown_individual_detection$sensitivity),
            pct(result$performance$unknown_individual_detection$specificity)),
    sprintf("Unknown cutoff: %s | simulated unknown samples: %s",
            if (is.null(unknown)) "see thresholds" else num(unknown$threshold),
            num(result$performance$unknown_simulation_n)),
    sprintf("Samples needing inspection: %d (showing up to 6)", length(flagged)),
    call_lines,
    "Internal simulation estimates depend on the simulated error mechanism;",
    "pre-swap labels are assumed, not verified truth. See predictions for all samples."
  )
  overview_data <- data.frame(x = 0, y = seq_along(overview), group = overview)
  overview_plot <- ggplot2::ggplot(overview_data, ggplot2::aes(x = x, y = -y, label = group)) +
    ggplot2::geom_text(hjust = 0, size = 3.2) + ggplot2::xlim(0, 1) +
    ggplot2::labs(title = "1. Annotation QC overview") +
    ggplot2::theme_void() + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14))
  draw_page(list(overview_plot, status_plot), 2, 1)

  # 2. First-pass calibration population and sample decisions.
  reference <- .report_holdout_reference(result)
  threshold <- result$calibration$exclusion_margin_threshold
  if (length(reference) >= 2L && length(unique(reference[is.finite(reference)])) > 1L) {
    d <- data.frame(x = reference[is.finite(reference)])
    holdout_density <- ggplot2::ggplot(d, ggplot2::aes(x = x)) +
      ggplot2::geom_density(fill = "#8bc2d7", colour = "#27658b", alpha = 0.5) +
      ggplot2::labs(title = "2. Holdout cutoff calibration",
                    subtitle = sprintf("Agreements with observed labels; all-alternative reference; n = %d margins",
                                       nrow(d)),
                    x = "Reference margin", y = "Density",
                    caption = paste("Dashed = cutoff", num(threshold))) + theme
    holdout_density <- cutoff_line(holdout_density, threshold)
  } else {
    holdout_density <- empty("2. Holdout cutoff calibration",
                             if (is.null(first)) "Requires diagnostics = TRUE" else
                               "Insufficient variation for a density plot")
  }
  if (!is.null(first) && nrow(first)) {
    d <- data.frame(group = ifelse(first$correct, "Agrees", "Disagrees"),
                    y = first$annotation_margin,
                    held = ifelse(first$held_out, "Held out", "Not held out"))
    holdout_points <- ggplot2::ggplot(d, ggplot2::aes(x = group, y = y)) +
      ggplot2::geom_boxplot(outlier.shape = NA, fill = "grey90") +
      ggplot2::geom_point(ggplot2::aes(colour = held),
                          position = ggplot2::position_jitter(width = 0.15, height = 0, seed = 1),
                          alpha = 0.65, size = 1.4) +
      ggplot2::scale_colour_manual(values = c("Held out" = "firebrick", "Not held out" = "#27658b")) +
      ggplot2::labs(title = "First-pass sample margins",
                    subtitle = "Agreement with observed labels, not verified truth",
                    x = NULL, y = "Annotation margin", colour = NULL,
                    caption = "The cutoff applies to disagreements only.") + theme
    holdout_points <- cutoff_line(holdout_points, threshold, horizontal = TRUE)
  } else {
    holdout_points <- empty("First-pass sample margins", "Requires diagnostics = TRUE")
  }
  draw_page(list(holdout_density, holdout_points), 2, 1)

  # 3. Mismatch simulation: threshold-only detection.
  cutoff <- result$thresholds$annotation_margin
  positives <- if (is.null(sim)) numeric() else sim$annotation_margin[sim$actually_misannotated]
  negatives <- if (is.null(sim)) numeric() else sim$annotation_margin[!sim$actually_misannotated]
  caption <- "Final cutoff uses both negative groups; curves are in-sample."
  mismatch_plot <- if (is.null(sim) || !length(positives) || !length(negatives)) {
    empty(paste("3. Misannotation detection |", verdict("misannotation_detection")),
          "Simulated swaps or unchanged samples unavailable")
  } else curve_plot(positives, negatives, cutoff, "upper",
                    paste("3. Misannotation detection |", verdict("misannotation_detection")),
                    caption, "Annotation margin cutoff")
  summary <- .report_detection_summary(result)
  summary_data <- data.frame(x = 0, y = -seq_along(summary), group = summary)
  summary_plot <- ggplot2::ggplot(summary_data, ggplot2::aes(x = x, y = y, label = group)) +
    ggplot2::geom_text(hjust = 0, size = 3.2) + ggplot2::xlim(0, 1) +
    ggplot2::theme_void()
  draw_page(list(mismatch_plot, summary_plot), 2, 1)

  # 4. Recovery coverage in the full and correction-calibration populations.
  eligible <- if (is.null(sim)) logical() else
    sim$actually_misannotated | !sim$agreement_simulated
  recovery_curves <- .report_recovery_curves(sim)
  if (nrow(recovery_curves)) {
    coverage <- result$calibration$correction_margin$coverage
    correction_cutoff <- result$thresholds$correction_margin
    has_cutoff <- length(coverage) == 1L && is.finite(coverage) &&
      length(correction_cutoff) == 1L && is.finite(correction_cutoff)
    singletons <- names(which(table(recovery_curves$group) == 1L))
    cutoff_caption <- if (has_cutoff) {
      paste0("Dashed = filtered coverage ", pct(coverage),
             " at filtered correction cutoff ", num(correction_cutoff))
    } else {
      "No filtered correction cutoff available"
    }
    recovery_plot <- ggplot2::ggplot(recovery_curves,
                                      ggplot2::aes(x = x, y = y, colour = group)) +
      ggplot2::geom_line(data = recovery_curves[!recovery_curves$group %in%
                                                  singletons, , drop = FALSE],
                         linewidth = 0.8) +
      ggplot2::geom_point(data = recovery_curves[recovery_curves$group %in%
                                                   singletons, , drop = FALSE],
                          size = 2) +
      ggplot2::scale_colour_manual(
        values = c("All simulated samples" = "#c46a24",
                   "Swaps and unchanged disagreements" = "#27658b"),
        name = NULL
      ) +
      ggplot2::scale_x_continuous(limits = c(0, 1)) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(title = paste("4. Identity recovery |", verdict("identity_recovery")),
                    subtitle = "Agreement with pre-swap labels",
                    x = "Coverage within each population", y = "Agreement",
                    caption = paste0(
                      cutoff_caption, "; filtered n=", sum(eligible), ", all n=", nrow(sim),
                      ".\nSame coverage can mean different cutoffs; curves are in-sample.",
                      "\nSuspicious, recovery-feasibility and unknown gates are not included.")) + theme
    if (has_cutoff) recovery_plot <- recovery_plot +
      ggplot2::geom_vline(xintercept = coverage, colour = "firebrick", linetype = "dashed")
  } else {
    recovery_plot <- empty(paste("4. Identity recovery |", verdict("identity_recovery")),
                           "Simulated recovery data unavailable")
  }
  draw_page(list(recovery_plot), 1, 1)

  # 5. Simulated unknowns versus reference samples.
  if (is.null(unknown)) {
    unknown_plot <- empty("5. Unknown detection", "Requires diagnostics = TRUE")
  } else if (!isTRUE(unknown$enabled)) {
    unknown_plot <- empty("5. Unknown detection", "Unknown detection was disabled")
  } else {
    direction <- if (identical(unknown$statistic, "top_margin")) "lower" else "upper"
    unknown_caption <- paste0(
      "Cutoff ", num(unknown$threshold), "; statistic: ", unknown$statistic,
      "; ", if (direction == "upper") "larger" else "smaller", " values are flagged.\n",
      "Simulated sensitivity ", pct(rate(unknown$simulated_unknown_values, unknown$threshold, direction)),
      " (n=", sum(is.finite(unknown$simulated_unknown_values)), "); reference specificity ",
      pct(1 - rate(unknown$known_values, unknown$threshold, direction)),
      " (n=", sum(is.finite(unknown$known_values)), ").\n",
      "The same simulated unknowns can be used to select the cutoff; estimates may be optimistic."
    )
    unknown_plot <- curve_plot(unknown$simulated_unknown_values, unknown$known_values,
                               unknown$threshold, direction,
                               paste("5. Unknown detection |", verdict("unknown_individual_detection")),
                               unknown_caption, if (direction == "upper") {
                                 "Nearest-class distance cutoff"
                               } else "Top-margin cutoff")
  }
  draw_page(list(unknown_plot), 1, 1)

  # 6. Actual sample distributions; historical results may include dropped-refit scores.
  final <- p[!p$excluded_from_refit, , drop = FALSE]
  real_holdout <- if (is.null(first)) {
    empty("6. First-pass holdout", "Requires diagnostics = TRUE")
  } else histogram(first$annotation_margin, result$calibration$exclusion_margin_threshold,
                   "6. First-pass holdout", "Cutoff applies only to disagreements.",
                   "First-pass annotation margin")
  missing_refit <- sum(p$excluded_from_refit)
  real_detection <- histogram(final$annotation_score_margin, result$thresholds$annotation_margin,
                              "Final annotation margins",
                              paste("Dashed = calibrated cutoff; excluded from refit:", missing_refit),
                              "Annotation margin")
  real_recovery <- histogram(final$top_score_margin, result$thresholds$correction_margin,
                             "Final top margins",
                             "Cutoff applies to disagreements; unknown flags block corrections, not clean predictions.",
                             "Top margin")
  real_unknown <- if (is.null(unknown)) {
    empty("Final unknown scores", "Requires diagnostics = TRUE")
  } else if (!isTRUE(unknown$enabled)) {
    empty("Final unknown scores", "Unknown detection was disabled")
  } else histogram(final$unknown_detection_score, unknown$threshold,
                   "Final unknown scores", paste("Unavailable scores:", sum(!is.finite(p$unknown_detection_score))),
                   "Nearest-class distance")
  draw_page(list(real_holdout, real_detection, real_recovery, real_unknown), 2, 2)
  invisible(file)
}
