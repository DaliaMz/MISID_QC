#' Plot identity agreement against prediction coverage
#'
#' Draws the risk-coverage curve from [check_annotations()]. Samples are
#' retained when their top-versus-second raw score margin is at or above a
#' cutoff; the curve shows agreement among retained samples as the cutoff
#' rises (coverage falls).
#'
#' Three curves are drawn where available: agreement with the observed labels
#' (a lower bound when some labels are wrong), agreement with the simulated
#' labels, and agreement with the true labels in the label-swap simulation
#' (`agreement_real`), which is the only curve with known ground truth.
#'
#' @param result Output of [check_annotations()].
#' @param ... Further arguments passed to [graphics::plot()].
#' @return The risk-coverage table, invisibly.
#' @export
plot_risk_coverage <- function(result, ...) {
  table <- result$risk_coverage
  if (!is.data.frame(table) || !nrow(table)) {
    stop("result does not contain a risk_coverage table.")
  }
  observed <- table[table$source == "observed", , drop = FALSE]
  simulated <- table[table$source == "simulated", , drop = FALSE]
  values <- c(observed$agreement, simulated$agreement,
              simulated$agreement_real)
  graphics::plot(
    NA, xlim = c(0, 1), ylim = range(values, na.rm = TRUE),
    xlab = "Proportion of samples retained",
    ylab = "Identity agreement",
    main = "Identity agreement vs prediction coverage", ...
  )
  graphics::lines(observed$coverage, observed$agreement, lwd = 2)
  labels <- "Observed labels"
  colors <- "black"
  if (nrow(simulated)) {
    graphics::lines(simulated$coverage, simulated$agreement, lwd = 2,
                    col = "orange")
    graphics::lines(simulated$coverage, simulated$agreement_real, lwd = 2,
                    col = "blue")
    labels <- c(labels, "Simulated labels", "True labels (simulation)")
    colors <- c(colors, "orange", "blue")
  }
  graphics::legend("bottomleft", legend = labels, col = colors, lwd = 2,
                   bty = "n")
  invisible(table)
}
