\name{annotation_report}
\alias{annotation_report}
\title{Write a printable annotation QC report}
\usage{
annotation_report(result, file = "annotation_report.pdf")
}
\arguments{
\item{result}{A result from \code{\link{check_annotations}}.}
\item{file}{Path to the PDF to create.}
}
\value{The PDF path, invisibly.}
\description{
Writes a six-page ggplot2 report covering the overview, first-pass holdout,
misannotation calibration, recovery coverage, unknown detection, and actual
sample distributions. Recovery coverage shows all simulated samples and the
filtered correction-calibration samples, each relative to its own population.
The marked correction cutoff is selected from the filtered samples. Run
\code{\link{check_annotations}} with
\code{diagnostics = TRUE} for the detailed plots; unavailable sections are labelled.
The misannotation page shows detection targets, in-sample and verdict
simulation rates, and exploratory estimated detection errors together.
}
\details{
Internal validation uses simulated errors and is not ground truth or a
sample-level probability. First-pass agreement is with observed labels.
Estimated detection errors are exploratory projections, not confidence
intervals; simulated swaps may not represent real errors or unknown people.
Margins and unknown scores are plotted on their original scales.
}
