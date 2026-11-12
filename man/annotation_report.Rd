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
sample distributions. Run \code{\link{check_annotations}} with
\code{diagnostics = TRUE} for the detailed plots; unavailable sections are labelled.
}
\details{
Internal validation uses simulated errors and is not ground truth or a
sample-level probability. First-pass agreement is with observed labels.
Margins and unknown scores are plotted on their original scales.
}
