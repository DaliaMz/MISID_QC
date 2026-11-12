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
Summarizes the first-pass holdout, detection and unknown thresholds, recovery
coverage, and final sample calls. Run \code{\link{check_annotations}} with
\code{diagnostics = TRUE} for the calibration plots. Without diagnostics, the
report explicitly marks those plots as unavailable.
}
\details{
Internal validation uses simulated errors and is not ground truth or a
sample-level probability. First-pass agreement is with observed labels.
}
