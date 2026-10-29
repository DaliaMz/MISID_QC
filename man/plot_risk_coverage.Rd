\name{plot_risk_coverage}
\alias{plot_risk_coverage}
\title{Plot identity agreement against prediction coverage}
\usage{
plot_risk_coverage(result, ...)
}
\arguments{
\item{result}{Output of \code{check_annotations()}.}

\item{...}{Further arguments passed to \code{graphics::plot()}.}
}
\value{
The risk-coverage table, invisibly.
}
\description{
Draws the risk-coverage curve from \code{check_annotations()}. Samples are
retained when their top-versus-second raw score margin is at or above a
cutoff; the curve shows agreement among retained samples as the cutoff rises.
}
\details{
Three curves are drawn where available: agreement with the observed labels (a
lower bound when some labels are wrong), agreement with the simulated labels,
and agreement with the true labels in the label-swap simulation, which is the
only curve with known ground truth.
}
