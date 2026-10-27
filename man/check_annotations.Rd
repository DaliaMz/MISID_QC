\name{check_annotations}
\alias{check_annotations}
\title{Assess repeated-measures RNA-seq sample annotations}
\usage{
check_annotations(expr, metadata, sample_col, id_col, shrinkage = 0.1,
  target_specificity = 0.95, assess_feasibility = TRUE)
}
\arguments{
\item{expr}{Numeric expression matrix, features by samples, with unique
sample names in \code{colnames(expr)}.}
\item{metadata}{Data frame containing sample and observed identity columns.}
\item{sample_col}{Name of the sample identifier column.}
\item{id_col}{Name of the observed individual identifier column.}
\item{shrinkage}{Shrinkage fraction for the diagonal covariance target.}
\item{target_specificity}{Target specificity used to select the lower-tail
thresholds for annotation support and unknown detection.}
\item{assess_feasibility}{Whether to return task-specific feasibility
assessments based on internal validation.}
}
\value{
A list with \code{predictions}, \code{performance}, \code{calibration},
\code{feasibility}, and \code{thresholds} components.
}
\description{
Uses leave-one-sample-out shrinkage linear discriminant scores and internal
artificial-label and held-out-individual experiments to describe dataset-
specific annotation evidence.
}
\details{
The expression matrix must contain transformed, finite values with features
in rows and samples in columns. Each enrolled individual needs at least two
samples. This implementation uses all supplied features and does not normalize
or select features. The evaluated sample is excluded from its own fitted
reference.

Temperature-scaled probabilities describe candidate identities, not the
probability that an observed annotation is correct. The latter is not reported
because simulated error prevalence does not identify real-data error
prevalence. Internal thresholds and performance are estimated from the same
out-of-fold predictions and should be treated as exploratory.
}
