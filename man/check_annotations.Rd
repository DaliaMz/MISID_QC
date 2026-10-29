\name{check_annotations}
\alias{check_annotations}
\title{Assess repeated-measures RNA-seq sample annotations}
\usage{
check_annotations(expr, metadata, sample_col, id_col, shrinkage = 0.1,
  target_specificity = 0.95, assess_feasibility = TRUE,
  progress = interactive(), cv_folds = NULL,
  min_agreement = 0.7, margin_coverage = 0.95, refit_clean = TRUE,
  mismatch_rate = 0.1, min_accuracy = 0.9, unknown_mad = 3)
}
\arguments{
\item{expr}{Numeric expression matrix, features by samples, with unique
sample names in \code{colnames(expr)}.}
\item{metadata}{Data frame containing sample and observed identity columns.}
\item{sample_col}{Name of the sample identifier column.}
\item{id_col}{Name of the observed individual identifier column.}
\item{shrinkage}{Shrinkage fraction for the diagonal covariance target.}
\item{target_specificity}{Target specificity used to select the lower-tail
annotation margin threshold.}
\item{assess_feasibility}{Whether to return task-specific feasibility
assessments based on internal validation.}
\item{progress}{Whether to report progress updates. Defaults to \code{TRUE}
in interactive sessions and \code{FALSE} otherwise.}
\item{cv_folds}{\code{NULL} (default) uses leave-one-sample-out validation;
otherwise, an integer from 2 to the smallest individual sample count selects
the number of stratified cross-validation folds. Fold assignment is
randomized; call \code{set.seed()} beforehand for reproducible results.}
\item{min_agreement}{Minimum fraction of samples whose out-of-fold top
predicted identity equals the observed identity. Below this, the function
stops before calibration.}
\item{margin_coverage}{Fraction of correctly predicted samples that must fall
at or below the top-versus-second raw score margin threshold. Samples
predicted incorrectly with an annotation margin (top versus observed score)
above this threshold are treated as confidently misassigned and excluded from
the threshold reference set.}
\item{refit_clean}{Whether to refit the out-of-fold models after removing
confidently misassigned samples from training, and use those second-pass
scores for calibration and all downstream metrics. Individuals with fewer than
two remaining clean samples are dropped (with all their samples) from the
second pass and listed in \code{calibration$dropped_individuals}.
\code{calibration$first_pass_agreement} and \code{calibration$oof_agreement}
report agreement before and after the refit; the minimum-agreement gate
applies to the first pass.}
\item{mismatch_rate}{Fraction of clean samples whose labels are swapped
pairwise in a simulation that gives known-truth agreement for the
risk-coverage curve. \code{NULL} or 0 disables the simulation.}
\item{min_accuracy}{Target agreement among samples that receive a prediction.
The abstention cutoff is the lowest top-versus-second margin at which simulated
true-label agreement (observed agreement when the simulation is off) among
retained samples reaches this value. Abstaining samples get \code{predicted_id =
NA}, are never flagged as suspicious and have status
\code{"no_confident_prediction"}; if the target is never reached, all samples
abstain.}
\item{unknown_mad}{Samples whose top margin is below
\code{median - unknown_mad * MAD} of the log top margin of enrolled samples
abstain and are marked as potential unknown
individuals. No share of samples is forced to be unknown.}
}
\value{
A list with sample-level \code{predictions}, raw \code{identity_scores},
\code{performance}, \code{calibration}, \code{risk_coverage},
\code{mismatch_simulation}, \code{held_out} (count, fraction and IDs of samples held
out as confidently misassigned, from the first pass), \code{feasibility}, and
\code{thresholds} components.
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

Evidence is expressed as raw score margins (top versus second, and top versus
observed label), not probabilities. The probability that an observed annotation
is correct is not reported
because simulated error prevalence does not identify real-data error
prevalence. Internal thresholds and performance are estimated from the same
out-of-fold predictions and should be treated as exploratory.
}
