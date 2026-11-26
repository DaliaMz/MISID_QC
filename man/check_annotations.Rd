\name{check_annotations}
\alias{check_annotations}
\title{Assess repeated-measures RNA-seq sample annotations}
\usage{
check_annotations(expr, metadata, sample_col, id_col, shrinkage = 0.1,
  target_specificity = 0.95, assess_feasibility = TRUE,
  progress = interactive(), cv_folds = NULL,
  min_agreement = 0.6, margin_coverage = 0.025,
  mismatch_rate = 0.1, min_accuracy = 0.9, detect_unknown = TRUE,
  grouped_unknown_validation = FALSE, unknown_stray_samples = 2L,
  mismatch_repeats = 3L, covariate_col = NULL,
  min_correction_samples = 20L, correction_repeat_tolerance = 0.05,
  diagnostics = FALSE)
}
\arguments{
\item{expr}{Numeric expression matrix, features by samples, with unique sample names in \code{colnames(expr)}.}
\item{metadata}{Data frame containing sample and observed identity columns.}
\item{sample_col}{Name of the sample identifier column.}
\item{id_col}{Name of the observed individual identifier column.}
\item{shrinkage}{Shrinkage fraction for the diagonal covariance target.}
\item{target_specificity}{Minimum specificity for reference and unchanged simulated samples.}
\item{assess_feasibility}{Whether to return task-specific feasibility assessments.}
\item{progress}{Whether to report progress updates.}
\item{cv_folds}{Leave-one-sample-out when \code{NULL}, otherwise a stratified fold count.}
\item{min_agreement}{Minimum out-of-fold agreement required before calibration.}
\item{margin_coverage}{Quantile used for the first-pass holdout cutoff. The reference pools top-minus-alternative margins.}
\item{mismatch_rate}{Positive numeric swap fraction or \code{"held_out"} to use the first-pass held-out fraction.}
\item{min_accuracy}{Target agreement for selecting the correction cutoff.}
\item{min_correction_samples}{Minimum eligible simulated rows retained at the correction cutoff before candidate corrections are enabled.}
\item{correction_repeat_tolerance}{Maximum permitted shortfall of the lowest repeat-level correction accuracy below \code{min_accuracy}.}
\item{detect_unknown}{Whether to run nearest-class unknown-individual validation.}
\item{grouped_unknown_validation}{Whether held-out individuals are processed in groups instead of individually.}
\item{unknown_stray_samples}{Maximum number of samples retained per simulated unknown individual.}
\item{mismatch_repeats}{Number of independent pairwise-swap simulation repeats.}
\item{covariate_col}{Optional metadata covariate column or columns to residualize.}
\item{diagnostics}{Whether to retain intermediate gate and threshold data.}
}
\value{
A list containing sample predictions, raw identity scores, performance,
calibration, covariate adjustment, mismatch simulation, unknown-individual
simulation, correction gate, held-out sample information, feasibility, thresholds, and
optional diagnostics.
}
\description{
Uses leave-one-sample-out or stratified cross-validation shrinkage LDA
scores and compulsory pairwise-swap simulations to describe dataset-specific
annotation quality. Unknown individuals are evaluated with nearest-class
Mahalanobis distance only.
}
