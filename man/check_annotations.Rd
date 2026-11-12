\name{check_annotations}
\alias{check_annotations}
\title{Assess repeated-measures RNA-seq sample annotations}
\usage{
check_annotations(expr, metadata, sample_col, id_col, shrinkage = 0.1,
  target_specificity = 0.95, assess_feasibility = TRUE,
  progress = interactive(), cv_folds = NULL,
  min_agreement = 0.7, margin_coverage = 0.025,
  holdout_reference = c("all_alternatives", "top"),
  mismatch_rate = 0.1, min_accuracy = 0.9, detect_unknown = TRUE,
  unknown_threshold_method = c("simulated", "mad"),
  unknown_mad = 3,
  grouped_unknown_validation = FALSE,
  unknown_stray_samples = 2L, mismatch_repeats = 3L,
  stop_on_low_agreement = TRUE, diagnostics = FALSE)
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
warns and returns every sample as
\code{"insufficient_evidence"} (no clean predictions, flags or corrections),
keeping raw scores, margins and the remaining outputs.}
\item{margin_coverage}{Quantile used to set the first-pass holdout cutoff.
The reference values are chosen by \code{holdout_reference}.}
\item{holdout_reference}{Whether to use the top-versus-second margin
(\code{"top"}) or pool the top-versus-each-alternative margins
(\code{"all_alternatives"}) from correctly predicted samples. Samples
predicted incorrectly with an annotation margin (top versus observed score)
above this threshold are treated as confidently misassigned and excluded from
the threshold reference set. The same choice is applied to simulated mismatches.
The original out-of-fold models are not refit after holdout.}
\item{mismatch_rate}{Fraction of clean samples whose labels are swapped
pairwise in a simulation that gives known-truth agreement for the
risk-coverage curve. \code{NULL} or 0 disables the simulation.}
\item{mismatch_repeats}{Number of independent label-swap simulations whose results
are pooled to tune the annotation cutoff. With at least two repeats, detection
sensitivity is cross-fitted: each repeat is scored with a cutoff tuned on the
other repeats. More repeats give a less extreme cutoff and a steadier
sensitivity estimate at the cost of roughly proportional simulation time. With 1,
sensitivity is in-sample.}
\item{stop_on_low_agreement}{If \code{TRUE} (default), the run stops with an
error when first-pass out-of-fold agreement is below 0.6, because the
calibration is then unreliable. Set to \code{FALSE} to proceed anyway. A
warning is always issued when first-pass agreement is below 0.75.}
\item{diagnostics}{Whether to include intermediate gate data in the result.
This is intended for understanding and plotting the validation steps; it can
add substantially to the result size.}
\item{min_accuracy}{Target agreement among samples that receive a prediction.
The abstention cutoff is the lowest top-versus-second margin at which simulated
true-label agreement (observed agreement when the simulation is off) among
retained samples reaches this value. Abstaining samples get \code{predicted_id_clean = NA}
(\code{predicted_id} keeps the raw top candidate). They can still be flagged as
suspicious when the annotation margin is large, but never receive a candidate
correction; if the target is never reached, all samples
abstain.}
\item{detect_unknown}{If \code{FALSE}, unknown-individual detection is switched
off: no stranger validation is run, no sample is marked as a potential
unknown individual, and abstention uses only the abstention cutoff.
\code{feasibility$unknown_individual_detection} is \code{"not_assessed"}.}

\item{unknown_threshold_method}{How the unknown-individual cutoff on the top
statistic is chosen. \code{"simulated"} (default) uses the Mahalanobis
distance (under the shrinkage covariance) from a sample to its nearest
enrolled individual's mean, with the cutoff chosen so at most
\code{1 - target_specificity} of enrolled reference samples exceed it while
catching the most simulated strangers; with no simulated strangers it is the
\code{target_specificity} quantile of the enrolled distances. \code{"mad"}
restores the earlier rule, a cutoff on the top-versus-second margin of
\code{median - unknown_mad * MAD} of the log top margin of enrolled samples,
which uses no simulated strangers. Samples beyond the cutoff abstain and are
marked as potential unknown individuals. No share
of samples is forced to be unknown.}
\item{unknown_mad}{Multiplier used only when
\code{unknown_threshold_method = "mad"}.}
\item{unknown_stray_samples}{Unknown-individual validation pretends each enrolled
individual is a stranger. Up to this many of their samples (never all of them)
are kept in the data relabelled as random other individuals, as a stranger's
samples would be in a real cohort, and the rest are dropped. Each kept sample
is scored by a model that never saw it, giving a margin for an unknown person.}
\item{grouped_unknown_validation}{By default the stranger simulation is repeated
for each individual separately. If \code{TRUE} and \code{cv_folds} is set, individuals are
processed in \code{cv_folds} groups instead, so far fewer models are fit, at the cost
of slightly optimistic (too large) unknown margins. It only affects the
reported unknown sensitivity and feasibility. With the default \code{unknown_threshold_method = "simulated"} it also affects the unknown cutoff.}
}
\value{
A list with sample-level \code{predictions}, raw \code{identity_scores},
\code{performance}, \code{calibration}, \code{risk_coverage},
\code{mismatch_simulation}, \code{held_out} (count, fraction and IDs of samples held
out as confidently misassigned, from the first pass), \code{feasibility}, and
\code{thresholds} components. \code{predictions$potential_unknown_individual} is a
separate flag, not an \code{annotation_status} value. Predictions also include
the nearest-class distance and the unknown-detection score and threshold.
\code{predictions$annotation_status} is descriptive (flag, abstention,
agreement) and does not depend on \code{feasibility};
\code{"insufficient_evidence"} is used when the minimum-agreement gate fails.
All samples, including held-out samples, are flagged as suspicious only when
their annotation margin exceeds the calibrated threshold. Every input sample has a row.
\code{performance$misannotation_detection} is validated on simulated swaps (see its
\code{source} column), \code{performance$all_alternatives_detection} keeps the
every-wrong-label diagnostic, and \code{performance$recovery_accuracy_misannotated}
and \code{performance$correction_validation} report recovery on simulated
misannotated samples and the precision and recall of proposed corrections.
Candidate corrections require identity recovery to appear feasible. If
\code{diagnostics = TRUE}, the result also
contains first-pass and final-pass sample margins, threshold populations,
abstention-curve data, and gate inputs under \code{diagnostics}.
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
