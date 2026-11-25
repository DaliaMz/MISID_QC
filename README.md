# MISID_QC

An R package for detecting and evaluating sample misannotations in
repeated-measures RNA-seq data using shrinkage LDA, internal validation,
probability calibration, and dataset-specific feasibility assessment.

## Initial implementation

The package accepts transformed expression data as a numeric matrix with
features in rows and sample identifiers in columns. Metadata must include one
row per sample and the observed individual identifier. Sample identifiers are
matched exactly; each individual needs at least two samples.

```r
results <- SampleAnnotationLikelihood::check_annotations(
  expr = expression_matrix,
  metadata = sample_metadata,
  sample_col = "sample",
  id_col = "individual"
)

results$predictions
results$identity_scores
results$performance
results$feasibility
plot_risk_coverage(results)
```

The initial implementation uses leave-one-sample-out shrinkage LDA and raw
score margins (top versus second, and top versus observed label). First-pass agreement
below 75% triggers a warning. The run stops if first-pass agreement is below
`min_agreement` (default 60%); set it to `0` to disable the gate when studying
low-agreement simulations. The former `stop_on_low_agreement` argument has
been removed. Confidently
misassigned samples are held out of the calibration reference set, but the
models are not refit. A sample is confidently misassigned when its prediction
disagrees with its label and its annotation margin (top versus observed label)
exceeds the `margin_coverage` quantile (default 0.025) of reference margins
from correctly predicted samples. `holdout_reference = "all_alternatives"`
(default) pools top-minus-each-alternative margins; `"top"` uses only
top-minus-second margins. The label-swap simulation applies the same rule.
Set `mismatch_rate = "held_out"` to use the first-pass held-out fraction as
the swap rate instead of the default 0.1. This is not an estimate of the true
error rate. Swaps are made in whole pairs, with at least one pair even when
the held-out fraction is zero; the actual fraction can differ and is recorded
as `calibration$mismatch_simulation$actual_rate`.
The former `refit_clean` argument has been removed; omit it from existing calls.
For compatibility, `calibration$dropped_individuals` is empty and
`predictions$excluded_from_refit` is always `FALSE`.
Flagging uses a margin threshold, and a label-swap simulation
supplies known-truth agreement for `plot_risk_coverage()`. Margins
are **not** annotation-correctness probabilities: no absolute probability
that a sample is correctly annotated is reported because the true error
prevalence is unknown. Internal artificial relabeling and leave-one-individual-
out experiments provide exploratory, dataset-specific performance summaries.
`results$unknown_simulation` contains the held-out-individual sample records,
including their original and randomly assigned labels, predictions, margins,
nearest-class distances, and unknown flags. The corresponding raw LDA scores
are in `results$unknown_simulation_scores`, in the same row order; held-out
classes have `NA` scores. Both are `NULL` when unknown detection is disabled
or no held-out-individual samples can be simulated. The scores come from
different models and are intended for class comparisons within each sample.
Relabeling metrics: detection sensitivity comes from the simulated pairwise
swaps and, with `mismatch_repeats >= 2`, is cross-fitted (each repeat is
scored with a cutoff tuned on simulated positives and unchanged samples from the
other repeats, plus the reference negatives; see
`performance$misannotation_detection$source`).
The final cutoff maximizes simulated-swap sensitivity while independently
requiring both reference and unchanged simulated specificity to meet
`target_specificity`. If simulations are unavailable, it falls back to
calibration against the reference and all wrong-label alternatives.
The detection verdict uses simulated sensitivity (target 80%) and specificity
(target `target_specificity`). It is `appears_feasible` when both meet their
targets, `unreliable` when either misses by more than 10 percentage points, and
`below_target` for smaller shortfalls (including exactly 10 points). The rates
are cross-fitted by repeat when possible, and in-sample with one repeat; without
simulated positives or unchanged samples the verdict is insufficient evidence.
This verdict does not change sample flags or candidate corrections.
`performance$misannotation_detection$simulated_specificity` is the verdict
estimate, while `specificity` is the original reference specificity.
`performance$estimated_detection_errors` uses the observed suspicious-flag
fraction and these simulated rates to estimate the misannotation fraction:
`(flag fraction - (1 - specificity)) / (sensitivity + specificity - 1)`.
It projects the numbers of false positives and missed misannotations at that
estimate and at +/- 2.5 percentage points, bounded to 0--100%. These are
expected counts under each assumed prevalence, not confidence intervals or
guaranteed counts; the alternative assumptions need not reproduce the observed
flag count. Estimates are unavailable when simulated rates are missing, do not
distinguish errors, or imply a prevalence outside 0--100%. When sensitivity
and specificity barely distinguish errors, the estimate can be highly unstable.
Cross-fitted rates use cutoffs from other repeats, which can differ from the
final cutoff used for the observed flags; a single repeat uses in-sample rates.
Swaps may not represent real misannotations, especially unknown individuals.
The report includes the estimates; `progress = TRUE` also displays them.
The report's cutoff curves show in-sample rates; cross-fitted specificity
can be below the target even when the calibration rates meet it.
The every-wrong-label and
all-pairs swap metrics are kept as diagnostics
(`performance$all_alternatives_detection`, `performance$pairwise_swap_detection`);
precision reflects those artificial scenarios, not a population error
prevalence. Recovery on simulated misannotated samples and the precision and
recall of proposed corrections are in `performance$recovery_accuracy_misannotated`
and `performance$correction_validation`.
Other summaries use the same out-of-fold predictions for threshold selection and
performance, so they are not independent estimates and should not be used as
confirmed error rates. `annotation_status` is descriptive (flag, confidence,
agreement) and is not overwritten by feasibility, which is reported separately
in `results$feasibility`. All samples, including held-out samples, are
`suspicious` only when their annotation margin exceeds the calibrated cutoff.
Every input sample has a row in `predictions`.
For correction calibration, simulations include swapped samples plus unchanged
samples predicted as a different identity. The lowest top-margin cutoff reaching
`min_accuracy` (default 90%) on that filtered set is recorded in
`thresholds$correction_margin`; the retained count, coverage and in-sample accuracy are in
`calibration$correction_margin`. Swaps predicted as their simulated label are
included, so this accuracy is not the precision of proposed corrections.
Recovery appears feasible if such a cutoff exists and there are at least 20
input samples; without swaps it is insufficient evidence. `predicted_id_clean`
returns agreeing predictions and disagreements at or above the correction cutoff,
even for unknown-flagged samples. Candidate corrections additionally require a
suspicious flag and must not be flagged as unknown. Without a usable cutoff,
disagreements have no clean prediction. The former general abstention threshold
and its outputs have been removed. Simulated correction validation does not
apply the default distance-based unknown gate or the feasibility verdict, so it
does not exactly reproduce final candidates. There is no
minimum retained simulation count: even one correct simulated row can meet the
target. The cutoff and accuracy use the same simulations, whose rows may be
dependent, so the estimate can be optimistic.

To inspect the individual stages and plot the margins and cutoffs, load the
development version and request diagnostics:

```r
devtools::load_all()
result <- check_annotations(
  expression_matrix, sample_metadata, "sample", "individual",
  diagnostics = TRUE
)
source("scripts/plot_check_annotation_gates.R")
diagnostic_tables <- plot_check_annotation_gates(
  result, file = "annotation-gates.pdf"
)

diagnostic_tables$first_pass
diagnostic_tables$final_pass
diagnostic_tables$annotation_threshold
```

For a printable six-page overview, holdout and calibration plots, recovery
coverage, unknown detection, and sample score distributions, write a PDF with:

```r
result <- check_annotations(expression_matrix, sample_metadata, "sample", "individual",
                            diagnostics = TRUE)
annotation_report(result, "annotation-report.pdf")
```

The report uses ggplot2, which must be installed before installing this package.
The misannotation page shows detection targets, in-sample and cross-simulation
rates, and estimated detection errors together. When cross-fitting is
unavailable, it labels the verdict rates as in-sample or unavailable.
Its recovery page compares all simulated samples with the filtered
correction-calibration samples. Coverage is relative to each population, so
equal coverage does not necessarily mean the same top-margin cutoff. The
marked correction cutoff is selected from the filtered samples.
Margins and unknown scores are plotted on their original, untransformed scales.
The report marks detailed calibration plots unavailable if diagnostics were
not retained. Its simulated validation estimates are exploratory, not true
error rates or sample-level probabilities.

The PDF plots boxplots of first-pass top and annotation margins for samples
whose prediction matches vs. differs from the observed label (agreement is in
the title), the first-pass exclusion margins,
the annotation-threshold populations, correction accuracy/coverage, unknown
individual margins, and feasibility gates. The first-pass table is kept
separate from the final decisions. Diagnostics are opt-in because they retain
additional intermediate score and margin data.

For one saved GTEx simulation scenario, run from the repository root:

```sh
Rscript scripts/run_gtex_scenario_diagnostics.R \
  scenario=12 \
  scenarios_file=simulations/results/gtex_scenarios.rds \
  expr=/path/to/GTEx_v11_VST.rds
```

The wrapper installs the current checkout into the default R library with
`devtools::install_local(force = TRUE)` before loading it. It then uses the scenario's
observed labels, selects the top variable features within that scenario
(default 2,000, matching the GTEx runner), disables unknown-individual
detection, and writes a gate PDF and an RDS bundle to
`simulations/results/gtex_diagnostics/`. Pass `out=...`, `n_features=...`,
`cv_folds=...` (or `cv_folds=loo`), and the other supported diagnostic settings
to change defaults. It also writes `scenario_<id>_first_pass_truth.pdf`, built
by `plot_first_pass_truth()` in `scripts/plot_first_pass_truth.R`: first-page
top/annotation margin distributions for correctly predicted versus
disagreeing samples, then first-pass top and annotation margins coloured by the true `is_misid` status, the
exclusion cutoff, which samples were held out, and holdout precision/recall.
The simulation seed is used by default; `seed=...` can
override it. The RDS includes the full result, diagnostic tables, scenario
truth and parameters.

This version does not normalize raw counts or select features. Supply
appropriately transformed expression values. For known tissue or batch
covariates, pass metadata column names with `covariate_col`:

```r
check_annotations(
  expr, metadata, "sample", "individual",
  covariate_col = "tissue"
)
```

`covariate_col` may be a character vector for multiple covariates. Categorical
covariates are fitted as fixed effects and numeric covariates as linear effects;
expression is residualized before individual, mislabel, and unknown validation.
Covariates are not modelled using individual identities or unknown-truth labels.
