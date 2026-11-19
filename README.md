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
Relabeling metrics: detection sensitivity comes from the simulated pairwise
swaps and, with `mismatch_repeats >= 2`, is cross-fitted (each repeat is
scored with a cutoff tuned on simulated positives and unchanged samples from the
other repeats, plus the reference negatives; see
`performance$misannotation_detection$source`).
The final cutoff maximizes simulated-swap sensitivity while independently
requiring both reference and unchanged simulated specificity to meet
`target_specificity`. If simulations are unavailable, it falls back to
calibration against the reference and all wrong-label alternatives.
The detection verdict also requires simulated specificity among unchanged
simulation samples to reach `target_specificity`. This specificity is
cross-fitted by repeat when possible, and in-sample with one repeat; without
simulated positives or unchanged samples the verdict is insufficient evidence.
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
confirmed error rates. `annotation_status` is descriptive (flag, abstention,
agreement) and is not overwritten by feasibility, which is reported separately
in `results$feasibility`. All samples, including held-out samples, are
`suspicious` only when their annotation margin exceeds the calibrated cutoff.
Every input sample has a row in `predictions`.
Candidate
corrections are emitted only when identity
recovery on simulated misannotated samples appears feasible.

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
Margins and unknown scores are plotted on their original, untransformed scales.
The report marks detailed calibration plots unavailable if diagnostics were
not retained. Its simulated validation estimates are exploratory, not true
error rates or sample-level probabilities.

The PDF plots boxplots of first-pass top and annotation margins for samples
whose prediction matches vs. differs from the observed label (agreement is in
the title), the first-pass exclusion margins,
the annotation-threshold populations, abstention accuracy/coverage, unknown
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

This version does not normalize raw counts, select features, or account for
tissue/batch covariates. Supply appropriately transformed expression values
and consider tissue-related expression variation in the interpretation.
