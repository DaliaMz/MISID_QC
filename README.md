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
score margins (top versus second, and top versus observed label). A first-pass
agreement gate (`min_agreement`, default 70%) must be met; confidently
misassigned samples are held out of the reference set and the models are refit
on clean samples. Flagging uses a margin threshold, and a label-swap simulation
supplies known-truth agreement for `plot_risk_coverage()`. Margins
are **not** annotation-correctness probabilities: no absolute probability
that a sample is correctly annotated is reported because the true error
prevalence is unknown. Internal artificial relabeling and leave-one-individual-
out experiments provide exploratory, dataset-specific performance summaries.
Relabeling metrics: detection sensitivity comes from the simulated pairwise
swaps and, with `mismatch_repeats >= 2`, is cross-fitted (each repeat is
scored with a cutoff tuned on the other repeats; see
`performance$misannotation_detection$source`). The every-wrong-label and
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
in `results$feasibility`. Samples held out in the first pass are always
`suspicious`, and every input sample has a row in `predictions`.
`insufficient_evidence` appears when the agreement gate fails and for unflagged
samples of individuals dropped at the refit (`excluded_from_refit`). Candidate
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

The PDF plots boxplots of first-pass top and annotation margins for samples
whose prediction matches vs. differs from the observed label (agreement is in
the title), the first-pass exclusion margins,
the annotation-threshold populations, abstention accuracy/coverage, unknown
individual margins, and feasibility gates. The first-pass table is kept
separate from the final/refit table so a post-refit disagreement is not
mistaken for a first-pass holdout. Diagnostics are opt-in because they retain
additional intermediate score and margin data.

For one saved GTEx simulation scenario, run from the repository root:

```sh
Rscript scripts/run_gtex_scenario_diagnostics.R \
  scenario=12 \
  scenarios_file=simulations/results/gtex_scenarios.rds \
  expr=/path/to/GTEx_v11_VST.rds
```

The wrapper installs the current checkout into
`/scratch/project/palpant_scratch/DMizikovsky/R_packages` with
`devtools::install_local(force = TRUE)` before loading it. Set
`MISID_QC_R_LIB` to use a different R library. It then uses the scenario's
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
