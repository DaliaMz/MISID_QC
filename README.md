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
Relabeling metrics include all wrong candidate identities and pairwise swaps;
precision therefore reflects those artificial scenarios, not a population error
prevalence.
Threshold selection and performance use the same out-of-fold predictions, so
these summaries are not independent estimates and should not be used as
confirmed error rates. Candidate corrections are emitted only when both
misannotation detection and identity recovery appear feasible.

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

The PDF plots the first-pass agreement gate, the first-pass exclusion margins,
the annotation-threshold populations, abstention accuracy/coverage, unknown
individual margins, and feasibility gates. The first-pass table is kept
separate from the final/refit table so a post-refit disagreement is not
mistaken for a first-pass holdout. Diagnostics are opt-in because they retain
additional intermediate score and margin data.

For one saved GTEx simulation scenario, install the current package with
`R CMD INSTALL .` and run from the repository root:

```sh
Rscript scripts/run_gtex_scenario_diagnostics.R \
  scenario=12 \
  scenarios_file=simulations/results/gtex_scenarios.rds \
  expr=/path/to/GTEx_v11_VST.rds
```

The wrapper uses the scenario's observed labels, selects the top variable
features within that scenario (default 2,000, matching the GTEx runner),
disables unknown-individual detection, and writes a gate PDF and an RDS bundle
to `simulations/results/gtex_diagnostics/`. Pass `out=...`, `n_features=...`,
`cv_folds=...` (or `cv_folds=loo`), and the other supported diagnostic settings
to change defaults. The simulation seed is used by default; `seed=...` can
override it. The RDS includes the full result, diagnostic tables, scenario
truth and parameters.

This version does not normalize raw counts, select features, or account for
tissue/batch covariates. Supply appropriately transformed expression values
and consider tissue-related expression variation in the interpretation.
