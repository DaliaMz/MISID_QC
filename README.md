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
results$performance
results$feasibility
```

The initial implementation uses leave-one-sample-out shrinkage LDA and
temperature-scaled candidate-identity probabilities. Candidate probabilities
are **not** annotation-correctness probabilities: no absolute probability
that a sample is correctly annotated is reported because the true error
prevalence is unknown. Internal artificial relabeling and leave-one-individual-
out experiments provide exploratory, dataset-specific performance summaries.
Threshold selection and performance use the same out-of-fold predictions, so
these summaries are not independent estimates and should not be used as
confirmed error rates. Candidate corrections are emitted only when both
misannotation detection and identity recovery appear feasible.

This version does not normalize raw counts, select features, or account for
tissue/batch covariates. Supply appropriately transformed expression values
and consider tissue-related expression variation in the interpretation.
