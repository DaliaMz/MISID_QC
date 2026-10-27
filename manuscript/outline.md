# Sample annotation likelihood manuscript outline

## Abstract

**Background:** [Complete after literature review.]

**Methods:** [Describe the validated workflow and resampling design.]

**Results:** [Insert only results from completed simulation and real-data
analyses.]

**Conclusions:** [Complete after results and limitations are established.]

## 1. Introduction

1. Sample misannotation in high-throughput omics.
2. Limitations of existing annotation quality-control procedures.
3. Repeated-measures molecular data as identity evidence.
4. Need for dataset-specific validation and confidence calibration.
5. Study objectives and accurately qualified methodological contribution.

## 2. Methods

### 2.1 Problem definition

Define observed annotations, enrolled candidate identities, and samples from
individuals absent from the reference set.

### 2.2 Shrinkage LDA identity scoring

Specify the within-class covariance estimator, shrinkage, class means, priors,
class-specific intercepts, and full discriminant scores.

### 2.3 Artificial annotation error generation

Describe correctly annotated controls, random reassignment, pairwise swaps,
and whole-individual holdout for unknown-individual experiments.

### 2.4 Internal resampling and leakage prevention

State the resampling unit, training-only preprocessing, reference exclusions,
and separation of calibration, threshold selection, and evaluation.

### 2.5 Probability calibration

Compare candidate-identity calibration methods and separately assess whether
annotation-correctness probabilities are identifiable under unknown error
prevalence.

### 2.6 Misannotation detection and identity recovery

Define detection outcomes, candidate correction, open-set decisions, and
discrimination and calibration metrics.

### 2.7 Threshold selection

Specify the operating criterion, threshold selection data, and uncertainty
estimation.

### 2.8 Dataset feasibility assessment

Define task-specific criteria for misannotation detection, identity recovery,
and unknown-individual detection, including an insufficient-evidence outcome.

### 2.9 Real-data application

Describe repeated-measures RNA-seq datasets, preprocessing, tissue effects,
and original annotation assessment. Treat unverified flags as candidates, not
confirmed errors.

### 2.10 Software implementation

Describe the R package, reproducibility, input requirements, and limitations.

## 3. Results

1. Shrinkage LDA identity-discriminant scores.
2. Internal experiments and annotation-detection performance.
3. Score calibration and threshold operating characteristics.
4. Dataset feasibility assessment.
5. Suspicious annotations and candidate alternatives in real data.
6. R package workflow.

All result sections remain placeholders until supported by analyses; no
simulation or GTEx results are currently available in this repository.

## 4. Discussion

Discuss within-dataset validation, dataset dependence, annotation-correctness
probability interpretation, open-set limitations, reference contamination,
expression and tissue variation, computational requirements, and extension
to longitudinal metadata.
