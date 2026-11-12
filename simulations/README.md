# Simulations

Development scripts for testing `check_annotations()`. Not part of the
installed package (excluded via `.Rbuildignore`). Install the package first
(`R CMD INSTALL .`), then run scripts from the repo root.

On the PBS cluster, run `bash simulations/reinstall_package.sh` before
submitting jobs. It activates the same `jupyter` Conda environment as the PBS
runners, reinstalls the package, and prints its installed location. Do not
reinstall while simulation jobs are using the package.

- `simulate.R`: `make_dataset()` builds Gaussian expression data with
  confusable pairs of individuals and random label swaps, and records which
  samples are truly mislabelled.
- `run_sweep.R`: runs `check_annotations()` over a grid of settings with fixed
  seeds and scores the calls against the true labels.
- `results/`: output, ignored by git.

## What the simulations are for

| # | Question | Evidence | Decision it drives |
|---|---|---|---|
| Q1 | How well does it detect mislabels, name the right identity and avoid false alarms? | True sensitivity, specificity, AUROC, correction accuracy against the truth | The operating range to document |
| Q2 | How far are its own estimates from the truth? | `sens_gap`, `spec_gap`, `recovery_gap`, `correction_precision_gap` (reported minus true) against `first_pass_agreement` | Whether the feasibility output can be trusted |
| Q3 | Is the verdict right, and where does it break down? | `verdict_correct`, `recovery_verdict_correct` by first-pass agreement | The threshold for an agreement flag |
| Q4 | Does the first-pass holdout work? | `held_precision`, `held_recall` | Whether `margin_coverage` can stay as it is |
| Q5 | Do the tuning parameters matter? | One-at-a-time sweeps of `mismatch_repeats`, `margin_coverage`, `target_specificity`, `cv_folds` on a subset | Which parameters to fix at a default |
| Q6 | Does it hold on real data? | The GTEx scenarios | Whether the simulated findings transfer |
| Q7 | False positives on clean data? | Mislabel rate 0 (specificity) | The false-alarm rate |

Q2 and Q3 are the main questions. The simulated data give controlled factors
(label signal, mislabel rate, replication) so you can see why something fails;
GTEx repeats Q2, Q3 and Q6 on real expression. The generators also include
unknown individuals; the benchmark runners default to `detect_unknown=FALSE`.
Not covered:
structured (non-random) mislabels or variation in the number of labels,
features or batches in the simulated data.

## Unknown individuals and tissues

By default, the synthetic scenarios and GTEx **individual** scenarios add one
unknown individual with two samples. GTEx **tissue** scenarios add one unknown
tissue class, sampled from an `SMTS` tissue absent from that scenario's enrolled
labels, with two samples from distinct donors per unknown tissue. Unknown samples
receive an observed label from an enrolled class and
have `is_unknown = TRUE` and `is_misid = TRUE`. Each base GTEx scenario also has
a paired version with no unknowns. `misid_rate` controls only the mislabels
among enrolled samples, so a zero mislabel rate is not a completely clean
scenario unless unknowns are disabled.

To change the defaults when sourcing the GTEx generator, set
`n_unknown_individuals`, `n_unknown_tissues` and
`unknown_samples_per_individual` first. Set the relevant count to `0L` to
disable unknowns for that task. The synthetic generator uses
`n_unknown_individuals` and `unknown_samples_per_individual`, and also accepts
these arguments in `simulate_dataset()`. GTEx stops with an error if too few
eligible unknown donors or tissues are available.

For example, from the repository root in an R session:

```r
n_unknown_individuals <- 1L
n_unknown_tissues <- 1L
unknown_samples_per_individual <- 2L
source("simulations/Simulate_full_dataset.R")
sims <- simulate_grid(grid)
```

For GTEx, set the same variables before sourcing
`simulations/Simulate_from_GTEX_grid.R`, then call `simulate_gtex_grid(grid)`.
Sourcing either script does not generate or save the full grid.

With `detect_unknown=TRUE`, both runners report `unknown_sens_true` and
`unknown_spec_true`, comparing `potential_unknown_individual` against
`is_unknown`. Sensitivity is the fraction of assessed unknown samples flagged;
specificity is the fraction of assessed enrolled samples not flagged, including
those with swapped labels. Samples with missing flags are excluded. The
`unknown_n_true`, `unknown_n_known`, `unknown_n_assessed_true` and
`unknown_n_assessed_known` columns show the total and assessed counts; an
absent class or disabled detection gives `NA` for the corresponding rate.
The existing misannotation and recovery metrics still group unknowns together
with enrolled-label swaps.

Existing saved results can be scored without refitting:

```sh
Rscript simulations/Run_across_simulated_data.R combine=TRUE rescore_unknown=TRUE
Rscript simulations/Aggregate_gtex_results.R rescore_unknown=TRUE
```

Add `out=...` for a non-default output directory. Rescoring reads each matching
file in `dirty` or `dirty_gtex` and updates the combined summary, not the
individual row files; it requires saved truth, flags and detection settings.
Regenerate the scenario files and use a new runner output directory (or
`force=TRUE`) when changing this design.

The runners no longer accept `refit_clean`: the annotation models are not
refit, and suspicious calls use the calibrated annotation-margin cutoff even
for held-out samples. New configuration names end in `_cutoffonly` to avoid
skipping older runs. Use a new output directory when comparing results, or
select the new configuration explicitly to avoid mixing historical results.

After the runs finish, `Rscript simulations/Summarize_by_question.R` writes
`summary_by_question.pdf` (one section per summary file, one page per question)
and prints the key numbers. Q1 to Q4 use one config (the most common, or
`cfg=...`); Q5 compares all configs.

## Simulated-dataset benchmark

`Simulate_full_dataset.R` builds 360 base settings as one crossed design with 3
repetitions per cell. Each setting has an unknown and a no-unknown version,
giving 720 scenarios total. Gaussian (324 base settings) is the full factorial of label variance
(0.05, 0.2, 0.5), nuisance variance (none, or discrete 0.3 plus continuous
0.2), samples per label (6, 10, 20) and mislabel rate (0, 0.02, 0.05, 0.1, 0.2,
0.4). RNA-seq and Olink (18 each) check that the findings carry over, at label
variance 0.2 and 10 samples per label, crossed with the same nuisance levels
and rates 0, 0.1, 0.4. The scenario numbering is new, so results from the
earlier grid are not comparable scenario by scenario. `introduce_misid()`
redraws when a label would keep fewer than 2 samples, which the package
requires. Use `subset="rep==1"` for a third of the grid.

`run_sims.sh` submits `Run_across_simulated_data.R` as a Slurm job array (60
tasks by default, 25 at a time). Each task writes one row file per scenario and
config under `<out>/rows`, and finished ones are skipped on resubmission. After
all tasks finish, run
`Rscript simulations/Run_across_simulated_data.R combine=TRUE` to write
`summary_v2.rds` and `summary_v2.csv`. Defaults match the GTEx runner (5-fold
CV, `mismatch_repeats=3`, `min_agreement=0`, new columns for the gaps between
reported and true sensitivity, recovery and correction precision). To evaluate
the other parameters (Q5), submit one-at-a-time sweeps on a subset, for example
`RUN_ARGS="cv_folds=5 mismatch_repeats=1,3,5 subset=rep==1"`. Use a new `out=`
directory or `force=TRUE` after changing the package.
Both runners also accept `detect_unknown=TRUE,FALSE` to run each scenario in
both detection modes; the default is `FALSE`. For PBS, submit with
`RUN_ARGS='cv_folds=5 detect_unknown=TRUE,FALSE' qsub simulations/run_sims.pbs`
or use `run_gtex.pbs` for the GTEx benchmark.
The holdout reference can also be compared with
`holdout_reference=all_alternatives,top`: the default pools margins to every
alternative identity, while `top` uses only the top-versus-second margin.
Both runners use the same selected rule for observed and simulated samples, and
include non-default choices in the config name so their output files do not
collide.

## GTEx benchmark scenarios and settings

`Simulate_from_GTEX_grid.R` builds 150 base settings: tissue (6 or 12 tissues, 10
or 20 samples each) and individual (20 or 40 individuals; 4, 10 or all
samples) tasks with random tissue sets only, at mislabel rates 0, 0.05, 0.1,
0.2 and 0.4, with 3 repetitions. Each base setting has an unknown-class and a
no-unknown version, for up to 300 scenarios. The numbering and seeds come from the original
grid, so scenarios shared with earlier runs keep their numbers and seeds. Set
`keep_sets`, `keep_rates` or `keep_individual_n` before sourcing to change the
subset. The runners default to 5-fold CV and
`mismatch_repeats=3`, so detection sensitivity is cross-fitted, and to
`min_agreement=0` (agreement gate off). `first_pass_agreement` is recorded per
scenario; apply any cutoff to it afterwards. Finished scenarios and configs are
skipped on rerun, so use a new `out=` directory or `force=TRUE` after changing
the package. Besides true sensitivity and specificity, the scored rows include
`sens_gap`, `recovery_gap` and `correction_precision_gap` (reported minus true)
and `detection_source`, to check how well the internal estimates match the
truth.

## GTEx benchmark performance

After GTEx scenarios finish, run `Rscript simulations/Aggregate_gtex_results.R`
to write `results/summary_gtex_performance_by_setting.pdf`. There is one page
per config, task and tissue set (for example tissue / random). Each page has
six panels: low-agreement rate; the per-sample annotation-status distribution
(stacked bars clustered by labels / samples per label, with mislabel rate on
the x-axis); and, for usable runs, true sensitivity, specificity,
annotation-margin AUROC and correction accuracy. Execution-error runs are
excluded from the annotation-status stacks because they have no per-sample
annotations. Execution errors
and usable-run fractions are not plotted; instead, the script prints an error
flag and count, also shown in each report-page header. Within a panel, each
line is a label count / samples-per-label setting (distinct colour, symbol and
line type, nudged sideways so overlapping lines stay visible); the x-axis is
the simulated mislabel rate. Points are means across repetitions, with
whiskers showing the observed range. Proportion axes zoom to the data, with a
minimum displayed range of 0.25. In an interactive R session the same pages
are also drawn on screen.

To estimate how much expression variance is associated with the labels in each
GTEx scenario, run
`Rscript simulations/Quantify_gtex_label_variance.R`. The script uses the same
top-variable-feature selection as the GTEx runner and reports pooled raw and
adjusted R-squared values for both the true labels and the simulated observed
labels. It also reports the expected null R-squared for each label set and the
number of distinct labels. The script writes per-scenario results to
`gtex_label_variance.csv` and `gtex_label_variance.rds`, plus grouped means and
standard deviations to `gtex_label_variance_by_setting.csv` under
`simulations/results` by default. When run interactively, it displays a ggplot
of raw and adjusted variance by mislabel rate, comparing true and observed
labels; colors distinguish settings and lines show means across
repetitions. Interactive plotting requires `ggplot2`. These empirical VST-scale
estimates are analogues of the simulation grid's `var_label`, not exact matches
to its latent generative variance component.
