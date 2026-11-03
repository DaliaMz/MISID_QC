# Simulations

Development scripts for testing `check_annotations()`. Not part of the
installed package (excluded via `.Rbuildignore`). Install the package first
(`R CMD INSTALL .`), then run scripts from the repo root.

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
GTEx repeats Q2, Q3 and Q6 on real expression. Not covered: structured
(non-random) mislabels, unknown individuals (`detect_unknown` is off) and
variation in the number of labels, features or batches in the simulated data.

After the runs finish, `Rscript simulations/Summarize_by_question.R` writes
`summary_by_question.pdf` (one section per summary file, one page per question)
and prints the key numbers. Q1 to Q4 use one config (the most common, or
`cfg=...`); Q5 compares all configs.

## Simulated-dataset benchmark

`Simulate_full_dataset.R` builds 360 scenarios as one crossed design with 3
repetitions per cell. Gaussian (324) is the full factorial of label variance
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

## GTEx benchmark scenarios and settings

`Simulate_from_GTEX_grid.R` builds 150 scenarios: tissue (6 or 12 tissues, 10
or 20 samples each) and individual (20 or 40 individuals; 4, 10 or all
samples) tasks with random tissue sets only, at mislabel rates 0, 0.05, 0.1,
0.2 and 0.4, with 3 repetitions. The numbering and seeds come from the original
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
