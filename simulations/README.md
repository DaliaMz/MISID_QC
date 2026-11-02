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

## GTEx benchmark performance

After GTEx scenarios finish, run `Rscript simulations/Aggregate_gtex_results.R`
to write `results/summary_gtex_performance_by_setting.pdf`. There is one page
per config, task and tissue set (for example tissue / random). Each page has
one panel per metric: fraction of runs that were usable, low-agreement rate,
error rate, and, for usable runs, true sensitivity, specificity,
annotation-margin AUROC and correction accuracy. Within a panel, each line is a
label count / samples-per-label setting (distinct colour, symbol and line type,
nudged sideways so overlapping lines stay visible); the x-axis is the simulated
mislabel rate. Points are means across repetitions, with whiskers showing the
observed range. In an interactive R session the same pages are also drawn on
screen. The status panels show whether the method cannot proceed (low
agreement) or failed to run (error), separately from the metric panels.
