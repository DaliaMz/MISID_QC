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
