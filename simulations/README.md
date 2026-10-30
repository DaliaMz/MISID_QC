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
