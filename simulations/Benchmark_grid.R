library(data.table)

# Shared benchmark dimensions for synthetic, GTEx tissue, and GTEx individual
# scenarios. Override these variables before sourcing a scenario grid.
if (!exists("benchmark_n_groups", inherits = FALSE)) {
  benchmark_n_groups <- c(3L, 5L, 10L, 20L)
}
if (!exists("benchmark_n_per_group", inherits = FALSE)) {
  benchmark_n_per_group <- c(3L, 6L, 10L, 20L, 30L)
}
if (!exists("benchmark_misid_rates", inherits = FALSE)) {
  benchmark_misid_rates <- c(0, 0.05, 0.1, 0.2, 0.4)
}
if (!exists("benchmark_reps", inherits = FALSE)) {
  benchmark_reps <- 1:3
}

if (!is.numeric(benchmark_n_groups) ||
    length(benchmark_n_groups) < 1L ||
    any(!is.finite(benchmark_n_groups)) ||
    any(benchmark_n_groups < 3L) ||
    any(benchmark_n_groups != floor(benchmark_n_groups))) {
  stop("benchmark_n_groups must contain integers of at least 3.")
}
if (!is.numeric(benchmark_n_per_group) ||
    length(benchmark_n_per_group) < 1L ||
    any(!is.finite(benchmark_n_per_group)) ||
    any(benchmark_n_per_group < 3L) ||
    any(benchmark_n_per_group != floor(benchmark_n_per_group))) {
  stop("benchmark_n_per_group must contain integers of at least 3.")
}
if (!is.numeric(benchmark_misid_rates) ||
    length(benchmark_misid_rates) < 1L ||
    any(!is.finite(benchmark_misid_rates)) ||
    any(benchmark_misid_rates < 0 | benchmark_misid_rates > 1)) {
  stop("benchmark_misid_rates must contain finite values from 0 to 1.")
}
if (!is.numeric(benchmark_reps) ||
    length(benchmark_reps) < 1L ||
    any(!is.finite(benchmark_reps)) ||
    any(benchmark_reps < 1L) ||
    any(benchmark_reps != floor(benchmark_reps))) {
  stop("benchmark_reps must contain positive integers.")
}

# Cross the shared dimensions and expose the size of the dimension product.
benchmark_grid <- function(n_groups = benchmark_n_groups,
                           n_per_group = benchmark_n_per_group,
                           misid_rates = benchmark_misid_rates,
                           reps = benchmark_reps) {
  CJ(n_groups = n_groups, n_per_group = n_per_group,
     misid_rate = misid_rates, rep = reps)
}
