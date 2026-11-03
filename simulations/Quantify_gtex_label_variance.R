custom_lib <- "/scratch/project/palpant_scratch/DMizikovsky/R_packages"
.libPaths(c(custom_lib, .libPaths()))

# Estimate the fraction of GTEx expression variance associated with the
# true and simulated observed labels in every GTEx scenario.
#
# Usage:
#   Rscript simulations/Quantify_gtex_label_variance.R \
#     [scenarios_file=simulations/results/gtex_scenarios.rds] \
#     [expr=/path/to/GTEx_v11_VST.rds] [out=simulations/results] [n_features=2000]
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
opt <- list(
  scenarios_file = "simulations/results/gtex_scenarios.rds",
  expr = "/scratch/project_mnt/S0138/DMizikovsky/misid/GTEX/GTEx_v11_VST.rds",
  out = "simulations/results",
  n_features = "2000"
)
for (arg in args[grepl("=", args)]) {
  key_value <- strsplit(arg, "=", fixed = TRUE)[[1]]
  opt[[key_value[1]]] <- paste(key_value[-1], collapse = "=")
}
n_features <- as.integer(opt$n_features)
if (length(n_features) != 1L || is.na(n_features) || n_features < 1L) {
  stop("n_features must be a positive integer.")
}

scenarios <- readRDS(opt$scenarios_file)
if (!length(scenarios)) stop("No scenarios found in ", opt$scenarios_file)
expr <- readRDS(opt$expr)
if (is.null(dim(expr)) || length(dim(expr)) != 2L || is.null(colnames(expr))) {
  stop("Expression input must be a genes-by-samples matrix with sample names.")
}

validate_labels <- function(label, n_samples = length(label)) {
  label <- as.character(label)

  if (!length(label) || length(label) != n_samples || anyNA(label)) {
    stop("Labels must contain one non-missing value per sample.")
  }

  label
}

pooled_label_variance <- function(x, label) {
  label <- validate_labels(label, ncol(x))
  n <- ncol(x)
  k <- length(unique(label))

  grand_mean <- rowMeans(x)
  total_ss <- sum(rowSums((x - grand_mean)^2))

  if (!is.finite(total_ss) || total_ss <= 0) {
    return(c(raw = NA_real_, adjusted = NA_real_))
  }

  between_ss <- 0

  for (indices in split(seq_along(label), label)) {
    group_mean <- rowMeans(x[, indices, drop = FALSE])

    between_ss <- between_ss +
      length(indices) * sum((group_mean - grand_mean)^2)
  }

  r2 <- between_ss / total_ss

  adjusted_r2 <- if (n > k) {
    1 - (1 - r2) * (n - 1) / (n - k)
  } else {
    NA_real_
  }

  c(raw = r2, adjusted = adjusted_r2)
}

expected_null_r2 <- function(label) {
  label <- validate_labels(label)
  n <- length(label)
  k <- length(unique(label))

  if (n <= 1L) return(NA_real_)

  (k - 1) / (n - 1)
}



rows <- vector("list", length(scenarios))
for (i in seq_along(scenarios)) {
  cat(sprintf("Quantifying scenario %d/%d\n", i, length(scenarios)))
  scenario <- scenarios[[i]]
  if (is.null(scenario$meta$sample) || is.null(scenario$meta$label_true) ||
      is.null(scenario$meta$label_observed)) {
    stop("Scenario ", i, " must contain sample, label_true, and label_observed.")
  }
  samples <- as.character(scenario$meta$sample)
  label_true <- validate_labels(scenario$meta$label_true, length(samples))
  label_observed <- validate_labels(scenario$meta$label_observed, length(samples))
  sample_index <- match(samples, colnames(expr))
  if (anyNA(sample_index)) {
    stop("Expression matrix is missing samples for scenario ", i, ".")
  }

  x <- expr[, sample_index, drop = FALSE]
  feature_variance <- apply(x, 1, stats::var)
  selected <- order(feature_variance, decreasing = TRUE, na.last = NA)
  selected <- head(selected, n_features)
  if (!length(selected)) stop("No finite-variance features in scenario ", i, ".")
  x <- x[selected, , drop = FALSE]
  if (any(!is.finite(x))) {
    stop("Expression contains non-finite values in scenario ", i, ".")
  }

  true_var <- pooled_label_variance(x, label_true)
  observed_var <- pooled_label_variance(x, label_observed)
  params <- as.data.table(scenario$params)
  rows[[i]] <- cbind(
    params,
    data.table(
      n_samples = length(samples),
      n_features = nrow(x),
      n_labels_true = uniqueN(label_true),
      n_labels_observed = uniqueN(label_observed),
      label_variance_true = true_var["raw"],
      label_variance_true_adjusted = true_var["adjusted"],
      label_variance_observed = observed_var["raw"],
      label_variance_observed_adjusted = observed_var["adjusted"],
      expected_null_r2_true = expected_null_r2(label_true),
      expected_null_r2_observed = expected_null_r2(label_observed)
    )
  )
}

results <- rbindlist(rows, fill = TRUE)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(opt$out, "gtex_label_variance.rds"))
fwrite(results, file.path(opt$out, "gtex_label_variance.csv"))

group_cols <- intersect(
  c("task", "set", "n_labels", "n_per_label", "misid_rate"),
  names(results)
)
if (length(group_cols)) {
  summary <- results[, .(
    scenarios = .N,
    label_variance_true_mean = mean(label_variance_true, na.rm = TRUE),
    label_variance_true_sd = stats::sd(label_variance_true, na.rm = TRUE),
    label_variance_true_adjusted_mean =
      mean(label_variance_true_adjusted, na.rm = TRUE),
    label_variance_true_adjusted_sd =
      stats::sd(label_variance_true_adjusted, na.rm = TRUE),
    label_variance_observed_mean = mean(label_variance_observed, na.rm = TRUE),
    label_variance_observed_sd = stats::sd(label_variance_observed, na.rm = TRUE),
    label_variance_observed_adjusted_mean =
      mean(label_variance_observed_adjusted, na.rm = TRUE),
    label_variance_observed_adjusted_sd =
      stats::sd(label_variance_observed_adjusted, na.rm = TRUE),
    expected_null_r2_true_mean = mean(expected_null_r2_true, na.rm = TRUE),
    expected_null_r2_observed_mean =
      mean(expected_null_r2_observed, na.rm = TRUE)
  ), by = group_cols]
  setorderv(summary, group_cols)
  fwrite(summary, file.path(opt$out, "gtex_label_variance_by_setting.csv"))
  print(summary)
}
cat(nrow(results), "scenarios quantified; results written to", opt$out, "\n")

if (interactive()) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 to display the interactive variance plot.")
  }

  plot_groups <- c("task", "set", "n_labels", "n_per_label", "misid_rate")
  plot_data <- melt(
    results,
    id.vars = plot_groups,
    measure.vars = c(
      "label_variance_true", "label_variance_observed",
      "label_variance_true_adjusted", "label_variance_observed_adjusted"
    ),
    variable.name = "metric",
    value.name = "r_squared"
  )
  plot_data[, `:=`(
    label_type = ifelse(grepl("_true", metric), "True labels", "Observed labels"),
    metric_type = ifelse(grepl("_adjusted$", metric), "Adjusted", "Raw"),
    setting = sprintf("%s labels / %s samples per label", n_labels, n_per_label)
  )]
  plot_summary <- plot_data[, .(
    mean = if (any(is.finite(r_squared))) mean(r_squared, na.rm = TRUE) else NA_real_
  ), by = c("task", "set", "setting", "misid_rate", "metric_type", "label_type")]

  variance_plot <- ggplot2::ggplot(
    plot_summary,
    ggplot2::aes(
      x = misid_rate, y = mean, color = setting, linetype = label_type,
      group = interaction(setting, label_type)
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::facet_grid(task + set ~ metric_type, scales = "free_y") +
    ggplot2::labs(
      title = "GTEx label-associated expression variance",
      x = "Mislabel rate",
      y = "Pooled R-squared",
      color = "Setting",
      linetype = "Labels"
    ) +
    ggplot2::theme_bw()
  print(variance_plot)
}
