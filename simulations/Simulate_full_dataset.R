custom_lib <- "/scratch/project/palpant_scratch/DMizikovsky/R_packages"
.libPaths(c(custom_lib, .libPaths()))

# Simulate repeated-sample data with label-driven and competing variance, plus mis-identified labels.
# Grid: data type x label variance share x discrete and continuous variance shares x samples/label.
library(tidyverse)
library(data.table)

# ---- Simulation -------------------------------------------------------------

# Latent log-scale = label + batch + residual; variance shares sum to 1 (times sigma2_total).
simulate_latent <- function(n_labels, n_per_label, n_features,
                            var_label, var_discrete, var_cont,
                            n_batches = 4, sigma2_total = 1, shared_time = FALSE) {
  stopifnot(var_label + var_discrete + var_cont <= 1)
  n     <- n_labels * n_per_label
  label <- rep(seq_len(n_labels), each = n_per_label)
  # Discrete batch and continuous time are assigned independently of label
  batch <- sample.int(n_batches, n, replace = TRUE)
  # shared_time = TRUE: every label has the same evenly spaced time points (sample j of each label at the j-th point)
  time  <- if (shared_time) rep(seq(0, 1, length.out = n_per_label), times = n_labels) else runif(n)

  fscale <- sqrt(rexp(n_features))
  fscale <- fscale / sqrt(mean(fscale^2))
  sd_mat <- function(k, v) matrix(rnorm(k * n_features), k, n_features) *
    rep(fscale * sqrt(sigma2_total * v), each = k)

  # Smooth time effect: linear + quadratic trend (orthogonal, unit-variance bases),
  # each with feature-specific random coefficients and half of the time variance share.
  time_effect <- function(v) {
    if (v == 0) return(matrix(0, n, n_features))
    deg <- min(2, length(unique(time)) - 1)
    B <- scale(poly(time, deg))
    # With only 2 unique time points a quadratic is unidentifiable, so the linear term takes all of v
    vk <- if (deg == 2) v / 2 else v
    Reduce(`+`, lapply(seq_len(deg), function(d) B[, d] %*% sd_mat(1, vk)))
  }

  z <- sd_mat(n_labels, var_label)[label, , drop = FALSE] +
    sd_mat(n_batches, var_discrete)[batch, , drop = FALSE] +
    time_effect(var_cont) +
    sd_mat(n, max(0, 1 - var_label - var_discrete - var_cont))
  list(z = z, label = label, batch = batch, time = time)
}

# Map latent values to the observed (log-scale) data type.
to_observed <- function(z, type) {
  n <- nrow(z); p <- ncol(z)
  switch(type,
    gaussian = z,
    rnaseq = {
      base <- rnorm(p, 5, 2)
      lib  <- exp(rnorm(n, 0, 0.3))
      mu   <- exp(z + rep(base, each = n)) * lib
      size <- 1 / (0.1 + 1 / sqrt(colMeans(mu) + 1))
      counts <- matrix(rnbinom(n * p, mu = mu, size = rep(size, each = n)), n, p)
      log2((counts + 0.5) / (rowSums(counts) + 1) * 1e6)
    },
    olink = {
      # NPX is log2; censor at a per-assay limit of detection
      npx <- z / log(2) + rep(rnorm(p, 4, 1.5), each = n) + rnorm(n * p, 0, 0.2)
      lod <- apply(npx, 2, quantile, probs = 0.05)
      pmax(npx, rep(lod, each = n))
    },
    stop("unknown type: ", type)
  )
}

# Relabel a fraction of samples to a different random label.
introduce_misid <- function(label, rate) {
  n <- length(label); K <- max(label)
  is_misid <- logical(n)
  is_misid[sample.int(n, round(n * rate))] <- TRUE
  new <- label
  # Offset in 1..(K-1) guarantees a different label
  new[is_misid] <- (label[is_misid] - 1L + sample.int(K - 1L, sum(is_misid), replace = TRUE)) %% K + 1L
  list(label = new, is_misid = is_misid)
}

simulate_dataset <- function(type, var_label, var_discrete, var_cont, n_per_label,
                             n_labels = 25, n_features = 300, misid_rate = 0.1,
                             n_batches = 4, shared_time = TRUE) {
  sim <- simulate_latent(n_labels, n_per_label, n_features,
                         var_label, var_discrete, var_cont, n_batches,
                         shared_time = shared_time)
  X   <- to_observed(sim$z, type)
  mi  <- introduce_misid(sim$label, misid_rate)
  ids <- sprintf("S%05d", seq_along(sim$label))
  dimnames(X) <- list(ids, sprintf("F%04d", seq_len(n_features)))
  list(
    X = X,
    meta = data.table(sample = ids, label_true = sim$label, label_observed = mi$label,
                      is_misid = mi$is_misid, batch = sim$batch, time = sim$time)
  )
}

# ---- Scenario grid ----------------------------------------------------------

# Easier than the first version: label variance below ~0.2 mostly fails the
# package's 70% agreement gate (nothing to evaluate), and n_per_label = 5 breaks
# cv_folds = 5 once misidentified samples shrink a label below 5 samples.
# Gaussian gets the full factorial; the other data types a reduced one.
grid <- rbind(
  CJ(type = "gaussian",
     var_label    = c(0.2, 0.35, 0.5),
     var_discrete = c(0, 0.2, 0.4),
     var_cont     = c(0, 0.2),
     n_per_label  = c(6, 10, 20),
     misid_rate   = c(0, 0.05, 0.1, 0.2),  # 0 measures false positives on clean data
     rep          = 1:3),
  CJ(type = c("rnaseq", "olink"),
     var_label    = c(0.3, 0.5),
     var_discrete = c(0, 0.3),
     var_cont     = c(0, 0.2),
     n_per_label  = c(10, 20, 30),
     misid_rate   = c(0, 0.1),
     rep          = 1:3),
  # Low label variance: many of these fail the 70% agreement gate, which is
  # the point (checks the package reports insufficient_evidence, not noise).
  # Appended last so earlier scenario numbers and seeds are unchanged.
  CJ(type = "gaussian",
     var_label    = c(0.05, 0.1),
     var_discrete = c(0, 0.2),
     var_cont     = 0,
     n_per_label  = c(10, 20),
     misid_rate   = c(0, 0.1),
     rep          = 1:3)
)[var_label + var_discrete + var_cont <= 1]
grid[, scenario := .I]
grid[, seed := 1000L + scenario]

# Each scenario has its own seed, so changing the grid does not change the
# data of the remaining scenarios.
simulate_grid <- function(grid, ...) {
  res <- lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i]
    set.seed(g$seed)
    d <- simulate_dataset(g$type, g$var_label, g$var_discrete, g$var_cont, g$n_per_label, misid_rate = g$misid_rate, ...)
    d$params <- g
    d
  })
  names(res) <- grid$scenario
  res
}

sims <- simulate_grid(grid)
# Run from the repo root so this lands where the runner looks for it.
dir.create("simulations/results", recursive = TRUE, showWarnings = FALSE)
saveRDS(sims, "simulations/results/simulated_datasets_v2.rds")
