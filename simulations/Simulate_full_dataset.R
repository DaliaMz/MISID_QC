custom_lib <- "/scratch/project/palpant_scratch/DMizikovsky/R_packages"
.libPaths(c(custom_lib, .libPaths()))

# Simulate repeated-sample data with label-driven and competing variance, plus mis-identified labels.
# Grid: data type x label variance share x discrete and continuous variance shares x samples/label.
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

# Relabel a fraction of samples to a different random label. The package needs
# at least 2 samples per observed label, so redraw (rarely needed at low rates)
# until every label keeps that many.
introduce_misid <- function(label, rate, min_count = 2L, max_tries = 200L) {
  n <- length(label); K <- max(label)
  for (try in seq_len(max_tries)) {
    is_misid <- logical(n)
    is_misid[sample.int(n, round(n * rate))] <- TRUE
    new <- label
    # Offset in 1..(K-1) guarantees a different label
    new[is_misid] <- (label[is_misid] - 1L + sample.int(K - 1L, sum(is_misid), replace = TRUE)) %% K + 1L
    if (all(tabulate(new, K) >= min_count)) break
  }
  list(label = new, is_misid = is_misid)
}

simulate_dataset <- function(type, var_label, var_discrete, var_cont, n_per_label,
                             n_labels = 25, n_features = 300, misid_rate = 0.1,
                             n_batches = 4, shared_time = TRUE,
                             n_unknown_individuals = 1L,
                             unknown_samples_per_individual = 2L) {
  if (length(n_unknown_individuals) != 1L || !is.numeric(n_unknown_individuals) ||
      !is.finite(n_unknown_individuals) || n_unknown_individuals < 0 ||
      n_unknown_individuals != floor(n_unknown_individuals) ||
      length(unknown_samples_per_individual) != 1L ||
      !is.numeric(unknown_samples_per_individual) ||
      !is.finite(unknown_samples_per_individual) ||
      unknown_samples_per_individual < 1 ||
      unknown_samples_per_individual != floor(unknown_samples_per_individual) ||
      unknown_samples_per_individual > n_per_label) {
    stop("n_unknown_individuals must be a non-negative integer and unknown_samples_per_individual must be an integer between 1 and n_per_label.")
  }
  sim <- simulate_latent(n_labels + n_unknown_individuals, n_per_label, n_features,
                         var_label, var_discrete, var_cont, n_batches,
                         shared_time = shared_time)
  enrolled <- which(sim$label <= n_labels)
  unknown <- if (n_unknown_individuals > 0L) {
    unlist(lapply(seq_len(n_unknown_individuals), function(i) {
      indices <- which(sim$label == n_labels + i)
      indices[sample.int(length(indices), unknown_samples_per_individual)]
    }), use.names = FALSE)
  } else {
    integer()
  }
  keep <- c(enrolled, unknown)
  X <- to_observed(sim$z[keep, , drop = FALSE], type)
  mi <- introduce_misid(sim$label[enrolled], misid_rate)
  unknown_labels <- if (length(unknown)) {
    sample.int(n_labels, length(unknown), replace = TRUE)
  } else {
    integer()
  }
  label_true <- sim$label[keep]
  is_unknown <- label_true > n_labels
  ids <- sprintf("S%05d", seq_along(keep))
  dimnames(X) <- list(ids, sprintf("F%04d", seq_len(n_features)))
  list(
    X = X,
    meta = data.table(sample = ids, label_true = label_true,
                      label_observed = c(mi$label, unknown_labels),
                      is_misid = c(mi$is_misid, rep(TRUE, length(unknown))),
                      is_unknown = is_unknown,
                      batch = sim$batch[keep], time = sim$time[keep])
  )
}

# ---- Scenario grid ----------------------------------------------------------

# One crossed design, so every effect is estimated on the same footing:
#   - gaussian: the full factorial of the factors below.
#   - rnaseq, olink: a small check that the Gaussian findings carry over, at
#     one label-variance level, crossed with the same nuisance levels.
# The questions this design answers (see simulations/README.md): how well
# detection and recovery work, and how far the package's internal estimates sit
# from the truth, mainly as a function of first-pass agreement. That is driven by
# label signal (var_label), mislabel rate and replication; the discrete and
# continuous nuisance variance are contrasted as one factor ("nuisance"), none
# versus both. rate 0 measures false positives on clean data. n_per_label >= 6
# keeps cv_folds = 5 usable. The agreement gate is off in the runner.
grid <- rbind(
  CJ(type = "gaussian",
     var_label    = c(0.05, 0.2, 0.5),
     nuisance     = c("none", "both"),
     n_per_label  = c(6, 10, 20),
     misid_rate   = c(0, 0.02, 0.05, 0.1, 0.2, 0.4),
     rep          = 1:3),
  CJ(type = c("rnaseq", "olink"),
     var_label    = 0.2,
     nuisance     = c("none", "both"),
     n_per_label  = 10,
     misid_rate   = c(0, 0.1, 0.4),
     rep          = 1:3)
)
grid[, `:=`(var_discrete = fifelse(nuisance == "both", 0.3, 0),
            var_cont = fifelse(nuisance == "both", 0.2, 0))]
stopifnot(all(grid$var_label + grid$var_discrete + grid$var_cont <= 1))
grid[, scenario := .I]
grid[, seed := 1000L + scenario]
if (!exists("n_unknown_individuals")) n_unknown_individuals <- 1L
if (!exists("unknown_samples_per_individual")) unknown_samples_per_individual <- 2L
grid[, `:=`(n_unknown_individuals = n_unknown_individuals,
            unknown_samples_per_individual = unknown_samples_per_individual)]

# Each scenario has its own seed, so changing the grid does not change the
# data of the remaining scenarios.
simulate_grid <- function(grid, ...) {
  res <- lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i]
    set.seed(g$seed)
    d <- simulate_dataset(g$type, g$var_label, g$var_discrete, g$var_cont,
                          g$n_per_label, misid_rate = g$misid_rate,
                          n_unknown_individuals = g$n_unknown_individuals,
                          unknown_samples_per_individual = g$unknown_samples_per_individual,
                          ...)
    d$params <- g
    d
  })
  names(res) <- grid$scenario
  res
}

if (sys.nframe() == 0L) {
  sims <- simulate_grid(grid)
  # Run from the repo root so this lands where the runner looks for it.
  dir.create("simulations/results", recursive = TRUE, showWarnings = FALSE)
  saveRDS(sims, "simulations/results/simulated_datasets_v2.rds")
}
