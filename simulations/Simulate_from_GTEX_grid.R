custom_lib <- "/scratch/project/palpant_scratch/DMizikovsky/R_packages"
.libPaths(c(custom_lib, .libPaths()))

# Subsampled GTEx datasets across a grid of scenarios, for running
# check_annotations() against known mis-identifications.
#
# Two tasks:
#   tissue     label = tissue (SMTS); one sample per individual per label, so
#              labels do not share donors.
#   individual label = donor; each donor contributes n_per_label samples from
#              a set of tissues. Tissue differences are then the competing
#              variance, so the choice of tissues matters:
#                similar     the k tissues closest to a random seed tissue
#                random      k random tissues
#                dissimilar  k tissues spread out (farthest-first)
#              "similar" should make donor discrimination easier because
#              tissue variance is smaller.
#
# Only sample IDs and labels are stored per scenario (small). Expression is
# sliced from the full matrix when a scenario is run: build_gtex_dataset().
#
# Expected inputs (as in the draft script): `expr_vst` genes x samples (VST),
# `sample_info` with SAMPID, SMTS, SMTSD.
suppressMessages(library(data.table))

if (!exists("expr_vst") || !exists("sample_info")) {
  ddir <- "/scratch/project_mnt/S0138/DMizikovsky/misid/GTEX/"
  sample_info <- fread(paste0(ddir, "data/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"))
  expr_vst <- readRDS(paste0(ddir, "GTEx_v11_VST.rds"))
}
if (!exists("min_tissue_samples")) min_tissue_samples <- 200  # tissues eligible for sets
if (!exists("n_features")) n_features <- 2000                  # most variable genes per dataset

meta_all <- data.table(sample = colnames(expr_vst))
meta_all[, ID := sub("^([^-]+-[^-]+).*", "\\1", sample)]
meta_all <- merge(meta_all, unique(sample_info[, .(sample = SAMPID, SMTS, SMTSD)]),
                  by = "sample")
meta_all <- meta_all[!is.na(SMTS) & SMTS != "" & !is.na(SMTSD) & SMTSD != ""]

# ---- Tissue similarity --------------------------------------------------------

# Correlation between tissue centroids over the most variable centroid genes.
tissue_similarity <- function(col, max_per_tissue = 100) {
  groups <- split(meta_all$sample, meta_all[[col]])
  groups <- groups[lengths(groups) >= 20]
  cent <- vapply(groups, function(s) {
    s <- if (length(s) > max_per_tissue) sample(s, max_per_tissue) else s
    rowMeans(expr_vst[, s, drop = FALSE])
  }, numeric(nrow(expr_vst)))
  top <- order(apply(cent, 1, var), decreasing = TRUE)[seq_len(min(5000, nrow(cent)))]
  cor(cent[top, , drop = FALSE])
}

set.seed(1)
sim_smtsd <- tissue_similarity("SMTSD")
sim_smts  <- tissue_similarity("SMTS")

# Choose k labels from `pool` by similarity to a random seed label.
pick_labels <- function(sim, pool, k, mode) {
  pool <- intersect(pool, rownames(sim))
  k <- min(k, length(pool))
  if (mode == "random") return(sample(pool, k))
  first <- sample(pool, 1)
  if (mode == "similar") {
    return(names(sort(sim[first, pool], decreasing = TRUE))[seq_len(k)])
  }
  chosen <- first  # dissimilar: farthest-first
  while (length(chosen) < k) {
    rest <- setdiff(pool, chosen)
    nearest <- vapply(rest, function(t) max(sim[t, chosen]), 1)
    chosen <- c(chosen, rest[which.min(nearest)])
  }
  chosen
}

# ---- Scenario construction ----------------------------------------------------

# Relabel a fraction of samples to a different label from the same set. The
# package needs at least 2 samples per observed label, so redraw until every
# label keeps that many.
introduce_misid <- function(label, rate, min_count = 2L, max_tries = 200L) {
  labs <- sort(unique(label))
  n_mis <- round(length(label) * rate)
  for (try in seq_len(max_tries)) {
    is_misid <- logical(length(label))
    is_misid[sample.int(length(label), n_mis)] <- TRUE
    new <- label
    new[is_misid] <- vapply(label[is_misid], function(l) sample(setdiff(labs, l), 1), "")
    if (all(table(factor(new, levels = labs)) >= min_count)) {
      return(list(label = new, is_misid = is_misid))
    }
  }
  NULL
}

make_tissue_scenario <- function(n_labels, n_per_label, set) {
  # Labels need enough distinct donors to sample from.
  donors <- meta_all[, .(n = uniqueN(ID)), by = SMTS][n >= 2 * n_per_label, SMTS]
  labs <- pick_labels(sim_smts, donors, n_labels, set)
  rbindlist(lapply(labs, function(l) {
    d <- meta_all[SMTS == l][sample(.N)]
    d <- d[!duplicated(ID)][seq_len(n_per_label)]  # one sample per donor
    d[, label_true := l]
  }))
}

make_individual_scenario <- function(n_labels, n_per_label, set, n_tissues = 8L) {
  pool <- meta_all[, .N, by = SMTSD][N >= min_tissue_samples, SMTSD]
  tis <- pick_labels(sim_smtsd, pool, n_tissues, set)
  d <- meta_all[SMTSD %in% tis]
  elig <- d[, .N, by = ID][N >= n_per_label, ID]
  if (length(elig) < 5L) return(NULL)
  keep <- sample(elig, min(n_labels, length(elig)))
  d <- d[ID %in% keep][sample(.N)]
  d <- d[, head(.SD, n_per_label), by = ID]
  d[, label_true := ID]
  d[, tissues := paste(tis, collapse = "|")]
  d
}

build_scenario <- function(g) {
  set.seed(g$seed)
  d <- if (g$task == "tissue") {
    make_tissue_scenario(g$n_labels, g$n_per_label, g$set)
  } else {
    make_individual_scenario(g$n_labels, g$n_per_label, g$set)
  }
  if (is.null(d) || uniqueN(d$label_true) < 3L) return(NULL)
  mi <- introduce_misid(d$label_true, g$misid_rate)
  if (is.null(mi)) return(NULL)
  d[, `:=`(label_observed = mi$label, is_misid = mi$is_misid)]
  d[, .(sample, ID, SMTS, SMTSD, label_true, label_observed, is_misid)]
}

# Samples x features matrix (most variable genes within the dataset).
build_gtex_dataset <- function(meta, n_features = 2000) {
  x <- expr_vst[, meta$sample, drop = FALSE]
  v <- apply(x, 1, var)
  x <- x[order(v, decreasing = TRUE)[seq_len(min(n_features, nrow(x)))], , drop = FALSE]
  t(x)
}

# ---- Scenario grid ------------------------------------------------------------

grid <- rbind(
  CJ(task = "tissue",
     set = c("random", "similar"),
     n_labels = c(6, 12),
     n_per_label = c(10, 20),
     misid_rate = c(0, 0.05, 0.1, 0.2),  # 0 measures false positives on clean data
     rep = 1:3),
  CJ(task = "individual",
     set = c("similar", "random", "dissimilar"),
     n_labels = c(20, 40),
     n_per_label = c(4, 6),
     misid_rate = c(0, 0.05, 0.1, 0.2),
     rep = 1:3)
)
grid[, scenario := .I]
grid[, seed := 5000L + scenario]

simulate_gtex_grid <- function(grid) {
  res <- lapply(seq_len(nrow(grid)), function(i) {
    m <- build_scenario(grid[i])
    if (is.null(m)) return(NULL)
    list(meta = m, params = grid[i])
  })
  names(res) <- grid$scenario
  res[!vapply(res, is.null, logical(1))]
}

if (sys.nframe() == 0L) {
  gtex_sims <- simulate_gtex_grid(grid)
  cat(length(gtex_sims), "of", nrow(grid), "scenarios built\n")
  # Run from the repo root so this lands where the runner looks for it.
  dir.create("simulations/results", recursive = TRUE, showWarnings = FALSE)
  saveRDS(gtex_sims, "simulations/results/gtex_scenarios.rds")
}

# Running one scenario:
#   s <- gtex_sims[["1"]]
#   X <- build_gtex_dataset(s$meta)                      # samples x features
#   check_annotations(t(X), as.data.frame(s$meta), "sample", "label_observed",
#                     cv_folds = min(5, min(table(s$meta$label_observed))),
#                     detect_unknown = FALSE)
