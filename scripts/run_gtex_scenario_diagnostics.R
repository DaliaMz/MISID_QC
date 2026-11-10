# Run one saved GTEx simulation scenario with intermediate gate diagnostics.
#
# From the repository root:
#   Rscript scripts/run_gtex_scenario_diagnostics.R \
#     scenario=12 scenarios_file=simulations/results/gtex_scenarios.rds \
#     expr=/path/to/GTEx_v11_VST.rds
#
# Optional arguments:
#   out=simulations/results/gtex_diagnostics
#   n_features=2000 cv_folds=5 mismatch_repeats=1
#   margin_coverage=0.025 target_specificity=0.95
#   holdout_reference=all_alternatives (or top)
#   min_agreement=0.7 min_accuracy=0.9 seed=123

args <- commandArgs(trailingOnly = TRUE)
opt <- list(
  scenario = NULL,
  scenarios_file = "simulations/results/gtex_scenarios.rds",
  expr = "/home/projects/ssi_gen1/people/dalmiz/Sample_Misassignment_0826/GTEX/data/GTEx_v11_VST.rds",
  out = "simulations/results/gtex_diagnostics",
  n_features = "2000",
  cv_folds = "5",
  mismatch_repeats = "1",
  margin_coverage = "0.025",
  holdout_reference = "all_alternatives",
  target_specificity = "0.95",
  min_agreement = "0.7",
  min_accuracy = "0.9",
  seed = ""
)
for (arg in args) {
  if (!grepl("=", arg, fixed = TRUE)) {
    stop("Arguments must be key=value: ", arg)
  }
  parts <- strsplit(arg, "=", fixed = TRUE)[[1L]]
  key <- parts[1L]
  if (!key %in% names(opt)) {
    stop("Unknown argument: ", key)
  }
  opt[[key]] <- paste(parts[-1L], collapse = "=")
}
if (is.null(opt$scenario) || !nzchar(opt$scenario)) {
  stop("Specify scenario=<scenario number>.")
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install devtools in this R environment before running this wrapper.")
}
devtools::install_local(
  path = normalizePath(".", mustWork = TRUE),
  force = TRUE,
  upgrade = "never",
  quiet = TRUE
)
suppressPackageStartupMessages(library(SampleAnnotationLikelihood))
check_annotations_fun <- getExportedValue(
  "SampleAnnotationLikelihood", "check_annotations"
)
if (!"diagnostics" %in% names(formals(check_annotations_fun))) {
  stop(
    "The loaded SampleAnnotationLikelihood package does not support ",
    "`diagnostics`. It was loaded from ",
    find.package("SampleAnnotationLikelihood"),
    ". Install the current repository into a library on .libPaths(), then ",
    "restart R and rerun this script."
  )
}

if (!file.exists(opt$scenarios_file)) {
  stop("Scenario file does not exist: ", opt$scenarios_file)
}
if (!file.exists(opt$expr)) {
  stop("Expression RDS does not exist: ", opt$expr)
}
scenario_id <- as.character(opt$scenario)
scenarios <- readRDS(opt$scenarios_file)
scenario_index <- if (is.null(names(scenarios))) integer() else
  match(scenario_id, names(scenarios))
if (length(scenario_index) != 1L || is.na(scenario_index)) {
  scenario_index <- which(vapply(scenarios, function(s) {
    as.character(s$params$scenario) == scenario_id
  }, logical(1)))
}
if (length(scenario_index) != 1L || is.na(scenario_index)) {
  stop("Could not uniquely find scenario ", scenario_id,
       " in ", opt$scenarios_file)
}
scenario <- scenarios[[scenario_index]]
meta <- as.data.frame(scenario$meta, stringsAsFactors = FALSE)
required <- c("sample", "label_observed", "label_true", "is_misid")
if (!all(required %in% names(meta))) {
  stop("Scenario metadata must contain: ", paste(required, collapse = ", "))
}

expr_vst <- readRDS(opt$expr)
if (!is.matrix(expr_vst) || is.null(colnames(expr_vst))) {
  stop("Expression RDS must be a genes-by-samples matrix with sample names.")
}
if (anyDuplicated(meta$sample) || !all(meta$sample %in% colnames(expr_vst))) {
  stop("Scenario sample IDs must be unique and present in the expression matrix.")
}
n_features <- as.integer(opt$n_features)
if (is.na(n_features) || n_features < 1L) {
  stop("n_features must be a positive integer.")
}
expr <- expr_vst[, meta$sample, drop = FALSE]
feature_variance <- apply(expr, 1L, stats::var)
selected <- order(feature_variance, decreasing = TRUE)
selected <- selected[seq_len(min(n_features, length(selected)))]
expr <- expr[selected, , drop = FALSE]

cv_folds <- if (tolower(opt$cv_folds) %in% c("loo", "null")) {
  NULL
} else {
  folds <- as.integer(opt$cv_folds)
  if (is.na(folds) || folds < 2L) {
    stop("cv_folds must be 'loo' or an integer of at least 2.")
  }
  min(folds, min(table(meta$label_observed)))
}
repeats <- as.integer(opt$mismatch_repeats)
if (is.na(repeats) || repeats < 1L) {
  stop("mismatch_repeats must be a positive integer.")
}
seed <- if (nzchar(opt$seed)) as.integer(opt$seed) else
  as.integer(scenario$params$seed)
if (length(seed) != 1L || is.na(seed)) {
  stop("Provide seed=<integer> or ensure the scenario has an integer seed.")
}
set.seed(seed)

result <- check_annotations_fun(
  expr = expr,
  metadata = data.frame(
    sample = meta$sample,
    label_observed = as.character(meta$label_observed),
    stringsAsFactors = FALSE
  ),
  sample_col = "sample",
  id_col = "label_observed",
  cv_folds = cv_folds,
  mismatch_repeats = repeats,
  margin_coverage = as.numeric(opt$margin_coverage),
  holdout_reference = opt$holdout_reference,
  target_specificity = as.numeric(opt$target_specificity),
  min_agreement = as.numeric(opt$min_agreement),
  min_accuracy = as.numeric(opt$min_accuracy),
  detect_unknown = FALSE,
  diagnostics = TRUE,
  progress = TRUE
)

source(file.path("scripts", "plot_check_annotation_gates.R"))
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(opt$out, paste0("scenario_", scenario_id))
pdf_file <- paste0(prefix, "_gates.pdf")
rds_file <- paste0(prefix, "_diagnostics.rds")
tables <- plot_check_annotation_gates(result, file = pdf_file)
source(file.path("scripts", "plot_first_pass_truth.R"))
truth_pdf_file <- paste0(prefix, "_first_pass_truth.pdf")
first_pass_truth <- plot_first_pass_truth(result, meta, file = truth_pdf_file)
saveRDS(
  list(
    result = result,
    diagnostics = tables,
    first_pass_truth = first_pass_truth,
    truth = meta,
    params = scenario$params,
    settings = list(
      n_features = length(selected),
      cv_folds = cv_folds,
      mismatch_repeats = repeats,
      margin_coverage = as.numeric(opt$margin_coverage),
      holdout_reference = opt$holdout_reference,
      seed = seed
    )
  ),
  rds_file
)
cat("Scenario", scenario_id, "diagnostics saved:\n",
    "  PDF:", pdf_file, "\n",
    "  Truth PDF:", truth_pdf_file, "\n",
    "  RDS:", rds_file, "\n")
