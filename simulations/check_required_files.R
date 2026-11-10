args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1L]] else "gtex"
if (!mode %in% c("gtex", "simulated", "all")) {
  stop("Usage: Rscript simulations/check_required_files.R [gtex|simulated|all]")
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
repo_root <- dirname(dirname(script_path))
gtex_data <- Sys.getenv(
  "GTEX_DATA_DIR",
  "/home/projects/ssi_gen1/people/dalmiz/Sample_Misassignment_0826/GTEX/data"
)

gtex_files <- c(
  "GTEx sample attributes" = file.path(gtex_data, "GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"),
  "GTEx expression matrix" = file.path(gtex_data, "GTEx_v11_VST.rds"),
  "Generated GTEx scenarios" = file.path(repo_root, "simulations/results/gtex_scenarios.rds")
)
simulated_files <- c(
  "Generated simulated datasets" = file.path(repo_root, "simulations/results/simulated_datasets_v2.rds")
)
required_files <- switch(
  mode,
  gtex = gtex_files,
  simulated = simulated_files,
  all = c(gtex_files, simulated_files)
)

present <- file.exists(required_files)
for (i in seq_along(required_files)) {
  cat(sprintf("%s  %s: %s\n", if (present[[i]]) "OK" else "MISSING",
              names(required_files)[[i]], required_files[[i]]))
}

if (any(!present)) {
  quit(save = "no", status = 1L)
}
cat("All required files exist.\n")