
library(data.table)
library(tidyverse)
library(ggplot2)
library(DESeq2)

ddir = "/scratch/project_mnt/S0138/DMizikovsky/misid/GTEX/data/"
sample_info <- fread(paste0(ddir, "GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"))
subject_info <- fread(paste0(ddir, "GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt"))
source("/scratch/project_mnt/S0138/DMizikovsky/misid/benchmarking/Benchmarking_Method_Functions_201026_v3.R")

expr_vst <- readRDS(paste0("/scratch/project_mnt/S0138/DMizikovsky/misid/GTEX/GTEx_v11_VST.rds"))
### make sample meta data:
meta <- data.frame(sample = colnames(expr_vst)) %>%
  mutate(ID = sub("^([^-]+-[^-]+).*", "\\1", sample)) %>%
  left_join(select(sample_info, SAMPID, SMTS, SMTSD), by = c("sample"="SAMPID"))

### add misids:
spike_misids <- function(meta, prop_misid = 0.1, seed = 123) {

  set.seed(seed)

  meta <- meta %>%
    mutate(
      label_true = as.character(label_true),
      label_observed = as.character(label_true)
    )

  # Number of samples to misannotate
  n_misid <- round(nrow(meta) * prop_misid)

  # Randomly select samples
  misid_idx <- sample(seq_len(nrow(meta)), n_misid)

  # Available labels
  labels <- unique(meta$label_true)

  # Assign each selected sample a different label
  meta$label_observed[misid_idx] <- vapply(
    misid_idx,
    function(i) {
      sample(
        labels[labels != meta$label_true[i]],
        size = 1
      )
    },
    character(1)
  )

  meta <- meta %>%
    mutate(
      is_misidentified = label_observed != label_true
    )

  return(meta)
}

# group column: 
group_col <- "ID"
n_groups <- 10
max_n_per_group <- 30

set.seed(123)

meta_sub <- meta %>%
  mutate(
    group_col = .data[[group_col]],
    label_true = as.character(group_col),
    label_observed = as.character(group_col),
    ID = sub("^([^-]+-[^-]+).*", "\\1", sample)
  ) %>%
  filter(
    group_col %in% sample(unique(group_col), n_groups)
  ) %>%
  group_by(group_col) %>%
  slice_sample(prop = 1) %>%
  slice_head(n = max_n_per_group) %>%
  ungroup()

expr_sub <- expr_vst[
  ,
  meta_sub$sample,
  drop = FALSE
]

stopifnot(
  identical(colnames(expr_sub), meta_sub$sample)
)

meta_spiked <- spike_misids(
  meta_sub,
  prop_misid = 0.05
)


results1 <- SampleAnnotationLikelihood::check_annotations(
  expr = expr_sub,
  metadata = meta_spiked,
  sample_col = "sample",
  id_col = "label_observed",
  cv_folds = 5,
  detect_unknown = FALSE,
  target_specificity = 0.95
)