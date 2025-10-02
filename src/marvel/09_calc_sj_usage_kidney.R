### Kidney splice junction usage calculations
# Author: Tabea M. Soelter (OG code by Emma Jones)
# Date: August 20th, 2025

## Goal: Calculate SJU for each cell type for each splice junction

## Reproducibility:
# * GitHub: lasseignelab/230926_EJ_Setbp1_AlternativeSplicing
# * Docker: emmafjones/setbp1_alternative_splicing
#       * Version: 1.1.0

## Data:
# MARVEL object
#  * Name: setbp1_kidney_marvel_aligned.rds
#  * Location: data/marvel/

## Analysis Plan:
#  * Load necessary packages
#  * Load data from Seurat & STAR Solo
#  * Prepare data
#  * Make MARVEL object
#  * Process MARVEL object
#  * Save MARVEL object

## Analysis:
# Setup
args <- R.utils::commandArgs(trailingOnly = TRUE)

wd <- args[1]
setwd(wd)

set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ComplexHeatmap)
  library(viridis)
  library(MARVEL)
  library(ggnewscale)
  library(ggrepel)
  library(reshape2)
  library(plyr)
  library(stringr)
  library(textclean)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(gtools)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(wiggleplotr)
  library(Matrix)
  library(data.table)
  library(gridExtra)
})

source(file.path(wd, "src", "marvel", "functions.R"))

# Load data
setbp1_marvel <- read_rds(file.path(wd, 
  "data", "marvel",
  "setbp1_kidney_marvel_aligned.rds"
))

# Prepare data
sample_metadata <- setbp1_marvel$sample.metadata

sj_metadata <- setbp1_marvel$sj.metadata

sj_counts <- setbp1_marvel[["sj.count.matrix"]]
sj_counts <- as(sj_counts, "CsparseMatrix")

gene_counts <- setbp1_marvel[["gene.count.matrix"]]
gene_counts <- as(gene_counts, "CsparseMatrix")

# Define cell groups
cell_group_list <- sample_metadata %>%
  group_split(cell_type, .keep = TRUE) %>%
  {set_names(map(., ~ set_names(.$cell.id, .$cell_type[1])), 
             map_chr(., ~ .$cell_type[1]))}

mutant_list <- sample_metadata %>%
  group_split(seq_folder, .keep = TRUE) %>%
  map(~ set_names(.$cell.id, .$seq_folder[[1]]))

mutant_list <- set_names(mutant_list, c(
  "Mutant", "Wildtype"
))

# Split matrices by cell type
split_matrices_list <- lapply(names(cell_group_list), function(x) {
  subset_cell_type_matrices(cell_type = x)
})
names(split_matrices_list) <- names(cell_group_list)

# Split matrices by condition
split_mutants_list <- lapply(
  split_matrices_list,
  function(x) subset_mutant_matrices(cell_type = x)
)

# Calculate cell type splice junction usage
sj_usage_cell_type_list <- lapply(split_matrices_list, get_sj_usage_cell_type)

sj_usage_cell_type_df <- do.call(rbind, sj_usage_cell_type_list)
sj_usage_cell_type_df <- t(sj_usage_cell_type_df)

# Save cell-type-specific SJU results
write_rds(
  sj_usage_cell_type_df,
  file.path(wd, "data", "marvel", "sj_usage_kidney_cell_type.rds")
)

# Calculate condition splice junction usage
sj_usage_condition_list <- lapply(split_mutants_list, get_sj_usage_condition)

sj_usage_condition_df <- do.call(rbind, sj_usage_condition_list)
sj_usage_condition_df <- t(sj_usage_condition_df)

# Save condition-specific SJU results
write_rds(
  sj_usage_condition_df,
  file.path(wd, "data", "marvel", "sj_usage_kidney_condition.rds")
)

# Save other results
saveRDS(sj_counts, file.path(wd, "data", "marvel", "sj_counts_kidney.rds"))
saveRDS(gene_counts, file.path(wd, "data", "marvel", "gene_counts_kidney.rds"))
saveRDS(split_matrices_list, file.path(wd, "data", "marvel", "split_matrices_list_kidney.rds"))

# Get session information
sessionInfo()
