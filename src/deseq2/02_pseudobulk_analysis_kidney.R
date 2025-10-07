### DEA using DESeq2 of S858R mouse kidney data
# Author: Tabea M. Soelter (Code adapted from Emma Jones)
# Date: September 30th, 2025

## Goal: Perform a differential gene expression analysis in the S858R kidneys.

## Reproducibility:
# * GitHub: lasseignelab/230926_EJ_Setbp1_AlternativeSplicing
# * Docker: emmafjones/setbp1_alternative_splicing
#       * Version: 1.0.8 

## Data:
# Annotated Seurat object
# * Consists of 2 conditions (S858R-mutants, wildtype) and 1 sex (male).
# * Name: annotated_kidney_samples.rds
# * Location: data/seurat/

## Analysis Plan:
# * Load necessary packages 
# * Load data 
# * Prepare pbject for DEA with DESeq2
# * DEA with DESeq2
# * Save all and significantly expressed genes

## Analysis:
# Set-up
args <- R.utils::commandArgs(trailingOnly = TRUE)

wd <- args[1]
setwd(wd)

set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ComplexHeatmap)
  library(viridis)
  library(Seurat)
  library(Matrix.utils)
  library(edgeR)
  library(Matrix)
  library(reshape2)
  library(S4Vectors)
  library(SingleCellExperiment)
  library(pheatmap)
  library(apeglm)
  library(png)
  library(DESeq2)
  library(RColorBrewer)
  library(data.table)
})

source(file.path(
  wd,
  "src",
  "functions_soelter.R"
))

# Load data
seurat <- readRDS(file.path(wd, "data", "seurat", "annotated_kidney_samples.rds"))

# Prepare SCE object
counts <- seurat@assays[["RNA"]]@layers[["counts"]]
rownames(counts) <- rownames(seurat)

metadata <- seurat@meta.data
metadata$cluster_id <- factor(metadata$cell_type)
metadata$sample_id <- gsub("_", "", metadata$sample_id)
metadata$sample_id <- factor(metadata$sample_id)

sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = metadata
)

# Aggregate counts
cluster_names <- levels(colData(sce)$cluster_id)

sample_names <- levels(colData(sce)$sample_id)

groups <- colData(sce)[, c("cluster_id", "sample_id")]

aggr_counts <- aggregate.Matrix(
  t(counts(sce)),
  groupings = groups,
  fun = "sum"
  )
aggr_counts <- t(aggr_counts)

# Create cell-type-specific counts matrices
counts_ls <- list()

for (i in seq_along(cluster_names)) {
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
  
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
  }

# Prepare cell-type-specific metadata aggregation
metadata <- colData(sce) %>%
  as.data.frame() %>%
  dplyr::select(seq_folder, sample_id)
metadata <- metadata[!duplicated(metadata), ]
colnames(metadata) <- c("group_id", "sample_id")
rownames(metadata) <- metadata$sample_id

t <- table(
  colData(sce)$sample_id,
  colData(sce)$cluster_id
)

# Create cell-type-specific metadata
metadata_ls <- list()

for (i in seq_along(counts_ls)) {
  df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
  
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  
  df$sample_id <- tstrsplit(df$cluster_sample_id, "_")[[2]]
  
  idx <- which(colnames(t) == unique(df$cluster_id))
  
  cell_counts <- t[, idx]
  cell_counts <- cell_counts[cell_counts > 0]
  sample_order <- match(df$sample_id, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  df$cell_count <- cell_counts
  df <- plyr::join(df, metadata,
                   by = intersect(names(df), names(metadata))
  )
  rownames(df) <- df$cluster_sample_id
  
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df$cluster_id)
}

# Run DESeq2 on all cell types
map(
  cluster_names,
  get_dds_resultsAvsB,
  A = "wildtype",
  B = "mutant",
  save_path = "results/deseq2/kidney/"
)

sessionInfo()
