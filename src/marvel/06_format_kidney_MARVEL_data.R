### MARVEL analysis on S858R kidney data
# Author: Tabea M. Soelter (OG code by Emma Jones)
# Date: August 7th, 2025

## Goal: Format Seurat expression data and SJ outputs from STARsolo for MARVEL.

## Reproducibility:
# * GitHub: lasseignelab/230926_EJ_Setbp1_AlternativeSplicing
# * Docker: emmafjones/setbp1_alternative_splicing
#       * Version: 1.1.0

## Data:
# STAR Solo outputs
#  * Name: N/A
#  * Location: data/star/
# Seurat object
#  * Name: annotated_kidney_samples.rds
#  * Location: data/seurat/

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
  library(Seurat)
  library(patchwork)
  library(harmony)
  library(presto)
  library(MARVEL)
  library(ggnewscale)
  library(ggrepel)
  library(reshape2)
  library(plyr)
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
  library(MatrixExtra)
  options("MatrixExtra.quick_show" = FALSE)
})

source(file.path(wd, "src", "marvel", "functions.R"))
source(file.path(wd, "src", "functions_soelter.R"))

# Load data
annotated_kidney_samples <- read_rds(file.path(wd,
                                               "data",
                                               "seurat",
                                               "annotated_kidney_samples.rds"
                                               )
                                     )

# Format data
gene_norm <- annotated_kidney_samples@assays[["RNA"]]@layers[["data"]]
gene_norm <- as(gene_norm, "TsparseMatrix")
gene_norm[gene_norm[, 1] != 0, ][1:5, 1:5]

gene_metadata <- annotated_kidney_samples@meta.data
gene_metadata <- tibble::rownames_to_column(gene_metadata, "cell.id")

gene_features <- rownames(annotated_kidney_samples)
gene_features <- as.data.frame(gene_features)
colnames(gene_features) <- "gene_short_name"

gene_counts <- annotated_kidney_samples@assays[["RNA"]]@layers[["counts"]]
gene_counts <- as(gene_counts, "TsparseMatrix")
gene_counts[gene_counts[, 1] != 0, ][1:5, 1:5]

# Modify K6 matrix file
# Remove line 5808747, which I identified using grep("4294967295", lines)
# The value 4294967295 was returned when splitting SJ info using the OG matrix.mtx file
K6_matrix <- readLines(
  file.path(
    wd,
    "data",
    "star",
    "K6",
    "Solo.out",
    "SJ",
    "raw",
    "matrix.mtx"
    )
  )
K6_matrix_fixed <- K6_matrix[-5808747]
writeLines(K6_matrix_fixed,
           file.path(
             wd, 
             "data", 
             "star",
             "K6",
             "Solo.out",
             "SJ",
             "raw",
             "matrix_fixed.mtx"
             )
           )

# Import splice junction data
samples <- c("K1", "K2", "K3", "K4", "K5", "K6")

for (sample in samples) {
  var_name <- paste0(sample, "_sj_info")
  
  assign(var_name, split_sj_info_ts(sample_id = sample))
}

# Combine data
matrices <- list(
  K2_sj_info, K4_sj_info, K6_sj_info, K1_sj_info, K3_sj_info,
  K5_sj_info
)

all_features <- unique(unlist(lapply(matrices, rownames)))

combined_matrix <- do.call(cbind, lapply(matrices, align_matrix, all_features))
colnames(combined_matrix) <- gene_metadata$cell.id
combined_matrix <- as(combined_matrix, "TsparseMatrix")

all_sj_counts <- combined_matrix

sj_metadata <- colnames(combined_matrix)
sj_metadata <- as.data.frame(sj_metadata)
colnames(sj_metadata) <- "cell.id"

sj_features <- rownames(combined_matrix)
sj_features <- as.data.frame(sj_features)
colnames(sj_features) <- "coord.intron"

# Get UMAP coordinates
umap_coords <- annotated_kidney_samples@reductions[["umap_harmony"]]@cell.embeddings
umap_coords <- tibble::rownames_to_column(as.data.frame(umap_coords), "cell.id")
colnames(umap_coords) <- c("cell.id", "x", "y")

# Get and prepare GTF file
genome_dir <- "/data/project/lasseigne_lab/GENOME_dir/GENCODE_mm39/release_M31/"
m31_gtf <- as.data.frame(
  fread(paste0(genome_dir, "GTF/gencode.vM31.primary_assembly.annotation.gtf")),
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE)
colnames(m31_gtf) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9")
m31_gtf$V1 <- substr(m31_gtf$V1, 4, nchar(m31_gtf$V1))

# Make MARVEL object
setbp1_marvel <- CreateMarvelObject.10x(
  gene.norm.matrix = gene_norm,
  gene.norm.pheno = gene_metadata,
  gene.norm.feature = gene_features,
  gene.count.matrix = gene_counts,
  gene.count.pheno = gene_metadata,
  gene.count.feature = gene_features,
  sj.count.matrix = all_sj_counts,
  sj.count.pheno = sj_metadata,
  sj.count.feature = sj_features,
  pca = umap_coords,
  gtf = m31_gtf
)

# Process MARVEL object
setbp1_marvel <- AnnotateGenes.10x(MarvelObject = setbp1_marvel)
setbp1_marvel <- AnnotateSJ.10x(MarvelObject = setbp1_marvel)
setbp1_marvel <- ValidateSJ.10x(MarvelObject = setbp1_marvel)
setbp1_marvel <- FilterGenes.10x(MarvelObject = setbp1_marvel,
                                 gene.type = "protein_coding")
setbp1_marvel <- CheckAlignment.10x(MarvelObject = setbp1_marvel)

# Save MARVEL object
write_rds(setbp1_marvel,
          file = file.path(
            wd,
            "data",
            "marvel",
            "setbp1_kidney_marvel_aligned.rds"
            )
          )

# Get session information
sessionInfo()
