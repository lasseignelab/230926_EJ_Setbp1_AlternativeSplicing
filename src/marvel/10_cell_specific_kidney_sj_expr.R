### Kidney splice junction usage calculations
# Author: Tabea M. Soelter (OG code by Emma Jones)
# Date: August 22nd, 2025

## Goal: Quantify the splice junction expression across cell types, both per
#        cell and aggregated for each cell type.

## Reproducibility:
# * GitHub: lasseignelab/230926_EJ_Setbp1_AlternativeSplicing
# * Docker: emmafjones/setbp1_alternative_splicing
#       * Version: 1.1.0

## Data:
# SJU by cell type object
#   * Name: sj_usage_kidney_cell_type.rds
#   * Location: data/marvel/
# SJU by condition object
#   * Name: sj_usage_kidney_condition.rds
#   * Location: data/marvel/
# SJ counts
#   * Name: sj_counts_kidney.rds
#   * Location: data/marvel/
# Gene counts
#   * Name: gene_counts_kidney.rds
#   * Location: data/marvel/
# MARVEL object
#  * Name: setbp1_kidney_marvel_aligned.rds
#  * Location: data/marvel/
# Split matrices
#   * Name: split_matrices_list_kidney.rds
#   * Location: data/marvel/

## Analysis Plan:
#  * Load necessary packages
#  * Load data
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
  library(cowplot)
  library(MARVEL)
  library(ggnewscale)
  library(ggrepel)
  library(reshape2)
  library(plyr)
  library(stringr)
  library(ggplot2)
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
  library(ggridges)
  library(GenomicFeatures)
  library(rtracklayer)
  library(ggpmisc)
  library(cowplot)
})

source(here::here("src", "marvel", "functions.R"))

# Set colors
cell_type_colors <- c(
  `Proximal tubule cells` = "#F8766D",
  `Thick ascending limb (LOH)` = "#AA937E",
  `PCT` = "#A85A5A",
  `Endothelial cells` = "#CF9400",
  `PST` = "#948802",
  `CDPC` = "#F03F00",
  `DCT` = "#6D0404",
  `Mesenchymal cells` = "#FFC8C4",
  `Thin ascending limb (LOH)` = "#AA6320",
  `CDIC-B` = "#EEDF37",
  `Dendritic cells` = "#C70F0F",
  `Thin descending limb (LOH)` = "#FFE196",
  `CDIC-A` = "#4E3801",
  `Podocytes` = "#BD085E",
  `B cells` = "#FFBE1D",
  `T regulatory cells` = "#00B0F6",
  `T cells` = "#FF7600",
  `Connecting tubule cells` = "#E76BF3"
)

# Load data
sj_usage_cell_type <- readRDS(here::here(
  "data", "marvel",
  "sj_usage_kidney_cell_type.rds"
))

sj_usage_condition <- readRDS(here::here(
  "data", "marvel",
  "sj_usage_kidney_condition.rds"
))

sj_counts <- readRDS(here::here("data", "marvel", "sj_counts_kidney.rds"))
sj_counts <- as(sj_counts, "CsparseMatrix")

gene_counts <- readRDS(here::here("data", "marvel", "gene_counts_kidney.rds"))
gene_counts <- as(gene_counts, "CsparseMatrix")

setbp1_marvel <- readRDS(here::here(
  "data", "marvel",
  "setbp1_kidney_marvel_aligned.rds"
))

split_matrices_list <- readRDS(here::here("data", "marvel", "split_matrices_list_kidney.rds"))

# Retrieve & prepare metadatas
sj_metadata <- setbp1_marvel$sj.metadata

sample_metadata <- setbp1_marvel$sample.metadata
sample_metadata$num_sjs <- diff(sj_counts@p)
sample_metadata$num_genes <- diff(gene_counts@p)
sample_metadata$num_sjs_genes <- diff(sj_counts@p) / diff(gene_counts@p)

# Plotting
genes_per_cell_violin <- ggplot(
  sample_metadata,
  aes(
    x = cell_type, y = num_genes,
    fill = cell_type, alpha = 0.7
  )
) +
  geom_violin() +
  scale_fill_manual(values = cell_type_colors) +
  guides(fill = "none", alpha = "none") +
  geom_boxplot(width = 0.1) +
  theme_minimal(base_size = 14) +
  xlab("Cell Type") +
  ylab("Genes per Cell") +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 0.95)
  )

sjs_per_cell_violin <- ggplot(
  sample_metadata,
  aes(
    x = cell_type, y = num_sjs,
    fill = cell_type, alpha = 0.7
  )
) +
  geom_violin() +
  scale_fill_manual(values = cell_type_colors) +
  guides(fill = "none", alpha = "none") +
  geom_boxplot(width = 0.1) +
  theme_minimal(base_size = 14) +
  xlab("Cell Type") +
  ylab("SJs per Cell") +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 0.95)
  )

sjs_genes_per_cell_violin <- ggplot(
  sample_metadata,
  aes(
    x = cell_type, y = num_sjs_genes,
    fill = cell_type, alpha = 0.7
  )
) +
  geom_violin() +
  scale_fill_manual(values = cell_type_colors) +
  guides(fill = "none", alpha = "none") +
  geom_boxplot(width = 0.1) +
  theme_minimal(base_size = 14) +
  xlab("Cell Type") +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 0.95)
  ) +
  ylab("SJs per Cell/\nGenes per Cell")

figure <- plot_grid(genes_per_cell_violin,
                       sjs_per_cell_violin,
                       sjs_genes_per_cell_violin,
                       ncol = 1,
                       labels = c("A", "B", "C")
)

# Save figure
png(here::here("results", "kidney_outputs", "genes_sjs_percell_violin.png"),
    width = 8, height = 12, units = "in", res = 300
)
figure
dev.off()

# Count total number of sjs detected in each cell type
all_expressed_genes <- rownames(gene_counts)[rowSums(gene_counts) > 0]

sj_per_gene <- matrix(NA, nrow = length(all_expressed_genes), ncol = 18)
rownames(sj_per_gene) <- sort(all_expressed_genes)

column <- 0

for (i in split_matrices_list) {
  column <- column + 1
  table <- get_sjs_per_gene(i)
  indices <- match(rownames(sj_per_gene), names(table))
  sj_per_gene[, column] <- ifelse(!is.na(indices), table[indices], 0)
}

colnames(sj_per_gene) <-
  names(split_matrices_list)

write.csv(sj_per_gene, here::here("results", "tables", "sj_per_gene_kidney.csv"))

# Get session information
sessionInfo()
