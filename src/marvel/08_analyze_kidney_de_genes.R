### Differential kidney MARVEL results
# Author: Tabea M. Soelter (OG code by Emma Jones)
# Date: August 20th, 2025

## Goal: Analyze kidney differential MARVEL analysis results

## Reproducibility:
# * GitHub: lasseignelab/230926_EJ_Setbp1_AlternativeSplicing
# * Docker: emmafjones/setbp1_alternative_splicing
#       * Version: 1.1.0

## Data:
# MARVEL object
#  * Name: setbp1_kidney_marvel_aligned.rds
#  * Location: data/marvel/
# MARVEL significant tables
#  * Name: significant_tables_kidney.RData
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
  library(here)
  library(styler)
  library(lintr)
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

source(here::here("src", "marvel", "functions.R"))
source(here::here("src", "marvel", "PlotSJPosition_modification.R"))

# Load data
setbp1_marvel <- read_rds(here::here(
  "data", "marvel",
  "setbp1_kidney_marvel_aligned.rds"
))

load(here::here("data", "marvel", "significant_tables_kidney.RData"))

# Prepare data for UpSet plotting
sample_metadata <- setbp1_marvel$sample.metadata

tables_list <- mget(ls(pattern = "_sig_table"))
tables_list <- tables_list[-1]

gene_lists <- map(tables_list, ~ select(.x, gene_short_name))
gene_lists <- map(gene_lists, ~ flatten(.))

gene_lists_mat <- ComplexHeatmap::list_to_matrix(gene_lists)
colnames(gene_lists_mat) <- gsub("_sig_table", "", colnames(gene_lists_mat))

cell_types <- setbp1_marvel[["sample.metadata"]][["cell_type"]] %>%
  unique()

to_col_format <- function(x) {
  tolower(gsub(" ", "_", x))
}

name_mapping <- setNames(cell_types, to_col_format(cell_types))

colnames(gene_lists_mat) <- name_mapping[colnames(gene_lists_mat)]

gene_lists_comb_mat <- make_comb_mat(gene_lists_mat)

# Draw & save UpSet plot
comb_degrees <- comb_degree(gene_lists_comb_mat)
color_indices <- rank(comb_degrees)
n_combinations <- length(comb_degrees)

# draw UpSet plot (ComplexHeatmap)
png(here::here("results", "upset_plots", "sig_gene_overlap_kidney.png"),
    width = 12, height = 6, units = "in", res = 300
)
gene_lists_upset_plot <- UpSet(
  gene_lists_comb_mat,
  comb_col = viridis(n_combinations)[color_indices],
  top_annotation = upset_top_annotation(
    gene_lists_comb_mat,
    add_numbers = TRUE, annotation_name_gp = gpar(fontface = "bold"),
    axis_param = list(gp = gpar(fontface = "bold")),
    numbers_gp = gpar(fontface = "bold"),
    numbers_rot = 0
  ),
  right_annotation = upset_right_annotation(
    gene_lists_comb_mat,
    add_numbers = TRUE, annotation_name_gp = gpar(fontface = "bold"),
    axis_param = list(gp = gpar(fontface = "bold")),
    numbers_gp = gpar(fontface = "bold")
  ),
  row_names_gp = gpar(fontface = "bold")
)
draw(gene_lists_upset_plot)
dev.off()

# Get session information
sessionInfo()
