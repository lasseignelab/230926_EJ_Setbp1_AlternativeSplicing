# calc_complexity_all_cells.R - Emma Jones
# This script is for calculating splicing complexity based on dividing number of
# splice junctions divided by number of genes detected for all cells.
# It is optimized to run as a submitted job on slurm

# load packages
suppressPackageStartupMessages({
  # Load Lasseigne Lab standard packages
  library(tidyverse)
  
  # Load MARVEL package
  library(MARVEL)
  
  # Load adjunct MARVEL packages
  library(ggnewscale)
  library(ggrepel)
  library(reshape2)
  library(plyr)
  library(stringr)
  library(textclean)
  library(Matrix)
  library(data.table)
  library(gridExtra)
})

print("loaded packages")

# source functions
source(here::here("src", "marvel", "functions.R"))

# Load MARVEL object
setbp1_marvel <- read_rds(here::here(
  "data", "marvel",
  "setbp1_marvel_aligned.rds"
))

print("loaded in setbp1 data")

# Retrieve sample metadata
sample_metadata <- setbp1_marvel$sample.metadata

# pull out dfs
sj_counts <- setbp1_marvel[["sj.count.matrix"]]
gene_counts <- setbp1_marvel[["gene.count.matrix"]]

calc_avg_sjs_per_gene <- function(cell_id) {
  sjs_detected <- sum(sj_counts[, cell_id] > 0)
  genes_detected <- sum(gene_counts[, cell_id] > 0)
  complexity_score <- sjs_detected / genes_detected
}


print("running function...")

splicing_complexity <- lapply(sample_metadata$cell.id, calc_avg_sjs_per_gene)

print("function finished, saving now...")

# save subset of splice junctions psi sparse matrix
write_csv(splicing_complexity, "data/marvel/complexity_per_cell.csv",
          quote = FALSE)

print("csv saved")
