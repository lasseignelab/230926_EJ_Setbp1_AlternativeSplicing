# calc_PSI_all_SJs.R - Emma Jones
# This script is for making a psi_matrix for all splice junctions.
# It is optimized to run as an array job on slurm

# enable usage of args
args <- R.utils::commandArgs(trailingOnly = TRUE)
print("enabled args usage")

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
  
  # Load this for data wrangling ease
  library(MatrixExtra)
  # Disable this option to use the following indexing method of checking values
  options("MatrixExtra.quick_show" = FALSE)
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

# get order of genes
gene_order <- setbp1_marvel[["sj.metadata"]][["gene_short_name.start"]]

# make empty sparse matrix for expanded gene counts
expanded_gene_counts <-
  emptySparse(nrow = nrow(sj_counts), ncol = ncol(sj_counts))

# repeat gene values for each splice junction
expanded_gene_counts <- gene_counts[gene_order, ]

print("expanded gene counts")

# rename so no rows have repeated names
rownames(expanded_gene_counts) <- rownames(sj_counts)

# select sj_list
chr <- args[1]
sj_names <- rownames(sj_counts)
sj_subset <- sj_names[grep(chr, sj_names)]

print("splice junctions subsetted")

print("running for loop...")

# run custom function on all sjs
subset_psi_matrix <- make_psi_matrix(sj_subset)

print("for loop finished")

# save subset of splice junctions psi sparse matrix
write_rds(subset_psi_matrix,
          paste0("data/marvel/",
                 paste0(sub(".$", "", chr), "_psi_matrix.Rds")))

print("psi matrix saved")