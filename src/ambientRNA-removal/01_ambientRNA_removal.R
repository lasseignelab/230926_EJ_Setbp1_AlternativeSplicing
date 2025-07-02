### Ambient RNA removal
# Author: Tabea M. Soelter
# Date: July 1st, 2025

## Goal: Remove ambient RNA for downstream analyses.

## Reproducibility:
# * GitHub: lasseignelab/230926_EJ_Setbp1_AlternativeSplicing
# * Docker: tsoelter/sn-ml-drug-repurposing
#       * Version: 1.0.0

## Data:
# STAR Solo outputs
#  * Name: N/A
#  * Location: data/star/

## Analysis Plan:
#  * Load necessary packages
#  * Load data from STAR Solo
#  * Make Soup Channel object
#  * Create & prepare Seurat object
#  * Profile the soup & remove ambient RNA
#  * Save filtered objects

## Analysis:
# Set-up
args <- R.utils::commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SoupX)
  library(DropletUtils)
})

wd <- args[1]
setwd(wd)

source(file.path(wd, "src", "functions_soelter.R"))

output_path <- file.path(
        wd,
        "data",
        "soupX"
      )

# Ambient RNA removal
if (dir.exists(output_path)) {
  print("SoupX output directory already exists.")
} else {
  remove_ambient_rna(
    inputs =
      file.path(
        wd,
        "data",
        "star"
      ),
    outputs =
      output_path,
    plots =
      file.path(
        wd,
        "results",
        "ambientRNA_removal"
      )
  )
}

# Session info (printed to .out file)
sessionInfo()
