### MARVEL analysis on S858R kidney data
# Author: Tabea M. Soelter (OG code by Emma Jones)
# Date: August 18th, 2025

## Goal: Run differential analysis on kidney S858R data with MARVEL.

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
#  * Load MARVEL objects
#  * Cell-type-specific MARVEL analysis between S858R and WT
#  * Save MARVEL results

## Analysis:
# Setup
args <- R.utils::commandArgs(trailingOnly = TRUE)

wd <- args[1]
setwd(wd)

set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
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
source(file.path(wd, "src", "marvel", "PlotSJPosition_modification.R"))
source(file.path(wd, "src", "functions_soelter.R"))

# Load MARVEL objects
setbp1_marvel <- read_rds(file.path(wd, 
  "data", "marvel",
  "setbp1_kidney_marvel_aligned.rds"
))

sample_metadata <- setbp1_marvel$sample.metadata

# Extract & save cell-type-specific MARVEL objects
cell_types <- setbp1_marvel[["sample.metadata"]][["cell_type"]] %>%
  unique()

for (cell_type in cell_types) {
  celltype_name <- gsub(" ", "_", cell_type) %>%
    tolower()
  
  object_name <- paste0(celltype_name, "_marvel")

  assign(object_name, run_marvel_cell_type_kidney(
    marvel_object = setbp1_marvel,
    cell_type = cell_type,
    results_path = file.path(wd, "data", "marvel")
  ))
}

# Set reference group & group oi for differential MARVEL analysis
all_types_marvel <- setbp1_marvel

index_1 <- which(sample_metadata$seq_folder == "wildtype")
cell_ids_1 <- sample_metadata[index_1, "cell.id"]

index_2 <- which(sample_metadata$seq_folder == "mutant")
cell_ids_2 <- sample_metadata[index_2, "cell.id"]

# Explore % of cells expressing genes
all_types_marvel <- PlotPctExprCells.Genes.10x(
  MarvelObject = all_types_marvel,
  cell.group.g1 = cell_ids_1,
  cell.group.g2 = cell_ids_2,
  min.pct.cells = 5
)

# Explore % of cells expressing junctions
all_types_marvel <- PlotPctExprCells.SJ.10x(
  MarvelObject = all_types_marvel,
  cell.group.g1 = cell_ids_1,
  cell.group.g2 = cell_ids_2,
  min.pct.cells.genes = 5,
  min.pct.cells.sj = 5,
  downsample = TRUE,
  downsample.pct.sj = 10
)

# Differential Splicing Analysis
all_types_marvel <- CompareValues.SJ.10x(
  MarvelObject = all_types_marvel,
  cell.group.g1 = cell_ids_1,
  cell.group.g2 = cell_ids_2,
  min.pct.cells.genes = 5,
  min.pct.cells.sj = 5,
  min.gene.norm = 1,
  seed = 1,
  n.iterations = 100,
  downsample = TRUE,
  show.progress = TRUE
)

# Differential Gene Analysis
all_types_marvel <- CompareValues.Genes.10x(
  MarvelObject = all_types_marvel,
  show.progress = TRUE
)

# Make volcano plot
all_types_marvel <- PlotDEValues.SJ.10x(
  MarvelObject = all_types_marvel,
  pval = 0.05,
  delta = 1,
  min.gene.norm = 1,
  anno = FALSE
)

# Assign kinds of iso-switching
all_types_marvel <- IsoSwitch.10x(
  MarvelObject = all_types_marvel,
  pval.sj = 0.05,
  delta.sj = 1,
  min.gene.norm = 1,
  pval.adj.gene = 0.05,
  log2fc.gene = 0.5
)

# Save all results
saveRDS(all_types_marvel,
        file = file.path(wd, 
          "data", "marvel",
          paste0("All_Cell_Types_kidney_marvel_object.rds")
        )
)

# Extract and save significant results
object_list <- append(cell_types, "all_types")

for (object in object_list) {
  modified_object <- gsub(" ", "_", object) %>% tolower()
  object_name <- paste0(modified_object, "_marvel")
  sig_table_name <- paste0(modified_object, "_sig_table")
  marvel_object <- get(object_name)
  assign(sig_table_name, marvel_object[["SJ.Gene.Cor"]][["Data"]])
}

save(
  list = ls(pattern = "_sig_table"),
  file = file.path(wd, "data", "marvel", "significant_tables_kidney.RData")
)

# Get session information
sessionInfo()
