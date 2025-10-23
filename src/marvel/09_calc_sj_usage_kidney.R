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
# R version 4.3.1 (2023-06-16)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 22.04.3 LTS
#
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0
#
# locale:
# [1] C
#
# time zone: Etc/UTC
# tzcode source: system (glibc)
#
# attached base packages:
# [1] stats4    grid      stats     graphics  grDevices utils     datasets 
# [8] methods   base     
#
# other attached packages:
#  [1] gridExtra_2.3         data.table_1.14.8     Matrix_1.6-1.1       
#  [4] wiggleplotr_1.24.0    GenomicRanges_1.52.1  GenomeInfoDb_1.36.4  
#  [7] gtools_3.9.4          org.Mm.eg.db_3.17.0   org.Hs.eg.db_3.17.0  
# [10] clusterProfiler_4.8.3 AnnotationDbi_1.62.2  IRanges_2.34.1       
# [13] S4Vectors_0.38.2      Biobase_2.60.0        BiocGenerics_0.46.0  
# [16] textclean_0.9.3       plyr_1.8.9            reshape2_1.4.4       
# [19] ggrepel_0.9.5         ggnewscale_0.4.9      MARVEL_2.0.5         
# [22] viridis_0.6.4         viridisLite_0.4.2     ComplexHeatmap_2.16.0
# [25] patchwork_1.1.3       lubridate_1.9.3       forcats_1.0.0        
# [28] stringr_1.5.1         dplyr_1.1.4           purrr_1.0.2          
# [31] readr_2.1.4           tidyr_1.3.0           tibble_3.2.1         
# [34] ggplot2_3.5.0         tidyverse_2.0.0      
#
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      jsonlite_1.8.7          shape_1.4.6            
#   [4] magrittr_2.0.3          farver_2.1.1            GlobalOptions_0.1.2    
#   [7] fs_1.6.3                zlibbioc_1.46.0         vctrs_0.6.5            
#  [10] memoise_2.0.1           RCurl_1.98-1.12         ggtree_3.8.2           
#  [13] qdapRegex_0.7.8         gridGraphics_0.5-1      cachem_1.0.8           
#  [16] igraph_1.5.1            lifecycle_1.0.4         iterators_1.0.14       
#  [19] pkgconfig_2.0.3         R6_2.5.1                fastmap_1.1.1          
#  [22] gson_0.1.0              GenomeInfoDbData_1.2.10 clue_0.3-65            
#  [25] digest_0.6.33           aplot_0.2.2             enrichplot_1.20.3      
#  [28] colorspace_2.1-0        RSQLite_2.3.1           fansi_1.0.6            
#  [31] timechange_0.2.0        httr_1.4.7              polyclip_1.10-6        
#  [34] compiler_4.3.1          bit64_4.0.5             withr_3.0.0            
#  [37] doParallel_1.0.17       downloader_0.4          BiocParallel_1.34.2    
#  [40] DBI_1.1.3               ggforce_0.4.1           R.utils_2.12.2         
#  [43] MASS_7.3-60             rjson_0.2.21            HDO.db_0.99.1          
#  [46] tools_4.3.1             ape_5.7-1               scatterpie_0.2.1       
#  [49] R.oo_1.25.0             glue_1.7.0              nlme_3.1-163           
#  [52] GOSemSim_2.26.1         shadowtext_0.1.2        cluster_2.1.4          
#  [55] fgsea_1.26.0            generics_0.1.3          gtable_0.3.4           
#  [58] tzdb_0.4.0              R.methodsS3_1.8.2       hms_1.1.3              
#  [61] tidygraph_1.2.3         utf8_1.2.4              XVector_0.40.0         
#  [64] foreach_1.5.2           pillar_1.9.0            yulab.utils_0.1.0      
#  [67] circlize_0.4.15         splines_4.3.1           tweenr_2.0.2           
#  [70] treeio_1.24.3           lattice_0.21-8          bit_4.0.5              
#  [73] tidyselect_1.2.1        GO.db_3.17.0            Biostrings_2.68.1      
#  [76] graphlayouts_1.0.1      matrixStats_1.0.0       stringi_1.8.1          
#  [79] lazyeval_0.2.2          ggfun_0.1.3             codetools_0.2-19       
#  [82] ggraph_2.1.0            qvalue_2.32.0           ggplotify_0.1.2        
#  [85] cli_3.6.2               munsell_0.5.1           Rcpp_1.0.12            
#  [88] png_0.1-8               parallel_4.3.1          blob_1.2.4             
#  [91] DOSE_3.26.2             bitops_1.0-7            tidytree_0.4.5         
#  [94] scales_1.3.0            crayon_1.5.2            GetoptLong_1.0.5       
#  [97] rlang_1.1.3             cowplot_1.1.1           fastmatch_1.1-4        
# [100] KEGGREST_1.40.1
