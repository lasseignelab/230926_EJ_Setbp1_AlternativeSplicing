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
#  [1] data.table_1.14.8           RColorBrewer_1.1-3         
#  [3] DESeq2_1.40.2               png_0.1-8                  
#  [5] apeglm_1.22.1               pheatmap_1.0.12            
#  [7] SingleCellExperiment_1.22.0 SummarizedExperiment_1.30.2
#  [9] Biobase_2.60.0              GenomicRanges_1.52.1       
# [11] GenomeInfoDb_1.36.4         IRanges_2.34.1             
# [13] MatrixGenerics_1.12.3       matrixStats_1.0.0          
# [15] S4Vectors_0.38.2            BiocGenerics_0.46.0        
# [17] reshape2_1.4.4              edgeR_3.42.4               
# [19] limma_3.56.2                Matrix.utils_0.9.7         
# [21] Matrix_1.6-1.1              Seurat_5.0.0               
# [23] SeuratObject_5.0.0          sp_2.1-1                   
# [25] viridis_0.6.4               viridisLite_0.4.2          
# [27] ComplexHeatmap_2.16.0       patchwork_1.1.3            
# [29] lubridate_1.9.3             forcats_1.0.0              
# [31] stringr_1.5.1               dplyr_1.1.3                
# [33] purrr_1.0.2                 readr_2.1.4                
# [35] tidyr_1.3.0                 tibble_3.2.1               
# [37] ggplot2_3.4.4               tidyverse_2.0.0            
#
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.21        splines_4.3.1           later_1.3.1            
#   [4] bitops_1.0-7            R.oo_1.25.0             polyclip_1.10-6        
#   [7] fastDummies_1.7.3       lifecycle_1.0.4         rprojroot_2.0.3        
#  [10] doParallel_1.0.17       globals_0.16.2          lattice_0.21-8         
#  [13] MASS_7.3-60             magrittr_2.0.3          plotly_4.10.3          
#  [16] httpuv_1.6.12           sctransform_0.4.1       spam_2.10-0            
#  [19] spatstat.sparse_3.0-3   reticulate_1.34.0       cowplot_1.1.1          
#  [22] pbapply_1.7-2           abind_1.4-5             zlibbioc_1.46.0        
#  [25] Rtsne_0.16              R.utils_2.12.2          RCurl_1.98-1.12        
#  [28] circlize_0.4.15         GenomeInfoDbData_1.2.10 ggrepel_0.9.4          
#  [31] irlba_2.3.5.1           listenv_0.9.0           spatstat.utils_3.0-4   
#  [34] goftest_1.2-3           RSpectra_0.16-1         spatstat.random_3.2-1  
#  [37] fitdistrplus_1.1-11     parallelly_1.36.0       leiden_0.4.3           
#  [40] codetools_0.2-19        DelayedArray_0.26.7     tidyselect_1.2.0       
#  [43] shape_1.4.6             farver_2.1.1            spatstat.explore_3.2-5 
#  [46] jsonlite_1.8.7          GetoptLong_1.0.5        ellipsis_0.3.2         
#  [49] progressr_0.14.0        ggridges_0.5.4          survival_3.5-7         
#  [52] iterators_1.0.14        systemfonts_1.0.4       bbmle_1.0.25.1         
#  [55] foreach_1.5.2           tools_4.3.1             ragg_1.2.5             
#  [58] ica_1.0-3               Rcpp_1.0.11             glue_1.6.2             
#  [61] gridExtra_2.3           here_1.0.1              withr_2.5.2            
#  [64] numDeriv_2016.8-1.1     fastmap_1.1.1           fansi_1.0.5            
#  [67] digest_0.6.33           timechange_0.2.0        R6_2.5.1               
#  [70] mime_0.12               textshaping_0.3.6       colorspace_2.1-0       
#  [73] scattermore_1.2         tensor_1.5              spatstat.data_3.0-3    
#  [76] R.methodsS3_1.8.2       utf8_1.2.4              generics_0.1.3         
#  [79] httr_1.4.7              htmlwidgets_1.6.2       S4Arrays_1.0.6         
#  [82] uwot_0.1.16             pkgconfig_2.0.3         gtable_0.3.4           
#  [85] lmtest_0.9-40           XVector_0.40.0          htmltools_0.5.6.1      
#  [88] dotCall64_1.1-0         clue_0.3-65             scales_1.2.1           
#  [91] tzdb_0.4.0              rjson_0.2.21            coda_0.19-4.1          
#  [94] nlme_3.1-163            bdsmatrix_1.3-7         zoo_1.8-12             
#  [97] GlobalOptions_0.1.2     KernSmooth_2.23-22      parallel_4.3.1         
# [100] miniUI_0.1.1.1          pillar_1.9.0            vctrs_0.6.4            
# [103] RANN_2.6.1              promises_1.2.1          xtable_1.8-4           
# [106] cluster_2.1.4           mvtnorm_1.2-3           cli_3.6.1              
# [109] locfit_1.5-9.8          compiler_4.3.1          rlang_1.1.2            
# [112] crayon_1.5.2            grr_0.9.5               future.apply_1.11.0    
# [115] labeling_0.4.3          emdbook_1.3.13          plyr_1.8.9             
# [118] stringi_1.8.1           BiocParallel_1.34.2     deldir_1.0-9           
# [121] munsell_0.5.0           lazyeval_0.2.2          spatstat.geom_3.2-7    
# [124] RcppHNSW_0.5.0          hms_1.1.3               future_1.33.0          
# [127] shiny_1.7.5.1           ROCR_1.0-11             igraph_1.5.1           
