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
# [1] stats4    stats     graphics  grDevices utils     datasets  methods  
# [8] base     
#
# other attached packages:
#  [1] MatrixExtra_0.1.14    gridExtra_2.3         Matrix_1.6-1.1       
#  [4] wiggleplotr_1.24.0    GenomicRanges_1.52.1  GenomeInfoDb_1.36.4  
#  [7] gtools_3.9.4          org.Mm.eg.db_3.17.0   org.Hs.eg.db_3.17.0  
# [10] clusterProfiler_4.8.3 AnnotationDbi_1.62.2  IRanges_2.34.1       
# [13] S4Vectors_0.38.2      Biobase_2.60.0        BiocGenerics_0.46.0  
# [16] textclean_0.9.3       plyr_1.8.9            reshape2_1.4.4       
# [19] ggrepel_0.9.5         ggnewscale_0.4.9      MARVEL_2.0.5         
# [22] presto_1.0.0          data.table_1.14.8     harmony_1.1.0        
# [25] Rcpp_1.0.12           patchwork_1.1.3       Seurat_5.0.0         
# [28] SeuratObject_5.0.0    sp_2.1-1              lubridate_1.9.3      
# [31] forcats_1.0.0         stringr_1.5.1         dplyr_1.1.4          
# [34] purrr_1.0.2           readr_2.1.4           tidyr_1.3.0          
# [37] tibble_3.2.1          ggplot2_3.5.0         tidyverse_2.0.0      
#
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.21        splines_4.3.1           later_1.3.1            
#   [4] ggplotify_0.1.2         bitops_1.0-7            R.oo_1.25.0            
#   [7] polyclip_1.10-6         fastDummies_1.7.3       lifecycle_1.0.4        
#  [10] rprojroot_2.0.3         globals_0.16.2          lattice_0.21-8         
#  [13] MASS_7.3-60             magrittr_2.0.3          plotly_4.10.3          
#  [16] httpuv_1.6.12           sctransform_0.4.1       spam_2.10-0            
#  [19] spatstat.sparse_3.0-3   reticulate_1.34.0       cowplot_1.1.1          
#  [22] pbapply_1.7-2           DBI_1.1.3               RColorBrewer_1.1-3     
#  [25] abind_1.4-5             zlibbioc_1.46.0         Rtsne_0.16             
#  [28] R.utils_2.12.2          ggraph_2.1.0            RCurl_1.98-1.12        
#  [31] yulab.utils_0.1.0       tweenr_2.0.2            float_0.3-2            
#  [34] GenomeInfoDbData_1.2.10 enrichplot_1.20.3       irlba_2.3.5.1          
#  [37] listenv_0.9.0           spatstat.utils_3.0-4    tidytree_0.4.5         
#  [40] qdapRegex_0.7.8         goftest_1.2-3           RSpectra_0.16-1        
#  [43] spatstat.random_3.2-1   fitdistrplus_1.1-11     parallelly_1.36.0      
#  [46] leiden_0.4.3            codetools_0.2-19        ggforce_0.4.1          
#  [49] DOSE_3.26.2             tidyselect_1.2.1        aplot_0.2.2            
#  [52] farver_2.1.1            viridis_0.6.4           matrixStats_1.0.0      
#  [55] spatstat.explore_3.2-5  jsonlite_1.8.7          tidygraph_1.2.3        
#  [58] ellipsis_0.3.2          progressr_0.14.0        ggridges_0.5.4         
#  [61] survival_3.5-7          tools_4.3.1             treeio_1.24.3          
#  [64] ica_1.0-3               glue_1.7.0              here_1.0.1             
#  [67] qvalue_2.32.0           withr_3.0.0             fastmap_1.1.1          
#  [70] fansi_1.0.6             digest_0.6.33           gridGraphics_0.5-1     
#  [73] timechange_0.2.0        R6_2.5.1                mime_0.12              
#  [76] colorspace_2.1-0        scattermore_1.2         GO.db_3.17.0           
#  [79] tensor_1.5              spatstat.data_3.0-3     RSQLite_2.3.1          
#  [82] R.methodsS3_1.8.2       RhpcBLASctl_0.23-42     utf8_1.2.4             
#  [85] generics_0.1.3          graphlayouts_1.0.1      httr_1.4.7             
#  [88] htmlwidgets_1.6.2       scatterpie_0.2.1        uwot_0.1.16            
#  [91] pkgconfig_2.0.3         gtable_0.3.4            blob_1.2.4             
#  [94] lmtest_0.9-40           XVector_0.40.0          shadowtext_0.1.2       
#  [97] htmltools_0.5.6.1       dotCall64_1.1-0         fgsea_1.26.0           
# [100] scales_1.3.0            png_0.1-8               ggfun_0.1.3            
# [103] tzdb_0.4.0              nlme_3.1-163            zoo_1.8-12             
# [106] cachem_1.0.8            KernSmooth_2.23-22      parallel_4.3.1         
# [109] miniUI_0.1.1.1          HDO.db_0.99.1           pillar_1.9.0           
# [112] grid_4.3.1              vctrs_0.6.5             RANN_2.6.1             
# [115] promises_1.2.1          xtable_1.8-4            cluster_2.1.4          
# [118] cli_3.6.2               compiler_4.3.1          rlang_1.1.3            
# [121] crayon_1.5.2            future.apply_1.11.0     fs_1.6.3               
# [124] stringi_1.8.1           viridisLite_0.4.2       deldir_1.0-9           
# [127] BiocParallel_1.34.2     munsell_0.5.1           Biostrings_2.68.1      
# [130] lazyeval_0.2.2          spatstat.geom_3.2-7     GOSemSim_2.26.1        
# [133] RcppHNSW_0.5.0          hms_1.1.3               bit64_4.0.5            
# [136] future_1.33.0           KEGGREST_1.40.1         shiny_1.7.5.1          
# [139] ROCR_1.0-11             igraph_1.5.1            memoise_2.0.1          
# [142] ggtree_3.8.2            fastmatch_1.1-4         bit_4.0.5              
# [145] downloader_0.4          gson_0.1.0              ape_5.7-1              
