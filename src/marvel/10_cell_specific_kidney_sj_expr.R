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
#  [1] ggpmisc_0.5.5          ggpp_0.5.6             rtracklayer_1.60.1    
#  [4] GenomicFeatures_1.52.2 ggridges_0.5.4         gridExtra_2.3         
#  [7] data.table_1.14.8      Matrix_1.6-1.1         wiggleplotr_1.24.0    
# [10] GenomicRanges_1.52.1   GenomeInfoDb_1.36.4    gtools_3.9.4          
# [13] org.Mm.eg.db_3.17.0    org.Hs.eg.db_3.17.0    clusterProfiler_4.8.3 
# [16] AnnotationDbi_1.62.2   IRanges_2.34.1         S4Vectors_0.38.2      
# [19] Biobase_2.60.0         BiocGenerics_0.46.0    textclean_0.9.3       
# [22] plyr_1.8.9             reshape2_1.4.4         ggrepel_0.9.5         
# [25] ggnewscale_0.4.9       MARVEL_2.0.5           cowplot_1.1.1         
# [28] viridis_0.6.4          viridisLite_0.4.2      ComplexHeatmap_2.16.0 
# [31] patchwork_1.1.3        lubridate_1.9.3        forcats_1.0.0         
# [34] stringr_1.5.1          dplyr_1.1.4            purrr_1.0.2           
# [37] readr_2.1.4            tidyr_1.3.0            tibble_3.2.1          
# [40] ggplot2_3.5.0          tidyverse_2.0.0       
#
# loaded via a namespace (and not attached):
#   [1] splines_4.3.1               BiocIO_1.10.0              
#   [3] bitops_1.0-7                ggplotify_0.1.2            
#   [5] filelock_1.0.2              R.oo_1.25.0                
#   [7] polyclip_1.10-6             XML_3.99-0.14              
#   [9] lifecycle_1.0.4             rprojroot_2.0.3            
#  [11] doParallel_1.0.17           lattice_0.21-8             
#  [13] MASS_7.3-60                 magrittr_2.0.3             
#  [15] yaml_2.3.7                  DBI_1.1.3                  
#  [17] RColorBrewer_1.1-3          abind_1.4-5                
#  [19] zlibbioc_1.46.0             R.utils_2.12.2             
#  [21] ggraph_2.1.0                RCurl_1.98-1.12            
#  [23] yulab.utils_0.1.0           tweenr_2.0.2               
#  [25] rappdirs_0.3.3              circlize_0.4.15            
#  [27] GenomeInfoDbData_1.2.10     enrichplot_1.20.3          
#  [29] tidytree_0.4.5              qdapRegex_0.7.8            
#  [31] MatrixModels_0.5-2          codetools_0.2-19           
#  [33] DelayedArray_0.26.7         DOSE_3.26.2                
#  [35] xml2_1.3.5                  ggforce_0.4.1              
#  [37] tidyselect_1.2.1            shape_1.4.6                
#  [39] aplot_0.2.2                 farver_2.1.1               
#  [41] matrixStats_1.0.0           BiocFileCache_2.8.0        
#  [43] GenomicAlignments_1.36.0    jsonlite_1.8.7             
#  [45] GetoptLong_1.0.5            tidygraph_1.2.3            
#  [47] survival_3.5-7              iterators_1.0.14           
#  [49] foreach_1.5.2               tools_4.3.1                
#  [51] progress_1.2.2              treeio_1.24.3              
#  [53] Rcpp_1.0.12                 glue_1.7.0                 
#  [55] here_1.0.1                  qvalue_2.32.0              
#  [57] MatrixGenerics_1.12.3       withr_3.0.0                
#  [59] fastmap_1.1.1               fansi_1.0.6                
#  [61] SparseM_1.81                digest_0.6.33              
#  [63] timechange_0.2.0            R6_2.5.1                   
#  [65] gridGraphics_0.5-1          colorspace_2.1-0           
#  [67] GO.db_3.17.0                biomaRt_2.56.1             
#  [69] RSQLite_2.3.1               R.methodsS3_1.8.2          
#  [71] utf8_1.2.4                  generics_0.1.3             
#  [73] prettyunits_1.2.0           graphlayouts_1.0.1         
#  [75] httr_1.4.7                  S4Arrays_1.0.6             
#  [77] scatterpie_0.2.1            pkgconfig_2.0.3            
#  [79] gtable_0.3.4                blob_1.2.4                 
#  [81] XVector_0.40.0              shadowtext_0.1.2           
#  [83] fgsea_1.26.0                clue_0.3-65                
#  [85] scales_1.3.0                png_0.1-8                  
#  [87] ggfun_0.1.3                 tzdb_0.4.0                 
#  [89] rjson_0.2.21                nlme_3.1-163               
#  [91] curl_5.1.0                  cachem_1.0.8               
#  [93] GlobalOptions_0.1.2         parallel_4.3.1             
#  [95] HDO.db_0.99.1               restfulr_0.0.15            
#  [97] pillar_1.9.0                vctrs_0.6.5                
#  [99] dbplyr_2.3.4                cluster_2.1.4              
# [101] cli_3.6.2                   compiler_4.3.1             
# [103] Rsamtools_2.16.0            rlang_1.1.3                
# [105] crayon_1.5.2                labeling_0.4.3             
# [107] fs_1.6.3                    stringi_1.8.1              
# [109] BiocParallel_1.34.2         munsell_0.5.1              
# [111] Biostrings_2.68.1           lazyeval_0.2.2             
# [113] GOSemSim_2.26.1             quantreg_5.97              
# [115] hms_1.1.3                   bit64_4.0.5                
# [117] KEGGREST_1.40.1             SummarizedExperiment_1.30.2
# [119] igraph_1.5.
