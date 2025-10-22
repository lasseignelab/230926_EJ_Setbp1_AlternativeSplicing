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

source(file.path(wd, "src", "marvel", "functions.R"))
source(file.path(wd, "src", "marvel", "PlotSJPosition_modification.R"))
source(file.path(wd, "src", "functions_soelter.R"))

# Load data
setbp1_marvel <- read_rds(file.path(wd, 
  "data", "marvel",
  "setbp1_kidney_marvel_aligned.rds"
))

load(file.path(wd, "data", "marvel", "significant_tables_kidney.RData"))

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

name_mapping <- setNames(cell_types, to_col_format(cell_types))

colnames(gene_lists_mat) <- name_mapping[colnames(gene_lists_mat)]

gene_lists_comb_mat <- make_comb_mat(gene_lists_mat)

# Draw & save UpSet plot
comb_degrees <- comb_degree(gene_lists_comb_mat)
color_indices <- rank(comb_degrees)
n_combinations <- length(comb_degrees)

# draw UpSet plot (ComplexHeatmap)
png(file.path(wd, "results", "upset_plots", "sig_gene_overlap_kidney.png"),
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
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 22.04.3 LTS

# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

# locale:
# [1] C

# time zone: Etc/UTC
# tzcode source: system (glibc)

# attached base packages:
# [1] stats4    grid      stats     graphics  grDevices utils     datasets 
# [8] methods   base     

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
# [34] ggplot2_3.5.0         tidyverse_2.0.0       lintr_3.1.0          
# [37] styler_1.10.2         here_1.0.1            

# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      jsonlite_1.8.7          shape_1.4.6            
#   [4] magrittr_2.0.3          farver_2.1.1            fs_1.6.3               
#   [7] GlobalOptions_0.1.2     zlibbioc_1.46.0         vctrs_0.6.5            
#  [10] memoise_2.0.1           RCurl_1.98-1.12         ggtree_3.8.2           
#  [13] qdapRegex_0.7.8         gridGraphics_0.5-1      desc_1.4.2             
#  [16] cachem_1.0.8            igraph_1.5.1            lifecycle_1.0.4        
#  [19] iterators_1.0.14        pkgconfig_2.0.3         gson_0.1.0             
#  [22] R6_2.5.1                fastmap_1.1.1           GenomeInfoDbData_1.2.10
#  [25] clue_0.3-65             digest_0.6.33           aplot_0.2.2            
#  [28] enrichplot_1.20.3       colorspace_2.1-0        ps_1.7.5               
#  [31] rprojroot_2.0.3         RSQLite_2.3.1           fansi_1.0.6            
#  [34] timechange_0.2.0        httr_1.4.7              polyclip_1.10-6        
#  [37] compiler_4.3.1          remotes_2.4.2.1         bit64_4.0.5            
#  [40] withr_3.0.0             doParallel_1.0.17       downloader_0.4         
#  [43] backports_1.4.1         BiocParallel_1.34.2     DBI_1.1.3              
#  [46] ggforce_0.4.1           R.utils_2.12.2          MASS_7.3-60            
#  [49] rjson_0.2.21            HDO.db_0.99.1           tools_4.3.1            
#  [52] scatterpie_0.2.1        ape_5.7-1               R.oo_1.25.0            
#  [55] glue_1.7.0              callr_3.7.3             nlme_3.1-163           
#  [58] R.cache_0.16.0          GOSemSim_2.26.1         shadowtext_0.1.2       
#  [61] cluster_2.1.4           fgsea_1.26.0            generics_0.1.3         
#  [64] gtable_0.3.4            tzdb_0.4.0              R.methodsS3_1.8.2      
#  [67] hms_1.1.3               tidygraph_1.2.3         xml2_1.3.5             
#  [70] utf8_1.2.4              XVector_0.40.0          foreach_1.5.2          
#  [73] pillar_1.9.0            yulab.utils_0.1.0       circlize_0.4.15        
#  [76] splines_4.3.1           tweenr_2.0.2            treeio_1.24.3          
#  [79] lattice_0.21-8          bit_4.0.5               tidyselect_1.2.1       
#  [82] GO.db_3.17.0            Biostrings_2.68.1       graphlayouts_1.0.1     
#  [85] matrixStats_1.0.0       rex_1.2.1               stringi_1.8.1          
#  [88] lazyeval_0.2.2          ggfun_0.1.3             codetools_0.2-19       
#  [91] ggraph_2.1.0            qvalue_2.32.0           ggplotify_0.1.2        
#  [94] cli_3.6.2               munsell_0.5.1           processx_3.8.2         
#  [97] Rcpp_1.0.12             png_0.1-8               parallel_4.3.1         
# [100] blob_1.2.4              DOSE_3.26.2             bitops_1.0-7           
# [103] tidytree_0.4.5          cyclocomp_1.1.1         scales_1.3.0           
# [106] crayon_1.5.2            GetoptLong_1.0.5        rlang_1.1.3            
# [109] cowplot_1.1.1           fastmatch_1.1-4         KEGGREST_1.40.1        
