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
# [1] stats4    stats     graphics  grDevices utils     datasets  methods  
# [8] base     

# other attached packages:
#  [1] gridExtra_2.3         data.table_1.14.8     Matrix_1.6-1.1       
#  [4] wiggleplotr_1.24.0    GenomicRanges_1.52.1  GenomeInfoDb_1.36.4  
#  [7] gtools_3.9.4          org.Mm.eg.db_3.17.0   org.Hs.eg.db_3.17.0  
# [10] clusterProfiler_4.8.3 AnnotationDbi_1.62.2  IRanges_2.34.1       
# [13] S4Vectors_0.38.2      Biobase_2.60.0        BiocGenerics_0.46.0  
# [16] textclean_0.9.3       plyr_1.8.9            reshape2_1.4.4       
# [19] ggrepel_0.9.5         ggnewscale_0.4.9      MARVEL_2.0.5         
# [22] patchwork_1.1.3       lubridate_1.9.3       forcats_1.0.0        
# [25] stringr_1.5.1         dplyr_1.1.4           purrr_1.0.2          
# [28] readr_2.1.4           tidyr_1.3.0           tibble_3.2.1         
# [31] ggplot2_3.5.0         tidyverse_2.0.0      

# loaded via a namespace (and not attached):
#  [1] DBI_1.1.3               qdapRegex_0.7.8         bitops_1.0-7           
#  [4] gson_0.1.0              shadowtext_0.1.2        rlang_1.1.3            
#  [7] magrittr_2.0.3          DOSE_3.26.2             compiler_4.3.1         
# [10] RSQLite_2.3.1           png_0.1-8               vctrs_0.6.5            
# [13] pkgconfig_2.0.3         crayon_1.5.2            fastmap_1.1.1          
# [16] XVector_0.40.0          ggraph_2.1.0            utf8_1.2.4             
# [19] HDO.db_0.99.1           tzdb_0.4.0              enrichplot_1.20.3      
# [22] bit_4.0.5               zlibbioc_1.46.0         cachem_1.0.8           
# [25] aplot_0.2.2             jsonlite_1.8.7          blob_1.2.4             
# [28] BiocParallel_1.34.2     tweenr_2.0.2            parallel_4.3.1         
# [31] R6_2.5.1                stringi_1.8.1           RColorBrewer_1.1-3     
# [34] GOSemSim_2.26.1         Rcpp_1.0.12             downloader_0.4         
# [37] R.utils_2.12.2          splines_4.3.1           igraph_1.5.1           
# [40] timechange_0.2.0        tidyselect_1.2.1        qvalue_2.32.0          
# [43] viridis_0.6.4           codetools_0.2-19        lattice_0.21-8         
# [46] treeio_1.24.3           withr_3.0.0             KEGGREST_1.40.1        
# [49] gridGraphics_0.5-1      scatterpie_0.2.1        polyclip_1.10-6        
# [52] Biostrings_2.68.1       ggtree_3.8.2            pillar_1.9.0           
# [55] ggfun_0.1.3             generics_0.1.3          RCurl_1.98-1.12        
# [58] hms_1.1.3               tidytree_0.4.5          munsell_0.5.1          
# [61] scales_1.3.0            glue_1.7.0              lazyeval_0.2.2         
# [64] tools_4.3.1             fgsea_1.26.0            fs_1.6.3               
# [67] graphlayouts_1.0.1      fastmatch_1.1-4         tidygraph_1.2.3        
# [70] cowplot_1.1.1           grid_4.3.1              ape_5.7-1              
# [73] colorspace_2.1-0        nlme_3.1-163            GenomeInfoDbData_1.2.10
# [76] ggforce_0.4.1           cli_3.6.2               fansi_1.0.6            
# [79] viridisLite_0.4.2       gtable_0.3.4            R.methodsS3_1.8.2      
# [82] yulab.utils_0.1.0       digest_0.6.33           ggplotify_0.1.2        
# [85] farver_2.1.1            memoise_2.0.1           R.oo_1.25.0            
# [88] lifecycle_1.0.4         httr_1.4.7              GO.db_3.17.0           
# [91] bit64_4.0.5             MASS_7.3-60            
