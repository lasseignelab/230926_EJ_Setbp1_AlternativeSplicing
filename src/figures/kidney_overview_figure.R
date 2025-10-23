### Plotting S858R mouse kidney data overview figure
# Author: Tabea M. Soelter (original code modified from Emma Jones)
# Date: October 15th, 2025

## Goal: Create S858R kidney data overview figure.

## Reproducibility:
# * GitHub: lasseignelab/230926_EJ_Setbp1_AlternativeSplicing
# * Docker: emmafjones/setbp1_alternative_splicing
#       * Version: 1.0.6 

## Data:
# Annotated Seurat object
# * Name: annotated_kidney_samples.rds
# * Location: data/seurat/

## Analysis Plan:
# * Load necessary packages 
# * Load data 
# * Plot individual panels
# * Compile overview figure
# * Save overview figure

## Set-up:
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
  library(cowplot)
})

## Analysis:
# Load data
annotated_kidney_samples <- read_rds(
  file = file.path(
    wd,
    "data",
    "seurat",
    "annotated_kidney_samples.rds"
  )
)

# Set color palettes
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

condition_colors <- c(
  `mutant` = "#666666",
  `wildtype` = "#C1C1C1"
)

other_colors <- c("#0346a3", "#076466", "#2a6607")

# Plot panel A
umap_1 <- DimPlot(annotated_kidney_samples,
                  cols = cell_type_colors,
                  label = TRUE,
                  label.box = TRUE,
                  label.color = "white",
                  label.size = 4.75,
                  repel = TRUE,
                  seed = 42
) +
  labs(color = "Cell Type") +
  theme(
    axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  NoLegend()

# Plot panel B
ordered_annotated_kidney_samples <- SetIdent(
  annotated_kidney_samples, 
  value = factor(Idents(annotated_kidney_samples),
                 levels = sort(levels(Idents(annotated_kidney_samples)),
                               decreasing = TRUE)
                 )
  )

dotplot <- DotPlot(ordered_annotated_kidney_samples,
                   cols = c("lightgray", other_colors[2]),
                   features = c(
                     "Aqp1", "Epha7", "Ptger3", "Ikzf2", "Cd247", "Slc22a7",
                     "Slc13a3", "Nphs1", "Slc5a2", "Gucy1a1",  "Flt1", "H2-Aa",
                     "Slc12a3", "Atp6v1b1", "Aqp2",  "Hmx2", "Aqp6", "Cd79a"
                   )
) +
  guides(
    size = guide_legend(title = "Percent\nExpressed"),
    color = guide_colorbar(title = "Average\nExpression")
  ) +
  theme(
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ylab("Cell Type")

# Plot panel C
umap_2 <- DimPlot(annotated_kidney_samples,
                  group.by = "seq_folder",
                  shuffle = TRUE,
                  seed = 123
) +
  theme(
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  scale_color_manual(
    labels = c("S858R+/-", "Wildtype"),
    values = condition_colors
  ) +
  labs(color = "Condition", title = NULL) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  NoLegend()

# Plot panel D
df <- table(
  Idents(annotated_kidney_samples),
  annotated_kidney_samples$seq_folder
) %>% as.data.frame()

df$Var1 <- as.character(df$Var1)
colnames(df) <- c("Cell_Type", "Condition", "Freq")

df3 <- data.frame(
  "Cell_Type" = as.numeric(c(NA, NA)),
  "Condition" = c("S858R+/-", "Wildtype"),
  "Freq" = as.numeric(c(NA, NA))
)


barplot <- ggplot(data = df, aes(
  x = fct_rev(Cell_Type),
  y = Freq, fill = Condition
)) +
  theme_minimal(base_size = 15) +
  geom_col(position = "fill", width = 0.75, show.legend = FALSE) +
  scale_fill_manual(
    labels = c("S858R+/-", "Wildtype"),
    values = condition_colors
  ) +
  theme(
    legend.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text()
  ) +
  coord_flip() +
  xlab("Cell Type") +
  ylab("Proportion") +
  geom_point(data = df3, aes(
    x = Cell_Type,
    y = Freq, color = Condition
  ), size = 4) +
  scale_color_manual(
    labels = c("S858R+/-", "Wildtype"),
    values = c("#666666", "#C1C1C1")
  ) +
  theme(legend.key = element_blank())

# Compile overview figure
figure <- plot_grid(umap_1, dotplot, umap_2, barplot,
                    labels = c("A", "B", "C", "D")
)

png(file.path(wd, "results", "kidney_figures", "kidney_overview_figure.png"),
    width = 17, height = 12, units = "in", res = 300
)
figure
dev.off()

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
# [1] stats     graphics  grDevices utils     datasets  methods   base     
#
# other attached packages:
#  [1] cowplot_1.1.1      presto_1.0.0       data.table_1.14.8  harmony_1.1.0     
#  [5] Rcpp_1.0.11        patchwork_1.1.3    Seurat_5.0.0       SeuratObject_5.0.0
#  [9] sp_2.1-1           lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1     
# [13] dplyr_1.1.3        purrr_1.0.2        readr_2.1.4        tidyr_1.3.0       
# [17] tibble_3.2.1       ggplot2_3.4.4      tidyverse_2.0.0   
#
# loaded via a namespace (and not attached):
#   [1] deldir_1.0-9           pbapply_1.7-2          gridExtra_2.3         
#   [4] rlang_1.1.2            magrittr_2.0.3         RcppAnnoy_0.0.21      
#   [7] spatstat.geom_3.2-7    matrixStats_1.0.0      ggridges_0.5.4        
#  [10] compiler_4.3.1         png_0.1-8              vctrs_0.6.4           
#  [13] reshape2_1.4.4         pkgconfig_2.0.3        fastmap_1.1.1         
#  [16] ellipsis_0.3.2         labeling_0.4.3         utf8_1.2.4            
#  [19] promises_1.2.1         tzdb_0.4.0             jsonlite_1.8.7        
#  [22] goftest_1.2-3          later_1.3.1            spatstat.utils_3.0-4  
#  [25] irlba_2.3.5.1          parallel_4.3.1         cluster_2.1.4         
#  [28] R6_2.5.1               ica_1.0-3              spatstat.data_3.0-3   
#  [31] stringi_1.8.1          RColorBrewer_1.1-3     reticulate_1.34.0     
#  [34] parallelly_1.36.0      lmtest_0.9-40          scattermore_1.2       
#  [37] tensor_1.5             future.apply_1.11.0    zoo_1.8-12            
#  [40] R.utils_2.12.2         sctransform_0.4.1      httpuv_1.6.12         
#  [43] Matrix_1.6-1.1         splines_4.3.1          igraph_1.5.1          
#  [46] timechange_0.2.0       tidyselect_1.2.0       abind_1.4-5           
#  [49] spatstat.random_3.2-1  codetools_0.2-19       miniUI_0.1.1.1        
#  [52] spatstat.explore_3.2-5 listenv_0.9.0          lattice_0.21-8        
#  [55] plyr_1.8.9             shiny_1.7.5.1          withr_2.5.2           
#  [58] ROCR_1.0-11            Rtsne_0.16             future_1.33.0         
#  [61] fastDummies_1.7.3      survival_3.5-7         polyclip_1.10-6       
#  [64] fitdistrplus_1.1-11    pillar_1.9.0           KernSmooth_2.23-22    
#  [67] plotly_4.10.3          generics_0.1.3         RcppHNSW_0.5.0        
#  [70] hms_1.1.3              munsell_0.5.0          scales_1.2.1          
#  [73] globals_0.16.2         xtable_1.8-4           glue_1.6.2            
#  [76] lazyeval_0.2.2         tools_4.3.1            RSpectra_0.16-1       
#  [79] RANN_2.6.1             leiden_0.4.3           dotCall64_1.1-0       
#  [82] grid_4.3.1             colorspace_2.1-0       nlme_3.1-163          
#  [85] cli_3.6.1              spatstat.sparse_3.0-3  spam_2.10-0           
#  [88] fansi_1.0.5            viridisLite_0.4.2      uwot_0.1.16           
#  [91] gtable_0.3.4           R.methodsS3_1.8.2      digest_0.6.33         
#  [94] progressr_0.14.0       ggrepel_0.9.4          farver_2.1.1          
#  [97] htmlwidgets_1.6.2      htmltools_0.5.6.1      R.oo_1.25.0           
# [100] lifecycle_1.0.4        httr_1.4.7             mime_0.12             
# [103] MASS_7.3-60
