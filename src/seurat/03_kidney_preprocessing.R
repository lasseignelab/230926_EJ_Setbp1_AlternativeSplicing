### Pre-processing of S858R mouse kidney data
# Author: Tabea M. Soelter
# Date: August 7th, 2025
  
## Goal: Process S858R (snRNA-seq) kidney data for downstream analyses.

## Reproducibility:
# * GitHub: lasseignelab/230926_EJ_Setbp1_AlternativeSplicing
# * Docker: tsoelter/setbp1_alternative_splicing
#       * Version: 2.0.0 

## Data:
# SoupX outputs
# * Consists of 2 conditions (S858R-mutants, wildtype) and 1 sex (male).
# * Name: N/A
# * Location: data/soupX/
  
## Analysis Plan:
# * Load necessary packages 
# * Load data 
# * Create seurat object
# * Perform quality control and filtering
# * Integration
# * Clustering
# * Cell type annotation
# * Save filtered, integrated, clustered, and annotated seurat objects

## Set-up:
args <- R.utils::commandArgs(trailingOnly = TRUE)
wd <- args[1]
setwd(wd)

set.seed(42)

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(harmony)
  library(stringr)
  library(gprofiler2)
  library(scCustomize)
})

source(file.path("src", "functions_soelter.R"))
source(file.path("src", "seurat", "functions.R"))

## Analysis:
# Object creation
merged_seurat <- make_seurat_object(file.path(
  "data",
  "soupX"
))
merged_seurat <- JoinLayers(merged_seurat)

# QC:
merged_seurat <- calc_qc(merged_seurat)

metadata <- format_metadata(merged_seurat)
merged_seurat@meta.data <- metadata

pdf(file = file.path("results", "seurat_outputs", "qc_plots_kidney.pdf"))
plot_qc(metadata)
dev.off()

# Filtering:
sub_seurat <- subset(merged_seurat,
                     subset =
                       nGene > 1000 &
                       nGene < 10000 &
                       mitoRatio < 0.2 &
                       log10GenesPerUMI > 0.80
)

counts <- GetAssayData(sub_seurat, layer = "counts")
nonzero <- counts > 0
keep_genes <- rowSums(nonzero) >= 10

filtered_counts <- counts[keep_genes, ]
filtered_seurat <- CreateSeuratObject(filtered_counts,
                                      meta.data = sub_seurat@meta.data
)

for (layer in names(filtered_seurat@commands)){
  filtered_seurat@commands[[layer]]@time.stamp <- as.POSIXct(character(0))
}

saveRDS(filtered_seurat,
        file = file.path("data",
                    "seurat",
                    "filtered_kidney_samples.rds"
                    )
        )

# Post-filtering QC:
metadata <- filtered_seurat@meta.data

pdf(file = file.path(
  "results",
  "seurat_outputs",
  "qc_plots_filtered_kidney.pdf"
))
plot_qc(metadata)
dev.off()

# Cell cycle scoring:
s.genes <- unlist(as.list(read_csv(
  file.path("doc/SGenes.csv"),
  col_names = FALSE,
  show_col_types = FALSE
)))
g2m.genes <- unlist(as.list(read_csv(
  file.path("doc/G2MGenes.csv"),
  col_names = FALSE,
  show_col_types = FALSE
)))

pdf(file.path("results", "seurat_outputs", "cellcycle_pca_kidney.pdf"), )
filtered_seurat <- cell_cycle_effects(filtered_seurat,
                                      s_genes = s.genes,
                                      g2m_genes = g2m.genes
                                      )
dev.off()

for (layer in names(filtered_seurat@commands)){
  filtered_seurat@commands[[layer]]@time.stamp <- as.POSIXct(character(0))
}

saveRDS(filtered_seurat,
        file = file.path("data", "seurat", "filtered_kidney_samples_pca.rds")
)

# Integration:
set.seed(42)
integrated_seurat <- RunHarmony(filtered_seurat,
                                group.by.vars = "sample_id"
                                )

integrated_seurat <- RunUMAP(integrated_seurat,
                             dims = 1:30,
                             reduction = "harmony",
                             reduction.name = "umap_harmony"
                             )

for (layer in names(integrated_seurat@commands)){
  integrated_seurat@commands[[layer]]@time.stamp <- as.POSIXct(character(0))
}

saveRDS(integrated_seurat,
        file = file.path("data", "seurat", "integrated_kidney_samples.rds")
        )

# Clustering:
clustered_seurat <- find_clusters(integrated_seurat,
                                  dims = 1:30,
                                  reduction = "harmony",
                                  resolutions = 1.25
                                  )

for (layer in names(clustered_seurat@commands)){
  clustered_seurat@commands[[layer]]@time.stamp <- as.POSIXct(character(0))
}

saveRDS(clustered_seurat,
        file = file.path("data", "seurat", "clustered_kidney_samples.rds")
        )

# Cell type annotation:
clustered_seurat <- SetIdent(clustered_seurat,
                             value = "RNA_snn_res.1.25"
)

marker_genes <- FindAllMarkers(clustered_seurat,
                               logfc.threshold = 0.2,
                               assay = "RNA",
                               only.pos = TRUE
                               )

saveRDS(marker_genes,
        file.path("results",
                  "seurat_outputs",
                  "kidney_marker_genes.rds"
                  )
        )

kidney_markers <- list(
  Podocytes = c("Nphs1", "Synpo"),
  Mesangial = c("Col1a2"),
  Proximal_Tubule = c("Slc34a1", "Slc13a3"),
  PCT = c("Slc5a2", "Slc6a19"),
  PST = c("Slc22a7", "Slc22a6"),
  LOH = c("Slc12a1"),
  Descending_thin_limb_LOH = c("Fst", "Aqp1"),
  Ascending_thin_limb_LOH = c("Epha7", "Mx2"),
  Thick_ascending_limb_LOH = c("Ptger3", "Tmem207"),
  DCT = c("Slc12a3"),
  CDPC = c("Aqp2", "Hsd11b2"),
  CDICA = c("Aqp6"),
  CDICB = c("Hmx2"),
  Connecting_Tubule = c("Atp6v1b1", "Fxyd4"),
  Fibroblasts = c("Col1a1", "Gucy1a1"),
  Pericytes = c("Pdgfrb"),
  Endothelial = c("Flt1", "Pecam1", "Cdh5"),
  Smooth_Muscle = c("Myh11"),
  Macrophages = c("Cd68", "C1qa"),
  Dendritic_Cells = c("Itgax", "H2-Aa"),
  T_Cells = c("Cxcr6", "Cd247"),
  T_regulatory_cells = c("Tnfrsf4", "Capg", "Ikzf2"),
  B_Cells = c("Cd79a", "Cd79b"),
  NK_Cells = c("Nkg7")
)

png(file.path("results", "seurat_outputs", "celltype_markers_kidney.png"),
    width = 48,
    height = 24,
    units = "cm",
    res = 300
    )
dotplot <- DotPlot(clustered_seurat, features = kidney_markers) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(angle = 90, hjust = 0))
print(dotplot)
dev.off()

annotated_seurat <- RenameIdents(clustered_seurat,
                                 `1` = "Endothelial cells",
                                 `2` = "PCT",
                                 `3` = "Mesenchymal cells",
                                 `4` = "PST",
                                 `5` = "Thick ascending limb (LOH)",
                                 `6` = "Proximal tubule cells",
                                 `7` = "PST",
                                 `8` = "Thick ascending limb (LOH)",
                                 `9` = "PCT",
                                 `10` = "Endothelial cells",
                                 `11` = "PST",
                                 `12` = "PST",
                                 `13` = "DCT",
                                 `14` = "CDPC",
                                 `15` = "DCT",
                                 `16` = "DCT",
                                 `17` = "Proximal tubule cells",
                                 `18` = "Proximal tubule cells",
                                 `19` = "Dendritic cells",
                                 `20` = "Mesenchymal cells",
                                 `21` = "CDPC",
                                 `22` = "Thin descending limb (LOH)",
                                 `23` = "CDIC-A",
                                 `24` = "CDIC-B",
                                 `25` = "Thin ascending limb (LOH)",
                                 `26` = "Mesenchymal cells",
                                 `27` = "Thick ascending limb (LOH)",
                                 `28` = "CDPC",
                                 `29` = "Proximal tubule cells",
                                 `30` = "Proximal tubule cells",
                                 `31` = "Dendritic cells",
                                 `32` = "Proximal tubule cells",
                                 `33` = "DCT",
                                 `34` = "Mesenchymal cells",
                                 `35` = "Thin ascending limb (LOH)",
                                 `36` = "Mesenchymal cells",
                                 `37` = "Proximal tubule cells",
                                 `38` = "T cells",
                                 `39` = "Podocytes",
                                 `40` = "Endothelial cells",
                                 `41` = "Podocytes",
                                 `42` = "B cells",
                                 `43` = "Mesenchymal cells",
                                 `44` = "Thick ascending limb (LOH)",
                                 `45` = "Thick ascending limb (LOH)",
                                 `46` = "Connecting tubule cells",
                                 `47` = "T regulatory cells",
                                 `48` = "PST",
                                 `49` = "Proximal tubule cells",
                                 `50` = "Proximal tubule cells"
                                 )

annotated_seurat <- AddMetaData(
  object = annotated_seurat,
  as.vector(annotated_seurat@active.ident),
  col.name = "cell_type"
  )

png(file.path("results", "seurat_outputs", "annotated_umap_kidney.png"))
DimPlot(annotated_seurat)
dev.off()

# Save final object for downstream analyses:
for (layer in names(annotated_seurat@commands)){
  annotated_seurat@commands[[layer]]@time.stamp <- as.POSIXct(character(0))
}

saveRDS(annotated_seurat,
        file = file.path("data", "seurat", "annotated_kidney_samples.rds")
        )

# Session Information:
sessionInfo()
#R version 4.3.1 (2023-06-16)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 22.04.3 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

#locale:
#  [1] C

#time zone: Etc/UTC
#tzcode source: system (glibc)

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] scCustomize_3.0.1  gprofiler2_0.2.3   harmony_1.1.0      Rcpp_1.0.12       
#[5] lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4       
#[9] purrr_1.0.2        readr_2.1.4        tidyr_1.3.0        tibble_3.2.1      
#[13] ggplot2_3.5.2      tidyverse_2.0.0    Seurat_5.0.0       SeuratObject_5.0.0
#[17] sp_2.1-1          

#loaded via a namespace (and not attached):
#  [1] RColorBrewer_1.1-3     shape_1.4.6            jsonlite_1.8.7        
#[4] magrittr_2.0.3         ggbeeswarm_0.7.2       spatstat.utils_3.0-4  
#[7] farver_2.1.1           GlobalOptions_0.1.2    vctrs_0.6.5           
#[10] ROCR_1.0-11            spatstat.explore_3.2-5 paletteer_1.6.0       
#[13] janitor_2.2.1          htmltools_0.5.6.1      sctransform_0.4.1     
#[16] parallelly_1.36.0      KernSmooth_2.23-22     htmlwidgets_1.6.2     
#[19] ica_1.0-3              plyr_1.8.9             plotly_4.10.3         
#[22] zoo_1.8-12             igraph_1.5.1           mime_0.12             
#[25] lifecycle_1.0.4        pkgconfig_2.0.3        Matrix_1.6-1.1        
#[28] R6_2.5.1               fastmap_1.1.1          snakecase_0.11.1      
#[31] fitdistrplus_1.1-11    future_1.33.0          shiny_1.7.5.1         
#[34] digest_0.6.33          colorspace_2.1-0       rematch2_2.1.2        
#[37] patchwork_1.3.1        rprojroot_2.0.3        tensor_1.5            
#[40] RSpectra_0.16-1        irlba_2.3.5.1          labeling_0.4.3        
#[43] progressr_0.14.0       fansi_1.0.6            spatstat.sparse_3.0-3 
#[46] timechange_0.2.0       mgcv_1.9-0             httr_1.4.7            
#[49] polyclip_1.10-6        abind_1.4-5            compiler_4.3.1        
#[52] here_1.0.1             bit64_4.0.5            withr_3.0.0           
#[55] fastDummies_1.7.3      R.utils_2.12.2         MASS_7.3-60           
#[58] rappdirs_0.3.3         tools_4.3.1            vipor_0.4.7           
#[61] lmtest_0.9-40          beeswarm_0.4.0         httpuv_1.6.12         
#[64] future.apply_1.11.0    goftest_1.2-3          R.oo_1.25.0           
#[67] glue_1.7.0             nlme_3.1-163           promises_1.2.1        
#[70] grid_4.3.1             Rtsne_0.16             cluster_2.1.4         
#[73] reshape2_1.4.4         generics_0.1.3         gtable_0.3.6          
#[76] spatstat.data_3.0-3    tzdb_0.4.0             R.methodsS3_1.8.2     
#[79] data.table_1.14.8      hms_1.1.3              utf8_1.2.4            
#[82] spatstat.geom_3.2-7    RcppAnnoy_0.0.21       ggrepel_0.9.5         
#[85] RANN_2.6.1             pillar_1.9.0           limma_3.56.2          
#[88] vroom_1.6.4            ggprism_1.0.6          spam_2.10-0           
#[91] RcppHNSW_0.5.0         later_1.3.1            circlize_0.4.15       
#[94] splines_4.3.1          lattice_0.21-8         bit_4.0.5             
#[97] survival_3.5-7         deldir_1.0-9           tidyselect_1.2.1      
#[100] miniUI_0.1.1.1         pbapply_1.7-2          gridExtra_2.3         
#[103] scattermore_1.2        RhpcBLASctl_0.23-42    matrixStats_1.0.0     
#[106] stringi_1.8.1          lazyeval_0.2.2         codetools_0.2-19      
#[109] cli_3.6.2              uwot_0.1.16            xtable_1.8-4          
#[112] reticulate_1.34.0      munsell_0.5.1          globals_0.16.2        
#[115] spatstat.random_3.2-1  png_0.1-8              ggrastr_1.0.2         
#[118] parallel_4.3.1         ellipsis_0.3.2         presto_1.0.0          
#[121] dotCall64_1.1-0        listenv_0.9.0          viridisLite_0.4.2     
#[124] scales_1.3.0           ggridges_0.5.4         crayon_1.5.2          
#[127] leiden_0.4.3           rlang_1.1.3            cowplot_1.1.1
