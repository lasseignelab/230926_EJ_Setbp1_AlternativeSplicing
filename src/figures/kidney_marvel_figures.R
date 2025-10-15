### Plotting S858R mouse kidney MARVEL results
# Author: Tabea M. Soelter (original code modified from Emma Jones)
# Date: October 9th, 2025

## Goal: Plot MARVEL results in the S858R kidney data.

## Reproducibility:
# * GitHub: lasseignelab/230926_EJ_Setbp1_AlternativeSplicing
# * Docker: emmafjones/setbp1_alternative_splicing
#       * Version: 1.0.9 

## Data:
# Annotated Seurat object
# * Name: annotated_kidney_samples.rds
# * Location: data/seurat/

## Analysis Plan:
# * Load necessary packages 
# * Load data 
# * Plot Setbp1 AS figure panels
# * Compile Setbp1 AS figure
# * Save Setbp1 AS figure
# * Plot & save AS splicing summary plots for Shiny app

## Set-up:
args <- R.utils::commandArgs(trailingOnly = TRUE)
wd <- args[1]
setwd(wd)

set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(viridis)
  library(MARVEL)
  library(ComplexHeatmap)
  library(Matrix)
  library(slam)
  library(circlize)
  library(ggtranscript)
  library(rtracklayer)
})

# source functions
source(file.path(wd, "src", "figures", "geom_split_violin.R"))
source(file.path(wd, "src", "figures", "functions.R"))
source(file.path(wd, "src", "functions_soelter.R"))

## Analysis:
# Load data
sj_usage_cell_type <- readRDS(file.path(
  wd, 
  "data", "marvel",
  "sj_usage_kidney_cell_type.rds"
))

sj_usage_condition <- readRDS(file.path(
  wd, 
  "data", "marvel",
  "sj_usage_kidney_condition.rds"
))

setbp1_marvel <- readRDS(file.path(
  wd, 
  "data", "marvel",
  "setbp1_kidney_marvel_aligned.rds"
))

sj_counts <- readRDS(file.path(wd, "data", "marvel", "sj_counts_kidney.rds"))
sj_counts <- as(sj_counts, "CsparseMatrix")

gene_counts <- readRDS(file.path(wd, "data", "marvel", "gene_counts_kidney.rds"))
gene_counts <- as(gene_counts, "CsparseMatrix")

# Set color palette
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

# Prepare data for plotting
normalized_sj_expression <- as(
  setbp1_marvel[["sj.count.matrix"]],
  "CsparseMatrix"
  )
normalized_sj_expression@x <- normalized_sj_expression@x /
  rep.int(colSums(normalized_sj_expression), diff(normalized_sj_expression@p))
normalized_sj_expression@x <- normalized_sj_expression@x * 1000
normalized_sj_expression@x <- log1p(normalized_sj_expression@x)

write_rds(
  normalized_sj_expression,
  file.path(wd, "data", "marvel", "normalized_kidney_sj_expression.Rds")
)

# Prepapre sample metadata
sample_metadata <- setbp1_marvel$sample.metadata
sample_metadata$num_sjs <- diff(sj_counts@p)
sample_metadata$num_genes <- diff(gene_counts@p)
sample_metadata$num_sjs_genes <- diff(sj_counts@p) / diff(gene_counts@p)

# Split data by cell type and condition
cell_group_list <- sample_metadata %>%
  group_split(cell_type, .keep = TRUE) %>%
  map(~ set_names(.$cell.id, .$cell_type[1]))

cell_types <- list()
for (i in seq_along(cell_group_list)) {
  cell_types[i] <- names(cell_group_list[[i]][1])
}

cell_group_list <- set_names(
  cell_group_list,
  cell_types
)

mutant_list <- sample_metadata %>%
  group_split(seq_folder, .keep = TRUE) %>%
  map(~ set_names(.$cell.id, .$seq_folder[[1]]))

condition_types <- list()
for (i in seq_along(mutant_list)) {
  condition_types[i] <- names(mutant_list[[i]][1])
}

mutant_list <- set_names(
  mutant_list,
  condition_types
)

mutant_ids <- mutant_list[["mutant"]]

# Grab mutant and wildtype SJE and gene information
mutant_sj <- normalized_sj_expression[
  ,
  colnames(normalized_sj_expression) == mutant_ids
  ]

mutant_gene <-
  setbp1_marvel$gene.norm.matrix[
    , colnames(setbp1_marvel$gene.norm.matrix)
    == mutant_ids
  ]

wildtype_sj <-
  normalized_sj_expression[, !colnames(normalized_sj_expression) == mutant_ids]

wildtype_gene <-
  setbp1_marvel$gene.norm.matrix[
    , !colnames(setbp1_marvel$gene.norm.matrix)
    == mutant_ids
  ]

# Reformat SJU information
mutant_sju <- sj_usage_condition[1, ]
mutant_sju <- do.call(cbind, mutant_sju)
colnames(mutant_sju) <- names(cell_group_list)

wt_sju <- sj_usage_condition[2, ]
wt_sju <- do.call(cbind, wt_sju)
colnames(wt_sju) <- names(cell_group_list)

# Reformat cell type level data
colnames(sj_usage_cell_type) <- names(cell_type_colors)

sj_usage_cell_type <- as.data.frame(sj_usage_cell_type)
sj_usage_cell_type$sj_position <- rownames(sj_usage_cell_type)
rownames(sj_usage_cell_type) <- NULL
sj_usage_cell_type <- pivot_longer(
  sj_usage_cell_type,
  cols = colnames(sj_usage_cell_type)[1:18],
  names_to = "cell_type",
  values_to = "sj_usage"
)

sj_usage_cell_type$sj_usage[is.infinite(sj_usage_cell_type$sj_usage)] <- NaN
sj_usage_cell_type$sj_usage[sj_usage_cell_type$sj_usage > 100] <- 100

mutant_sju[is.infinite(mutant_sju)] <- NaN
mutant_sju[mutant_sju > 100] <- 100

colnames(mutant_sju) <- paste0(colnames(mutant_sju), "_mutant")

wt_sju[is.infinite(wt_sju)] <- NaN
wt_sju[wt_sju > 100] <- 100

colnames(wt_sju) <- paste0(colnames(wt_sju), "_wildtype")

# Bind conditions together
conditions_sju <- cbind(mutant_sju, wt_sju)

# Create delta SJU matrix
delta_sju <- mutant_sju - wt_sju

base_df <- sample_metadata[, c("cell.id", "seq_folder", "cell_type")]

setbp1_gene_expr_df <- make_norm_gene_expr_df("Setbp1")

# Plot Setbp1 expression violin
setbp1_split_violin <- ggplot(
  setbp1_gene_expr_df,
  aes(
    x = cell_type, y = norm_expr,
    fill = cell_type, alpha = seq_folder
  )
) +
  geom_split_violin(show.legend = TRUE) +
  scale_alpha_manual(
    name = "Condition",
    labels = c("S858R+/-", "Wildtype"),
    values = c(1, 0.4),
    guide = guide_legend(override.aes = list(fill = "gray40", color = "black"))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    labels = names(cell_type_colors),
    values = cell_type_colors,
    guide = "none"
  ) +
  xlab("Cell Type") +
  ylab(expression(bold(atop(
    paste("Normalized ", bolditalic("Setbp1")),
    "Gene Expression"
  ))))

# Split SJE by cell type
# mutants
mutant_sj_expr_means <- list()

for (i in seq_along(cell_group_list)) {

  cell_group_mutants <- intersect(cell_group_list[[i]], mutant_list[["mutant"]])

  cell_type_sj_expr <-
    normalized_sj_expression[
      ,
      colnames(normalized_sj_expression) %in% cell_group_mutants
    ]

  mutant_sj_expr_means[[i]] <- rowMeans(cell_type_sj_expr)
}
names(mutant_sj_expr_means) <- names(cell_group_list)
mutant_sj_expr_means <- do.call(cbind, mutant_sj_expr_means)
mutant_sj_expr_means <- as.data.frame(mutant_sj_expr_means)
colnames(mutant_sj_expr_means) <- paste0(colnames(mutant_sj_expr_means),
                                         "_mutant"
                                         )

# wildtypes
sj_expr_means <- rowMeans(normalized_sj_expression)

wildtype_sj_expr_means <- list()

for (i in seq_along(cell_group_list)) {

  cell_group_wildtypes <-
    intersect(cell_group_list[[i]], mutant_list[["wildtype"]])

  cell_type_sj_expr <-
    normalized_sj_expression[
      ,
      colnames(normalized_sj_expression) %in% cell_group_wildtypes
    ]

  wildtype_sj_expr_means[[i]] <- rowMeans(cell_type_sj_expr)
}
names(wildtype_sj_expr_means) <- names(cell_group_list)
wildtype_sj_expr_means <- do.call(cbind, wildtype_sj_expr_means)
wildtype_sj_expr_means <- as.data.frame(wildtype_sj_expr_means)
colnames(wildtype_sj_expr_means) <- paste0(colnames(wildtype_sj_expr_means),
                                           "_wildtype")

# Combine conditions into one matrix
sj_expr_for_heatmap <- cbind(mutant_sj_expr_means, wildtype_sj_expr_means)
sj_expr_for_heatmap <- as.matrix(sj_expr_for_heatmap)

# Make delta SJE matrix
delta_expr_for_heatmap <-
  as.matrix(mutant_sj_expr_means) - as.matrix(wildtype_sj_expr_means)

# Create heatmap color annotation
annot_colors_list <- list(
  `Cell Type` = cell_type_colors,
  `Condition` = c("S858R+/-" = "#666666", "Wildtype" = "#C1C1C1")
)

# Get setbp1 SJs & SJE
setbp1_sjs <-
  setbp1_marvel[["sj.metadata"]]$coord.intron[
    setbp1_marvel[["sj.metadata"]]$gene_short_name.start == "Setbp1"
  ]

setbp1_sj_mean_expr <- sj_expr_means[setbp1_sjs]

# Set heatmap annotation & colors
col_annot_conditions <- HeatmapAnnotation(
  `Cell Type` = rep(names(cell_group_list), each = 2),
  `Condition` = rep(c("S858R+/-", "Wildtype"), times = 18),
  col = annot_colors_list,
  annotation_name_gp = gpar(fontface = "bold")
)

gene_expr_cols <- c("#FFFFFF", RColorBrewer::brewer.pal(name = "GnBu", n = 9))

sju_cols <- c("#FFFFFF", RColorBrewer::brewer.pal(name = "YlOrBr", n = 9))

# Extract cell type names to create the interleaved column order
cell_types <- names(cell_group_list)

desired_order <- as.vector(sapply(cell_types, function(ct) {
  c(paste0(ct, "_mutant"), paste0(ct, "_wildtype"))
}))

# Make SJE heatmap
setbp1_sje_heatmap_df <-
  sj_expr_for_heatmap[rownames(sj_expr_for_heatmap) %in% setbp1_sjs, ]

rownames(setbp1_sje_heatmap_df) <- c(paste0("SJ-", rev(1:6)))

setbp1_sje_heatmap_df <- setbp1_sje_heatmap_df[, desired_order]

gene_heatmap <- Heatmap(setbp1_sje_heatmap_df,
                        name = "Mean Norm\nSJ Expression",
                        col = gene_expr_cols,
                        top_annotation = col_annot_conditions, cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        show_column_names = FALSE,
                        row_title = "SJ Expression", border = TRUE,
                        row_title_gp = gpar(fontface = "bold"),
                        row_names_gp = gpar(fontface = "bold")
)

# make SJU heatmap
setbp1_sju_heatmap_df <-
  conditions_sju[rownames(conditions_sju) %in% setbp1_sjs, ]

rownames(setbp1_sju_heatmap_df) <- c(paste0("SJ-", rev(1:6)))

setbp1_sju_heatmap_df <- setbp1_sju_heatmap_df[, desired_order]

sj_heatmap <- Heatmap(setbp1_sju_heatmap_df,
                      name = "SJ Usage",
                      col = sju_cols,
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      show_column_names = FALSE,
                      row_title = "SJ Usage", border = TRUE,
                      row_title_gp = gpar(fontface = "bold"),
                      row_names_gp = gpar(fontface = "bold")
)

# Combine heatmaps
both_heatmaps <- gene_heatmap %v% sj_heatmap

sj_expr_usage_heatmaps <- grid.grabExpr(draw(both_heatmaps))

# Save Setbp1 SJE & SJU heatmaps
png(file.path(wd, "results", "kidney_figures", "setbp1_sj_expr_usage.png"),
    width = 8, height = 6, units = "in", res = 300
)
draw(both_heatmaps)
dev.off()

# Make annotation for delta heatmaps
col_annot <- HeatmapAnnotation(
  `Cell Type` = names(cell_group_list),
  col = annot_colors_list, annotation_name_gp = gpar(fontface = "bold")
)

# Make delta SJE heatmap
setbp1_delta_sje_heatmap_df <-
  delta_expr_for_heatmap[rownames(delta_expr_for_heatmap) %in% setbp1_sjs, ]

rownames(setbp1_delta_sje_heatmap_df) <- c(paste0("SJ-", rev(1:6)))

# Set colors
delta_gene_expr_cols <-
  colorRamp2(
    c(
      min(setbp1_delta_sje_heatmap_df), 0,
      max(setbp1_delta_sje_heatmap_df)
    ),
    c(RColorBrewer::brewer.pal(name = "BrBG", n = 3))
  )

delta_gene_heatmap <- Heatmap(setbp1_delta_sje_heatmap_df,
                              name = "\u0394 Mean Norm\nSJ Expression",
                              column_order = sort(colnames(setbp1_delta_sje_heatmap_df)),
                              col = delta_gene_expr_cols,
                              top_annotation = col_annot, cluster_rows = FALSE,
                              show_column_names = FALSE,
                              row_title = expression(bold(paste(
                                bolditalic("Setbp1 "),
                                "SJ Expression"
                              ))), border = TRUE,
                              row_title_gp = gpar(fontface = "bold"),
                              row_names_gp = gpar(fontface = "bold")
)

# Make delta SJU heatmap
setbp1_delta_sju_heatmap_df <-
  delta_sju[rownames(delta_sju) %in% setbp1_sjs, ]

rownames(setbp1_delta_sju_heatmap_df) <- c(paste0("SJ-", rev(1:6)))

# Set colors
delta_sju_cols <-
  colorRamp2(
    c(
      min(setbp1_delta_sju_heatmap_df), 0,
      max(setbp1_delta_sju_heatmap_df)
    ),
    c(RColorBrewer::brewer.pal(name = "PRGn", n = 3))
  )

delta_sj_heatmap <- Heatmap(setbp1_delta_sju_heatmap_df,
                            name = "\u0394 SJ Usage",
                            column_order = sort(colnames(setbp1_delta_sju_heatmap_df)),
                            col = delta_sju_cols,
                            cluster_rows = FALSE,
                            show_column_names = FALSE,
                            row_title = expression(bold(paste(
                              bolditalic("Setbp1 "),
                              "SJ Usage"
                            ))), border = TRUE,
                            row_title_gp = gpar(fontface = "bold"),
                            row_names_gp = gpar(fontface = "bold")
)

# Combine delta heatmaps for Setbp1 figure
both_delta_heatmaps <- delta_gene_heatmap %v% delta_sj_heatmap

delta_sj_expr_usage_heatmaps <- grid.grabExpr(draw(both_delta_heatmaps))

# Visualize Setbp1 transcripts
gtf_file_path <- "/data/project/lasseigne_lab/GENOME_dir/GENCODE_mm39/release_M31/GTF/gencode.vM31.primary_assembly.annotation.gtf"

gtf <- rtracklayer::import(gtf_file_path)
gtf <- gtf %>% dplyr::as_tibble()

setbp1_annotation_from_gtf <- gtf %>%
  dplyr::filter(
    !is.na(gene_name),
    gene_name == "Setbp1"
  )
setbp1_annotation_from_gtf <- setbp1_annotation_from_gtf %>%
  dplyr::select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_name,
    transcript_name,
    transcript_type
  )

setbp1_exons <- setbp1_annotation_from_gtf %>% dplyr::filter(type == "exon")
setbp1_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_name
  )) +
  geom_range(
    aes(fill = transcript_type)
  ) +
  geom_intron(
    data = to_intron(setbp1_exons, "transcript_name"),
    aes(strand = strand)
  )

setbp1_junctions <- strsplit(setbp1_sjs, ":")
setbp1_junctions <- as_tibble(do.call(rbind, setbp1_junctions))
colnames(setbp1_junctions) <- c("seqnames", "start", "end")
setbp1_junctions$start <- as.numeric(setbp1_junctions$start)
setbp1_junctions$end <- as.numeric(setbp1_junctions$end)
setbp1_junctions$strand <- "-"
setbp1_junctions$mean_count <- setbp1_sj_mean_expr
setbp1_junctions <- setbp1_junctions %>%
  dplyr::mutate(transcript_name = "Setbp1-201")
setbp1_junctions$sj_label <- paste0("SJ-", rev(1:6))
 
setbp1_transcript_label <- setbp1_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_name
  )) +
  geom_range(
    aes(fill = transcript_type)
  ) +
  geom_intron(
    data = to_intron(setbp1_exons, "transcript_name"),
    aes(strand = strand)
  ) +
  geom_junction(
    data = setbp1_junctions, junction.y.max = 0.5,
    aes(linewidth = mean_count)
  ) +
  geom_junction_label_repel(
    data = setbp1_junctions,
    aes(label = sj_label), junction.y.max = 0.5,
    segment.color = NA
  ) +
  scale_linewidth(range = c(0.1, 1)) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title.y = element_blank()
  ) +
  guides(fill = FALSE) +
  labs(linewidth = "Normalized SJ\nExpression") +
  xlab("Genomic Location")

# Combine Setbp1 panels
gene_expr_vln <- plot_grid(setbp1_transcript_label, setbp1_split_violin,
                           ncol = 1,
                           labels = c("A", "B")
                           )

setbp1_figure <- plot_grid(gene_expr_vln, delta_sj_expr_usage_heatmaps,
                           labels = c("", "C"), nrow = 1
                           )



# Save Setbp1 figure
png(file.path(wd, "results", "kidney_figures", "setbp1.png"),
    width = 14, height = 9, units = "in", res = 300
)
setbp1_figure
dev.off()

# Get AS summary plots for all significant genes
load(file.path(wd, "data", "marvel", "significant_tables_kidney.RData"))

matching_dfs <- mget(ls(pattern = "sig_table"))

sig_sj_genes <- unique(unlist(lapply(matching_dfs, function(df_name) {
  df_name$gene_short_name
})))

filepath <- file.path(wd, "results", "as_gene_summaries", "kidney")
dir.create(filepath, recursive = TRUE, showWarnings = FALSE)

lapply(sig_sj_genes, function(gene) {
  make_gene_sj_expr_usage_plots_ts(
    gene_of_interest = gene,
    save_path = filepath
  )
})

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
#  [1] rtracklayer_1.60.1    GenomicRanges_1.52.1  GenomeInfoDb_1.36.4  
#  [4] IRanges_2.34.1        S4Vectors_0.38.2      BiocGenerics_0.46.0  
#  [7] ggtranscript_0.99.9   circlize_0.4.15       slam_0.1-50          
# [10] Matrix_1.6-1.1        ComplexHeatmap_2.16.0 MARVEL_2.0.5         
# [13] viridis_0.6.4         viridisLite_0.4.2     patchwork_1.1.3      
# [16] cowplot_1.1.1         lubridate_1.9.3       forcats_1.0.0        
# [19] stringr_1.5.1         dplyr_1.1.4           purrr_1.0.2          
# [22] readr_2.1.4           tidyr_1.3.0           tibble_3.2.1         
# [25] ggplot2_3.5.0         tidyverse_2.0.0      
#
# loaded via a namespace (and not attached):
#  [1] tidyselect_1.2.1            farver_2.1.1               
#  [3] R.utils_2.12.2              Biostrings_2.68.1          
#  [5] bitops_1.0-7                RCurl_1.98-1.12            
#  [7] GenomicAlignments_1.36.0    XML_3.99-0.14              
#  [9] digest_0.6.33               timechange_0.2.0           
# [11] lifecycle_1.0.4             cluster_2.1.4              
# [13] magrittr_2.0.3              compiler_4.3.1             
# [15] rlang_1.1.3                 tools_4.3.1                
# [17] utf8_1.2.4                  yaml_2.3.7                 
# [19] labeling_0.4.3              S4Arrays_1.0.6             
# [21] DelayedArray_0.26.7         plyr_1.8.9                 
# [23] RColorBrewer_1.1-3          abind_1.4-5                
# [25] BiocParallel_1.34.2         withr_3.0.0                
# [27] R.oo_1.25.0                 fansi_1.0.6                
# [29] colorspace_2.1-0            scales_1.3.0               
# [31] iterators_1.0.14            SummarizedExperiment_1.30.2
# [33] cli_3.6.2                   crayon_1.5.2               
# [35] generics_0.1.3              tzdb_0.4.0                 
# [37] rjson_0.2.21                zlibbioc_1.46.0            
# [39] parallel_4.3.1              XVector_0.40.0             
# [41] restfulr_0.0.15             matrixStats_1.0.0          
# [43] vctrs_0.6.5                 hms_1.1.3                  
# [45] GetoptLong_1.0.5            ggrepel_0.9.5              
# [47] clue_0.3-65                 foreach_1.5.2              
# [49] glue_1.7.0                  codetools_0.2-19           
# [51] stringi_1.8.1               shape_1.4.6                
# [53] gtable_0.3.4                BiocIO_1.10.0              
# [55] munsell_0.5.1               pillar_1.9.0               
# [57] GenomeInfoDbData_1.2.10     R6_2.5.1                   
# [59] doParallel_1.0.17           Biobase_2.60.0             
# [61] lattice_0.21-8              R.methodsS3_1.8.2          
# [63] png_0.1-8                   Rsamtools_2.16.0           
# [65] Rcpp_1.0.12                 gridExtra_2.3              
# [67] MatrixGenerics_1.12.3       pkgconfig_2.0.3            
# [69] GlobalOptions_0.1.2
