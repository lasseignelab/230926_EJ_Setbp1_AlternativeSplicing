### Functions for the Setbp1 alternative splicing project
# Tabea M. Soelter

## find_file
# Helper function to find STAR Solo files with or without .gz extension to
# accommodate EFJ's existing file naming conventions.
find_file <- function(dir, basename) {
  gz_file <- file.path(dir, paste0(basename, ".gz"))
  unzipped_file <- file.path(dir, basename)
  
  if (file.exists(gz_file)) {
    return(gz_file)
  } else if (file.exists(unzipped_file)) {
    return(unzipped_file)
  } else {
    stop(paste("Neither", gz_file, "nor", unzipped_file, "exists"))
  }
}


## remove_ambient_rna
# A function which removes ambient RNA from matrix, barcodes, and feature files
# from STAR Solo (inputs) and generates filtered barcodes, features, and matrix
# tsv files (outputs), which serve as input for downstream preprocessing.
remove_ambient_rna <- function(inputs, outputs, plots) {
  print("Making list of objects")
  sample_dir_paths <- list.dirs(inputs, full.names = TRUE, recursive = FALSE)
  sample_dir_paths <- sample_dir_paths[grepl("^[JK]",
                                       basename(sample_dir_paths)
                                            )
                                      ]

  pdf(file.path(plots, "rho_density_plots.pdf"))
  for (sample_dir_path in sample_dir_paths) {
    sample_name <- basename(sample_dir_path)
    print(sample_name)
    set.seed(42)

    # Adjust path to match STARsolo structure
    base_dir <- file.path(sample_dir_path, "Solo.out", "GeneFull_Ex50pAS")
    filt_dir <- file.path(base_dir, "filtered")
    raw_dir  <- file.path(base_dir, "raw")

    # Load filtered matrix
    print("Loading filtered STARsolo matrix")
    filt_matrix <- Matrix::readMM(find_file(filt_dir, "matrix.mtx"))
    filt_features <- read.delim(find_file(filt_dir, "features.tsv"),
                                header = FALSE
                                )
    filt_barcodes <- read.delim(find_file(filt_dir, "barcodes.tsv"),
                                header = FALSE
                                )
    rownames(filt_matrix) <- make.unique(filt_features$V2)
    colnames(filt_matrix) <- filt_barcodes$V1

    # Load raw matrix
    print("Loading raw STARsolo matrix")
    raw_matrix <- Matrix::readMM(find_file(raw_dir, "matrix.mtx"))
    raw_features <- read.delim(find_file(raw_dir, "features.tsv"),
                               header = FALSE
                               )
    raw_barcodes <- read.delim(find_file(raw_dir, "barcodes.tsv"),
                               header = FALSE
                               )
    rownames(raw_matrix) <- make.unique(raw_features$V2)
    colnames(raw_matrix) <- raw_barcodes$V1

    # Create Seurat object from filtered matrix
    print("Creating Seurat object for clustering")
    object <- CreateSeuratObject(counts = filt_matrix)
    object <- NormalizeData(object, verbose = FALSE)
    object <- FindVariableFeatures(object,
                                   selection.method = "vst",
                                   nfeatures = 5000,
                                   verbose = FALSE
                                   )
    object <- ScaleData(object, verbose = FALSE)
    object <- RunPCA(object, approx = FALSE, verbose = FALSE)
    object <- FindNeighbors(object, dims = 1:30, verbose = FALSE)
    object <- FindClusters(object, verbose = TRUE, graph.name = "RNA_snn")
    object <- RunUMAP(object, dims = 1:30, verbose = FALSE)

    # Prepare metadata and UMAP for SoupX
    meta <- object@meta.data
    umap <- Embeddings(object, "umap")

    # Create SoupX channel
    print("Creating SoupChannel object")
    sco <- SoupChannel(tod = raw_matrix, toc = filt_matrix)
    sco <- setClusters(sco, setNames(meta$seurat_clusters, rownames(meta)))
    sco <- setDR(sco, umap)

    print("Profiling the soup")
    sco <- autoEstCont(sco)

    print("Adjusting counts")
    adjusted_matrix <- adjustCounts(sco, roundToInt = TRUE)

    # Save output in 10X format
    if (!dir.exists(outputs)) dir.create(outputs)
    output_path <- file.path(outputs, sample_name)
    print(paste0("Saving filtered objects to: ", output_path))
    DropletUtils::write10xCounts(output_path, adjusted_matrix)
  }
  dev.off()
}

## make_seurat_object
# A function which takes a path to sample folders with the three CellRanger
# output files and creates a merged seurat object
make_seurat_object <- function(path) {
  counts_list <- list.dirs(path, 
                           full.names = TRUE, 
                           recursive = FALSE)
  counts_list <- grep("K", counts_list, value = TRUE)
  object_list <- vector('list')
  for (i in counts_list){
    counts <- Read10X(i)
    print(i)
    sample_name <- basename(i)
    seurat_object <- CreateSeuratObject(counts = counts,
                                        project = sample_name,
                                        min.features = 200
                                        )
    if (substring(sample_name, nchar(sample_name)) %in% c("1", "3", "5")) {
      seurat_object$condition <- "wildtype"
    } else {
      seurat_object$condition <- "mutant"
    }
    seurat_object$sample_id <- sample_name
    object_list[[i]] <- seurat_object
  } 
  print("Making wildtype Seurat Object")
  WT_list <- object_list[grepl("K[135]", names(object_list))]
  for (i in names(WT_list)) {
    sample_name <- basename(i)
    WT_list[[i]] <- RenameCells(WT_list[[i]],
                                add.cell.id = sample_name)
  }
  wildtypes <- Merge_Seurat_List(WT_list)
  
  print("Making mutant Seurat Object")
  mutant_list <- object_list[grepl("K[246]", names(object_list))]
  for (i in names(mutant_list)) {
    sample_name <- basename(i)
    mutant_list[[i]] <- RenameCells(mutant_list[[i]],
                                    add.cell.id = sample_name)
  }
  mutants <- Merge_Seurat_List(mutant_list)
  
  print("Making merged Seurat Object")
  merged_seurat <- merge(x = wildtypes,
                         y = mutants,
                         add.cell.id = c("WT", "S858R"))
  colnames(merged_seurat) <- gsub("__", "_", colnames(merged_seurat))
  return(merged_seurat)
}


## calc_qc
# A function which calculates quality control metrics for a merged seurat object
calc_qc <- function(seurat_object) {
  seurat_object$log10GenesPerUMI <- log10(seurat_object$nFeature_RNA) / log10(seurat_object$nCount_RNA)
  seurat_object$mitoRatio <- PercentageFeatureSet(object = seurat_object, 
                                                  pattern = "^mt-")
  seurat_object$mitoRatio <- seurat_object@meta.data$mitoRatio / 100
  return(seurat_object)
}


## cell_cycle_effects
# A function which calculates and plots the effect of cell cycle on the data
# using a filtered seurat object as input. It also performs log normalization,
# scaling, and dimension reduction using PCA
cell_cycle_effects <- function(filtered_seurat, g2m_genes, s_genes) {
  # log normalize -----
  filtered_seurat <- NormalizeData(filtered_seurat)
  # score cells based in gex of genes -----
  filtered_seurat <- CellCycleScoring(filtered_seurat,
                                      g2m.features = g2m_genes,
                                      s.features = s_genes)
  set.seed(42)
  filtered_seurat <- FindVariableFeatures(filtered_seurat,
                                          selection.method = "vst",
                                          verbose = FALSE)
  # scale data -----
  filtered_seurat <- ScaleData(filtered_seurat)
  # run pca -----
  set.seed(42)
  filtered_seurat <- RunPCA(filtered_seurat, approx = FALSE)
  # plot pca -----
  elbow <- ElbowPlot(filtered_seurat, reduction = "pca", ndims = 50)
  # plot cell cycle scoring -----
  cell_cycle_plot <- DimPlot(filtered_seurat,
                             reduction = "pca",
                             group.by = "Phase",
                             split.by = "Phase")
  plot(cell_cycle_plot)
  plot(elbow)
  return(filtered_seurat)
}

## make_gene_sj_expr_usage_plots
# Original code by Emma Jones as found in src/figures/functions.R
# This function has additional functionality to account for differing numbers
# of SJs for each gene causing different needs for plot size. Therefore,
# my adjustments allow for dynamic plot sizing.
make_gene_sj_expr_usage_plots_ts <- function(gene_of_interest, save_path) {
  #print gene name for error handling
  print(gene_of_interest)
  
  # run own function to get normalized expression
  gene_expr_df <- make_norm_gene_expr_df(gene_of_interest)
  
  # get sjs
  gene_sjs <-
    setbp1_marvel[["sj.metadata"]]$coord.intron[
      setbp1_marvel[["sj.metadata"]]$gene_short_name.start == gene_of_interest
    ]
  
  # get mean expression
  sj_mean_expr <- sj_expr_means[gene_sjs]
  
  # get gene splice junctions
  gene_junctions <- strsplit(gene_sjs, ":")
  
  gene_junctions <- as_tibble(do.call(rbind, gene_junctions))
  
  colnames(gene_junctions) <- c("seqnames", "start", "end")
  
  gene_junctions$start <- as.numeric(gene_junctions$start)
  
  gene_junctions$end <- as.numeric(gene_junctions$end)
  
  gene_junctions$strand <- gtf$strand[gtf$gene_name == gene_of_interest][1]
  
  gene_junctions$mean_count <- sj_mean_expr
  
  gene_junctions <- gene_junctions %>%
    dplyr::mutate(transcript_name = paste0(gene_of_interest, "-201"))
  
  
  # make split violin plot of gene expression
  split_violin <- ggplot(
    gene_expr_df,
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
      guide = guide_legend(
        override.aes =
          list(fill = "gray40", color = "black")
      )
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
    ylab("Normalized\nGene Expression")
  
  # get sj expr
  sj_mean_expr <- sj_expr_means[gene_sjs]
  
  # make splice junction expression heatmap
  delta_sje_heatmap_df <-
    delta_expr_for_heatmap[rownames(delta_expr_for_heatmap) %in% gene_sjs, ]
  
  print(gene_junctions$strand[1])
  
  ifelse(gene_junctions$strand[1] == "+", rownames(delta_sje_heatmap_df) <- c(paste0("SJ-", seq_len(nrow(gene_junctions)))),
         rownames(delta_sje_heatmap_df) <- c(paste0("SJ-", rev(seq_len(nrow(gene_junctions))))))
  
  # set colors
  delta_gene_expr_cols <-
    colorRamp2(
      c(
        min(delta_sje_heatmap_df), 0,
        max(delta_sje_heatmap_df)
      ),
      c(RColorBrewer::brewer.pal(name = "BrBG", n = 3))
    )
  
  # create heatmap
  delta_gene_heatmap <- Heatmap(delta_sje_heatmap_df,
                                name = "\u0394 Mean Norm\nSJ Expression",
                                column_order = sort(colnames(delta_sje_heatmap_df)),
                                col = delta_gene_expr_cols,
                                top_annotation = col_annot, cluster_rows = FALSE,
                                show_column_names = FALSE,
                                row_title = "SJ Expression", border = TRUE,
                                row_title_gp = gpar(fontface = "bold"),
                                row_names_gp = gpar(fontface = "bold")
  )
  
  # make splice junction usage heatmap
  delta_sju_heatmap_df <-
    delta_sju[rownames(delta_sju) %in% gene_sjs, ]
  
  ifelse(gene_junctions$strand[1] == "+", rownames(delta_sju_heatmap_df) <- c(paste0("SJ-", seq_len(nrow(gene_junctions)))),
         rownames(delta_sju_heatmap_df) <- c(paste0("SJ-", rev(seq_len(nrow(gene_junctions))))))
  
  # set colors
  delta_sju_cols <-
    colorRamp2(
      c(
        min(delta_sju_heatmap_df), 0,
        max(delta_sju_heatmap_df)
      ),
      c(RColorBrewer::brewer.pal(name = "PRGn", n = 3))
    )
  
  # make heatmap
  delta_sj_heatmap <- Heatmap(delta_sju_heatmap_df,
                              name = "\u0394 SJ Usage",
                              column_order = sort(colnames(delta_sju_heatmap_df)),
                              col = delta_sju_cols,
                              cluster_rows = FALSE,
                              show_column_names = FALSE,
                              row_title = "SJ Usage", border = TRUE,
                              row_title_gp = gpar(fontface = "bold"),
                              row_names_gp = gpar(fontface = "bold")
  )
  
  # Get number of splice junctions for conditional layout
  n_sjs <- nrow(gene_junctions)
  
  # Combine heatmaps - horizontally if >50 SJs, vertically otherwise
  if (n_sjs > 50) {
    col_annot_50sjs <- HeatmapAnnotation(
      `Cell Type` = names(cell_group_list),
      col = annot_colors_list, annotation_name_gp = gpar(fontface = "bold"),
      show_annotation_name = FALSE
    )
    
    # Recreate heatmaps with adjustments for horizontal layout
    delta_gene_heatmap <- Heatmap(delta_sje_heatmap_df,
                                  name = "\u0394 Mean Norm\nSJ Expression",
                                  column_order = sort(colnames(delta_sje_heatmap_df)),
                                  col = delta_gene_expr_cols,
                                  top_annotation = col_annot_50sjs, 
                                  cluster_rows = FALSE,
                                  show_column_names = FALSE,
                                  column_title = "SJ Expression",
                                  column_title_side = "bottom",
                                  border = TRUE,
                                  column_title_gp = gpar(fontface = "bold"),
                                  row_names_gp = gpar(fontface = "bold"),
                                  row_title = " ",
                                  row_title_side = "left",
                                  row_title_gp = gpar(fontsize = 12)
    )
    
    delta_sj_heatmap <- Heatmap(delta_sju_heatmap_df,
                                name = "\u0394 SJ Usage",
                                column_order = sort(colnames(delta_sju_heatmap_df)),
                                col = delta_sju_cols,
                                top_annotation = col_annot,
                                cluster_rows = FALSE,
                                show_column_names = FALSE,
                                column_title = "SJ Usage",
                                column_title_side = "bottom",
                                border = TRUE,
                                column_title_gp = gpar(fontface = "bold"),
                                row_names_gp = gpar(fontface = "bold")
    )
    
    both_delta_heatmaps <- delta_gene_heatmap + delta_sj_heatmap  # side-by-side
    
    # Adjust dimensions for horizontal layout
    dynamic_height <- max(7, 5 + (n_sjs * 0.15))
    dynamic_height <- min(dynamic_height, 20)
    dynamic_width <- max(18, 16 + (ncol(delta_sje_heatmap_df) * 2 * 0.25))
    dynamic_width <- min(dynamic_width, 35)
  } else {
    both_delta_heatmaps <- delta_gene_heatmap %v% delta_sj_heatmap  # stacked
    
    # Original dimensions for vertical layout
    dynamic_height <- max(10, 5 + (n_sjs * 0.25))
    dynamic_width <- max(14, 12 + (ncol(delta_sje_heatmap_df) * 0.25))
  }
  
  # make heatmaps compatible with cowplot
  delta_sj_expr_usage_heatmaps <- grid.grabExpr(draw(both_delta_heatmaps))
  
  # filter annotation for specific gene
  gene_annotation_from_gtf <- gtf %>%
    dplyr::filter(
      !is.na(gene_name),
      gene_name == gene_of_interest
    )
  
  gene_annotation_from_gtf <- gene_annotation_from_gtf %>%
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
  
  # extract exons
  gene_exons <- gene_annotation_from_gtf %>% dplyr::filter(type == "exon")
  
  ifelse(gene_junctions$strand[1] == "+", gene_junctions$sj_label <- c(paste0("SJ-", seq_len(nrow(gene_junctions)))),
         gene_junctions$sj_label <- c(paste0("SJ-", rev(seq_len(nrow(gene_junctions))))))
  
  # plot labeled transcripts with splice junctions
  gene_transcript_label <- gene_exons %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = transcript_name
    )) +
    geom_range(
      aes(fill = transcript_type)
    ) +
    geom_intron(
      data = to_intron(gene_exons, "transcript_name"),
      aes(strand = strand)
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold", color = "black"),
      legend.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      axis.title.y = element_blank()
    ) +
    labs(linewidth = "Normalized SJ\nExpression", fill = "Transcript Type") +
    scale_fill_manual(labels = c("NMD", "Protein Coding", "Protein Coding\nCDS Not Defined", "Retained Intron"),
                      values = c(viridis(6)[2:5])) +
    xlab("Genomic Location")
  
  # arrange vln plot
  gene_expr_vln <- plot_grid(gene_transcript_label, split_violin,
                             ncol = 1, rel_heights = c(0.66, 0.33),
                             labels = c("A", "B")
  )
  
  # arrange panels
  paneled_figure <- plot_grid(gene_expr_vln, delta_sj_expr_usage_heatmaps,
                              labels = c("", "C"), nrow = 1
  )
  
  # save
  png(paste0(save_path, "/", gene_of_interest, ".png"),
      width = dynamic_width, height = dynamic_height, units = "in", res = 300
  )
  print(paneled_figure)
  dev.off()
}
  
## split_sj_info_ts
# This function was originally written by Emma Jones and named split_sj_info. It can
# be found in src/marvel/functions.R
# Since the function was originally written for the brain data, it was not applicable
# to the kidney data and modifications were necessary. I removed code that edited
# barcodes, as it led to mismatch problems. Additionally, I added an if statement
# to treat K6 sample differently due to an outlier being manually removed.
split_sj_info_ts <- function(sample_id) {
  if (sample_id %in% c("K1", "K3", "K5")) {
    condition <- "WT"
  } else {
    condition <- "S858R"
  }
  # split the sample ID string into individual characters
  char_vector <- strsplit(sample_id, "")[[1]]
  # get cell barcodes
  sj_barcodes <- read.table(
    here::here(
      "data", "star", sample_id, "Solo.out", "SJ", "raw",
      "barcodes.tsv"
    )
  )
  colnames(sj_barcodes) <- "cell.id"
  sj_barcodes$cell.id <- paste0(
    condition, "_", sample_id, "_",
    sj_barcodes$cell.id
  )
  # subset barcodes and get order
  subset_barcodes <-
    gene_metadata[gene_metadata$sample_id == sample_id, "cell.id"]
  barcode_order <- match(subset_barcodes, sj_barcodes$cell.id)
  # subset the sj counts matrix
  if (sample_id == "K6") {
    sj_matrix <- readMM(
      here::here("data", "star", sample_id, "Solo.out", "SJ", "raw", "matrix_fixed.mtx")
    )
    sj_matrix <- sj_matrix[, barcode_order]
  } else {
    sj_matrix <- readMM(
      here::here("data", "star", sample_id, "Solo.out", "SJ", "raw", "matrix.mtx")
    )
    sj_matrix <- sj_matrix[, barcode_order]
  }
  # import sj features
  sj_features <- read.table(
    here::here(
      "data", "star", sample_id, "Solo.out", "SJ", "raw",
      "features.tsv"
    )
  )
  sj_features <- paste(sj_features$V1, sj_features$V2,
                       sj_features$V3,
                       sep = ":"
  )
  # make everything into a single dataframe
  colnames(sj_matrix) <- sj_barcodes[barcode_order, ]
  rownames(sj_matrix) <- sj_features
  
  # export data in long dataframe format for that sample
  return(sj_matrix)
}


## run_marvel_cell_type_kidney
# Adapated from the run_marvel_cell_type function by Emma Jones
# Original function and description can be found in src/marvel/functions.R
# Due to differing cell type names in the kidney, I changed the way the results
# are saved by creating a modified cell type name to avoid spaces in final objects
run_marvel_cell_type_kidney <- function(marvel_object, cell_type, min_pct_cells = 5, 
                                        min_pct_cells_gene = 5, min_pct_cells_sj = 5,
                                        min_gene_norm = 1, results_path) {
  
  # Assign MARVEL object to start from
  marvel_object <- setbp1_marvel
  
  # Group 1 (reference)
  index_1 <- which(sample_metadata$cell_type == cell_type & 
                     sample_metadata$seq_folder == "wildtype")
  cell_ids_1 <- sample_metadata[index_1, "cell.id"]
  
  # Group 2
  index_2 <- which(sample_metadata$cell_type == cell_type & 
                     sample_metadata$seq_folder == "mutant")
  cell_ids_2 <- sample_metadata[index_2, "cell.id"]
  
  # Explore % of cells expressing genes
  marvel_object <- PlotPctExprCells.Genes.10x(
    MarvelObject = marvel_object,
    cell.group.g1 = cell_ids_1,
    cell.group.g2 = cell_ids_2,
    min.pct.cells = min_pct_cells
  )
  
  # Explore % of cells expressing junctions
  marvel_object <- PlotPctExprCells.SJ.10x(
    MarvelObject = marvel_object,
    cell.group.g1 = cell_ids_1,
    cell.group.g2 = cell_ids_2,
    min.pct.cells.genes = min_pct_cells_gene,
    min.pct.cells.sj = min_pct_cells_sj,
    downsample = TRUE,
    downsample.pct.sj = 10
  )
  
  # Differential Splicing Analysis
  marvel_object <- CompareValues.SJ.10x(
    MarvelObject = marvel_object,
    cell.group.g1 = cell_ids_1,
    cell.group.g2 = cell_ids_2,
    min.pct.cells.genes = min_pct_cells_gene,
    min.pct.cells.sj = min_pct_cells_sj,
    min.gene.norm = min_gene_norm,
    seed = 1,
    n.iterations = 100,
    downsample = TRUE,
    show.progress = TRUE
  )
  
  # Differential Gene Analysis
  marvel_object <- CompareValues.Genes.10x(
    MarvelObject = marvel_object,
    show.progress = TRUE
  )
  
  # Make volcano plot
  marvel_object <- PlotDEValues.SJ.10x(
    MarvelObject = marvel_object,
    pval = 0.05,
    delta = 1,
    min.gene.norm = min_gene_norm,
    anno = FALSE
  )
  # Assign kinds of iso-switching
  marvel_object <- IsoSwitch.10x(
    MarvelObject = marvel_object,
    pval.sj = 0.05,
    delta.sj = 1,
    min.gene.norm = min_gene_norm,
    pval.adj.gene = 0.05,
    log2fc.gene = 0.5
  )
  
  # Pull significant genes
  significant_genes <- marvel_object[["SJ.Gene.Cor"]][["Data"]]$gene_short_name
  
  modified_cell_type <- gsub(" ", "-", cell_type)
  
  # Save MARVEL object
  write_rds(marvel_object, file = paste0(results_path, "/", modified_cell_type,
                                         "_marvel_object.rds")
  )
  
  # Return list
  return(marvel_object)
}

## to_col_format
# Function that helps formatting cell type names
to_col_format <- function(x) {
  tolower(gsub(" ", "_", x))
}
