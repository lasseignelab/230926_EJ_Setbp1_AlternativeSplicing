### MARVEL functions for Setbp1 Alternative Splicing Project
# Emma F. Jones (EJ)
## split_sj_info - Emma Jones
# This function gets the splice junction information for each sample
split_sj_info <- function(sample_id, condition) {
  # split the sample ID string into individual characters
  char_vector <- strsplit(sample_id, "")[[1]]
  # add an underscore after the first character
  sample_id_mod <- paste0(char_vector[1], "_", paste0(char_vector[-1], collapse = ""))
  # get cell barcodes
  sj_barcodes <- read.table(
    here::here(
      "data", "star", sample_id, "Solo.out", "SJ", "raw",
      "barcodes.tsv"
    )
  )
  colnames(sj_barcodes) <- "cell.id"
  sj_barcodes$cell.id <- paste0(
    sample_id, "_", condition, "_",
    sj_barcodes$cell.id
  )
  # subset barcodes and get order
  subset_barcodes <-
    gene_metadata[gene_metadata$sample_id == sample_id_mod, "cell.id"]
  barcode_order <- match(subset_barcodes, sj_barcodes$cell.id)
  # subset the sj counts matrix
  sj_matrix <- readMM(
    here::here("data", "star", sample_id, "Solo.out", "SJ", "raw", "matrix.mtx")
  )
  sj_matrix <- sj_matrix[, barcode_order]
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

## align_matrix - Emma Jones
# This function is made to align sparse matrix rows so they can be combined
align_matrix <- function(mat, all_features) {
  # Create an empty sparse matrix with all features
  aligned_mat <- emptySparse(nrow = length(all_features), ncol = ncol(mat))
  rownames(aligned_mat) <- all_features

  # Fill in existing values
  intersecting_features <- intersect(rownames(mat), all_features)
  aligned_mat[intersecting_features, ] <- mat[intersecting_features, ]
  return(aligned_mat)
}

## run_marvel_cell_type - Emma Jones
# This function is a wrapper function for comparing mutant and wildtype splicing
# within a given cell type. The thresholds may need adjusting, but have defaults.
run_marvel_cell_type <- function(marvel_object, cell_type, min_pct_cells = 5, 
                       min_pct_cells_gene = 5, min_pct_cells_sj = 5,
                       min_gene_norm = 1) {
  
  # Assign MARVEL object to start from
  marvel_object <- setbp1_marvel
  
  # Group 1 (reference)
  index_1 <- which(sample_metadata$cell_type == cell_type & 
                     sample_metadata$seq_folder == "mutant")
  cell_ids_1 <- sample_metadata[index_1, "cell.id"]
  
  # Group 2
  index_2 <- which(sample_metadata$cell_type == cell_type & 
                     sample_metadata$seq_folder == "wildtype")
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
  significant_genes <- marvel_object[["SJ.Gene.Cor"]][["Data"]]
  
  # Return list
  return(list(marvel_object, significant_genes))
}

## run_marvel_cell_type - Emma Jones
# This function is a wrapper function for comparing mutant and wildtype splicing
# within a given cell type. The thresholds may need adjusting, but have defaults.
run_marvel_cell_type <- function(marvel_object, cell_type, min_pct_cells = 5, 
                       min_pct_cells_gene = 5, min_pct_cells_sj = 5,
                       min_gene_norm = 1) {
  
  # Assign MARVEL object to start from
  marvel_object <- setbp1_marvel
  
  # Group 1 (reference)
  index_1 <- which(sample_metadata$cell_type == cell_type & 
                     sample_metadata$seq_folder == "mutant")
  cell_ids_1 <- sample_metadata[index_1, "cell.id"]
  
  # Group 2
  index_2 <- which(sample_metadata$cell_type == cell_type & 
                     sample_metadata$seq_folder == "wildtype")
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
  
  # Save MARVEL object
  write_rds(marvel_object,
          file = here::here("data", "marvel",
                            paste0(cell_type, "_marvel_object.rds"))
  )
  
  # Return list
  return(marvel_object)
}

## plot_marvel_umap - Emma Jones
# This function a wrapper function for plotting marvel umaps for a cell type,
# gene, and splice junction locus. User can also customize the colors used.
plot_marvel_umap <- function(marvel_object, prefix, gene, sj_loc,
                             color_grad_gene = c("grey", "cyan", "green",
                                                 "yellow", "red"),
                             color_grad_psi = c("grey", "cyan", "green",
                                                 "yellow", "red")) {
  # Group 1 (reference)
  index_1 <- which(sample_metadata$seq_folder == "mutant")
  cell_ids_1 <- sample_metadata[index_1, "cell.id"]
  
  # Group 2
  index_2 <- which(sample_metadata$seq_folder == "wildtype")
  cell_ids_2 <- sample_metadata[index_2, "cell.id"]
  
  # Save into list
  cell_group_list <- list(
    "Mutant" = cell_ids_1,
    "Wildtype" = cell_ids_2
  )
  
  # Print generation message
  message("Generating Plots...")
  
  # Plot cell groups
  marvel_object <- PlotValues.PCA.CellGroup.10x(
    MarvelObject = marvel_object,
    cell.group.list = cell_group_list,
    legendtitle = "Cell group",
    type = "umap"
  )
  
  plot_group <- marvel_object$adhocPlot$PCA$CellGroup
  
  # Plot gene expression
  marvel_object <- PlotValues.PCA.Gene.10x(
    MarvelObject = marvel_object,
    gene_short_name = gene,
    color.gradient = color_grad_gene,
    type = "umap"
  )
  
  
  plot_gene <- excitatory_marvel$adhocPlot$PCA$Gene
  
  # Plot PSI
  marvel_object <- PlotValues.PCA.PSI.10x(
    MarvelObject = marvel_object,
    coord.intron = sj_loc,
    min.gene.count = 3,
    log2.transform = FALSE,
    color.gradient = color_grad_psi,
    type = "umap"
  )
  
  plot_sj <- marvel_object$adhocPlot$PCA$PSI
  
  # Print plotting message
  message("Saving Plot...")
  
  # Save UMAP
  png(here::here("results", "marvel_outputs",
                 paste0(prefix, "_", gene, ".png")),
      width = 1200, height = 600)
  grid.arrange(plot_group, plot_gene, plot_sj, nrow = 1)
  dev.off()
  
  # Return object
  return(marvel_object)
}

