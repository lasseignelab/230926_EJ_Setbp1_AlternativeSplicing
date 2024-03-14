### MARVEL functions for Setbp1 Alternative Splicing Project
# Emma F. Jones (EJ)

## split_sj_info - Emma Jones
# This function gets the splice junction information for each sample
split_sj_info <- function(sample_id, condition) {
  # split the sample ID string into individual characters
  char_vector <- strsplit(sample_id, "")[[1]]
  # add an underscore after the first character
  sample_id_mod <- paste0(char_vector[1], "_",
                          paste0(char_vector[-1], collapse = ""))
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
# within a given cell type.
# The thresholds may need adjusting, but have defaults.
run_marvel_cell_type <- function(marvel_object, cell_type, min_pct_cells = 5, 
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
  
  # Save MARVEL object
  write_rds(marvel_object, file = paste0(results_path, "/", cell_type,
                                         "_marvel_object.rds")
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
                                                 "yellow", "red"),
                             results_path) {
  # Group 1 (reference)
  index_1 <- which(sample_metadata$seq_folder == "wildtype")
  cell_ids_1 <- sample_metadata[index_1, "cell.id"]
  
  # Group 2
  index_2 <- which(sample_metadata$seq_folder == "mutant")
  cell_ids_2 <- sample_metadata[index_2, "cell.id"]
  
  # Save into list
  cell_group_list <- list(
    "Wildtype" = cell_ids_1,
    "Mutant" = cell_ids_2
  )
  
  # Print generation message
  message("Generating Plots...")
  
  # Plot cell groups
  marvel_object <- PlotValues_PCA_CellGroup_EJ(
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
  
  plot_gene <- marvel_object$adhocPlot$PCA$Gene
  
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
  png(paste0(results_path, "/", prefix, "_", gene, ".png"),
      width = 1200, height = 600)
  grid.arrange(plot_group, plot_gene, plot_sj, nrow = 1)
  dev.off()
  
  # Return object
  return(marvel_object)
}

### EJ EDIT OF MARVEL PLOTTING FUNCTION TO SHUFFLE POINTS
# Most code is copied directly from PlotValues_PCA_CellGroup
PlotValues_PCA_CellGroup_EJ <-
function (MarvelObject, cell.group.list, legendtitle = "Cell group", 
          alpha = 0.75, point.size = 1, point.stroke = 0.1, point.colors = NULL, 
          point.size.legend = 2, type, shuffle = TRUE) 
{
  MarvelObject <- MarvelObject
  df <- MarvelObject$pca
  cell.group.list <- cell.group.list
  legendtitle <- legendtitle
  alpha <- alpha
  point.size <- point.size
  point.stroke <- point.stroke
  point.colors <- point.colors
  point.size.legend <- point.size.legend
  type <- type
  .list <- list()
  for (i in 1:length(cell.group.list)) {
    . <- data.frame(cell.id = cell.group.list[[i]],
                    group = names(cell.group.list)[i])
    .list[[i]] <- .
  }
  ref <- do.call(rbind.data.frame, .list)
  df <- join(df, ref, by = "cell.id", type = "left")
  df$group <- factor(df$group, levels = names(cell.group.list))
  n.anno.missing <- sum(is.na(df$group))
  if (n.anno.missing == 0) {
    message("All cells defined with coordinates found")
  }
  else {
    message(paste(n.anno.missing, " cells defined with no coordinates found", 
                  sep = ""))
  }
  data <- df
  # EJ ADDITION TO FUNCTION HERE
  if (isTRUE(x = shuffle)) {
    data <- data[sample(x = 1:nrow(x = data)), ]
  }
  x <- data$x
  y <- data$y
  z <- data$group
  if (type == "umap") {
    xtitle <- "UMAP-1"
    ytitle <- "UMAP-2"
  }
  else if (type == "tsne") {
    xtitle <- "tSNE-1"
    ytitle <- "tSNE-2"
  }
  if (is.null(point.colors[1])) {
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    n = length(levels(z))
    point.colors = gg_color_hue(n)
  }
  if (length(unique(z)) > 20) {
    ncol.legend <- 3
  }
  else if (length(unique(z)) > 10) {
    ncol.legend <- 2
  }
  else if (length(unique(z)) <= 10) {
    ncol.legend <- 1
  }
  plot <- ggplot() + 
    geom_point(data, mapping = aes(x = x, y = y, fill = z),
               color = "black", pch = 21, size = point.size, 
                                alpha = alpha, stroke = point.stroke) + 
    scale_fill_manual(values = point.colors) + 
    labs(title = NULL, x = xtitle, y = ytitle, fill = legendtitle) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          plot.title = element_text(size = 12, hjust = 0.5),
          axis.line = element_line(colour = "black"), 
          axis.title = element_text(size = 12),
          axis.text.x = element_text(size = 10, colour = "black"),
          axis.text.y = element_text(size = 10, colour = "black"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8)) +
    guides(fill = guide_legend(override.aes = 
                                 list(size = point.size.legend, alpha = alpha,
                                      stroke = point.stroke),
                               ncol = ncol.legend))
  MarvelObject$adhocPlot$PCA$CellGroup <- plot
  return(MarvelObject)
}

## subset_sparse_matrix - Emma Jones
# The purpose of this function is to subset a sparse matrix by cell type.
# It also drop rows that are all zeroes.
subset_sparse_matrix <- function(matrix, cell_type){
  # subset matrix based on cell type
  subset_matrix <- matrix[, colnames(matrix) %in% cell_group_list[[cell_type]]]
  # return matrix subset
  return(subset_matrix)
}

## subset_cell_type - Emma Jones
# This function uses subset_sparse_matrix 2 times for gene, splice junction,
# for a given cell type. The default args are already in the
# environment and are the full sized versions of the counts matrices.
subset_cell_type_matrices <- function(gene_matrix = gene_counts,
                                      sj_matrix = sj_counts, cell_type) {
  # run function
  gene_subset <- subset_sparse_matrix(gene_matrix, cell_type)
  sj_subset <- subset_sparse_matrix(sj_matrix, cell_type)
  # make into single list
  list <- list(gene_subset, sj_subset)
  # name list
  names(list) <- c("Gene Counts", "Splice Junction Counts")
  # return object
  return(list)
}

## subset_mutant_matrices - Emma Jones
# The purpose of this function is to subset a list of sparse matrices by
# condition. It also drop rows that are all zeroes.
subset_mutant_matrices <- function(cell_type_subset_list){
  ## start with mutant
  # initialize empty list
  subset_matrices_list_1 <- list()
  # loop through each list
  for(i in 1:length(cell_type_subset_list)){
    # set matrix
    matrix <- cell_type_subset_list[[i]]
    # subset matrix based on cell type
    subset_matrix <- matrix[, colnames(matrix) %in% mutant_list[["Mutant"]]]
    # assign element in list
    subset_matrices_list_1[[i]] <- subset_matrix
  }
  # rename list elements
  names(subset_matrices_list_1) <- c("Gene Counts", "Splice Junction Counts")
  ## now for wildtype
  # initialize empty list
  subset_matrices_list_2 <- list()
  # loop through each list
  for(i in 1:length(cell_type_subset_list)){
    # set matrix
    matrix <- cell_type_subset_list[[i]]
    # subset matrix based on cell type
    subset_matrix <- matrix[, colnames(matrix) %in% mutant_list[["Wildtype"]]]
    # assign element in list
    subset_matrices_list_2[[i]] <- subset_matrix
  }
  # rename list elements
  names(subset_matrices_list_2) <- c("Gene Counts", "Splice Junction Counts")
  ## combine lists
  subset_matrices_lists <- list(subset_matrices_list_1, subset_matrices_list_2)
  # rename list elements
  names(subset_matrices_lists) <- c("Mutant", "Wildtype")
  # return matrix subset
  return(subset_matrices_lists)
}

## get_sjs_per_gene - Emma Jones
# The purpose of this function is to get the number of splice junctions
# detected for each gene for a certain cell type
get_sjs_per_gene <- function(cell_type_matrices) {
  # assign matrix
  sj_matrix <- cell_type_matrices[["Splice Junction Counts"]]
  # drop rows with zero
  sj_matrix_drop0 <- sj_matrix[rowSums(sj_matrix) > 0, ]
  # pull names of splice junctions
  detected_sjs <- rownames(sj_matrix_drop0)
  # get gene names
  detected_genes <- sj_metadata$gene_short_name.start[sj_metadata$coord.intron %in% detected_sjs]
  # make a table of gene name values
  sjs_per_gene <- table(detected_genes)
  # return
  return(sjs_per_gene)
}


## remove_0_variance - quick method for remove zero variance 
# pulled from https://www.ttested.com/removing-zero-variance-columns/
remove_0_variance <- function(df){
  df[, !sapply(df, function(x) min(x) == max(x))]
}

## get_sj_usage_cell_type - Emma Jones
# The purpose of this function is to calculate splice junction usage for each
# cell type.
get_sj_usage_cell_type <- function(counts_matrix_list){
  # get sums of all gene counts for a given cell type
  gene_sums <- rowSums(counts_matrix_list[["Gene Counts"]])
  # get sums of all splice junction counts for a given cell type
  sj_sums <- rowSums(counts_matrix_list[["Splice Junction Counts"]])
  # inititalize list
  sj_usage_cell_type <- c()
  # run for loop for calculating splice junction usage
  for(sj in names(sj_sums)){
    sj_count_sum <- sj_sums[sj]
    gene_name <- sj_metadata$gene_short_name.start[sj == sj_metadata$coord.intron]
    gene_count_sum <- gene_sums[gene_name]
    sj_usage <- (sj_count_sum / gene_count_sum) * 100
    sj_usage_cell_type[sj] <- sj_usage
  }
  # return value
  return(sj_usage_cell_type)
}

## get_sj_usage_condition - Emma Jones
# The purpose of this function is to calculate splice junction usage for each
# condition, mutant or wildtype, for every cell type.
get_sj_usage_condition <- function(counts_matrix_list){
  ## mutants ##
  counts_matrix_list_1 <- counts_matrix_list[["Mutant"]]
  # get sums of all gene counts for a given cell type
  gene_sums_1 <- rowSums(counts_matrix_list_1[["Gene Counts"]])
  # get sums of all splice junction counts for a given cell type
  sj_sums_1 <- rowSums(counts_matrix_list_1[["Splice Junction Counts"]])
  # inititalize list
  sj_usage_condition_1 <- c()
  # run for loop for calculating splice junction usage
  for(sj in names(sj_sums_1)){
    sj_count_sum <- sj_sums_1[sj]
    gene_name <- sj_metadata$gene_short_name.start[sj == sj_metadata$coord.intron]
    gene_count_sum <- gene_sums_1[gene_name]
    sj_usage <- (sj_count_sum / gene_count_sum) * 100
    sj_usage_condition_1[sj] <- sj_usage
  }
  ## wildtypes ##
  counts_matrix_list_2 <- counts_matrix_list[["Wildtype"]]
  # get sums of all gene counts for a given cell type
  gene_sums_2 <- rowSums(counts_matrix_list_2[["Gene Counts"]])
  # get sums of all splice junction counts for a given cell type
  sj_sums_2 <- rowSums(counts_matrix_list_2[["Splice Junction Counts"]])
  # inititalize list
  sj_usage_condition_2 <- c()
  # run for loop for calculating splice junction usage
  for(sj in names(sj_sums_2)){
    sj_count_sum <- sj_sums_2[sj]
    gene_name <- sj_metadata$gene_short_name.start[sj == sj_metadata$coord.intron]
    gene_count_sum <- gene_sums_2[gene_name]
    sj_usage <- (sj_count_sum / gene_count_sum) * 100
    sj_usage_condition_2[sj] <- sj_usage
  }
  # final result
  sj_usage_result <- list(sj_usage_condition_1, sj_usage_condition_2)
  # name lists
  names(sj_usage_result) <- c("Mutant", "Wildtype")
  # return value
  return(sj_usage_result)
}
