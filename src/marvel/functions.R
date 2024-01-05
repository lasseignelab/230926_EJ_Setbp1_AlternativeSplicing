### MARVEL functions for Setbp1 Alternative Splicing Project
# Emma F. Jones (EJ)
## split_sj_info - Emma Jones
# This function gets the splice junction information for each sample
split_sj_info <- function(sample_id, condition, sample_id_mod) {
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
# Function to align matrix rows
align_matrix <- function(mat, all_features) {
  # Create an empty sparse matrix with all features
  aligned_mat <- emptySparse(nrow = length(all_features), ncol = ncol(mat))
  rownames(aligned_mat) <- all_features

  # Fill in existing values
  intersecting_features <- intersect(rownames(mat), all_features)
  aligned_mat[intersecting_features, ] <- mat[intersecting_features, ]
  return(aligned_mat)
}
