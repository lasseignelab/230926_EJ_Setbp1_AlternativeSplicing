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

