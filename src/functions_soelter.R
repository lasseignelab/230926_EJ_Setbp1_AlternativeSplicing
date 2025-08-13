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
make_seurat_object <- function(path){
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
  return(merged_seurat)
}


## calc_qc
# A function which calculates quality control metrics for a merged seurat object
calc_qc <- function(seurat_object){
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
cell_cycle_effects <- function(filtered_seurat, g2m_genes, s_genes){
  # log normalize -----
  filtered_seurat <- NormalizeData(filtered_seurat)
  # score cells based in gex of genes -----
  filtered_seurat <- CellCycleScoring(filtered_seurat,
                                      g2m.features = g2m_genes,
                                      s.features = s_genes)
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
