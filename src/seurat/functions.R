### Seurat functions for Setbp1 Alternative Splicing Project
# Emma F. Jones (EJ)
# Note: These functions are borrowed/modified from Tabea Soelter
# Those modified functions will be marked accordingly.

## calculate_qc - c/o Tabea Soelter, modified by EJ
# A function which calculates quality control metrics for a merged seurat object
calculate_qc <- function(seurat_object){
  seurat_object$log10GenesPerUMI <- log10(seurat_object$nFeature_RNA) / log10(seurat_object$nCount_RNA)
  seurat_object$mitoRatio <- seurat_object[["percent_mt"]] / 100
  return(seurat_object)
}

## format_metadata - c/o Tabea Soelter, modified by EJ
# A function which extracts and formats metadata of a seurat object
format_metadata <- function(seurat_object){
  metadata <- seurat_object@meta.data
  metadata$cells <- rownames(metadata)
  # Rename columns -----
  metadata <- metadata %>% dplyr::rename(seq_folder = condition,
                                         nUMI = nCount_RNA,
                                         nGene = nFeature_RNA)
  return(metadata)
}

## plot_qc - c/o Tabea Soelter
# A function which takes seurat metadata and plots quality control metrics for filtering purposes
plot_qc <- function(metadata) {
  # Visualize the number of cell counts per condition
  number_of_cells <- metadata %>% 
    ggplot(aes(x = seq_folder, fill = seq_folder)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle("NCells") 
  # Visualize the number UMIs/transcripts per cell
  number_of_umis <- metadata %>% 
    ggplot(aes(color = seq_folder, x = nUMI, fill = seq_folder)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)
  # Visualize the distribution of genes detected per cell
  dist_genes_per_cell <- metadata %>% 
    ggplot(aes(color = seq_folder, x = nGene, fill = seq_folder)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = 300)
  # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
  novelty_score <- metadata %>%
    ggplot(aes(x = log10GenesPerUMI, color = seq_folder, fill = seq_folder)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8)
  # Visualize the distribution of mitochondrial gene expression detected per cell
  dist_mito_gex <- metadata %>% 
    ggplot(aes(color = seq_folder, x = mitoRatio, fill = seq_folder)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = 0.2)
  # Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells   with low numbers of genes/UMIs
  cor <- metadata %>% 
    ggplot(aes(x = nUMI, y = nGene, color = mitoRatio)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method = lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 250) +
    facet_wrap(~seq_folder)
  # Plot QC metrics
  plot(number_of_cells) 
  plot(number_of_umis)
  plot(dist_genes_per_cell)
  plot(novelty_score)
  plot(dist_mito_gex)
  plot(cor)
}

## convert_human_gene_list - c/o Tabea Soelter
# Source: https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/ 
# I adapted the function to use an archived version of ensembl, as there were mirror issues due to update in Feb 2023
# A function to convert human to mouse gene names 
convert_human_gene_list <- function(x) {
  require("biomaRt")
  human = useMart("ensembl",
                  dataset = "hsapiens_gene_ensembl",
                  host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl",
                  dataset = "mmusculus_gene_ensembl",
                  host = "https://dec2021.archive.ensembl.org/")
  genesV2 = getLDS(attributes = c("hgnc_symbol"),
                   filters = "hgnc_symbol",
                   values = x ,
                   mart = human,
                   attributesL = c("mgi_symbol"),
                   martL = mouse,
                   uniqueRows = T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}

## evaluate_cell_cycle is a modified version of 
## cell_cycle_effects() - c/o Tabea Soelter, modified by EJ
# A function which calculates and plots the effect of cell cycle on the data using a filtered seurat object as input. It also performs log normalization, scaling, and dimension reduction using PCA
evaluate_cell_cycle <- function(filtered_seurat){
  # log normalize -----
  filtered_seurat <- NormalizeData(filtered_seurat)
  # convert human cell cycle markers to mouse -----
  s.genes <- convert_human_gene_list(cc.genes.updated.2019$s.genes)
  g2m.genes <- convert_human_gene_list(cc.genes.updated.2019$g2m.genes)
  # EJ V5 BUG FIX: JOIN LAYERS
  filtered_seurat <- JoinLayers(filtered_seurat)
  # score cells based in gex of genes -----
  filtered_seurat <- CellCycleScoring(filtered_seurat,
                                      g2m.features = g2m.genes,
                                      s.features = s.genes)
  filtered_seurat <- FindVariableFeatures(filtered_seurat,
                                          selection.method = "vst",
                                          verbose = FALSE)
  # scale data -----
  filtered_seurat <- ScaleData(filtered_seurat)
  # run pca -----
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

## find_clusters - c/o Tabea Soelter
# A function which finds clusters for a Seurat Object at user determined resolutions
find_clusters <- function(object, dims, reduction, resolutions) {
  # set reduction method to harmony -----
  object <- FindNeighbors(object, 
                          dims = dims, 
                          reduction = reduction)
  # clustering (Leiden aka algorithm 4)
  for (res in resolutions) {
    object <- FindClusters(object,
                           graph.name = "RNA_snn",
                           resolution = res,
                           algorithm = 4,
                           method = "igraph")
  } 
  return(object)
}