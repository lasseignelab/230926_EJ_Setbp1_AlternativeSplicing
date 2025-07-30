# The following function is from https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html
# Edited by Emma Jones
# Function to run DESeq2 Wald Test and get results for any cluster:
## clustx is the name of the cluster (cell type) on which to run the function
## A is the sample group to compare (e.g. stimulated condition)
## B is the sample group to compare against (base/control level)
## padj_cutoff defines the adjusted p-value cutoff for significance (set to 0.05 by default)

## This function assumes the counts matrices and metadata for all clusters have been prepared
## and arranged in matching named lists (as illustrated in tutorial above)
## This function assumes the contrast (e.g. stim vs. control) is stored in a variable named "group_id"

get_dds_resultsAvsB <- function(clustx, A, B, padj_cutoff = 0.05, save_path) {
  
  print(clustx) # useful for debugging
  
  # Extract counts matrix and metadata for cluster x
  idx <- which(names(counts_ls) == clustx)
  cluster_counts <- counts_ls[[idx]]
  cluster_metadata <- metadata_ls[[idx]]
  
  # Print error message if sample names do not match
  if ( all(colnames(cluster_counts) != rownames(cluster_metadata)) ) {
    print("ERROR: sample names in counts matrix columns and metadata rows do not match!")
  }
  
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ group_id)
  
  # Transform counts for data visualization
  rld <- rlog(dds, blind = TRUE)
  
  # Modify cluster name for saving
  clustx_modified <- gsub(" ", "_", clustx)
  
  # Generate QC plots
  
  ## Plot and save PCA plot
  DESeq2::plotPCA(rld, intgroup = "group_id")
  ggsave(here::here(paste0(save_path, clustx_modified, "_specific_PCAplot.png")))
  
  ## Extract rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  
  ## Plot and save heatmap
  png(here::here(paste0(save_path, clustx_modified, "_specific_heatmap.png")),
      height = 6, width = 7.5, units = "in", res = 300)
  pheatmap(rld_cor, annotation = cluster_metadata[, c("group_id"), drop = FALSE])
  dev.off()
  
  
  # Run DESeq2 differential expression analysis
  dds <- DESeq(dds)
  
  ## Plot dispersion estimates
  png(here::here(paste0(save_path, clustx_modified, "_dispersion_plot.png")),
      height = 5, width = 6, units = "in", res = 300)
  plotDispEsts(dds)
  dev.off()
  
  ## Output and shrink results of Wald test for contrast A vs B
  contrast <- paste(c("group_id", A, "vs", B), collapse = "_")
  
  res <- results(dds, name = contrast, alpha = 0.05)
  res <- lfcShrink(dds, coef = contrast, res = res)
  
  ## Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()
  
  write.csv(res_tbl,
            here::here(paste0(save_path, clustx_modified, "_", contrast, "_all_genes.csv")),
            quote = FALSE, 
            row.names = FALSE)
  
  ## Subset the significant results
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
    dplyr::arrange(padj)
  
  write.csv(sig_res,
            here::here(paste0(save_path, clustx_modified, "_", contrast, "_signif_genes.csv")),
            quote = FALSE, 
            row.names = FALSE)
  
  
  # Generate results visualization plots
  
  ## Extract normalized counts from dds object
  normalized_counts <- counts(dds, normalized = TRUE)
  
  ## Extract top 20 DEG from resLFC (make sure to order by padj)
  top20_sig_genes <- sig_res %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(gene) %>%
    head(n = 20)
  
  # added step to only proceed with following steps if you have 20 genes
  if(length(top20_sig_genes) == 20)
  
  ## Extract matching normalized count values from matrix
  {top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]
  
  ## Convert wide matrix to long data frame for ggplot2
  top20_sig_df <- data.frame(top20_sig_counts)
  top20_sig_df$gene <- rownames(top20_sig_counts)
  
  top20_sig_df <- melt(setDT(top20_sig_df), 
                       id.vars = c("gene"),
                       variable.name = "cluster_sample_id") %>% 
    data.frame()
  
  ## Replace "." by " " in cluster_sample_id variable (melt() introduced the ".")
  top20_sig_df$cluster_sample_id <- gsub("\\.", " ", top20_sig_df$cluster_sample_id)
  top20_sig_df$cluster_sample_id <- gsub("\\  ", "+ ", top20_sig_df$cluster_sample_id)
  
  ## Join counts data frame with metadata
  top20_sig_df <- plyr::join(top20_sig_df, as.data.frame(colData(dds)),
                             by = "cluster_sample_id")
  
  ## Generate plot
  ggplot(top20_sig_df, aes(y = value, x = group_id, col = group_id)) +
    geom_jitter(height = 0, width = 0.15) +
    scale_y_continuous(trans = 'log10') +
    ylab("log10 of normalized expression level") +
    xlab("condition") +
    ggtitle("Top 20 Significant DE Genes") +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~ gene)
  
  ggsave(here::here(paste0(save_path, clustx_modified, "_", contrast, "_top20_DE_genes.png")))}
  
  # how many genes are there?
  number <- length(top20_sig_genes)
  
  # added step to only plot if you have between 1 and 19 genes
  if(between(number, 1, 19))
    
    ## Extract matching normalized count values from matrix
  {top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]
  
  ## Convert wide matrix to long data frame for ggplot2
  top20_sig_df <- data.frame(top20_sig_counts)
  top20_sig_df$gene <- rownames(top20_sig_counts)
  
  top20_sig_df <- melt(setDT(top20_sig_df), 
                       id.vars = c("gene"),
                       variable.name = "cluster_sample_id") %>% 
    data.frame()
  
  ## Replace "." by " " in cluster_sample_id variable (melt() introduced the ".")
  top20_sig_df$cluster_sample_id <- gsub("\\.", " ", top20_sig_df$cluster_sample_id)
  top20_sig_df$cluster_sample_id <- gsub("\\  ", "+ ", top20_sig_df$cluster_sample_id)
  
  ## Join counts data frame with metadata
  top20_sig_df <- plyr::join(top20_sig_df, as.data.frame(colData(dds)),
                             by = "cluster_sample_id")
  
  ## Generate plot
  ggplot(top20_sig_df, aes(y = value, x = group_id, col = group_id)) +
    geom_jitter(height = 0, width = 0.15) +
    scale_y_continuous(trans = 'log10') +
    ylab("log10 of normalized expression level") +
    xlab("condition") +
    ggtitle(paste0("Top", number, "Significant DE Genes")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~ gene)
  
  ggsave(here::here(paste0(save_path, clustx_modified, "_", contrast, "_top", number, "_DE_genes.png")))}
  
}