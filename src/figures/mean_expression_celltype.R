# Install and load Seurat if not already installed

library(Seurat)

# Assuming seurat_obj is your Seurat object and cell types are in the metadata column named 'cell_type'
# Extract normalized expression data
expression_data <- GetAssayData(annotated_brain_samples, slot = "data")

# Extract cell type annotations
cell_types <- annotated_brain_samples@meta.data$cell_type

# Compute mean expression by cell type
mean_expression <- aggregate(t(as.matrix(expression_data)), by = list(cell_type = cell_types), FUN = mean)

# Clean up the resulting data frame
rownames(mean_expression) <- mean_expression$cell_type
mean_expression <- mean_expression[, -1]
mean_expression <- t(mean_expression)
