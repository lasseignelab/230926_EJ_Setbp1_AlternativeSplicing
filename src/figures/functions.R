# The purpose of this script is to provide any functions needed for generating
# figures
# Emma Jones
# Created March 4, 2024

wilcoxon_test_celltype <- function(cell_type) {
  result <- wilcox.test(mutant_metadata$num_sjs_genes[mutant_metadata$cell_type == cell_type], wildtype_metadata$num_sjs_genes[wildtype_metadata$cell_type == cell_type])
  
  return(result[["p.value"]])
}