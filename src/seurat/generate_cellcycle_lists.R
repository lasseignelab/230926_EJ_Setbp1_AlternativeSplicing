convert_human_gene_list <- function(gene_list) {
  require("biomaRt")
  # get human biomart results
  human <- useMart("ensembl",
                   dataset = "hsapiens_gene_ensembl",
                   host = "https://dec2021.archive.ensembl.org/"
  )
  
  # get mouse biomart results
  mouse <- useMart("ensembl",
                   dataset = "mmusculus_gene_ensembl",
                   host = "https://dec2021.archive.ensembl.org/"
  )
  
  # get query biomart for matching genes
  genes_lds <- getLDS(
    attributes = c("hgnc_symbol"),
    filters = "hgnc_symbol",
    values = gene_list,
    mart = human,
    attributesL = c("mgi_symbol"),
    martL = mouse,
    uniqueRows = TRUE
  )
  # pull only converted genes of interest
  converted_genes <- unique(genes_lds[, 2])
  # return converted genes of interest
  return(converted_genes)
}

s.genes <- convert_human_gene_list(cc.genes.updated.2019$s.genes)
g2m.genes <- convert_human_gene_list(cc.genes.updated.2019$g2m.genes)

wd <- "/data/user/tchowton/codeReview/230926_EJ_Setbp1_AlternativeSplicing/doc/"

write.table(
  data.frame(s.genes),
  file = paste0(wd, 'SGenes.csv'),
  append= TRUE,
  row.names = FALSE,
  col.names = FALSE,
  sep=','
)

write.table(
  data.frame(g2m.genes),
  file = paste0(wd, 'G2MGenes.csv'),
  append= TRUE,
  row.names = FALSE,
  col.names = FALSE,
  sep=','
)
