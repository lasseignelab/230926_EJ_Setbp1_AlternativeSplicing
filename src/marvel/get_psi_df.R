# This code is a chopped-up version of MARVEL's PlotValues.PCA.PSI.10x to break
# down how they get PSI values for a single gene.
# Edited by Emma Jones

get_psi_df <- function (MarvelObject, cell.ids = NULL, coord.intron, 
          log2.transform = FALSE) 
{
  MarvelObject <- MarvelObject
  df.coord <- MarvelObject$pca
  sj.metadata <- MarvelObject$sj.metadata
  df.gene.count <- MarvelObject$gene.count.matrix
  df.sj.count <- MarvelObject$sj.count.matrix
  cell.ids <- cell.ids
  coord.intron <- coord.intron
  log2.transform <- log2.transform
  df.sj.count <- df.sj.count[coord.intron, ]
  df.sj.count <- data.frame(cell.id = names(df.sj.count), sj.count = as.numeric(df.sj.count), 
                            stringsAsFactors = FALSE)
  df.coord <- join(df.coord, df.sj.count, by = "cell.id", type = "left")
  gene_short_name <- sj.metadata[which(sj.metadata$coord.intron == 
                                         coord.intron), "gene_short_name.start"]
  df.gene.count <- df.gene.count[gene_short_name, ]
  df.gene.count <- data.frame(cell.id = names(df.gene.count), 
                              gene.count = as.numeric(df.gene.count), stringsAsFactors = FALSE)
  df.coord <- join(df.coord, df.gene.count, by = "cell.id", 
                   type = "left")
  df.coord$psi <- df.coord$sj.count/df.coord$gene.count * 100
  if (!is.null(cell.ids[1])) {
    df.coord <- df.coord[which(df.coord$cell.id %in% cell.ids), 
    ]
  }
  df.coord$psi[which(df.coord$psi > 100)] <- 100
  df.coord <- df.coord[order(df.coord$psi), ]
  data <- df.coord

  return(data)
}