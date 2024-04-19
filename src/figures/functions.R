# The purpose of this script is to provide any functions needed for generating
# figures
# Emma Jones
# Created March 4, 2024

# wilcoxon_test_celltype - function for running a wilcoxon rank sum test
wilcoxon_test_celltype <- function(cell_type) {
  result <- wilcox.test(mutant_metadata$num_sjs_genes[mutant_metadata$cell_type == cell_type], wildtype_metadata$num_sjs_genes[wildtype_metadata$cell_type == cell_type])
  
  return(result[["p.value"]])
}

# make_norm_gene_expr_df - function to make normalized gene expression dataframe
# for split violin plots
make_norm_gene_expr_df <- function(gene_name) {
  # pull out gene values
  gene_vals <- setbp1_marvel[["gene.norm.matrix"]][gene_name,]
  gene_vals <- as.data.frame(gene_vals)
  colnames(gene_vals) <- "norm_expr"
  # create new dataframe
  new_df <- cbind(base_df, gene_vals)
  # return new dataframe
  return(new_df)
}

# make_norm_sj_expr_df - function to make normalized sj expression dataframe
# for split violin plots
make_norm_sj_expr_df <- function(sj_coords) {
  # pull out gene values
  sj_vals <- normalized_sj_expression[sj_coords,]
  sj_vals <- as.data.frame(sj_vals)
  colnames(sj_vals) <- "norm_expr"
  # create new dataframe
  new_df <- cbind(base_df, sj_vals)
  # return new dataframe
  return(new_df)
}

# make_sj_split_violin - make sj split violins
make_sj_split_violin <- function(sj_expr_df, color_var) {sj_split_violin <- 
  ggplot(sj_expr_df,
  aes(
    x = cell_type, y = norm_expr, alpha = seq_folder, fill = color_var
  )
) +
  geom_split_violin(show.legend = TRUE) +
  scale_alpha_manual(
    name = "Condition",
    labels = c("Mutant", "Wildtype"),
    values = c(1, 0.4),
    guide = guide_legend(override.aes = list(fill = color_var))
  ) +
  scale_fill_manual(values = color_var) +
  theme_minimal(base_size = 15) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(
    fill = "none"
  ) +
  xlab("Cell Type") +
  ylab("Normalized SJ Expression")

return(sj_split_violin)}

## make_gene_sj_expr_usage_plots - Emma Jones
# the purpose of this function is for plotting gene expression for a given gene,
# splice junction expression, splice junction usage, and transcript structure.
make_gene_sj_expr_usage_plots <- function(gene_of_interest, save_path) {
  #print gene name for error handling
  print(gene_of_interest)
  
  # run own function to get normalized expression
  gene_expr_df <- make_norm_gene_expr_df(gene_of_interest)
  
  # get sjs
  gene_sjs <-
    setbp1_marvel[["sj.metadata"]]$coord.intron[
      setbp1_marvel[["sj.metadata"]]$gene_short_name.start == gene_of_interest
    ]
  
  # get mean expression
  sj_mean_expr <- sj_expr_means[gene_sjs]
  
  # get gene splice junctions
  gene_junctions <- strsplit(gene_sjs, ":")
  
  gene_junctions <- as_tibble(do.call(rbind, gene_junctions))
  
  colnames(gene_junctions) <- c("seqnames", "start", "end")
  
  gene_junctions$start <- as.numeric(gene_junctions$start)
  
  gene_junctions$end <- as.numeric(gene_junctions$end)
  
  gene_junctions$strand <- gtf$strand[gtf$gene_name == gene_of_interest][1]
  
  gene_junctions$mean_count <- sj_mean_expr
  
  gene_junctions <- gene_junctions %>%
    dplyr::mutate(transcript_name = paste0(gene_of_interest, "-201"))
  
  
  # make split violin plot of gene expression
  split_violin <- ggplot(
    gene_expr_df,
    aes(
      x = cell_type, y = norm_expr,
      fill = cell_type, alpha = seq_folder
    )
  ) +
    geom_split_violin(show.legend = TRUE) +
    scale_alpha_manual(
      name = "Condition",
      labels = c("S858R+/-", "Wildtype"),
      values = c(1, 0.4),
      guide = guide_legend(
        override.aes =
          list(fill = "gray40", color = "black")
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold", color = "black"),
      legend.text = element_text(face = "bold", size = 10),
      legend.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_fill_manual(
      labels = names(cell_type_colors),
      values = cell_type_colors,
      guide = "none"
    ) +
    xlab("Cell Type") +
    ylab("Normalized\nGene Expression")
  
  # get sj expr
  sj_mean_expr <- sj_expr_means[gene_sjs]
  
  # make splice junction expression heatmap
  delta_sje_heatmap_df <-
    delta_expr_for_heatmap[rownames(delta_expr_for_heatmap) %in% gene_sjs, ]
  
  print(gene_junctions$strand[1])
  
  ifelse(gene_junctions$strand[1] == "+", rownames(delta_sje_heatmap_df) <- c(paste0("SJ-", seq_len(nrow(gene_junctions)))),
         rownames(delta_sje_heatmap_df) <- c(paste0("SJ-", rev(seq_len(nrow(gene_junctions))))))

  # set colors
  delta_gene_expr_cols <-
    colorRamp2(
      c(
        min(delta_sje_heatmap_df), 0,
        max(delta_sje_heatmap_df)
      ),
      c(RColorBrewer::brewer.pal(name = "BrBG", n = 3))
    )
  
  # create heatmap
  delta_gene_heatmap <- Heatmap(delta_sje_heatmap_df,
                                name = "\u0394 Mean Norm\nSJ Expression",
                                column_order = sort(colnames(delta_sje_heatmap_df)),
                                col = delta_gene_expr_cols,
                                top_annotation = col_annot, cluster_rows = FALSE,
                                show_column_names = FALSE,
                                row_title = "SJ Expression", border = TRUE,
                                row_title_gp = gpar(fontface = "bold"),
                                row_names_gp = gpar(fontface = "bold")
  )
  
  # make splice junction usage heatmap
  delta_sju_heatmap_df <-
    delta_sju[rownames(delta_sju) %in% gene_sjs, ]
  
  if(gene_junctions$strand[1] == "+") {
    rownames(delta_sju_heatmap_df) <- paste0("SJ-", seq_len(nrow(gene_junctions)))  
  } 
  {rownames(delta_sju_heatmap_df) <- paste0("SJ-", rev(seq_len(nrow(gene_junctions))))
  }
  
  # set colors
  delta_sju_cols <-
    colorRamp2(
      c(
        min(delta_sju_heatmap_df), 0,
        max(delta_sju_heatmap_df)
      ),
      c(RColorBrewer::brewer.pal(name = "PRGn", n = 3))
    )
  
  # make heatmap
  delta_sj_heatmap <- Heatmap(delta_sju_heatmap_df,
                              name = "\u0394 SJ Usage",
                              column_order = sort(colnames(delta_sju_heatmap_df)),
                              col = delta_sju_cols,
                              cluster_rows = FALSE,
                              show_column_names = FALSE,
                              row_title = "SJ Usage", border = TRUE,
                              row_title_gp = gpar(fontface = "bold"),
                              row_names_gp = gpar(fontface = "bold")
  )
  
  # combine heatmaps
  both_delta_heatmaps <- delta_gene_heatmap %v% delta_sj_heatmap
  
  # make heatmaps compatible with cowplot
  delta_sj_expr_usage_heatmaps <- grid.grabExpr(draw(both_delta_heatmaps))
  
  # filter annotation for specific gene
  gene_annotation_from_gtf <- gtf %>%
    dplyr::filter(
      !is.na(gene_name),
      gene_name == gene_of_interest
    )
  
  gene_annotation_from_gtf <- gene_annotation_from_gtf %>%
    dplyr::select(
      seqnames,
      start,
      end,
      strand,
      type,
      gene_name,
      transcript_name,
      transcript_type
    )
  
  # extract exons
  gene_exons <- gene_annotation_from_gtf %>% dplyr::filter(type == "exon")
  
  gene_exons %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = transcript_name
    )) +
    geom_range(
      aes(fill = transcript_type)
    ) +
    geom_intron(
      data = to_intron(gene_exons, "transcript_name"),
      aes(strand = strand)
    )
  
  if(gene_junctions$strand[1] == "+") {
    gene_junctions$sj_label <- paste0("SJ-", seq_len(nrow(gene_junctions)))  
  } 
  { gene_junctions$sj_label <- paste0("SJ-", rev(seq_len(nrow(gene_junctions))))
  }
  
  # plot labeled transcripts with splice junctions
  gene_transcript_label <- gene_exons %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = transcript_name
    )) +
    geom_range(
      aes(fill = transcript_type)
    ) +
    geom_intron(
      data = to_intron(gene_exons, "transcript_name"),
      aes(strand = strand)
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold", color = "black"),
      legend.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      axis.title.y = element_blank()
    ) +
    labs(linewidth = "Normalized SJ\nExpression", fill = "Transcript Type") +
    scale_fill_manual(labels = c("NMD", "Protein Coding", "Protein Coding\nCDS Not Defined", "Retained Intron"),
                      values = c(viridis(6)[2:5])) +
    xlab("Genomic Location")
  
  # arrange vln plot
  gene_expr_vln <- plot_grid(gene_transcript_label, split_violin,
                             ncol = 1, rel_heights = c(0.66, 0.33),
                             labels = c("A", "B")
  )
  
  # arrange panels
  paneled_figure <- plot_grid(gene_expr_vln, delta_sj_expr_usage_heatmaps,
                              labels = c("", "C"), nrow = 1
  )
  
  # save
  png(paste0(save_path, "/", gene_of_interest, ".png"),
      width = 12, height = 7, units = "in", res = 300
  )
  print(paneled_figure)
  dev.off()
}


