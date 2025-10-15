### Plotting S858R mouse kidney data overview figure
# Author: Tabea M. Soelter (original code modified from Emma Jones)
# Date: October 15th, 2025

## Goal: Create S858R kidney data overview figure.

## Reproducibility:
# * GitHub: lasseignelab/230926_EJ_Setbp1_AlternativeSplicing
# * Docker: emmafjones/setbp1_alternative_splicing
#       * Version: 1.0.6 

## Data:
# Annotated Seurat object
# * Name: annotated_kidney_samples.rds
# * Location: data/seurat/

## Analysis Plan:
# * Load necessary packages 
# * Load data 
# * Plot individual panels
# * Compile overview figure
# * Save overview figure

## Set-up:
args <- R.utils::commandArgs(trailingOnly = TRUE)
wd <- args[1]
setwd(wd)

set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(patchwork)
  library(harmony)
  library(presto)
  library(cowplot)
})

## Analysis:
# Load data
annotated_kidney_samples <- read_rds(
  file = file.path(
    wd,
    "data",
    "seurat",
    "annotated_kidney_samples.rds"
  )
)

# Set color palettes
cell_type_colors <- c(
  `Proximal tubule cells` = "#F8766D",
  `Thick ascending limb (LOH)` = "#AA937E",
  `PCT` = "#A85A5A",
  `Endothelial cells` = "#CF9400",
  `PST` = "#948802",
  `CDPC` = "#F03F00",
  `DCT` = "#6D0404",
  `Mesenchymal cells` = "#FFC8C4",
  `Thin ascending limb (LOH)` = "#AA6320",
  `CDIC-B` = "#EEDF37",
  `Dendritic cells` = "#C70F0F",
  `Thin descending limb (LOH)` = "#FFE196",
  `CDIC-A` = "#4E3801",
  `Podocytes` = "#BD085E",
  `B cells` = "#FFBE1D",
  `T regulatory cells` = "#00B0F6",
  `T cells` = "#FF7600",
  `Connecting tubule cells` = "#E76BF3"
)

condition_colors <- c(
  `mutant` = "#666666",
  `wildtype` = "#C1C1C1"
)

other_colors <- c("#0346a3", "#076466", "#2a6607")

# Plot panel A
umap_1 <- DimPlot(annotated_kidney_samples,
                  cols = cell_type_colors,
                  label = TRUE,
                  label.box = TRUE,
                  label.color = "white",
                  label.size = 3.75,
                  repel = TRUE,
                  seed = 42
) +
  labs(color = "Cell Type") +
  theme(
    axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  NoLegend()

# Plot panel B
ordered_annotated_kidney_samples <- SetIdent(
  annotated_kidney_samples,
  value = factor(Idents(annotated_kidney_samples),
                 levels = sort(levels(Idents(annotated_kidney_samples)))
  )
)

dotplot <- DotPlot(ordered_annotated_kidney_samples,
                   cols = c("lightgray", other_colors[2]),
                   features = c(
                     "Aqp1", "Epha7", "Ptger3", "Ikzf2", "Cd247", "Slc22a7",
                     "Slc13a3", "Nphs1", "Slc5a2", "Gucy1a1",  "Flt1", "H2-Aa",
                     "Slc12a3", "Atp6v1b1", "Aqp2",  "Hmx2", "Aqp6", "Cd79a"
                   )
) +
  theme(
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ylab("Cell Type")

# Plot panel C
umap_2 <- DimPlot(annotated_kidney_samples,
                  group.by = "seq_folder",
                  shuffle = TRUE,
                  seed = 123
) +
  theme(
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  scale_color_manual(
    labels = c("S858R+/-", "Wildtype"),
    values = condition_colors
  ) +
  labs(color = "Condition", title = NULL) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  NoLegend()

# Plot panel D
df <- table(
  Idents(annotated_kidney_samples),
  annotated_kidney_samples$seq_folder
) %>% as.data.frame()

df$Var1 <- as.character(df$Var1)
colnames(df) <- c("Cell_Type", "Condition", "Freq")

df3 <- data.frame(
  "Cell_Type" = as.numeric(c(NA, NA)),
  "Condition" = c("S858R+/-", "Wildtype"),
  "Freq" = as.numeric(c(NA, NA))
)


barplot <- ggplot(data = df, aes(
  x = fct_rev(Cell_Type),
  y = Freq, fill = Condition
)) +
  theme_minimal(base_size = 15) +
  geom_col(position = "fill", width = 0.75, show.legend = FALSE) +
  scale_fill_manual(
    labels = c("S858R+/-", "Wildtype"),
    values = condition_colors
  ) +
  theme(
    legend.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text()
  ) +
  coord_flip() +
  xlab("Cell Type") +
  ylab("Proportion") +
  geom_point(data = df3, aes(
    x = Cell_Type,
    y = Freq, color = Condition
  ), size = 4) +
  scale_color_manual(
    labels = c("S858R+/-", "Wildtype"),
    values = c("#666666", "#C1C1C1")
  ) +
  theme(legend.key = element_blank())

# Compile overview figure
figure <- plot_grid(umap_1, dotplot, umap_2, barplot,
                    labels = c("A", "B", "C", "D")
)

png(file.path(wd, "results", "kidney_figures", "kidney_overview_figure.png"),
    width = 14, height = 10, units = "in", res = 300
)
figure
dev.off()

# Get session information
sessionInfo()
