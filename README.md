Cell-type-specific alternative splicing in the cerebral cortex of a
Schinzel-Giedion Syndrome patient variant mouse model
================
2024-06-25

## Authors

Emma F. Jones, Timothy C. Howton, Tabea M. Soelter, Anthony B. Crumley,
Brittany N. Lasseigne

## Project Overview

Schinzel-Giedion Syndrome (SGS) is an ultra-rare Mendelian disorder
caused by gain-of-function mutations in the SETBP1 gene. While previous
studies determined multiple roles for how SETBP1 and associated pathways
may cause disease manifestation, they have not assessed whether
cell-type-specific alternative splicing (AS) plays a role in SGS. We
used STARsolo to quantify gene and splice junction (SJ) expression for
51,465 nuclei previously generated from the cerebral cortex of atypical
Setbp1S858R SGS patient variant mice (n = 3) and wild-type control mice
(n = 3). After cell type annotation, we performed pseudobulk
differential gene expression and SJ usage (SJU) analyses across cell
types and conditions. We identified 34 genes with statistically
significant alterations in SJU. Oligodendrocytes had the most genes with
changes in SJU, followed by astrocytes, excitatory, and inhibitory
neurons. One gene, Son, a splicing cofactor known to cause the
neurodevelopmental disorder ZTTK Syndrome, had SJU changes in all six
non-vascular cell types we measured in Setbp1S858R compared to controls.
This is the first research to report AS changes in the cerebral cortex
of an SGS model and the first study to link SGS to perturbations in Son.

![Graphical Abstract](https://github.com/lasseignelab/230926_EJ_Setbp1_AlternativeSplicing/assets/85246122/a089c747-468b-4ad4-ae5c-3a696d72334f)
(A) Schematic overview of our processing and analysis pipeline. (B) We analyzed pseudobulk gene expression and calculated SJU for each cell type and condition. (C) We compared SJU values for each cell type using a permutation test (Methods) to identify cell-type-specific differences in AS between *Setbp1*<sup>S858R</sup> and wild-type mouse brain tissue. (D) Next, we visualized all annotated transcripts and splice junction locations for each significant SJU gene. (E) Finally, we compared the genes and pathways identified through functional enrichment analysis that overlap between cell types and predict their biological relevance.

## Scripts

#### STARsolo

The following bash scripts are for running STARsolo on the raw fastq
files. These scripts are not required to run if you download processed data
from zenodo.

    ## src/starsolo_conda
    ## ├── 01_build_STAR_genome.sh
    ## └── 02_run_STARsolo.sh

- 01_build_STAR_genome.sh - The purpose of this script is to build a
  genome to run STAR.
- 02_run_STARsolo.sh - The purpose of this script is to run STARsolo
  with custom filtering parameters to quantify gene and splice junction
  expression.

#### Seurat

The following scripts are for processing and filtering the gene
expression data.

    ## src/seurat
    ## ├── 01_import_filter_data.Rmd
    ## ├── 02_annotate_cell_types.Rmd
    ## ├── functions.R
    ## └── generate_cellcycle_lists.R

- 01_import_filter_data.Rmd - The purpose of this script is to import
  and filter data using seurat. It is dependent on running the scripts
  in STARsolo_conda or downloading the data to start. Please use docker
  image setbp1_alternative_splicing:1.0.4.
- 02_annotate_cell_types.Rmd - The purpose of this script is to annotate
  cell types using seurat. It is dependent on seurat script 01. Please
  use docker image setbp1_alternative_splicing:1.0.4.
- functions.R - The purpose of this script is to provide function
  scripts necessary for analysis.
- generate_cellcycle_lists.R

#### MARVEL

The following scripts are for formatting, importing, and analyzing
splice junction and expression data using the MARVEL package.

    ## src/marvel
    ## ├── 01_format_MARVEL_data.Rmd
    ## ├── 02_MARVEL_differential_analysis.Rmd
    ## ├── 03_analyze_de_genes.Rmd
    ## ├── 04_calc_sj_usage.Rmd
    ## ├── 05_cell_specific_sj_expr.Rmd
    ## └── functions.R

- 01_format_MARVEL_data.Rmd - The purpose of this script is to format
  expression data from seurat and SJ outputs from STARsolo for MARVEL.
  It is dependent on running the seurat scripts 01 and 02. Please use
  docker image setbp1_alternative_splicing:1.0.5.
- 02_MARVEL_differential_analysis.Rmd - The purpose of this script is to
  run differential analysis with MARVEL. It is dependent on running all
  seurat scripts, and marvel script 01. Please use docker image
  setbp1_alternative_splicing:1.0.5.
- 03_analyze_de_genes.Rmd - The purpose of this script is to analyze
  differential analysis results MARVEL. It is dependent on running all
  seurat scripts, and marvel scripts 01 and 02. Please use docker image
  setbp1_alternative_splicing:1.0.6, since you need ComplexHeatmap so
  make the UpSet plots.
- 04_calc_sj_usage.Rmd - The purpose of this script is to calculate
  splice junction usage for each cell type for each splice junction. The
  data is too sparse to calculate splice junction usage for single
  cells. It is dependent on running all seurat scripts, and marvel
  scripts 01 through 03. Please use docker image
  setbp1_alternative_splicing:1.0.6.
- 05_cell_specific_sj_expr.Rmd - The purpose of this script is to
  quantify the splice junction expression across cell types, both per
  cell and aggregated for each cell type. It creates a new supplemental
  table that is referred to in the manuscript. It is dependent on
  running all seurat scripts, and marvel scripts 01 through 04. Please
  use docker image setbp1_alternative_splicing:1.1.0.
- functions.R - The purpose of this script is to provide function
  scripts necessary for analysis.

#### DESeq2

The following scripts are for running a pseudobulk differential gene
expression analysis

    ## src/deseq2
    ## ├── 01_pseudobulk_analysis.Rmd
    ## └── functions.R

- 01_pseudobulk_analysis.Rmd - The purpose of this script is to run a
  pseudobulk gene expression analysis using the cell type data from
  seurat. It is dependent on seurat scripts 01 and 02. It should be run
  in docker 1.0.8.
- functions.R - The purpose of this script is to provide function
  scripts necessary for analysis.

#### Figures

The following script are for creating finalized figures for the
manuscript

    ## src/figures
    ## ├── figure_2.Rmd
    ## ├── figure_3-4.Rmd
    ## ├── figure_5.Rmd
    ## ├── functions.R
    ## ├── geom_split_violin.R
    ## ├── mean_expression_celltype.R
    ## └── supp_figure_2.Rmd

- figure_2.Rmd - The purpose of this script is to create a finalized
  version of figure 2. It is dependent on Marvel scripts 01 through 06.
  Run in docker 1.0.9.
- figure_3-4.Rmd - The purpose of this script is to make and complete
  figures 3 and 4. Figure 3 has an upset plot of all alternatively
  spliced genes, and Figure 4 has functional enrichment analysis of
  cell-type specific AS genes.
- figure_5.Rmd - The point of this figure is to convey the importance of
  the Son gene and also add an additional panel to show where the splice
  junctions lit on transcript structures.
- functions.R - The purpose of this script is to provide function
  scripts necessary for analysis.
- geom_split_violin.R - This is a split violin ggplot function borrowed
  from jan-glx on stackoverflow
- mean_expression_celltype.R - This is a short R script to get average
  expression for each cell type of a Seurat object
- supp_figure_2.Rmd - The purpose of this script is to create a
  finalized version of supporting information figure 2. It is dependent
  on Seurat scripts 01 and 02. Run in docker 1.0.6.

## Lasseigne Lab

[What is Happening in the Lasseigne Lab?](https://www.lasseigne.org/)

<img src="https://www.lasseigne.org/img/main/lablogo.png" width="75" height="75">

## Funding

This work was funded by the UAB Lasseigne Lab funds and the UAB Pilot
Center for Precision Animal Modeling (C-PAM) (1U54OD030167).

## Acknowledgements

We would like to acknowledge all current and former members of the
Lasseigne Lab for their thoughtful feedback, especially Elizabeth J.
Wilk, Amanda D. Clark, and Vishal H. Oza. The graphical abstract was
created using BioRender.

## License

This repository is licensed under the MIT License, see LICENSE
documentation within this repository for more details.
