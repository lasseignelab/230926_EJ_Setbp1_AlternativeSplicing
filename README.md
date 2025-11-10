Cell-type-specific alternative splicing in the cerebral cortex and
kidney of a <i>Setbp1<sup>S858R</sup></i> Schinzel-Giedion Syndrome
patient variant mouse
================
2025-11-06

## Authors

**Tabea M. Soelter<sup>\#</sup>, Emma F. Jones<sup>\#</sup>, Timothy C.
Howton, Anthony B. Crumley, Elizabeth J. Wilk, Brittany N. Lasseigne**\*

\#Equal contribution  
\*Corresponding author

## Project Overview

Schinzel-Giedion Syndrome (SGS) is an ultra-rare Mendelian disorder
caused by gain-of-function mutations in the *SETBP1* gene. While
previous studies determined multiple roles for how *SETBP1* and
associated pathways may cause disease manifestation, they have not
assessed whether cell-type-specific alternative splicing (AS) plays a
role in SGS. We quantified gene and splice junction (SJ) expression from
snRNA-seq data we previously generated from the cerebral cortex and the
kidney of an atypical <i>Setbp1<sup>S858R</sup></i> SGS patient variant
(n = 3) and wild-type (n = 3) mice. We performed pseudobulk differential
gene expression and SJ usage (SJU) analyses across cell types and
conditions. We identified 33 and 62 genes with statistically significant
alterations in SJU in the brain and the kidney, respectively. Astrocytes
and T cells had the most genes with cell-type-specific changes in SJU (n
= 6 each) in the brain and kidney, respectively. We identified
significant SJU in a member of the heterogeneous nuclear
ribonucleoprotein family, *Hnrnpa2b1.* These findings were
cell-type-specific for inhibitory neurons in the cerebral cortex and
cell-type-agnostic in the kidney, suggesting tissue-specificity of AS in
<i>Setbp1<sup>S858R</sup></i> mice. To broaden the impact of our results
for the rare disease community, we developed a point-and-click web
application as a resource for users to explore single-cell resolution
changes in the presence of <i>Setbp1<sup>S858R</sup></i> at the gene and
splice junction level. Overall, we find that AS may be implicated in a
tissue- and cell-type-specific manner in the cerebral cortex and kidney
of <i>Setbp1<sup>S858R</sup></i> mice.

![](results/final_outputs/figure1.png)

## Scripts

    ## src
    ## ├── README
    ## ├── ambientRNA-removal
    ## │   ├── 01_ambientRNA_removal.R
    ## │   └── 01_ambientRNA_removal.sh
    ## ├── deseq2
    ## │   ├── 01_pseudobulk_analysis.Rmd
    ## │   ├── 02_pseudobulk_analysis_kidney.R
    ## │   ├── 02_pseudobulk_analysis_kidney.sh
    ## │   └── functions.R
    ## ├── figures
    ## │   ├── figure_2.Rmd
    ## │   ├── figure_3-4.Rmd
    ## │   ├── figure_5.Rmd
    ## │   ├── functions.R
    ## │   ├── geom_split_violin.R
    ## │   ├── kidney_marvel_figures.R
    ## │   ├── kidney_marvel_figures.sh
    ## │   ├── kidney_overview_figure.R
    ## │   ├── kidney_overview_figure.sh
    ## │   ├── mean_expression_celltype.R
    ## │   └── supp_figure_2.Rmd
    ## ├── functions_soelter.R
    ## ├── marvel
    ## │   ├── 01_format_MARVEL_data.Rmd
    ## │   ├── 02_MARVEL_differential_analysis.Rmd
    ## │   ├── 03_analyze_de_genes.Rmd
    ## │   ├── 04_calc_sj_usage.Rmd
    ## │   ├── 05_cell_specific_sj_expr.Rmd
    ## │   ├── 06_format_kidney_MARVEL_data.R
    ## │   ├── 06_format_kidney_MARVEL_data.sh
    ## │   ├── 07_kidney_MARVEL_differential_analysis.R
    ## │   ├── 07_kidney_MARVEL_differential_analysis.sh
    ## │   ├── 08_analyze_kidney_de_genes.R
    ## │   ├── 08_analyze_kidney_de_genes.sh
    ## │   ├── 09_calc_sj_usage_kidney.R
    ## │   ├── 09_calc_sj_usage_kidney.sh
    ## │   ├── 10_cell_specific_kidney_sj_expr.R
    ## │   ├── 10_cell_specific_kidney_sj_expr.sh
    ## │   ├── PlotSJPosition_modification.R
    ## │   └── functions.R
    ## ├── samtools
    ## │   └── sort_index_bam.sh
    ## ├── seurat
    ## │   ├── 01_import_filter_data.Rmd
    ## │   ├── 02_annotate_cell_types.Rmd
    ## │   ├── 03_kidney_preprocessing.R
    ## │   ├── 03_kidney_preprocessing.sh
    ## │   ├── functions.R
    ## │   └── generate_cellcycle_lists.R
    ## └── starsolo_conda
    ##     ├── 01_build_STAR_genome.sh
    ##     ├── 02_run_STARsolo.sh
    ##     └── README

## Lasseigne Lab

[What is Happening in the Lasseigne Lab?](https://www.lasseigne.org/)

<!-- markdownlint-disable-next-line MD013 -->

<img src="https://static.wixstatic.com/media/29dcd0_945b269e399141df9792e4ce78f832ff~mv2.png" width="75" height="75">

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
