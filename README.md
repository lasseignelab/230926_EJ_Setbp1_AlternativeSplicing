Cell-type-specific alternative splicing in the brain and kidney of a
<i>Setbp1</i><sup>S858R</sup> Schinzel-Giedion syndrome mouse
================
2026-02-27

Publication:  
[![DOI:10.1242/dmm.052402](https://img.shields.io/badge/DOI-10.1242/dmm/052402-blue)](https://doi.org/10.1242/dmm.052402)  
Featured First Author Interview:  
[![DOI:10.1242/dmm.052868](https://img.shields.io/badge/DOI-10.1242/dmm/052868-blue)](https://doi.org/10.1242/dmm.052868)

## Authors

**Tabea M. Soelter<sup>\#</sup>, Emma F. Jones<sup>\#</sup>, Timothy C.
Howton, Anthony B. Crumley, Elizabeth J. Wilk, Brittany N. Lasseigne**\*

\#Equal contribution  
\*Corresponding author

## Project Overview

Schinzel-Giedion syndrome (SGS) is an ultra-rare Mendelian disorder
caused by gain-of-function variants in the *SETBP1* gene. Although
previous studies determined multiple roles for *SETBP1* and its
associated pathways in disease manifestation, they did not assess
whether cell-type-specific alternative splicing (AS) plays a role in
SGS. We quantified gene and splice junction expression from
single-nuclei RNA-sequencing data from the cerebral cortex and the
kidney of atypical <i>Setbp1</i><sup>S858R</sup> SGS patient variant and
wild-type mice. We identified 33 and 62 genes with statistically
significant alterations in splice junction usage in the brain and the
kidney, respectively. We identified significant splice junction usage in
a member of the heterogeneous nuclear ribonucleoprotein family,
*Hnrnpa2b1.* These findings were cell-type-specific in the cerebral
cortex and cell-type-agnostic in the kidney, suggesting
tissue-specificity of AS in <i>Setbp1</i><sup>S858R</sup> mice. To
broaden the impact of our results for the rare disease community, we
developed a point-and-click web application that enables users to
explore single-cell-resolution changes at the gene and splice junction
levels. Overall, our findings implicate AS in a tissue- and
cell-type-specific manner in the cerebral cortex and kidney of
<i>Setbp1</i><sup>S858R</sup> mice.

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
Lasseigne Lab for their thoughtful feedback, especially Amanda D. Clark
and Vishal H. Oza. The graphical abstract was created using BioRender.

## License

This repository is licensed under the MIT License, see LICENSE
documentation within this repository for more details.
