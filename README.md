Cell-type-specific alternative splicing in the cerebral cortex of a Schinzel-Giedion Syndrome patient variant mouse model
================
2024-06-25

## Authors

Emma F. Jones, Timothy C. Howton, Tabea M. Soelter, Anthony B. Crumley, Brittany N.
Lasseigne

## Project Overview

Schinzel-Giedion Syndrome (SGS) is an ultra-rare Mendelian disorder caused by gain-of-function mutations in the *SETBP1* gene. While previous studies determined multiple roles for how *SETBP1* and associated pathways may cause disease manifestation, they have not assessed whether cell-type-specific alternative splicing (AS) plays a role in SGS. We used STARsolo to quantify gene and splice junction (SJ) expression for 51,465 nuclei previously generated from the cerebral cortex of atypical *Setbp1*<sup>S858R</sup> SGS patient variant mice (n = 3) and wild-type control mice (n = 3). After cell type annotation, we performed pseudobulk differential gene expression and SJ usage (SJU) analyses across cell types and conditions. We identified 34 genes with statistically significant alterations in SJU. Oligodendrocytes had the most genes with changes in SJU, followed by astrocytes, excitatory, and inhibitory neurons. One gene, *Son*, a splicing cofactor known to cause the neurodevelopmental disorder ZTTK Syndrome, had SJU changes in all six non-vascular cell types we measured in *Setbp1*<sup>S858R</sup> compared to controls. This is the first research to report AS changes in the cerebral cortex of an SGS model and the first study to link SGS to perturbations in *Son*.

![Project Overview Jun 24 (2)](https://github.com/lasseignelab/230926_EJ_Setbp1_AlternativeSplicing/assets/85246122/a089c747-468b-4ad4-ae5c-3a696d72334f)
(A) Schematic overview of our processing and analysis pipeline. (B) We analyzed pseudobulk gene expression and calculated SJU for each cell type and condition. (C) We compared SJU values for each cell type using a permutation test (Methods) to identify cell-type-specific differences in AS between Setbp1S858R and wild-type mouse brain tissue. (D) Next, we visualized all annotated transcripts and splice junction locations for each significant SJU gene. (E) Finally, we compared the genes and pathways identified through functional enrichment analysis that overlap between cell types and predict their biological relevance.

## Scripts

#### STARsolo

The following bash scripts are for running STARsolo on the raw fastq
files. The J samples are all brain samples, while K samples are kidney.
These scripts are not required to run if you download processed data
from zenodo.

    ## src/starsolo_conda
    ## ├── 01_build_STAR_genome.sh
    ## └── 02_run_STARsolo.sh

- 01_build_STAR_genome.sh
- 02_run_STARsolo.sh

#### Seurat

The following scripts are for processing and filtering the gene
expression data.

    ## src/seurat
    ## ├── 01_import_filter_data.Rmd
    ## ├── 02_annotate_cell_types.Rmd
    ## ├── functions.R
    ## └── generate_cellcycle_lists.R

- 01_import_filter_data.Rmd
- 02_annotate_cell_types.Rmd
- functions.R

#### MARVEL

The following scripts are for formatting and importing splice junction
and expression data for MARVEL analyses.

    ## src/marvel
    ## ├── 01_format_MARVEL_data.Rmd
    ## ├── 02_MARVEL_differential_analysis.Rmd
    ## ├── 03_analyze_de_genes.Rmd
    ## ├── 04_calc_sj_usage.Rmd
    ## ├── PlotSJPosition_modification.R
    ## ├── cv_sj_celltypes.Rmd
    ## ├── explore_sj_usage.Rmd
    ## ├── functions.R
    ## └── setbp1_targets_properties.Rmd

- 01_format_MARVEL_data.Rmd
- 02_MARVEL_differential_analysis.Rmd
- 03_analyze_de_genes.Rmd
- 04_calc_sj_usage.Rmd
- setbp1_targets_properties.Rmd
- functions.R

#### DESeq2

The following scripts are for running a pseudobulk differential gene
expression analysis

    ## src/deseq2
    ## ├── 01_pseudobulk_analysis.Rmd
    ## └── functions.R

- 01_pseudobulk_analysis.Rmd
- functions.R

## Lasseigne Lab

[What is Happening in the Lasseigne Lab?](https://www.lasseigne.org/)

<img src="https://www.lasseigne.org/img/main/lablogo.png" width="75" height="75">

## Funding

This work was funded in part by the UAB Lasseigne Lab funds and the UAB Pilot Center for Precision Animal Modeling (C-PAM) (1U54OD030167).

## Acknowledgements

We would like to acknowledge all current and prior members of the Lasseigne Lab for their thoughtful feedback, especially Elizabeth J. Wilk, Amanda D. Clark, and Vishal H. Oza. All figures and cartoons were assembled using BioRender.

## License

This repository is licensed under the MIT License, see LICENSE
documentation within this repository for more details.
