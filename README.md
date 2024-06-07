Cerebral cortex cell-type-specific alternative splicing in a Schinzel-Giedion Syndrome patient variant mouse model
================
2024-05-02

## Authors

Emma F. Jones, Timothy C. Howton, Tabea M. Soelter, Brittany N.
Lasseigne

## Purpose

Schinzel-Giedion Syndrome (SGS) is an ultra-rare Mendelian disorder caused by gain-of-function mutations in the SETBP1 gene. While previous studies determined multiple roles for how SETBP1 and associated pathways may cause disease manifestation, they have not assessed whether alternative splicing (AS) plays a role in SGS. We used STARsolo to quantify gene and splice junction (SJ) expression for 51,465 nuclei previously generated from Setbp1S858R+/- SGS patient variant (n=3) and wild-type control (n=3) mouse snRNA-Seq cerebral cortex tissue. We annotated cell types with Seurat. We then performed pseudobulk differential gene expression and SJ usage (SJU) analyses across patient variant mice and wild-type controls. We identified 34 cell-type-specific and shared genes with statistically significant alterations in SJU. Oligodendrocytes had the most genes with changes in SJU, followed by astrocytes, excitatory, and inhibitory neurons. One gene, Son, a splicing cofactor known to cause the neurodevelopmental disorder ZTTK Syndrome, had SJU changes in all six non-vascular cell types in the Setbp1S858R+/- variant mice compared to controls. This is the first research to report neural changes in AS in SGS and the first study to link SGS to perturbations in Son, which may help explain SGS’s severe neurological phenotype.

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
