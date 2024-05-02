Examining cell-type specific patterns of alternative splicing of a
Schinzel-Giedion Syndrome mouse model using snRNA-Seq data
================
2024-05-02

## Authors

Emma F. Jones, Timothy C. Howton, Tabea M. Soelter, Brittany N.
Lasseigne

## Purpose

The purpose of this project is to analyze cell-type specific patterns of
alternative splicing in a rare neurodevelopmental disease model.

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
    ## ├── 02_compare_degs.Rmd
    ## └── functions.R

- 01_pseudobulk_analysis.Rmd
- 02_compare_degs.Rmd
- functions.R

## Lasseigne Lab

[What is Happening in the Lasseigne Lab?](https://www.lasseigne.org/)

<img src="https://www.lasseigne.org/img/main/lablogo.png" width="75" height="75">

## Funding

List project funding sources.

## Acknowledgements

We would like to thank all members of the Lasseigne Lab for their
thoughtful feedback.

## License

This repository is licensed under the MIT License, see LICENSE
documentation within this repository for more details.
