Examining cell-type specific patterns of alternative splicing of a
Schinzel-Giedion Syndrome mouse model using snRNA-Seq data
================
2024-01-08

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
    ## ├── run_STARsolo_J1.sh
    ## ├── run_STARsolo_J13.sh
    ## ├── run_STARsolo_J15.sh
    ## ├── run_STARsolo_J2.sh
    ## ├── run_STARsolo_J3.sh
    ## ├── run_STARsolo_J4.sh
    ## ├── run_STARsolo_K1.sh
    ## ├── run_STARsolo_K2.sh
    ## ├── run_STARsolo_K3.sh
    ## ├── run_STARsolo_K4.sh
    ## ├── run_STARsolo_K5.sh
    ## └── run_STARsolo_K6.sh

- 01_build_STAR_genome.sh
- run_STARsolo\_\*

#### Seurat

The following scripts are for processing and filtering the gene
expression data.

    ## src/seurat
    ## ├── 01_import_filter_data.Rmd
    ## ├── 02_annotate_cell_types.Rmd
    ## └── functions.R

- 01_import_filter_data.Rmd
- 02_annotate_cell_types.Rmd
- functions.R

#### MARVEL

The following scripts are for formatting and importing splice junction
and expression data for MARVEL analyses.

    ## src/marvel
    ## ├── 01_format_MARVEL_data.Rmd
    ## ├── MARVEL_tutorial.Rmd
    ## └── functions.R

- 01_format_MARVEL_data.Rmd
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
