#!/bin/bash

#SBATCH --ntasks=4          
#SBATCH --time=12:00:00                    
#SBATCH --mem=30G
#SBATCH --partition=short
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################
module load Singularity/3.5.2-GCC-5.4.0-2.26

PROJECT_PATH="$(cd "$SLURM_SUBMIT_DIR/../.." && pwd)"

singularity exec --cleanenv \
                 --containall \
                 -B "${PROJECT_PATH}" \
                 "${PROJECT_PATH}"/bin/docker/setbp1_alternative_splicing_1.1.0.sif \
                 Rscript --vanilla "${PROJECT_PATH}"/src/seurat/03_kidney_preprocessing.R \
                 "${PROJECT_PATH}"