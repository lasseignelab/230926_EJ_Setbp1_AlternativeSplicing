#!/bin/bash

#SBATCH --ntasks=4          
#SBATCH --time=02:00:00                    
#SBATCH --mem=30G
#SBATCH --partition=express
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################
module load Singularity/3.5.2-GCC-5.4.0-2.26

PROJECT_PATH="$(cd "$SLURM_SUBMIT_DIR/../.." && pwd)"

mkdir -p "$PROJECT_PATH/results/ambientRNA_removal/"

singularity exec --cleanenv \
                 --containall \
                 -B "${PROJECT_PATH}" \
                 "${PROJECT_PATH}"/bin/docker/sn-ml-drug-repurposing_1.0.0.sif \
                 Rscript --vanilla "${PROJECT_PATH}"/src/ambientRNA-removal/01_ambientRNA_removal.R \
                 "${PROJECT_PATH}"
