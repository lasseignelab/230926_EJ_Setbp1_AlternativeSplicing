#!/bin/bash

#SBATCH --ntasks=1          
#SBATCH --time=00:30:00                    
#SBATCH --mem=16G
#SBATCH --partition=express
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################
PROJECT_PATH="$(cd "$SLURM_SUBMIT_DIR/../.." && pwd)"

module load Singularity/3.5.2-GCC-5.4.0-2.26

singularity exec --cleanenv \
                 --containall \
                 -B "${PROJECT_PATH}" \
                 "${PROJECT_PATH}"/bin/docker/setbp1_alternative_splicing_1.0.6.sif \
                 Rscript --vanilla "${PROJECT_PATH}"/src/figures/kidney_overview_figure.R \
                 "${PROJECT_PATH}"