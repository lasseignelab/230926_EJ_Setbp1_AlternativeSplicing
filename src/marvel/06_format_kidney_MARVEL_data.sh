#!/bin/bash

#SBATCH --ntasks=4          
#SBATCH --time=01:00:00                    
#SBATCH --mem=64G
#SBATCH --partition=express
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
                 -B /data/project/lasseigne_lab/ \
                 "${PROJECT_PATH}"/bin/docker/setbp1_alternative_splicing_1.1.0.sif \
                 Rscript --vanilla "${PROJECT_PATH}"/src/marvel/06_format_kidney_MARVEL_data.R \
                 "${PROJECT_PATH}"
