#!/bin/bash

#SBATCH --ntasks=4          
#SBATCH --time=2:00:00                    
#SBATCH --mem=64G
#SBATCH --partition=express
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################
PROJECT_PATH="$(cd "$SLURM_SUBMIT_DIR/../.." && pwd)"

TMPDIR="$PROJECT_PATH/tmp_dir"
mkdir -p "${TMPDIR}"
export SINGULARITYENV_TMPDIR="${TMPDIR}"

module load Singularity/3.5.2-GCC-5.4.0-2.26

singularity exec --cleanenv \
                 --containall \
                 -B "${TMPDIR}:$TMPDIR" \
                 -B "${PROJECT_PATH}" \
                 "${PROJECT_PATH}"/bin/docker/setbp1_alternative_splicing_2.0.0.sif \
                 Rscript --vanilla "${PROJECT_PATH}"/src/seurat/03_kidney_preprocessing.R \
                 "${PROJECT_PATH}"
                 
rm -rf "${TMPDIR}"