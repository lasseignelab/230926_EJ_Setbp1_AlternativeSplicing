#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${USER}@uab.edu
#SBATCH --job-name=calc_complexity
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=128G
#SBATCH --partition=short
#SBATCH --error=logs/%x_%j_%a.err.txt
#SBATCH --output=logs/%x_%j_%a.out.txt

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################

# Load modules
ml Singularity

# Paths
wd="$USER_DATA/230926_EJ_Setbp1_AlternativeSplicing"

# Code
cd ${wd}
singularity exec -B ${wd} ${wd}/bin/docker/setbp1_alternative_splicing_1.0.6.sif Rscript --vanilla ${wd}/src/marvel/07_calc_complexity_all_cells.R

