#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${USER}@uab.edu
#SBATCH --job-name=calc_psi
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=128G
#SBATCH --partition=short
#SBATCH --error=%x_%j_%a.err.txt
#SBATCH --output=%x_%j_%a.out.txt
#SBATCH --array=0-20

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################

# Load modules
ml Singularity

# Paths
wd="$USER_DATA/230926_EJ_Setbp1_AlternativeSplicing"

# Create variable for list of chromosome numbers
CHR_LIST="${wd}/bin/slurm/chromosome_prefixes.txt"
CHR_ARRAY=(`cat ${CHR_LIST}`)
INPUT=`echo ${CHR_ARRAY[$SLURM_ARRAY_TASK_ID]}`

# Code
cd ${wd}
singularity exec -B ${wd} ${wd}/bin/docker/setbp1_alternative_splicing_1.0.6.sif Rscript --vanilla ${wd}/src/marvel/04_calc_psi_all_sjs.R ${INPUT}
