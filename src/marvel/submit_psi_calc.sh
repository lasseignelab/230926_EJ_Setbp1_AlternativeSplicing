#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${USER}@uab.edu
#SBATCH --job-name=calc_psi
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=250G
#SBATCH --partition=short
#SBATCH --error=%x_%j.err.txt
#SBATCH --output=%x_%j.out.txt

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################

# Load modules
ml Singularity

# Paths
wd="$USER_DATA/230926_EJ_Setbp1_AlternativeSplicing"

# Create variable for list of chromosome numbers
SAMPLE_LIST="${wd}/bin/slurm/chromosome_prefixes.txt"
SAMPLE_ARRAY=(`cat ${SAMPLE_LIST}`)
INPUT=`echo ${SAMPLE_ARRAY[$SLURM_ARRAY_TASK_ID]}`

# Code
cd ${wd}
singularity exec -B ${wd} ${wd}/bin/docker/setbp1_alternative_splicing_1.0.6.sif Rscript --vanilla ${wd}/src/marvel/calc_PSI_all_sjs.R