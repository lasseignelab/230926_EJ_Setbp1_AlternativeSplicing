#!/bin/bash
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=efjones@uab.edu
#SBATCH --job-name=samtools_sort_index
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=80G
#SBATCH --nodes=1
#SABTCH --cpus-per-task=12
#SBATCH --time=12:00:00
#SBATCH --partition=short
#SBATCH --error=%x_%j.err.txt
#SBATCH --output=%x_%j.out.txt

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################
# Run this script with sbatch sort_index_bam.sh <sample id>
# For example: sbatch sort_index_bam.sh J1

cd /data/user/efjones/230926_EJ_Setbp1_AlternativeSplicing/data/igv/

module load SAMtools

samtools sort ${1}.bam -o sorted/${1}.sorted.bam

cd sorted/

samtools index ${1}.sorted.bam