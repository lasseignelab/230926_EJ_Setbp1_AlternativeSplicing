#!/bin/bash
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=efjones@uab.edu
#SBATCH --job-name=STAR_genome
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

cd /data/user/efjones/230926_EJ_Setbp1_AlternativeSplicing/data/

module load Anaconda3
conda activate star

mdkir ./star/genome/

STAR --runThreadN 12 \
     --runMode genomeGenerate \
     --genomeDir ./star/genome/ \
     --genomeFastaFiles /data/project/lasseigne_lab/GENOME_dir/GENCODE_mm39/release_M31/GRCm39.primary_assembly.genome.fa \
     --sjdbGTFfile /data/project/lasseigne_lab/GENOME_dir/GENCODE_mm39/release_M31/GTF/gencode.vM31.primary_assembly.annotation.gtf