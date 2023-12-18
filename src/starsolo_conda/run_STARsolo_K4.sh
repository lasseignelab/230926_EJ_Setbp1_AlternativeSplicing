#!/bin/bash
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=efjones@uab.edu
#SBATCH --job-name=STARsolo
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

cd /data/project/lasseigne_lab/DATASET_dir/Setbp1S858R_snRNAseq/

module load Anaconda3
conda activate star

STAR --runThreadN 12 \
     --genomeDir /data/user/efjones/230926_EJ_Setbp1_AlternativeSplicing/data/star/genome/ \
     --soloType CB_UMI_Simple \
     --readFilesIn K4_S4_L001_R2_001.fastq.gz,K4_S4_L002_R2_001.fastq.gz,K4_S4_L003_R2_001.fastq.gz,K4_S4_L004_R2_001.fastq.gz K4_S4_L001_R1_001.fastq.gz,K4_S4_L002_R1_001.fastq.gz,K4_S4_L003_R1_001.fastq.gz,K4_S4_L004_R1_001.fastq.gz \
     --readFilesCommand zcat \
     --soloCBwhitelist /data/user/efjones/230926_EJ_Setbp1_AlternativeSplicing/bin/3M-february-2018.txt \
     --soloFeatures GeneFull_Ex50pAS SJ \
     --soloCellFilter EmptyDrops_CR \
     --sjdbGTFfile /data/project/lasseigne_lab/GENOME_dir/GENCODE_mm39/release_M31/GTF/gencode.vM31.primary_assembly.annotation.gtf \
     --soloUMIlen 12 \
     --clipAdapterType CellRanger4 \
     --outFilterScoreMin 30 \
     --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
     --soloUMIfiltering MultiGeneUMI_CR \
     --soloUMIdedup 1MM_CR \
     --outFileNamePrefix /data/user/efjones/230926_EJ_Setbp1_AlternativeSplicing/data/star/K4/
