#!/bin/bash
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${USER}@uab.edu
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
# Run this with sbatch run_STARsolo.sh <output dir> <fastq prefix> for example sbatch run_STARsolo.sh J1 J-1_S7


# Define variables
GENOME="${USER_DATA}/codeReview"
OUTDIR="${USER_SCRATCH}"

cd /data/project/lasseigne_lab/DATASET_dir/Setbp1S858R_snRNAseq/

module load Anaconda3
conda activate star

STAR --runThreadN 12 \
     --genomeDir ${GENOME}/230926_EJ_Setbp1_AlternativeSplicing/data/star/genome/ \
     --soloType CB_UMI_Simple \
     --readFilesIn ${2}_L001_R2_001.fastq.gz,${2}_L002_R2_001.fastq.gz,${2}_L003_R2_001.fastq.gz,${2}_L004_R2_001.fastq.gz ${2}_L001_R1_001.fastq.gz,${2}_L002_R1_001.fastq.gz,${2}_L003_R1_001.fastq.gz,${2}_L004_R1_001.fastq.gz \
     --readFilesCommand zcat \
     --soloCBwhitelist /data/project/lasseigne_lab/EmmaJones/3M-february-2018.txt \
     --soloFeatures GeneFull_Ex50pAS SJ \
     --soloCellFilter EmptyDrops_CR \
     --sjdbGTFfile /data/project/lasseigne_lab/GENOME_dir/GENCODE_mm39/release_M31/GTF/gencode.vM31.primary_assembly.annotation.gtf \
     --soloUMIlen 12 \
     --clipAdapterType CellRanger4 \
     --outFilterScoreMin 30 \
     --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
     --soloUMIfiltering MultiGeneUMI_CR \
     --soloUMIdedup 1MM_CR \
     --outFileNamePrefix ${OUTDIR}/230926_EJ_Setbp1_AlternativeSplicing/data/star/${1}/

