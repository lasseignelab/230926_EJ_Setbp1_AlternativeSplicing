#!/bin/bash

cd ../bin/nextflow

module load Singularity
module load Anaconda3
conda activate nextflow

nextflow run nf-core/scrnaseq -r 2.4.1 \
   -profile cheaha \
   --input samplesheet.csv \
   --fasta /data/project/lasseigne_lab/GENOME_dir/GENCODE_mm39/release_M31/GRCm39.primary_assembly.genome.fa \
   --gtf /data/project/lasseigne_lab/GENOME_dir/GENCODE_mm39/release_M31/GTF/gencode.vM31.primary_assembly.annotation.gtf \
   --protocol 10XV3 \
   --aligner star \
   --outdir ../../data/nextflow/
