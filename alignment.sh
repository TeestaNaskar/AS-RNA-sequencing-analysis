---
title: "Identification of Alternative Splicing and polyadenylation in bulk RNA-seq data"
author: "Teesta"
date: "07/22/2024"
output: html_notebook
---

# Introduction
NOTE: The context of this R notebook includes all the R code requied for the analysis of Differential splicing analysis of bulk RNA-Seq data. Refer to the PART 1- Bulk RNA-Seq analysis of the Protocol section.

# 1. Data downloading and pre-processing 
NOTE: Download the raw fastq.gz file or transfer from the source.

#BSUB -J alternative_splicing
#BSUB -P acc_DADisorders
#BSUB -q premium
#BSUB -n 12
#BSUB -R rusage[mem=8000]
#BSUB -R span[hosts=1]
#BSUB -W 48:00
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

RAW_DATA="/sc/arion/scratch/naskat01/AC_RNAseq/00_fastq"
ODIR="/sc/arion/scratch/naskat01/AC_RNAseq/00_fastq"  # Replace with your actual output directory path

### Ensure the output directory exists and is writable
mkdir -p $ODIR
if [ ! -w $ODIR ]; then
  echo "ERROR: Cannot write to output directory $ODIR"
  exit 1
fi

### LOAD MODULES HERE
ml star
ml samtools

# Align reads to the genome 
for fq1 in $RAW_DATA/*R1_001.fastq.gz; do
  fq2=$(echo $fq1 | sed 's/R1_001.fastq.gz/R2_001.fastq.gz/g')
  OUTPUT=$(basename ${fq1} | sed 's/R1_001.fastq.gz//g')
  STAR --genomeDir /sc/arion/projects/MetaDope/genomes/rat/starIndex \
    --runThreadN 12 \
    --readFilesCommand zcat \
    --readFilesIn ${fq1} ${fq2} \
    --outFileNamePrefix ${ODIR}/${OUTPUT} \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard
done
