#!/bin/bash
#BSUB -P acc_DADisorders
#BSUB -q premium
#BSUB -n 8
#BSUB -R span[ptile=24]
#BSUB -R rusage[mem=60000]
#BSUB -W 120:00
#BSUB -o stdout.%J.%I
#BSUB -e stderr.%J.%I

module load R

# Path to the GTF file
GTF_FILE="/sc/arion/projects/MetaDope/genomes/rat/ensembl_rat_genome_jan182023/Rattus_norvegicus.mRatBN7.2.108.gtf"

# Run the R script with the GTF file path as an argument
Rscript /sc/arion/scratch/naskat01/AC_RNAseq/00_fastq/script/prepare.annotation.R $GTF_FILE
