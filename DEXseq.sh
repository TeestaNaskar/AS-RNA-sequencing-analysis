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

R CMD BATCH DEXseq.R
