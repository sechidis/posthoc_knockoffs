#!/bin/bash
#BSUB -q short
#BSUB -J knockoff_sim[1-2000]    
#BSUB -o log/knockoff-%I.out
#BSUB -e log/knockoff-%I.err
#BSUB -W 480
#BSUB -M 10000

mkdir -p log ./results/gaussian ./results/logistic

echo "LSB_JOBINDEX=${LSB_JOBINDEX}"
echo "Running on host: $(hostname)"

Rscript simulations_PH_vs_original_lsf.R "${LSB_JOBINDEX}"
