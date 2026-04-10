#!/bin/bash
#BSUB -q short
#BSUB -J knockoff_sim[1-5400]
#BSUB -o log/knockoff-%I.out
#BSUB -e log/knockoff-%I.err
#BSUB -W 480
#BSUB -M 10000

# Create directories
mkdir -p log
mkdir -p ./results/gaussian
mkdir -p ./results/logistic

# Debug LSF environment
echo "LSB_JOBINDEX=${LSB_JOBINDEX}"
echo "Running on host: $(hostname)"

# Run R script
Rscript simulations_derandomized_FDR_lsf.R ${LSB_JOBINDEX}
