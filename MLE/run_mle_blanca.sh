#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --qos=blanca-chbe-rdi
#SBATCH --job-name=MLE
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --output=mle.%j.out

source ~/.bashrc
conda activate ml_env
python3 MLE_poisson_CI_barcodes.py