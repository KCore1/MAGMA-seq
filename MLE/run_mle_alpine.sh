#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --job-name=MLE
#SBATCH --partition=amilan
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --output=mle.%j.out

source ~/.bashrc
conda activate ml_env
python3 MLE_poisson_CI_barcodes.py