#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --job-name=remVH10_2
#SBATCH --partition=amilan
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=rem10_2.%j.out

source ~/.bashrc
conda activate ml_env
python3 remove_muts.py /projects/brpe7306/fab-library-barcoding/maps_B15G10A5/VH10/merged.csv /scratch/alpine/moki5314/MLE_25Feb2023/fab-library-barcoding-main-3/MBK_test_25Feb2023/17Apr/MBK-4_S4_L001_R1_001.fastq.gz VH10_rep2_rewrite.csv