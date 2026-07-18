#!/bin/bash

# Author: Temi
# Date: Thursday October 16 2025
# Description: cluster wrapper for collect_enpact_weights.R
# Usage: sbatch collect_enpact_weights.sbatch models_metadata weights_file_basename

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=collect_enpact_weights
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/collect_enpact_weights.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/collect_enpact_weights.err
#SBATCH --time=01:00:00
#SBATCH --partition=caslake

date

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools

Rscript /beagle3/haky/users/temi/projects/TFXcan/src/collect_enpact_weights.R \
    --models_metadata ${1} \
    --weights_file_basename ${2}

status=$?
printf "INFO - All jobs are done. Exit status is $status.\n"

# /beagle3/haky/users/temi/software/gcta-1.94.4-linux-kernel-3-x86_64/gcta64 --reml --grm /project2/haky/Data/baca_cwas/vcfs/hg38/grm_files/chr22.dose --pheno /beagle3/haky/users/temi/projects/Enpact/experiments/heritability/phenotypes/chr10_100185843_100186293.peaks.txt --out hey


