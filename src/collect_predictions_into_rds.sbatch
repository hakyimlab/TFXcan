#!/bin/bash

# Author: Temi
# Date: Tuesday February 24 2026
# Description: cluster wrapper for collect_predictions_into_rds.R
# Usage: sbatch collect_predictions_into_rds.sbatch individuals gene_expression predictions_directory output_file

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=24G
#SBATCH --job-name=collect_predictions_into_rds
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/collect_predictions_into_rds.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/collect_predictions_into_rds.err
#SBATCH --time=06:00:00
#SBATCH --partition=bigmem

date

source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/renviron

individuals=${1}
gene_expression=${2}
predictions_directory=${3}
output_file=${4}

Rscript /beagle3/haky/users/temi/projects/TFXcan/src/collect_predictions_into_rds.R \
    --individuals ${individuals} \
    --gene_expression=${gene_expression} \
    --predictions_directory=${predictions_directory} \
    --output_file=${output_file}


status=$?
echo "Exit status of job is: $status"