#!/bin/bash

# Author: Temi
# Date: Thursday October 16 2025
# Description: apply an Enpact weights matrix to a feature matrix, via the external enpact_predict.py (TFXcan-snakemake)
# Usage: sbatch predict_enpact_scores.sbatch matrix weights metadata output_basename

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=predict_enpact_scores
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/predict_enpact_scores.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/predict_enpact_scores.err
#SBATCH --time=01:00:00
#SBATCH --partition=caslake

date

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools

# UNSURE: --split's exact behavior lives in the external enpact_predict.py, not in this repo
python3 /beagle3/haky/users/temi/projects/TFXcan-snakemake/workflow/process/enpact_predict.py --matrix ${1} --weights ${2} --metadata ${3} --split --output_basename ${4}

status=$?
echo "Exit status of job is: $status"
