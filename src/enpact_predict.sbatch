#!/bin/bash

# Author: Temi
# Description: apply an Enpact weights matrix to a feature matrix, via the external enpact_predict.py (TFXcan-snakemake)
# Usage: sbatch enpact_predict.sbatch matrix weights metadata output_basename

#SBATCH --job-name=enpact_predict
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/enpact_predict.err
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/enpact_predict.out
#SBATCH --account=pi-haky

date

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/beagle3/haky/users/shared_software/TFXcan-pipeline-tools/lib

# make sure that in your config, provider is set to "local" before running this script
# UNSURE: --split's exact behavior lives in the external enpact_predict.py, not in this repo
python3 /beagle3/haky/users/temi/projects/TFXcan-snakemake/workflow/process/enpact_predict.py --matrix ${1} --weights ${2} --metadata ${3} --output_basename ${4} --split

