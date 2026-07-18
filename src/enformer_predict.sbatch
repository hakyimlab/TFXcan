#!/bin/bash

# Author: Temi
# Description: run raw Enformer predictions for a custom set of regions (e.g. gene-module ACR experiments), via the TFXcan-snakemake repo's enformer_predict.py
# Usage: sbatch enformer_predict.sbatch parameters_yaml

#SBATCH --job-name=enformer_predict
#SBATCH --partition=beagle3
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=8
#SBATCH --time=36:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/enformer_predict.err
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/enformer_predict.out
#SBATCH --account=beagle3-exusers

date

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/beagle3/haky/users/shared_software/TFXcan-pipeline-tools/lib

nvidia-smi

# make sure that in your config, provider is set to "local" before running this script
python3 /beagle3/haky/users/temi/projects/TFXcan-snakemake/workflow/enformer/enformer_predict.py --parameters ${1}