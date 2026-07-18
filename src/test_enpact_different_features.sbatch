#!/bin/bash

# Author: Temi
# Date: Monday January 26 2026
# Description: cluster wrapper for test_enpact_different_features.R
# Usage: sbatch test_enpact_different_features.sbatch test_data_file features_to_remove metadata weights_file predictions_file

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=test_enpact_features
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/test_enpact_features.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/test_enpact_features.err
#SBATCH --time=00:10:00
#SBATCH --partition=caslake

date

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools

Rscript /beagle3/haky/users/temi/projects/TFXcan/src/test_enpact_different_features.R \
    --test_data_file ${1} \
    --features_to_remove ${2} \
    --metadata ${3} \
    --weights_file ${4} \
    --predictions_file ${5}

status=$?
echo "Exit status of job is: $status"