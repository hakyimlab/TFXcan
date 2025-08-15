#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=12G    
#SBATCH --job-name=train_enpact_features
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/Enpact/logs/train_enpact_features.out
#SBATCH --error=/beagle3/haky/users/temi/projects/Enpact/logs/train_enpact_features.err
#SBATCH --time=01:00:00
#SBATCH --partition=caslake

date

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools

Rscript /beagle3/haky/users/temi/projects/TFXcan/src/train_enpact_different_features.R \
    --train_data_file ${1} \
    --test_data_file ${2} \
    --features_to_remove ${3} \
    --nfolds ${4} \
    --metadata ${5} \
    --weights_file ${6} \
    --predictions_file ${7} \
    --model_rds_file ${8}

status=$?
echo "Exit status of job is: $status"