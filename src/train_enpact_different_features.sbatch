#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=12G    
#SBATCH --job-name=train_enpact_features
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/train_enpact_features.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/train_enpact_features.err
#SBATCH --time=02:00:00
#SBATCH --partition=bigmem

date

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools

training_data=${1}
output_basename=${2}
metadata=${3}

# Rscript /beagle3/haky/users/temi/projects/TFXcan/src/train_enpact_different_features.R \
#     --train_data_file ${1} \
#     --features_to_remove ${2} \
#     --nfolds ${3} \
#     --metadata ${4} \
#     --weights_file ${5} \
#     --model_rds_file ${6}

Rscript /beagle3/haky/users/temi/projects/TFXcan/src/train_enpact.R \
    --train_data_file ${training_data} \
    --metadata ${metadata} \
    --weights_file ${output_basename}.weights.tsv.gz \
    --model_rds_file ${output_basename}.glmnet.rds

status=$?
echo "Exit status of job is: $status"