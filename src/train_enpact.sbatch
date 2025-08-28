#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=12G    
#SBATCH --job-name=train_enpact
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/train_enpact.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/train_enpact.err
#SBATCH --time=01:00:00
#SBATCH --partition=caslake

date

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools

Rscript /beagle3/haky/users/temi/projects/TFXcan/src/trainEnpactModel.R \
    --train_data_file ${1} \
    --model_rds_file ${2} \
    --nfolds ${3} \
    --weights_file ${4} \
    --metadata ${5}

status=$?
echo "Exit status of job is: $status"