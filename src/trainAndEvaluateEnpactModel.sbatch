#!/bin/bash

# Author: Temi
# Description: S-PrediXcan on EnpactScores
# Usage: sbatch s-predixcan.sbatch
# Date: Wed Apr 3 2024
# Dependencies: 

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=train_evaluate_enpact
#SBATCH --account=pi-haky
#SBATCH --output=/project/haky/users/temi/projects/Enpact/logs/train_evaluate_enpact.out
#SBATCH --error=/project/haky/users/temi/projects/Enpact/logs/train_evaluate_enpact.err
#SBATCH --time=02:00:00	
#SBATCH --partition=caslake

# module load openmpi
# module load parallel

date
slurm_workdir=${SLURM_SUBMIT_DIR}
SLURM_O_WORKDIR=${slurm_workdir}/run
mkdir -p ${SLURM_O_WORKDIR}
echo Working directory is $SLURM_O_WORKDIR
cd $SLURM_O_WORKDIR

echo Jobid: $SLURM_JOBID
echo Running on host `hostname`

printf "Starting to run\n"

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools

# variables
train_script=${1}
evaluate_script=${2}
train_data=${3}
eval_data=${4}
model_rds=${5}
basename_eval_output=${6}

# train
if [[ ! -f ${model_rds} ]]; then
    printf "Training model\n"
    Rscript ${train_script} --train_data_file ${train_data} --rds_file ${model_rds}
fi

# test
printf "Evaluating model\n"
if [[ ! -f ${model_rds} ]]; then
    echo "Model not found"
    exit 1
elif [[ -f ${model_rds} ]]; then
    echo "Model found"
    Rscript ${evaluate_script} --logistic_model ${model_rds} --train_data_file ${train_data} --test_data_file ${eval_data} --eval_output ${basename_eval_output}
fi
printf "Done\n"