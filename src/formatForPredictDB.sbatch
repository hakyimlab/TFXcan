#!/bin/bash

# Author: Temi
# Description: Predict using cwas in PrediXcan framework
# Usage: sbatch predict_baca_with_cwas.sbatch
# Date: Mon July 24 2023
# Dependencies: 

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=formatForPredictDB
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/Enpact/logs/formatForPredictDB.out
#SBATCH --error=/beagle3/haky/users/temi/projects/Enpact/logs/formatForPredictDB.err
#SBATCH --time=00:30:00	
#SBATCH --partition=caslake

date 
echo $SLURM_O_WORKDIR

SLURM_O_WORKDIR=${SLURM_SUBMIT_DIR}/run
mkdir -p ${SLURM_O_WORKDIR}
echo Working directory is $SLURM_O_WORKDIR
cd $SLURM_O_WORKDIR

echo Jobid: $SLURM_JOBID
echo Running on host `hostname`

printf "Starting to run\n"

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools

format_script=${1}
Rscript ${format_script} --enpact_scores_file ${2} --formatted_escores_file ${3} --formatted_annot_file ${4} --include_model ${5}

printf "Done\n"