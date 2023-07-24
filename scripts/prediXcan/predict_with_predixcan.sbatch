#!/bin/bash

# Author: Temi
# Description: Predict using PrediXcan
# Usage: sbatch predict_with_predixcan.sbatch
# Date: Mon July 24 2023
# Dependencies: 

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=predict_with_predixcan
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/predict_with_predixcan.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/predict_with_predixcan.err
#SBATCH --time=01:00:00	
#SBATCH --partition=caslake

module load openmpi
module load parallel

slurm_workdir=${SLURM_SUBMIT_DIR}
SLURM_O_WORKDIR=${slurm_workdir}/run
mkdir -p ${SLURM_O_WORKDIR}
echo Working directory is $SLURM_O_WORKDIR
cd $SLURM_O_WORKDIR

echo Jobid: $SLURM_JOBID
echo Running on host `hostname`

printf "Starting to run\n"

source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/imlabtools 

# variables
db_name="EUR_AR_Prostate"
output_folder="/beagle3/haky/users/temi/projects/TFXcan/experiments/compare_predictors/l-tfpred_predictions"
mkdir -p ${output_folder}
model_db=/beagle3/haky/users/temi/projects/l-TFPred/results/filtered_db/gtex_v7_EUR_AR_Prostate_filtered_signif.db
exec_file=/beagle3/haky/users/temi/software/MetaXcan/software/Predict.py
txt_genotypes=/project2/haky/Data/baca_cwas/formatted_geno/all_chrs.text_dosages.txt.gz
txt_samples=/project2/haky/Data/baca_cwas/formatted_geno/samples.text_dosages.txt


${exec_file} \
    --model_db_path ${model_db} \
    --text_genotypes ${txt_genotypes} \
    --text_sample_ids ${txt_samples} \
    --prediction_output "${output_folder}/${db_name}/baca_cwas_predict.txt" \
    --prediction_summary_output "${output_folder}/${db_name}/baca_cwas_summary.txt" \
    --verbosity 9 \
    --model_db_snp_key varID \
    --throw    