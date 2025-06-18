#!/bin/bash

# Author: Temi
# Description: Predict using ltfpred in PrediXcan framework
# Usage: sbatch predict_baca_with_ltfpred.sbatch
# Date: Mon July 24 2023
# Dependencies: 

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=predict_baca_with_lenpact
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/Enpact/logs/predict_lenpact_scores.out
#SBATCH --error=/beagle3/haky/users/temi/projects/Enpact/logs/predict_lenpact_scores.err
#SBATCH --time=01:00:00	
#SBATCH --partition=caslake

date 
slurm_workdir=${SLURM_SUBMIT_DIR}
SLURM_O_WORKDIR=${slurm_workdir}/run
mkdir -p ${SLURM_O_WORKDIR}
echo "INFO - Working directory is $SLURM_O_WORKDIR"
cd $SLURM_O_WORKDIR

echo "INFO - Jobid: $SLURM_JOBID"
echo "INFO - Running on host `hostname`"

source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/imlabtools 

# variables
db_name=${1}
model_db=${2}
txt_genotypes=${3}
txt_samples=${4}
exec_file=${5}
output_folder=${6}
mkdir -p ${output_folder}

python ${exec_file} \
    --model_db_path ${model_db} \
    --text_genotypes ${txt_genotypes} \
    --text_sample_ids ${txt_samples} \
    --prediction_output "${output_folder}/${db_name}/baca_lEnpact_predict.txt" \
    --prediction_summary_output "${output_folder}/${db_name}/baca_lEnpact_summary.txt" \
    --verbosity 9 \
    --model_db_snp_key varID \
    --throw






# python ${exec_file} \
#     --model_db_path "${model_db}" \
#     --vcf_genotypes ${vcf_folder}/chr*.dose.vcf.gz \
#     --vcf_mode genotyped \
#     --on_the_fly_mapping METADATA "{}_{}_{}_{}_b38" \
#     --variant_mapping "/project2/haky/temi/projects/TFXcan/prediXcan/tutorial/data/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz" id rsid \
#     --prediction_output "${output_folder}/${db_name}/baca_ltfpred_predict.txt" \
#     --prediction_summary_output "${output_folder}/${db_name}/baca_ltfpred_predict.txt" \
#     --verbosity 9 \
#     --text_sample_ids ${txt_samples} \
#     --throw