#!/bin/bash

# Author: Temi
# Description: S-PrediXcan on EnpactScores
# Usage: sbatch s-predixcan.sbatch
# Date: Wed Apr 3 2024
# Dependencies: 

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=12G
#SBATCH --job-name=tfxcan_pcr
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/Enpact/logs/tfxcan_pcr.out
#SBATCH --error=/beagle3/haky/users/temi/projects/Enpact/logs/tfxcan_pcr.err
#SBATCH --time=24:00:00	
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
conda activate /beagle3/haky/users/temi/software/conda_envs/imlabtools 

# variables
db_name=${1}
model_db=${2}
covariances=${3}
gwas_folder=${4}
gwas_file_pattern=${5}
exec_file=${6}
output_file=${7}

# note that the ...
python ${exec_file} \
    --model_db_path ${model_db} \
    --model_db_snp_key varID \
    --covariance ${covariances} \
    --gwas_folder ${gwas_folder} \
    --gwas_file_pattern ${gwas_file_pattern} \
    --chromosome_column chrom \
    --snp_column variant_id \
    --effect_allele_column A1 \
    --non_effect_allele_column A2 \
    --position_column pos \
    --freq_column maf \
    --beta_column beta \
    --se_column se \
    --pvalue_column p_value \
    --zscore_column zscore \
    --keep_non_rsid \
    --output_file ${output_file} \
    --verbosity 9 \
    --throw

    #     --gwas_folder ${base_dir}/data/GWAS \
    # --gwas_file_pattern ".*_Schumacher.gwas_ss.txt.gz" \


    # awk '{gsub(/:/, "_"); print}' Covariances.varID.txt > Covariances.varID.txt

    # sed -i -e 's/:/_/g' Covariances.varID.txt

    #sed -i -e 's/:/_/g' Covariances.txt > Covariances.varID.txt
    # sed -e 's/:/_/g' predict_db_EUR_AR_Prostate_logistic.txt > Covariances.varID.txt