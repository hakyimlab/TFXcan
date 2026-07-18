#!/bin/bash

# Author: Temi
# Date: Monday January 26 2026
# Description: run S-PrediXcan ("TFXcan") across many TF-tissue model DBs (from models_list) against one phenotype's GWAS
# Usage: sbatch run_summary_tfxcan.sbatch phenotype models_list output_directory gwas_folder gwas_pattern exec_file conda_environment

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=summary_TFXcan
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/summary_TFXcan.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/summary_TFXcan.err
#SBATCH --time=02:00:00
#SBATCH --partition=caslake

date

slurm_workdir=${SLURM_SUBMIT_DIR}
# SLURM_O_WORKDIR=${slurm_workdir}/run
# mkdir -p ${SLURM_O_WORKDIR}
# cd $SLURM_O_WORKDIR
echo "INFO - Working directory is $SLURM_SUBMIT_DIR"
echo "INFO - Jobid is $SLURM_JOBID"
echo "INFO - Running on host `hostname`"

source ~/.bashrc
module load parallel
job_log='/beagle3/haky/users/temi/projects/TFXcan/logs/summary_TFXcan.log'

# variables
phenotype=${1}
models_list=${2}
output_directory=${3}
gwas_folder=${4}
gwas_pattern=${5}
exec_file=${6}
conda_environment=${7}

conda activate ${conda_environment} 

printf "INFO - Running S-TFXcan for %s\n" ${phenotype}

function tfxcan_association(){

    model_db=${1}
    covariances=${2}
    gwas_folder=${3}
    gwas_pattern=${4}
    outputfile=${5}
    exec_file=${6}

    python3 ${exec_file} \
        --model_db_path ${model_db} \
        --model_db_snp_key varID \
        --covariance ${covariances} \
        --gwas_folder ${gwas_folder} \
        --gwas_file_pattern ${gwas_pattern} \
        --chromosome_column chrom \
        --snp_column SNP \
        --effect_allele_column alt \
        --non_effect_allele_column ref \
        --zscore_column zscore \
        --keep_non_rsid \
        --output_file ${outputfile} \
        --verbosity 9 \
        --throw
}

export -f tfxcan_association

parallel -j 1 --joblog ${job_log} --skip-first-line --colsep '\t' -a ${models_list} "tfxcan_association {2} {3} ${gwas_folder} ${gwas_pattern} ${output_directory}/{1}.${phenotype}.summaryTFXcan.csv ${exec_file}"

status=$?

printf "\nINFO - Finished running summary TFXcan for\n"
echo "Exit status of summary TFXcan is: $status"





    #     --gwas_folder ${base_dir}/data/GWAS \
    # --gwas_file_pattern ".*_Schumacher.gwas_ss.txt.gz" \


    # awk '{gsub(/:/, "_"); print}' Covariances.varID.txt > Covariances.varID.txt

    # sed -i -e 's/:/_/g' Covariances.varID.txt

    #sed -e 's/:/_/g' Covariances.txt > Covariances.varID.txt

    # gzip -d predict_db_EUR_AR_Prostate_logistic_filtered.txt.gz
    # sed -e 's/:/_/g' predict_db_EUR_AR_Prostate_logistic_filtered.txt > Covariances.varID.txt