#!/bin/bash

# Author: Temi
# Description: Predict using cwas in PrediXcan framework
# Usage: sbatch predict_baca_with_cwas.sbatch
# Date: Mon July 24 2023
# Dependencies: 

#SBATCH --nodes=6
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=predict_baca_with_cwas
#SBATCH --account=pi-haky
#SBATCH --output=/project/haky/users/temi/projects/Enpact/logs/predict_cwas_scores.out
#SBATCH --error=/project/haky/users/temi/projects/Enpact/logs/predict_cwas_scores.err
#SBATCH --time=01:30:00	
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
conda activate /beagle3/haky/users/temi/software/conda_envs/imlabtools 
module load openmpi
echo "INFO - using this mpirun: `which mpirun`"

# paths
db_folder=${1}
txt_genotypes=${2}
txt_samples=${3}
exec_file=${4}
output_folder=${5}

# parallelization set up 
printf '%s\n' ${db_folder}/*.db > ./dbs.txt
dbs_list=./dbs.txt
x=`wc -l ${dbs_list} | awk '{print $1}'`
y=${SLURM_NNODES} 
ll=$(( ($x + $y - 1) / $y ))

split --lines="${ll}" --numeric-suffixes=1 --suffix-length=1 "${dbs_list}" "db."

# prepare nodes
nodelist=$(scontrol show hostname $SLURM_NODELIST)
printf "%s\n" "${nodelist[@]}" > local_hostfiles.txt
# Increase value of suffix-length if more than 99 jobs
split --lines=1 --numeric-suffixes=1 --suffix-length=1 local_hostfiles.txt "local_hostfile."

# the predict function is defined here
function predict () {
    dblist=${1}
    execfile=${2}
    txt_geno=${3}
    txt_ids=${4}
    outfolder=${5}
    readarray -t arrr < ${dblist}

    for model_db in ${arrr[@]}; do
        m_name=$( echo ${model_db} | rev | cut -d '/' -f 1 | cut -d '_' -f 1 | rev )
        m_name=${m_name%.*}
        echo "INFO - Running ${m_name}"

        # if [[ ${m_name} != 'top1' ]]; then
        #     continue
        # else
        python ${execfile} \
        --model_db_path ${model_db} \
        --text_genotypes ${txt_geno} \
        --text_sample_ids ${txt_ids} \
        --prediction_output "${outfolder}/${m_name}/baca_cwas_predict.txt" \
        --prediction_summary_output "${outfolder}/${m_name}/baca_cwas_summary.txt" \
        --verbosity 9 \
        --model_db_snp_key varID \
        --throw  
    done
    wait
}

export -f predict

# now, run and godspeed
for suf in `seq 1 ${SLURM_NNODES}`; do
  (
    echo "INFO - Job ${suf} running on `cat local_hostfile.${suf}` using db.${suf}"
    mpirun -np 1 --hostfile "local_hostfile.${suf}" bash -c "predict db.${suf} ${exec_file} ${txt_genotypes} ${txt_samples} ${output_folder}"

  ) & sleep 1
done
wait