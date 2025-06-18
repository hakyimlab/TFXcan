#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=12G    
#SBATCH --job-name=consensusFlashier
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/Enpact/logs/consensusFlashier.out
#SBATCH --error=/beagle3/haky/users/temi/projects/Enpact/logs/consensusFlashier.err
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem

date

source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/tfxcan-env
module load parallel

datafile=${1}
output_basename=${2}
priorL=${3}
priorF=${4}
greedy_Kmax=${5}
batch_list=${6}

splits_directory='/beagle3/haky/users/temi/projects/Enpact/data/tenerife/repeat_flash'

parallel -a ${batch_list} --jobs 20 --joblog /beagle3/haky/users/temi/projects/Enpact/logs/consensusFlashier.log \
    Rscript /beagle3/haky/users/temi/projects/Enpact/src/repeat_flash.R \
    --data ${datafile} \
    --splits ${splits_directory}/prca_risk.random_subsets.{}.rds \
    --batch {} \
    --output_basename ${output_basename} \
    --priorL ${priorL} \
    --priorF ${priorF} \
    --greedy_Kmax ${greedy_Kmax}

status=$?

if [ $status -eq 0 ]; then
    echo "INFO - Compiling job"
    Rscript /beagle3/haky/users/temi/projects/Enpact/src/gather_repeats.R \
        --output_basename ${output_basename} \
        --batch_list ${batch_list} \
        --output_file_basename ${output_basename}.repeat_flashier.rds.gz
fi

status=$?
echo "Exit status of job is: $status"




# if [[ ${run_in_parallel} == 'TRUE' ]]; then
#     printf "[INFO] Running consensus_flashier.R in parallel mode\n"
#     Rscript /beagle3/haky/users/temi/projects/Enpact/src/consensus_flashier.R \
#         --data ${datafile} \
#         --output_basename ${output_basename} \
#         --priorL ${priorL} \
#         --priorF ${priorF} \
#         --greedy_Kmax ${greedy_Kmax} \
#         --subset 0.8 \
#         --niterations 500
# elif [[ ${run_in_parallel} == 'FALSE' ]]; then
#     printf "[INFO] Running repeat_flash.R in serial mode\n"
#     Rscript /beagle3/haky/users/temi/projects/Enpact/src/repeat_flash.R \
#         --data ${datafile} \
#         --output_basename ${output_basename} \
#         --priorL ${priorL} \
#         --priorF ${priorF} \
#         --greedy_Kmax ${greedy_Kmax} \
#         --subset 0.8 \
#         --niterations 1000
#     else
#         echo "Transpose argument must be TRUE or FALSE"
#         exit 1
#     fi
# fi

