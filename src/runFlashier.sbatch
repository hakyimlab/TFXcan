#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G    
#SBATCH --job-name=flashier
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/Enpact/logs/flashier.out
#SBATCH --error=/beagle3/haky/users/temi/projects/Enpact/logs/flashier.err
#SBATCH --open-mode=append
#SBATCH --time=06:00:00	
#SBATCH --partition=caslake

date
module load parallel

source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/tfxcan-env

datafile=${1}
column_for_rownames=${2}
output_basename=${3}
priorL=${4}
priorF=${5}
transpose=${6}
greedy_Kmax=${7}

if [[ ${transpose} == "TRUE" ]]; then
    Rscript /beagle3/haky/users/temi/projects/Enpact/src/runFlashier.R \
        --data ${datafile} \
        --column_for_rownames ${column_for_rownames} \
        --output_basename ${output_basename} \
        --priorL ${priorL} \
        --priorF ${priorF} \
        --transpose \
        --greedy_Kmax ${greedy_Kmax}
elif [[ ${transpose} == "FALSE" ]]; then
    Rscript /beagle3/haky/users/temi/projects/Enpact/src/runFlashier.R \
        --data ${datafile} \
        --column_for_rownames ${column_for_rownames} \
        --output_basename ${output_basename} \
        --priorL ${priorL} \
        --priorF ${priorF} \
        --greedy_Kmax ${greedy_Kmax}
else
    echo "Transpose argument must be TRUE or FALSE"
    exit 1
fi

status=$?
echo "Exit status of lifting over is: $status"