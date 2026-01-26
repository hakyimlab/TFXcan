#!/bin/bash

#SBATCH --job-name=estimate_null_ppi
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/estimate_null_ppi.err
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/estimate_null_ppi.out
#SBATCH --account=pi-haky

date
source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/renviron

Rscript /beagle3/haky/users/temi/projects/TFXcan/src/estimate_null_ppi.R --ppi_db ${1} --query_db ${2} --programs ${3} --outputfile ${4}

status=$?
echo "INFO - Exit status is: $status"