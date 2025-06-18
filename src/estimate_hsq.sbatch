#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4G    
#SBATCH --job-name=estimate_hsq
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/Enpact/logs/estimate_hsq.out
#SBATCH --error=/beagle3/haky/users/temi/projects/Enpact/logs/estimate_hsq.err
#SBATCH --time=02:00:00
#SBATCH --partition=caslake

date
source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/renviron

Rscript /beagle3/haky/users/temi/projects/Enpact/src/estimate_hsq.R

status=$?
printf "INFO - All jobs are done. Exit status is $status.\n"