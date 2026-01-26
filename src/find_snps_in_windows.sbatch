#!/bin/bash

#SBATCH --job-name=find_snps_in_windows
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/find_snps_in_windows.err
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/find_snps_in_windows.out
#SBATCH --account=pi-haky

date
source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/renviron

Rscript /beagle3/haky/users/temi/projects/TFXcan/src/find_snps_in_windows.R

status=$?
echo "INFO - Exit status is: $status"