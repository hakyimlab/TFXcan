#!/bin/bash
# Author: Temi
# Date: Wednesday June 18 2025
# Description: cluster wrapper for estimate_hsq_one_grm.R (single merged-GRM heritability estimates). No CLI
#   args; all paths are hardcoded inside the R script.
# Usage: sbatch estimate_hsq_one_grm.sbatch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G    
#SBATCH --job-name=estimate_hsq
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/Enpact/logs/estimate_hsq_one_grm.out
#SBATCH --error=/beagle3/haky/users/temi/projects/Enpact/logs/estimate_hsq_one_grm.err
#SBATCH --time=02:00:00
#SBATCH --partition=caslake

date
source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/renviron

Rscript /beagle3/haky/users/temi/projects/Enpact/src/estimate_hsq_one_grm.R

status=$?
printf "INFO - All jobs are done. Exit status is $status.\n"