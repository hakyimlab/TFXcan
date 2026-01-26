#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G    
#SBATCH --job-name=train_enpact
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/accessibility.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/accessibility.err
#SBATCH --time=02:00:00
#SBATCH --partition=caslake

date

source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/renviron

Rscript /beagle3/haky/users/temi/projects/TFXcan/src/accessibility_effect_on_binding.ttest.R

status=$?
echo "Exit status of job is: $status"