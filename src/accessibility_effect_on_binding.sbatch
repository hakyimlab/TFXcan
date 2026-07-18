#!/bin/bash
# Author: Temi
# Date: Monday January 26 2026
# Description: wrapper for accessibility_effect_on_binding.ttest.R (the Wilcoxon-test variant); no args, all
#   inputs/outputs are hardcoded in the R script itself
# Usage: sbatch accessibility_effect_on_binding.sbatch
# NOTE: job-name below says "train_enpact" -- stale, leftover from a copy-pasted sbatch template

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