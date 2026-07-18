#!/bin/bash
# Author: Temi
# Date: Friday August 15 2025
# Description: cluster wrapper for calculate_wilcoxon_test.R (bound vs. unbound Wilcoxon test per TF-tissue model).
# Usage: sbatch calculate_wilcoxon_test.sbatch {models_list} {path_pattern} {output_file}

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G    
#SBATCH --job-name=calculate_wilcoxon_test
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/calculate_wilcoxon_test.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/calculate_wilcoxon_test.err
#SBATCH --time=02:00:00
#SBATCH --partition=caslake

date

source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/renviron
module load parallel

Rscript /beagle3/haky/users/temi/projects/TFXcan/src/calculate_wilcoxon_test.R \
    --models_list ${1} \
    --path_pattern ${2} \
    --output_file ${3}
    
status=$?
echo "Exit status of job is: $status"