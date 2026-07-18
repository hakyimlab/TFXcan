#!/bin/bash
# Author: Temi
# Description: wrapper for correlate_binding_expression.acrs.R
# Usage: sbatch correlate_binding_expression.acrs.sbatch <binding_matrix.tsv> <expression_matrix.tsv> <temp_output_basename> <final_output_file>
# NOTE: takes 4 positional args (matches the R script's 4 flags); some notebook call sites use only 3 args for this
#   script -- those look stale relative to the current on-disk version.

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=correlate_acrs
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/correlate_acrs.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/correlate_acrs.err
#SBATCH --time=06:00:00
#SBATCH --partition=caslake

date

source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/renviron

binding_matrix=${1}
gene_expression=${2}
temp_output_basename=${3}
final_output_file=${4}


Rscript /beagle3/haky/users/temi/projects/TFXcan/src/correlate_binding_expression.acrs.R \
    --binding_matrix ${binding_matrix} \
    --expression_matrix ${gene_expression} \
    --temp_output_basename ${temp_output_basename} \
    --final_output_file ${final_output_file}

status=$?
echo "Exit status of job is: $status"