#!/bin/bash
# Author: Temi
# Date: Tuesday February 24 2026
# Description: wrapper for correlate_binding_expression.R
# Usage: sbatch correlate_binding_expression.sbatch <binding_matrix.rds> <gene_expression.tsv.gz> <output_file_basename>

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=24G
#SBATCH --job-name=correlate_binding_expression
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/correlate_binding_expression.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/correlate_binding_expression.err
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem

date

source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/renviron

binding_matrix=${1}
gene_expression=${2}
output_file_basename=${3}

Rscript /beagle3/haky/users/temi/projects/TFXcan/src/correlate_binding_expression.R \
    --binding_matrix ${binding_matrix} \
    --gene_expression ${gene_expression} \
    --output_file_basename ${output_file_basename}


status=$?
echo "Exit status of job is: $status"