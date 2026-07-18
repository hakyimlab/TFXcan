#!/bin/bash

# Author: Temi
# Date: Monday January 26 2026
# Description: merge per-model S-PrediXcan (TFXcan) output CSVs from run_summary_tfxcan.sbatch into one results file
# Usage: sbatch gather_sTFXcan.sbatch file_glob_pattern phenotype output_file

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=gather_sTFXcan
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/gather_sTFXcan.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/gather_sTFXcan.err
#SBATCH --time=01:00:00
#SBATCH --partition=caslake

date

slurm_workdir=${SLURM_SUBMIT_DIR}
# SLURM_O_WORKDIR=${slurm_workdir}/run
# mkdir -p ${SLURM_O_WORKDIR}
# cd $SLURM_O_WORKDIR
echo "INFO - Working directory is $SLURM_SUBMIT_DIR"
echo "INFO - Jobid is $SLURM_JOBID"
echo "INFO - Running on host `hostname`"

source ~/.bashrc

# variables
fpattern="${1}"
phenotype=${2}
output_file=${3}

echo "${fpattern}"
echo ${phenotype}
echo ${output_file}

conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools

# calls the merge script from the sibling TFXcan-snakemake repo, not anything local to this src/
Rscript /beagle3/haky/users/temi/projects/TFXcan-snakemake/workflow/process/collect_summaryTFXcan_results.R \
    --input_files_pattern "${fpattern}" \
    --phenotype ${phenotype} \
    --output_file ${output_file}

status=$?

printf "\nINFO - Finished collecting summary TFXcan results\n"
echo "Exit status of collecting summary TFXcan results is: $status"