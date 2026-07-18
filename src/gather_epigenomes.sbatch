#!/bin/bash

# Author: Temi
# Date: Friday August 15 2025
# Description: aggregate per-locus Enformer predictions into a feature CSV; calls the TFPred-snakemake repo's copy of aggregate_epigenomes.py, NOT the local src/aggregate_epigenomes.py
# Usage: sbatch gather_epigenomes.sbatch loci_file reference_epigenome_dir output_file

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G    
#SBATCH --job-name=gather_epigenomes
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/gather_epigenomes.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/gather_epigenomes.err
#SBATCH --time=01:00:00
#SBATCH --partition=caslake

date

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools

python3 /beagle3/haky/users/temi/projects/TFPred-snakemake/workflow/src/aggregate_epigenomes.py --loci_file ${1} --reference_epigenome_dir ${2} --output_file ${3} --use_multiprocessing --pad_bins 0 # pad with zeros here

status=$?
echo "Exit status of job is: $status"