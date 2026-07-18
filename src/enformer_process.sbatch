#!/bin/bash

# Author: Temi
# Description: post-process a merged Enformer predictions hdf5 (sum by haplotype) into a final matrix + metadata, via the TFXcan-snakemake repo's enformer_process.py
# Usage: sbatch enformer_process.sbatch merged_h5_file output_basename

#SBATCH --job-name=enformer_process
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/enformer_process.err
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/enformer_process.out
#SBATCH --account=pi-haky

date

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/beagle3/haky/users/shared_software/TFXcan-pipeline-tools/lib

# make sure that in your config, provider is set to "local" before running this script
python3 /beagle3/haky/users/temi/projects/TFXcan-snakemake/workflow/enformer/enformer_process.py --merged_h5_file ${1} --process_by_haplotype --process_function 'sum' --output_basename ${2}

