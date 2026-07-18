#!/bin/bash

# Author: Temi
# Description: merge raw per-region Enformer predictions (from enformer_predict.sbatch) into one hdf5, via the TFXcan-snakemake repo's enformer_merge.py
# Usage: sbatch enformer_merge.sbatch aggregation_config_yaml epigenomes_input_dir output_h5_name

#SBATCH --job-name=enformer_merge
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/enformer_merge.err
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/enformer_merge.out
#SBATCH --account=pi-haky

date

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/beagle3/haky/users/shared_software/TFXcan-pipeline-tools/lib

# make sure that in your config, provider is set to "local" before running this script
python3 /beagle3/haky/users/temi/projects/TFXcan-snakemake/workflow/enformer/enformer_merge.py --config ${1} --output_directory ${2} --output_filename ${3}