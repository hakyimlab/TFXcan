#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=aggregate_epigenomes
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/Enpact/logs/aggregate_epigenomes.out
#SBATCH --error=/beagle3/haky/users/temi/projects/Enpact/logs/aggregate_epigenomes.err
#SBATCH --time=04:00:00	
#SBATCH --partition=caslake

# module load openmpi
# module load parallel

loci_file=${1}
reference_epigenome_dir=${2}
output_file=${3}
pad_bins=${4}

date

printf "Starting to run\n"

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools

agg_exec="/beagle3/haky/users/temi/projects/Enpact/src/aggregate_epigenomes.py"

python3 ${agg_exec} \
    --loci_file ${loci_file} \
    --reference_epigenome_dir ${reference_epigenome_dir} \
    --output_file ${output_file} \
    --use_multiprocessing \
    --pad_bins 0


# sbatch 