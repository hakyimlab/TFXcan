#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=12G    
#SBATCH --job-name=count_training_test_overlaps
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/Enpact/logs/count_training_test_overlaps.out
#SBATCH --error=/beagle3/haky/users/temi/projects/Enpact/logs/count_training_test_overlaps.err
#SBATCH --time=01:00:00
#SBATCH --partition=caslake

date

source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/renviron
module load parallel

Rscript /beagle3/haky/users/temi/projects/TFXcan/src/count_training_test_overlaps.R \
    --enpact_metadata_file ${1} \
    --files_directory ${2} \
    --output_file ${3}
    
status=$?
echo "Exit status of job is: $status"