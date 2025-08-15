#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G    
#SBATCH --job-name=calculate_auprc
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/Enpact/logs/calculate_auprc.out
#SBATCH --error=/beagle3/haky/users/temi/projects/Enpact/logs/calculate_auprc.err
#SBATCH --time=02:00:00
#SBATCH --partition=caslake

date

source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/renviron
module load parallel



Rscript /beagle3/haky/users/temi/projects/TFXcan/src/calculate_auprc.R \
    --transcription_factors ${1} \
    --tissues ${2} \
    --path_pattern ${3} \
    --output_file ${4}
    
status=$?
echo "Exit status of job is: $status"