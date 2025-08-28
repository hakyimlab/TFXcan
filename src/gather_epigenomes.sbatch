#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=12G    
#SBATCH --job-name=gather_epigenomes
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/Enpact/logs/gather_epigenomes.out
#SBATCH --error=/beagle3/haky/users/temi/projects/Enpact/logs/gather_epigenomes.err
#SBATCH --time=04:00:00
#SBATCH --partition=caslake

date

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools

python3 /beagle3/haky/users/temi/projects/TFPred-snakemake/workflow/src/aggregate_epigenomes.py --loci_file ${1} --reference_epigenome_dir ${2} --output_file ${3} --use_multiprocessing --pad_bins 0 # pas with zeros here

status=$?
echo "Exit status of job is: $status"