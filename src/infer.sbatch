#!/bin/bash

# Author: Temi
# Description: S-PrediXcan on EnpactScores
# Usage: sbatch s-predixcan.sbatch
# Date: Wed Apr 3 2024
# Dependencies: 

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=infer_epigenome
#SBATCH --account=pi-haky
#SBATCH --output=/project/haky/users/temi/projects/Enpact/logs/infer_epigenome.out
#SBATCH --error=/project/haky/users/temi/projects/Enpact/logs/infer_epigenome.err
#SBATCH --time=03:00:00	
#SBATCH --partition=caslake

# module load openmpi
# module load parallel

date
slurm_workdir=${SLURM_SUBMIT_DIR}
SLURM_O_WORKDIR=${slurm_workdir}/run
mkdir -p ${SLURM_O_WORKDIR}
echo Working directory is $SLURM_O_WORKDIR
cd $SLURM_O_WORKDIR

echo Jobid: $SLURM_JOBID
echo Running on host `hostname`

printf "Starting to run\n"

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools

# variables
exec_py=${1}
input_file=${2}
ref_epigenome=${3}
output_file=${4}

python3 ${exec_py} --loci_file ${input_file} --reference_epigenome_dir ${ref_epigenome} --output_file ${output_file} --use_multiprocessing

printf "Done\n"