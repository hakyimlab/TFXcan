#!/bin/bash

# Author: Temi
# Description: S-PrediXcan on EnpactScores
# Usage: sbatch s-predixcan.sbatch
# Date: Wed Apr 3 2024
# Dependencies: 

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=aggregate_enpact
#SBATCH --account=pi-haky
#SBATCH --output=/project/haky/users/temi/projects/Enpact/logs/aggregate_enpact.out
#SBATCH --error=/project/haky/users/temi/projects/Enpact/logs/aggregate_enpact.err
#SBATCH --time=02:00:00	
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

agg_exec="/beagle3/haky/users/shared_pipelines/enformer_pipeline_aggregate/scripts/aggregate/aggregate.py"

python3 ${agg_exec} --metadata_file ${1} --agg_types "aggByCollect" --output_directory ${2} --hpc 'caslake' --parsl_executor 'local'



# sbatch aggregate.sbatch /project/haky/users/temi/projects/Enpact/data/enpact/files/aggregation.FOXA1_Prostate.config.json /project/haky/users/temi/projects/Enpact/data/enpact/files

# sbatch aggregate.sbatch /project/haky/users/temi/projects/Enpact/data/enpact/files/aggregation.HOXB13_Prostate.config.json /project/haky/users/temi/projects/Enpact/data/enpact/files

# sbatch aggregate.sbatch /project/haky/users/temi/projects/Enpact/data/enpact/files/aggregation.GATA2_Prostate.config.json /project/haky/users/temi/projects/Enpact/data/enpact/files