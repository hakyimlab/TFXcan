#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=16G    
#SBATCH --job-name=create_h5database_from_csvs
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/create_h5database_from_csvs.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/create_h5database_from_csvs.err
#SBATCH --time=06:00:00
#SBATCH --partition=bigmem

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
metadata=${1}
output_basename=${2}
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools
python3 /beagle3/haky/users/temi/projects/TFXcan/src/create_h5database_from_csvs.py \
    --metadata ${metadata} \
    --output_basename ${output_basename}

status=$?

printf "\nINFO - Finished collecting summary TFXcan results\n"
echo "Exit status of collecting summary TFXcan results is: $status"