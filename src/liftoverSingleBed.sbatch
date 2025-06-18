#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G    
#SBATCH --job-name=liftover
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/Enpact/logs/liftover.out
#SBATCH --error=/beagle3/haky/users/temi/projects/Enpact/logs/liftover.err
#SBATCH --time=00:30:00	
#SBATCH --partition=caslake

date

liftover_exe=${1}
chain_file=${2}
input_bed=${3}
output_prefix=${4}
unmapped_prefix=${5}

${liftover_exe} ${input_bed} ${chain_file} ${output_prefix} ${unmapped_prefix}

status=$?

printf "\nINFO - Finished lifting over all files\n"
echo "Exit status of lifting over is: $status"