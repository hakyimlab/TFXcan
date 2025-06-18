#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4G    
#SBATCH --job-name=make_grm_single
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/Enpact/logs/make_grm_single.out
#SBATCH --error=/beagle3/haky/users/temi/projects/Enpact/logs/make_grm_single.err
#SBATCH --time=01:00:00
#SBATCH --partition=caslake

date

module load parallel

input_directory=${1}
output_directory=${2}
input_basenames=${3}

IFS=',' read -r -a input_array <<< "${input_basenames}"

gcta_executable="/beagle3/haky/users/temi/software/gcta-1.94.4-linux-kernel-3-x86_64/gcta64"

parallel -j 5 "
  ${gcta_executable} --bfile ${input_directory}/{1} --make-grm --out ${output_directory}/{1} --reml-no-constrain" ::: ${input_array[@]}

status=$?
printf "INFO - All jobs are done. Exit status is $status.\n"