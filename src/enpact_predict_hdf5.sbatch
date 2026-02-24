#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=12G    
#SBATCH --job-name=enpact_predict_hdf5
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/enpact_predict_hdf5.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/enpact_predict_hdf5.err
#SBATCH --time=06:00:00
#SBATCH --partition=bigmem

date

source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/renviron
module load parallel

loci_list=${1}
weights=${2}
individuals_list=${3}
input_directory=${4}
output_directory=${5}

echo ${individuals_list}

readarray -t myarray < ${individuals_list}

parallel --colsep "\t" --jobs 50 --joblog /beagle3/haky/users/temi/projects/TFXcan/logs/enpact_predict_hdf5.log \
    Rscript /beagle3/haky/users/temi/projects/TFXcan/src/enpact_predict_hdf5.R \
    --loci_list ${loci_list} \
    --weights ${weights} \
    --individual {} \
    --input_directory ${input_directory} \
    --output_basename ${output_directory}/{} ::: ${myarray[@]}

status=$?
echo "Exit status of job is: $status"