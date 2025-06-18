#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4G    
#SBATCH --job-name=calculate_distance
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/Enpact/logs/calculate_distance.out
#SBATCH --error=/beagle3/haky/users/temi/projects/Enpact/logs/calculate_distance.err
#SBATCH --time=01:00:00
#SBATCH --partition=caslake

date

source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/renviron

# use parallel
module load parallel

output_basename=${1}
list_of_transcription_factors=${2}

start=`date +%s`

# run the R script
parallel -j 50 "Rscript /beagle3/haky/users/temi/projects/Enpact/src/calculateDistanceToMotifs.R \
    --transcription_factor {} --output_basename ${output_basename}" :::: ${list_of_transcription_factors}

# Rscript /beagle3/haky/users/temi/projects/Enpact/src/calculateDistanceToMotifs.R \
#     --transcription_factor ${transcription_factor} --output_basename ${output_basename}

end=`date +%s`

runtime=$((end-start))
   
status=$?

echo "Exit status of lifting over is: $status"
echo "Runtime is: $runtime seconds"