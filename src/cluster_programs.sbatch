#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4G    
#SBATCH --job-name=cluster_flash
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/Enpact/logs/cluster_flash.out
#SBATCH --error=/beagle3/haky/users/temi/projects/Enpact/logs/cluster_flash.err
#SBATCH --time=06:00:00	
#SBATCH --partition=caslake

date

source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/renviron

datafile=${1}
outputbasename=${2}
normalization=${3}

start=`date +%s`

Rscript /beagle3/haky/users/temi/projects/Enpact/src/cluster_programs.R \
    --flash_results ${datafile} --output_basename ${outputbasename}\
    --normalization ${normalization}

end=`date +%s`

runtime=$((end-start))
   
status=$?

echo "Exit status of consensus clustering is: $status"
echo "Runtime is: $runtime seconds"

# example usage 
# sbatch /beagle3/haky/users/temi/projects/Enpact/src/repeat_flash.sbatch /beagle3/haky/users/temi/projects/Enpact/data/tenerife/PrCa.tfxcan.zratios.matrix.rds 1 /beagle3/haky/users/temi/projects/Enpact/data/tenerife/prca_flashier_repeats ebnm_point_exponential ebnm_point_exponential FALSE 1000

# sbatch /beagle3/haky/users/temi/projects/Enpact/src/cluster_programs.sbatch /beagle3/haky/users/temi/projects/Enpact/data/tenerife/prca_flashier_repeats.ebnm_point_exponential-ebnm_point_exponential.1000Iters.rds.gz /beagle3/haky/users/temi/projects/Enpact/misc/prca_program_clusters.1000.txt.gz /beagle3/haky/users/temi/projects/Enpact/misc/prca_program_matrix.1000.rds.gz /beagle3/haky/users/temi/projects/Enpact/misc/prca_program_loci.1000.txt.gz












# '/beagle3/haky/users/temi/projects/Enpact/data/tenerife/prca_flashier_repeats.ebnm_point_exponential-ebnm_point_exponential.1000Iters.rds.gz'