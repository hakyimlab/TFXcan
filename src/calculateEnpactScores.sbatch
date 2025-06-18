#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G    
#SBATCH --job-name=calculate_enpact_scores
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/Enpact/logs/calculate_enpact_scores.out
#SBATCH --error=/beagle3/haky/users/temi/projects/Enpact/logs/calculate_enpact_scores.err
#SBATCH --time=04:00:00	
#SBATCH --partition=caslake

date

source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/renviron

data_directory=${1}
individuals_list=${2}
enpact_weights=${3}
output_file=${4}
files_pattern=${5}

start=`date +%s`

Rscript /beagle3/haky/users/temi/projects/Enpact/src/calculateEnpactScores.R \
    --data_directory ${data_directory} --individuals_list ${individuals_list} \
    --enpact_weights ${enpact_weights} --output_file ${output_file} \
    --files_pattern ${files_pattern}


end=`date +%s`

runtime=$((end-start))
   
status=$?

echo "Exit status of predicing binding with Enpact is: $status"
echo "Runtime is: $runtime seconds"

# example usage  for 1kg individuals
# sbatch /beagle3/haky/users/temi/projects/Enpact/src/calculateEnpactScores.sbatch /beagle3/haky/users/temi/data/1KG_AR_prostate /beagle3/haky/users/temi/projects/TFXcan-snakemake/metadata/EUR_individuals.1KG.txt /beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2025-04-24/statistics/ENPACT_734_2025-04-24.compiled_weights.lambda.1se.txt.gz /beagle3/haky/users/temi/projects/Enpact/misc/reruns/enpact_predictions/ENPACT_48.1KG.CWAS_ARBS.predictions.2025-04-28.rds.gz _aggByCollect_AR_Prostate.csv

# sbatch /beagle3/haky/users/temi/projects/Enpact/src/cluster_programs.sbatch /beagle3/haky/users/temi/projects/Enpact/data/tenerife/prca_flashier_repeats.ebnm_point_exponential-ebnm_point_exponential.1000Iters.rds.gz /beagle3/haky/users/temi/projects/Enpact/misc/prca_program_clusters.1000.txt.gz /beagle3/haky/users/temi/projects/Enpact/misc/prca_program_matrix.1000.rds.gz /beagle3/haky/users/temi/projects/Enpact/misc/prca_program_loci.1000.txt.gz












# '/beagle3/haky/users/temi/projects/Enpact/data/tenerife/prca_flashier_repeats.ebnm_point_exponential-ebnm_point_exponential.1000Iters.rds.gz'