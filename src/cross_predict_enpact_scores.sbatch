#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G    
#SBATCH --job-name=cross_predict
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/cross_predict.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/cross_predict.err
#SBATCH --time=06:00:00
#SBATCH --partition=caslake
#SBATCH --open-mode=append

date
source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools
module load parallel/

# input_metadata="/beagle3/haky/users/temi/projects/TFXcan/experiments/validation/enpact_test_models.chippeaks.tsv" #"/beagle3/haky/users/temi/projects/TFXcan/experiments/validation/enpact_test_models.tsv"
# weights_file="/beagle3/haky/users/temi/projects/TFXcan/data/enpact/weights/ENPACT_734_2025-04-24.compiled_weights.with_intercept.lambda.1se.txt.gz" #"/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_CHIPPEAKS_2026-01-02/statistics/ENPACT_CHIPPEAKS_2026-01-02.compiled_weights.lambda.1se.txt.gz"
# output_basename="/scratch/midway3/temi/scores/motifs/predictions"
input_metadata=${1}
weights_file=${2}
output_basename=${3}

# create data
parallel -j 20 --joblog /beagle3/haky/users/temi/projects/TFXcan/logs/cross_predict.job --colsep '\t' -a ${input_metadata} "Rscript /beagle3/haky/users/temi/projects/TFXcan/src/cross_predict_enpact_scores.R --tf_tissue {1} --input_file {2} --output_basename ${output_basename} --weights_file ${weights_file}" 

status=$?