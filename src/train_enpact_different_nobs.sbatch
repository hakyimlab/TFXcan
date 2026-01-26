#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=12G    
#SBATCH --job-name=train_different_nobs
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/train_different_nobs.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/train_different_nobs.err
#SBATCH --time=06:00:00
#SBATCH --partition=bigmem

date
source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools
module load parallel/

full_data=${1}
output_basename=${2}
reference_epigenome=${3} # /beagle3/haky/data/enformer-reference-epigenome

array_of_nobs=(50 100 200 500 1000 3000 5000 8000 10000 12500 15000 18000) # (50 100 200)   #
n_array="${#array_of_nobs[@]}"

# create data
parallel -j ${n_array} --joblog /beagle3/haky/users/temi/projects/TFXcan/logs/train_different_nobs.job "Rscript /beagle3/haky/users/temi/projects/TFXcan/src/create_training_data.R --full_data ${full_data} --output_basename ${output_basename}_{} --num_bound {}" ::: "${array_of_nobs[@]}" 

status=$?

echo "Exit status of creating data is: $status"

# collect epigenomes
if [ ${status} -eq 0 ]; then
    parallel -j ${n_array} --joblog /beagle3/haky/users/temi/projects/TFXcan/logs/train_different_nobs.job "python3 /beagle3/haky/users/temi/projects/TFPred-snakemake/workflow/src/aggregate_epigenomes.py --loci_file ${output_basename}_{}.loci.tsv --reference_epigenome_dir ${reference_epigenome} --output_file ${output_basename}_{}.epigenome.csv.gz --use_multiprocessing --pad_bins 0" ::: "${array_of_nobs[@]}"
fi

echo "Exit status of extracting epigenome is: $status"

# merge with ground truth 
if [ ${status} -eq 0 ]; then
    parallel -j ${n_array} --joblog /beagle3/haky/users/temi/projects/TFXcan/logs/train_different_nobs.job "Rscript /beagle3/haky/users/temi/projects/TFXcan/src/merge_epigenome_ground_truth.R --epigenome ${output_basename}_{}.epigenome.csv.gz --ground_truth ${output_basename}_{}.ground_truth.tsv --output_filename ${output_basename}_{}.train.tsv" ::: "${array_of_nobs[@]}"
fi

echo "Exit status of merging is: $status"


#train enpact

if [ ${status} -eq 0 ]; then
    parallel -j ${n_array} --joblog /beagle3/haky/users/temi/projects/TFXcan/logs/train_different_nobs.job "Rscript /beagle3/haky/users/temi/projects/TFXcan/src/train_enpact.R --train_data_file ${output_basename}_{}.train.tsv --metadata score_{} --weights_file ${output_basename}_{}.weights.tsv.gz --model_rds_file ${output_basename}_{}.glmnet.rds" ::: "${array_of_nobs[@]}"
fi

echo "Exit status of training is: $status"

