#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G    
#SBATCH --job-name=model_grm
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/Enpact/logs/model_grm.out
#SBATCH --error=/beagle3/haky/users/temi/projects/Enpact/logs/model_grm.err
#SBATCH --time=04:00:00
#SBATCH --partition=caslake

date

module load parallel

output_basename=${1}
grm_list=${2}
phenotype=${3}
number_of_phenotypes=${4}

# IFS=',' read -r -a input_array <<< "${input_basenames}"

gcta_executable="/beagle3/haky/users/temi/software/gcta-1.94.4-linux-kernel-3-x86_64/gcta64"


parallel -j 20 \
    --joblog /beagle3/haky/users/temi/projects/Enpact/logs/model_grm.log \
    "${gcta_executable} --reml --grm "/project2/haky/Data/baca_cwas/vcfs/hg38/merged_files/grm/baca_individuals.plink.merged" --pheno ${phenotype} --out ${output_basename}_{} --mpheno {}" ::: `seq ${number_of_phenotypes}`

# parallel -j 20 \
#     --joblog /beagle3/haky/users/temi/projects/Enpact/logs/model_grm.log \
#     "${gcta_executable} --reml --mgrm ${grm_list} --pheno ${phenotype} --out ${output_basename}_{} --mpheno {} --reml-no-constrain" ::: `seq ${number_of_phenotypes}`

status=$?
printf "INFO - All jobs are done. Exit status is $status.\n"

# /beagle3/haky/users/temi/software/gcta-1.94.4-linux-kernel-3-x86_64/gcta64 --reml --grm /project2/haky/Data/baca_cwas/vcfs/hg38/grm_files/chr22.dose --pheno /beagle3/haky/users/temi/projects/Enpact/experiments/heritability/phenotypes/chr10_100185843_100186293.peaks.txt --out hey


