#!/bin/bash

# Author: Temi
# Description: Convert vcfs to dosages.txt file
# Usage: sbatch convert_vcf_to_dosages.sbatch
# Date: Fri July 21 2023
# Dependencies: 

#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=convert_vcf_to_dosages
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/convert_vcf_to_dosages.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/convert_vcf_to_dosages.err
#SBATCH --time=00:30:00	
#SBATCH --partition=caslake

module load openmpi
module load parallel

slurm_workdir=${SLURM_SUBMIT_DIR}
SLURM_O_WORKDIR=${slurm_workdir}/run
mkdir -p ${SLURM_O_WORKDIR}
echo Working directory is $SLURM_O_WORKDIR
cd $SLURM_O_WORKDIR

echo Jobid: $SLURM_JOBID
echo Running on host `hostname`

printf "Starting to run\n"

source ~/.bashrc
conda activate compbio-tools

# variables
mpiexec="/opt/cray/pe/pals/1.1.7/bin/mpiexec"

# 
parent_folder=/project2/haky/Data/baca_cwas
vcfs_path=${parent_folder}/imputed_vcf_hg38_snps_only
plink_results_folder=${parent_folder}/plink_geno
formatted_files_folder=${parent_folder}/formatted_geno

# exec files
exec_file=/beagle3/haky/users/temi/projects/TFXcan/scripts/manipulate_vcfs/format_files.sh
rscript=/home/temi/miniconda3/envs/r-env/bin/Rscript
utility_functions=/beagle3/haky/users/temi/projects/TFXcan/scripts/manipulate_vcfs/vcf_utils.R

mkdir -p ${plink_results_folder}
mkdir -p ${formatted_files_folder}

# create a vcf_list and split based on how you want to parallelize
printf '%s\n' ${vcfs_path}/*.vcf.gz > ./chr_vcfs.txt
vcf_list=./chr_vcfs.txt # created by: printf '%s\n' /lus/grand/projects/TFXcan/imlab/data/GEUVADIS/vcf_snps_only/*.gz > chr_vcfs.txt
#split -n ${NNODES} ${vcf_list} vcf_ --numeric-suffixes=1 --suffix-length=2
x=`wc -l ${vcf_list} | awk '{print $1}'`
y=${SLURM_NNODES} 
ll=$(( ($x + $y - 1) / $y ))

split --lines="${ll}" --numeric-suffixes=1 --suffix-length=1 "${vcf_list}" "vcf."

# prepare nodes
nodelist=$(scontrol show hostname $SLURM_NODELIST)
printf "%s\n" "${nodelist[@]}" > local_hostfiles.txt
# Increase value of suffix-length if more than 99 jobs
split --lines=1 --numeric-suffixes=1 --suffix-length=1 local_hostfiles.txt "local_hostfile."


for suf in `seq 1 ${SLURM_NNODES}`; do
  #IFS=,$'\n' read -d '' -r -a vcf_arr < "vcf_${suf}"
  (
    echo "INFO - Job ${suf} running on `cat local_hostfile.${suf}` using vcf.${suf}"
    srun --nodelist=`cat local_hostfile.${suf}` -n 1 -N 1 ${exec_file} ${rscript} ${utility_functions} "vcf.${suf}" "${plink_results_folder}" "${formatted_files_folder}"

  ) & sleep 1
done
wait

printf '\n%s\n' "INFO - merging all files into different snp and geno files."

${rscript} ${utility_functions} --command "merge_genotype_dosage_and_snp_annot_files" --files_folder ${formatted_files_folder} --output_folder ${formatted_files_folder} --combine 'yes'

# printf '\n%s\n' "INFO - combining snp and geno files into one."

# ${rscript} ${utility_functions} --command "create_merged_text_genotypes_for_prediction" --snp_annot_file ${formatted_files_folder}/all_chrs.snp_annot.txt --geno_file ${formatted_files_folder}/all_chrs.geno.txt --combine 'yes'

printf '\n%s\n' "INFO - Finished preparing all VCF files\n"