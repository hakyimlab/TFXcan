#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=merge_vcf_files
#SBATCH --account=pi-haky		# Account to be charged
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/merge_vcf_files.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/merge_vcf_files.err	# Error file name
#SBATCH --time=02:00:00			# Job duration (wall time)
#SBATCH --partition=caslake

source ~/.bashrc
module load parallel
conda activate compbio-tools #/home/temi/miniconda3/envs/compbio-tools

vcfs_folder=/beagle3/haky/data/CWAS/vcf
vcf_files=$( ls ${vcfs_folder}/*.vcf.gz )
#echo ${vcf_files}

#parallel -j 16 bcftools index {} ::: ${vcf_files}
# create the index files files
#bcftools index /beagle3/haky/data/CWAS/vcf/*.vcf.gz
bcftools merge --output /beagle3/haky/data/CWAS/merged/merged_cwas_genotypes.vcf.gz --output-type z /beagle3/haky/data/CWAS/vcf/*.vcf.gz #`echo ${vcfs_folder}` ##