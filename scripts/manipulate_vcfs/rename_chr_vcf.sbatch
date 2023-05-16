#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=rename_chr_vcf
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/rename_chr_vcf.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/rename_chr_vcf.err
#SBATCH --time=01:00:00	
#SBATCH --partition=caslake

slurm_workdir=${SLURM_SUBMIT_DIR}
SLURM_O_WORKDIR=${slurm_workdir}

mkdir -p ${SLURM_O_WORKDIR}
echo Working directory is $SLURM_O_WORKDIR
cd $SLURM_O_WORKDIR

echo Jobid: $SLURM_JOBID
echo Running on host `hostname`
#echo Running on nodes `cat $SLURM_JOB_NODELIST`

printf "Starting to run\n"
source ~/.bashrc
conda activate compbio-tools

input_vcf="/beagle3/haky/data/CWAS/merged/merged_cwas_genotypes.vcf.gz"
output_vcf="/beagle3/haky/data/CWAS/merged/merged_cwas_genotypes_chrRenamed.vcf.gz"
mtdata_dir="/beagle3/haky/users/temi/projects/TFXcan/metadata"
chr_list="${mtdata_dir}/num_chrom.txt"

bash /beagle3/haky/users/temi/projects/TFXcan/scripts/manipulate_vcfs/vcf_utils.sh rename_chromosomes ${input_vcf} ${output_vcf} ${chr_list}