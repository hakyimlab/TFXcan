#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=split_vcf_files
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/split_vcf_files.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/split_vcf_files.err
#SBATCH --time=02:00:00	
#SBATCH --partition=caslake

module load openmpi
module load parallel

slurm_workdir=${SLURM_SUBMIT_DIR}
SLURM_O_WORKDIR=${slurm_workdir}

mkdir -p ${SLURM_O_WORKDIR}
echo Working directory is $SLURM_O_WORKDIR
cd $SLURM_O_WORKDIR

echo Jobid: $SLURM_JOBID
echo Running on host `hostname`

printf "Starting to run\n"
source ~/.bashrc
conda activate compbio-tools

output_folder="/beagle3/haky/data/CWAS/split_vcf"
mkdir -p ${output_folder}
vcf_file="/beagle3/haky/data/CWAS/merged/merged_cwas_genotypes_chrRenamed.vcf.gz"

mtdata_dir="/beagle3/haky/users/temi/projects/TFXcan/metadata"
chr_list="${mtdata_dir}/chromosomes.txt"

# Instead of $SLURM_NTASKS, use $SLURM_NNODES to determine how
# many jobs should be run simultaneously.
## define srun directives
srun="srun --exclusive -N1 -n1 --cpus-per-task $SLURM_CPUS_PER_TASK"
parallel="parallel --delay 0.2 -j $SLURM_NNODES --joblog /beagle3/haky/users/temi/projects/TFXcan/logs/split_vcf_task.log"

echo ${srun}
echo ${parallel}

$parallel "$srun /beagle3/haky/users/temi/projects/TFXcan/scripts/manipulate_vcfs/vcf_utils.sh split_vcf_per_chromosome {1} ${vcf_file} ${output_folder}" ::: < ${chr_list}

printf "Finished separating VCF file by chromosomes\n"