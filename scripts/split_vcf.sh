#!/bin/bash
#PBS -l select=10:system=polaris
#PBS -l walltime=01:00:00,filesystems=home:grand
#PBS -A covid-ct
#PBS -q preemptable   
#PBS -N split_variants
#PBS -k doe
#PBS -o /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/logs/split_variants.out
#PBS -e /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/logs/split_variants.err


pbs_workdir=${PBS_O_WORKDIR}/${PBS_JOBNAME}
PBS_O_WORKDIR=${pbs_workdir}

mkdir ${PBS_O_WORKDIR}
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

echo Jobid: $PBS_JOBID
echo Running on host `hostname`
echo Running on nodes `cat $PBS_NODEFILE`

NNODES=`wc -l < $PBS_NODEFILE`
NUM_NODES_PER_MPI=1
NRANKS_PER_NODE=10
NDEPTH=8
NTHREADS=1

NTOTRANKS=$(( NUM_NODES_PER_MPI * NRANKS_PER_NODE ))

echo "NUM_OF_NODES=${NNODES}  TOTAL_NUM_RANKS=${NTOTRANKS}  RANKS_PER_NODE=${NRANKS_PER_NODE}  THREADS_PER_RANK=${NTHREADS}"

printf "Starting to run\n"
source ~/.bashrc
conda activate compbio-tools

echo "PBS_JOBID = " $PBS_JOBID

output_folder="/lus/grand/projects/covid-ct/imlab/data/freedman_vcfs/vcfs_2023-01-25"

vcf_file="/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/genotypes/prj6_genotypes/merged_phased_SNPs.vcf.gz"

mpiexec="/opt/cray/pe/pals/1.1.7/bin/mpiexec"

mtdata_dir="/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/metadata"
chr_list="${mtdata_dir}/chromosomes.txt"

split -n 10 ${chr_list} chr_ --numeric-suffixes=1 --suffix-length=2

# Increase value of suffix-length if more than 99 jobs
split --lines=${NUM_NODES_PER_MPI} --numeric-suffixes=1 --suffix-length=2 $PBS_NODEFILE local_hostfile.

for suf in `echo {01..10}`; do
    #echo ${suf}
    IFS=,$'\n' read -d '' -r -a chr_arr < "chr_${suf}"
    (
        for j in "${chr_arr[@]}"; do
            echo "${suf}(`cat local_hostfile.${suf}`)-chr${j}"

            ${mpiexec} -n 1 --ppn 1 --hostfile "local_hostfile.${suf}" --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS="${NTHREADS}" bcftools filter "${vcf_file}" -r "chr${j}" --output-type z --output "${output_folder}/freedman_chr${j}_shapeit4_SNPs_2023-01-25_GRCh38_phased.vcf.gz" && tabix -p vcf "${output_folder}/freedman_chr${j}_shapeit4_SNPs_2023-01-25_GRCh38_phased.vcf.gz" \
            & sleep 1
        done
        wait
    ) & sleep 1
done
wait

printf "Finished separating all chromosome from VCF file\n"