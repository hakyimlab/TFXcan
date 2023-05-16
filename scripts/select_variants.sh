#!/bin/bash
#PBS -l select=10:system=polaris
#PBS -l walltime=02:00:00,filesystems=home:grand
#PBS -A covid-ct
#PBS -q preemptable   
#PBS -N select_SNPS
#PBS -k doe
#PBS -o /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/logs/select_SNPS.out
#PBS -e /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/logs/select_SNPS.err

echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

echo Jobid: $PBS_JOBID
echo Running on host `hostname`
echo Running on nodes `cat $PBS_NODEFILE`

NNODES=`wc -l < $PBS_NODEFILE`
# NRANKS=1
# NDEPTH=24
# NTHREADS=2

NUM_NODES_PER_MPI=1
NRANKS_PER_NODE=10
NDEPTH=8
NTHREADS=1

NTOTRANKS=$(( NUM_NODES_PER_MPI * NRANKS_PER_NODE ))

#NTOTRANKS=$(( NNODES * NRANKS ))
#NTOTRANKS=$(( NNODES * NRANKS_PER_NODE ))

echo "NUM_OF_NODES=${NNODES}  TOTAL_NUM_RANKS=${NTOTRANKS}  RANKS_PER_NODE=${NRANKS_PER_NODE}  THREADS_PER_RANK=${NTHREADS}"

printf "Starting to run\n"
source ~/.bashrc
conda activate compbio-tools

echo "PBS_JOBID = " $PBS_JOBID

reference_genome="/lus/grand/projects/covid-ct/imlab/data/hg_sequences/hg38/Homo_sapiens_assembly38.fasta"
output_folder="/lus/grand/projects/covid-ct/imlab/data/GEUVADIS/vcf_snps_only"
vcf_folder="/lus/grand/projects/covid-ct/imlab/data/GEUVADIS/VCFs_chr_annotated"

mpiexec="/opt/cray/pe/pals/1.1.7/bin/mpiexec"

mtdata_dir="/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/metadata"
chr_list="${mtdata_dir}/chromosomes.txt"

split -n 10 ${chr_list} chr_ --numeric-suffixes=1 --suffix-length=2

# Increase value of suffix-length if more than 99 jobs
split --lines=${NUM_NODES_PER_MPI} --numeric-suffixes=1 --suffix-length=2 $PBS_NODEFILE local_hostfile.

for suf in `echo {01..10}`; do
    #echo ${suf}
    IFS=,$'\n' read -d '' -r -a chr_arr < "${mtdata_dir}/chr_${suf}"
    (
        for j in "${chr_arr[@]}"; do
            echo ${suf}-chr${j}

            vcf_file="${vcf_folder}/ALL.chr${j}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.vcf.gz"

            ${mpiexec} -n 1 --ppn 1 --hostfile "local_hostfile.${suf}" --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS="${NTHREADS}" gatk SelectVariants -R "${reference_genome}" -V "${vcf_file}" --select-type-to-include SNP -O "${output_folder}/ALL.chr${j}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.vcf.gz" \
            & sleep 1
        done
        wait
    ) & sleep 1
done
wait

printf "Finished selecting all SNPs\n"

# for lh in local_hostfile*
#     IFS=,$'\n' read -d '' -r -a chr_arr < "${grid_parameters}"
#     for i in `seq 1 22` X; do

#         echo "chr${i}"
#         vcf_file="${vcf_folder}/ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"

#         ${mpiexec} -n ${NTOTRANKS} --ppn ${NRANKS_PER_NODE} --hostfile ${lh} --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS="${NTHREADS}" gatk SelectVariants -R "${reference_genome}" -V "${vcf_file}" --select-type-to-include SNP -O "${output_folder}/ALL.chr${i}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.vcf.gz" \
#         & sleep 1
#     done
# wait


# #project_dir="/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/impact_pipeline"
# #homer_dir="~/miniconda3/envs/homer-env/share/homer"
# mpiexec="/opt/cray/pe/pals/1.1.7/bin/mpiexec"
# perl=`which perl`
# homer_cmd=`ls ${homer_dir}/bin/scanMotifGenomeWide.pl`
# output_file="${output_basename}.txt"

# if [[ -f ${homer_cmd} ]]; then 
#     printf "${homer_cmd} exists."
# fi

# ${mpiexec} -n ${NRANKS} --ppn ${NRANKS} --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS="${NTHREADS}" ${perl} ${homer_cmd} ${motif_file} ${genome} > ${output_file}

# #sort -t $'\t' -k6,6rn ${output_file} > ${output_file}

# status=$?
# echo "Exit status of job is: $status"


#split --lines=1 --numeric-suffixes=1 --suffix-length=1 $PBS_NODEFILE local_hostfile.