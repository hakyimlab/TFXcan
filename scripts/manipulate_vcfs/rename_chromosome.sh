#!/bin/bash

#PBS -l select=10:system=polaris
#PBS -l walltime=02:00:00,filesystems=home:grand
#PBS -A covid-ct
#PBS -q preemptable   
#PBS -N rename_VCF
#PBS -k doe
#PBS -o /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/logs/rename_VCF.out
#PBS -e /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/logs/rename_VCF.err

echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

echo Jobid: $PBS_JOBID
echo Running on host `hostname`
echo Running on nodes `cat $PBS_NODEFILE`

NNODES=`wc -l < $PBS_NODEFILE`
NRANKS=1
NDEPTH=24
NTHREADS=2

#NTOTRANKS=$(( NNODES * NRANKS ))
NTOTRANKS=$(( NNODES * NRANKS_PER_NODE ))

#cat $PE_HOSTFILE | awk ’{print $1" slots=12"}’ > mf

echo "NUM_OF_NODES=${NNODES}  TOTAL_NUM_RANKS=${NTOTRANKS}  RANKS_PER_NODE=${NRANKS}  THREADS_PER_RANK=${NTHREADS}"

printf "Starting to run\n"
source ~/.bashrc
conda activate compbio-tools

echo "PBS_JOBID = " $PBS_JOBID

output_folder="/lus/grand/projects/covid-ct/imlab/data/GEUVADIS/VCFs_chr_annotated"
vcf_folder="/lus/grand/projects/covid-ct/imlab/data/GEUVADIS/VCFs"
chr_annots="/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/metadata/chromosomes_annotation.txt"

mpiexec="/opt/cray/pe/pals/1.1.7/bin/mpiexec"

for j in `seq 1 22` X; do 
    echo "chr${j}"
    vcf_file="${vcf_folder}/ALL.chr${j}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
    # ${mpiexec} -n ${NRANKS} --ppn ${NRANKS} --hostfile ${PBS_NODEFILE} --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS="${NTHREADS}" bcftools annotate --rename-chrs "${chr_annots}" "${vcf_file}" -Oz -o "${output_folder}/ALL.chr${j}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.vcf.gz" --threads 24 \
    # & sleep 1
    ${mpiexec} -n ${NRANKS} --ppn ${NRANKS} --hostfile ${PBS_NODEFILE} --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS="${NTHREADS}" tabix -p vcf "${output_folder}/ALL.chr${j}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.vcf.gz" \
    & sleep 1
done 
wait

# N=10

# (
# for j in `seq 1 22` X; do 
#     ((i=i%N)); ((i++==0)) && wait

#     echo "chr${j}"
#     vcf_file="${vcf_folder}/ALL.chr${j}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
#     bcftools annotate --rename-chrs "${chr_annots}" "${vcf_file}" -Oz -o "${output_folder}/ALL.chr${j}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.vcf.gz" --threads 32
#     tabix -p vcf "${output_folder}/ALL.chr${j}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.vcf.gz" &
# done
# wait
# )

printf "Finished renaming all vcf files\n"