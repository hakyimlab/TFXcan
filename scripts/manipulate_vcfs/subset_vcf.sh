#!/bin/bash
#PBS -l select=10:system=polaris
#PBS -l walltime=01:00:00,filesystems=home:grand
#PBS -A covid-ct
#PBS -q preemptable   
#PBS -N select_individuals_vcf
#PBS -k doe
#PBS -o /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/geuvadis_individuals/log/select_individuals_vcf.out
#PBS -e /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/geuvadis_individuals/log/select_individuals_vcf.err

mkdir $PBS_O_WORKDIR/rundir
PBS_O_WORKDIR=$PBS_O_WORKDIR/rundir
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

#vcf_folder="/lus/grand/projects/covid-ct/imlab/data/GEUVADIS/VCFs_chr_annotated"

mpiexec="/opt/cray/pe/pals/1.1.7/bin/mpiexec"

mtdata_dir="/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/metadata"
chr_list="${mtdata_dir}/chromosomes.txt"

split -n 10 ${chr_list} chr_ --numeric-suffixes=1 --suffix-length=2

# Increase value of suffix-length if more than 99 jobs
split --lines=${NUM_NODES_PER_MPI} --numeric-suffixes=1 --suffix-length=2 $PBS_NODEFILE local_hostfile.

echo "vcf dir is ${vcf_directory}"
echo "output dir is ${output_directory}"
echo "individuals file dir is ${individuals_file}"


if [[ ! -f "${individuals_file}" ]]; then
    echo "${individuals_file} does not exist" && exit 1
fi

if [[ ! -d "${output_directory}" ]]; then
    echo "${output_directory} does not exist" && exit 1
fi

if [[ ! -d "${vcf_directory}" ]]; then
    echo "${vcf_directory} does not exist" && exit 1
fi


for suf in `echo {01..10}`; do
    #echo ${suf}
    IFS=,$'\n' read -d '' -r -a chr_arr < "${PBS_O_WORKDIR}/chr_${suf}"
    (
        for j in "${chr_arr[@]}"; do
            echo ${suf}-chr${j}

            vcf_file="${vcf_directory}/ALL.chr${j}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.vcf.gz"
            if [[ ! -f "${vcf_file}" ]]; then
                echo "${vcf_file} does not exist" 
            else
                echo "${vcf_file} exists. All good."
                ${mpiexec} -n 1 --ppn 1 --hostfile "${PBS_O_WORKDIR}/local_hostfile.${suf}" --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS="${NTHREADS}" bcftools view --output-type z --output-file "${output_directory}/chr${j}_99_individuals.vcf.gz" --regions "chr${j}" --samples-file "${individuals_file}" "${vcf_file}" \
                & sleep 1
            fi
        done
        wait
    ) & sleep 1
done
wait

rm `echo local_hostfile.{01..10}`
rm `echo chr_{01..10}`

printf "Finished subsetting all individuals\n"

#qsub -v vcf_directory=/lus/grand/projects/covid-ct/imlab/data/GEUVADIS/vcf_snps_only,individuals_file=/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/metadata/99_individuals.txt,output_directory=/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/geuvadis_individuals/vcf_files ./subset_vcf.sh