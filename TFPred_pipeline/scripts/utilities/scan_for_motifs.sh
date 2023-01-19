#!/bin/bash
#PBS -l select=1:system=polaris
#PBS -l walltime=00:59:00,filesystems=home:grand
#PBS -A covid-ct
#PBS -q preemptable   
#PBS -N scan_genomewide_motifs
#PBS -k doe
#PBS -o /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/logs/scan_genomewide_motifs.out
#PBS -e /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/logs/scan_genomewide_motifs.err

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

echo "NUM_OF_NODES=${NNODES}  TOTAL_NUM_RANKS=${NTOTRANKS}  RANKS_PER_NODE=${NRANKS}  THREADS_PER_RANK=${NTHREADS}"

printf "Starting to run\n"
source ~/.bashrc
conda activate r-env

echo "PBS_JOBID = " $PBS_JOBID

#project_dir="/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/impact_pipeline"
#homer_dir="~/miniconda3/envs/homer-env/share/homer"
mpiexec="/opt/cray/pe/pals/1.1.7/bin/mpiexec"
perl=`which perl`
homer_cmd=`ls ${homer_dir}/bin/scanMotifGenomeWide.pl`
output_file="${output_basename}.txt"

if [[ -f ${homer_cmd} ]]; then 
    printf "${homer_cmd} exists."
fi

${mpiexec} -n ${NRANKS} --ppn ${NRANKS} --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS="${NTHREADS}" ${perl} ${homer_cmd} ${motif_file} ${genome} > ${output_file}

#sort -t $'\t' -k6,6rn ${output_file} > ${output_file}

status=$?
echo "Exit status of job is: $status"