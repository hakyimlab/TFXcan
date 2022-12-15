#!/bin/bash
#PBS -l walltime=00:59:00,filesystems=grand
#PBS -A covid-ct
#PBS -q debug-scaling    
#PBS -N find_motifs_genome
#PBS -k doe
#PBS -o /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/impact_pipeline/log/find_motifs_genome.out
#PBS -e /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/impact_pipeline/log/find_motifs_genome.err

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

project_dir="/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/impact_pipeline"
homer_dir="~/miniconda3/envs/homer-env/share/homer"
mpiexec="/opt/cray/pe/pals/1.1.7/bin/mpiexec"

${mpiexec} -n ${NRANKS} --ppn ${NRANKS} --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS="${NTHREADS}" perl ${homer_dir}/bin/findMotifsGenome.pl ${peaks} ${genome} ${homer_output} -size ${size} -find ${motif_file} > ${output}

status=$?
echo "Exit status of job is: $status"