#!/bin/bash
#PBS -l select=10:system=polaris
#PBS -l walltime=00:30:00,filesystems=grand
#PBS -A covid-ct
#PBS -q preemptable  
#PBS -N enformer_predict
#PBS -k doe
#PBS -o /lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFpred_pipeline/logs/enformer_predict.out
#PBS -e /lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFpred_pipeline/logs/enformer_predict.err

echo "Working directory is ${PBS_O_WORKDIR}"
cd "${PBS_O_WORKDIR}"

echo Jobid: "${PBS_JOBID}"
echo Running on host `hostname`
echo Running on nodes `cat $PBS_NODEFILE`

NNODES=`wc -l < $PBS_NODEFILE`
NRANKS=1
NDEPTH=24
NTHREADS=2

NTOTRANKS=$(( NNODES * NRANKS ))

echo "NUM_OF_NODES=${NNODES}  TOTAL_NUM_RANKS=${NTOTRANKS}  RANKS_PER_NODE=${NRANKS}  THREADS_PER_RANK=${NTHREADS}"
#NTOTRANKS=$(( NNODES * NRANKS_PER_NODE ))

printf "Starting to run\n"
source ~/.bashrc
conda activate dl-tools

echo "PBS_JOBID = " $PBS_JOBID

mpiexec -n 1 --ppn ${NRANKS} --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS="${NTHREADS}" python3 "${enformer_predict_py}" --param_config "${param_file}"

status=$?
echo "Exit status of training run is: $status"