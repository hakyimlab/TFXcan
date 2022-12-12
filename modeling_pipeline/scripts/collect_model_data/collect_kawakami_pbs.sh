#!/bin/bash
#PBS -l select=10:system=polaris
#PBS -l walltime=02:00:00,filesystems=grand
#PBS -A covid-ct
#PBS -q prod  
#PBS -N create_model_data_for_kawakami
#PBS -k doe
#PBS -o /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/log/create_model_data_for_kawakami.out
#PBS -e /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/log/create_model_data_for_kawakami.err

echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

echo Jobid: $PBS_JOBID
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

project_dir="/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline"
id="kawakami"

mpiexec -n 1 --ppn ${NRANKS} --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS="${NTHREADS}" python3 "${project_dir}/scripts/collect_model_data/collect_kawakami.py" "${id}"

status=$?
echo "Exit status of training run is: $status"

# aprun -n ${NTOTRANKS} -N ${NRANKS_PER_NODE} -d ${NDEPTH} -j ${NTHREADS_PER_CORE} -cc depth -e OMP_NUM_THREADS=2 ~/miniconda3/envs/r-env/bin/Rscript /projects/covid-ct/imlab/users/temi/projects/running-knl/scripts/xgboost-minimal-2.R

# status=$?

# echo "Exit status of training run is: $status"
# exit $status