#!/bin/bash
#PBS -l select=10:system=polaris
#PBS -l walltime=03:00:00,filesystems=grand
#PBS -A covid-ct
#PBS -q prod    
#PBS -N train_xgboost
#PBS -k doe
#PBS -o /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/log/train_xgboost.out
#PBS -e /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/log/train_xgboost.err

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

project_dir="/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline"
mpiexec="/opt/cray/pe/pals/1.1.7/bin/mpiexec"
grid_parameters="${project_dir}/metadata/xgb_gridcv_parameters.txt"

IFS=,$'\n' read -d '' -r -a grid_array < "${grid_parameters}" 

for params in "${grid_array[@]}"; do
    printf "I am line: ${params}\n"
    ${mpiexec} -n 1 --ppn ${NRANKS} --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS="${NTHREADS}" ~/miniconda3/envs/r-env/bin/Rscript "${project_dir}/scripts/xgboost/cv_xgboost.R" ${params} \
    & sleep 1
done
wait

status=$?
echo "Exit status of training run is: $status"
exit $status