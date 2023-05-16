#!/bin/bash
#PBS -l select=1:system=polaris
#PBS -l walltime=04:59:00,filesystems=grand
#PBS -A covid-ct
#PBS -q preemptable    
#PBS -N train_enet_attempt_1
#PBS -k doe
#PBS -o /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/log/train_enet_attempt_1.out
#PBS -e /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/log/train_enet_attempt_2.err

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

project_dir="/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline"
#data_file=${1} # path to the data file
#"${project_dir}/data/enet_data/data_2022-12-13/balanced_enet.csv.gz"
#Rscript="~/miniconda3/envs/r-env/bin/Rscript"
mpiexec="/opt/cray/pe/pals/1.1.7/bin/mpiexec"

${mpiexec} -n ${NRANKS} --ppn ${NRANKS} --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS="${NTHREADS}" ~/miniconda3/envs/r-env/bin/Rscript "${project_dir}/scripts/enet/train_enet_model.R" "${data_file}" "${metainfo}"

status=$?
echo "Exit status of training run is: $status"

# qsub -v 'data_file=/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/data/train-test-val/kawakami/data_2022-12-12/kawakami_aggByCenter_FOXA1_old.csv.gz,metainfo=old' train_enet_model_pbs.sh