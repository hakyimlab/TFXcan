#!/bin/bash
#PBS -l walltime=00:59:00,filesystems=grand
#PBS -A covid-ct
#PBS -q debug-scaling    
#PBS -N create_model_data_for_freedman
#PBS -k doe
#PBS -o /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/log/create_model_data_for_freedman.out
#PBS -e /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/log/create_model_data_for_freedman.err

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
#echo "PBS_JOBSIZE (nodes) =" $PBS_JOBSIZE
#echo "PBS_PARTNAME = " $PBS_PARTNAME

project_dir="/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline"
individuals_list="${project_dir}/metadata/individuals.txt"
IFS=,$'\n' read -d '' -r -a individuals_array < "${individuals_list}" 

for each_individual in "${individuals_array[@]}"; do
    #printf "${each_individual}"
    mpiexec -n 1 --ppn ${NRANKS} --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS="${NTHREADS}" python3 "${project_dir}/scripts/collect_model_data/collect_freedman.py" "${each_individual}" &
    sleep 1
done

wait

status=$?
echo "Exit status of training run is: $status"

# aprun -n ${NTOTRANKS} -N ${NRANKS_PER_NODE} -d ${NDEPTH} -j ${NTHREADS_PER_CORE} -cc depth -e OMP_NUM_THREADS=2 ~/miniconda3/envs/r-env/bin/Rscript /projects/covid-ct/imlab/users/temi/projects/running-knl/scripts/xgboost-minimal-2.R

# status=$?

# echo "Exit status of training run is: $status"
# exit $statusls
