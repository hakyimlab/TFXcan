#!/bin/bash
#COBALT -n 128
#COBALT -t 03:00:00
#COBALT -A covid-ct
#COBALT -q default      
#COBALT --attrs filesystems=theta-fs0
#COBALT --jobname=create_model_data
#COBALT -o /projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/log/create_model_data.out
#COBALT -e /projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/log/create_model_data.err
#COBALT --debuglog /projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/log/create_model_data.cobalt

NNODES=${COBALT_JOBSIZE}
#NRANKS_PER_NODE=1
NTHREADS_PER_CORE=2
NDEPTH=2

#NTOTRANKS=$(( NNODES * NRANKS_PER_NODE ))

printf "Starting to run\n"
source ~/.bashrc
conda activate dl-tools

echo "COBALT_JOBID = " $COBALT_JOBID
echo "COBALT_JOBSIZE (nodes) =" $COBALT_JOBSIZE
echo "COBALT_PARTNAME = " $COBALT_PARTNAME

project_dir="/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline"
individuals_list="${project_dir}/metadata/individuals.txt"
IFS=,$'\n' read -d '' -r -a individuals_array < "${individuals_list}" 

for each_individual in "${individuals_array[@]}"; do
    aprun -n 1 -N 1 -d ${NDEPTH} -j ${NTHREADS_PER_CORE} -cc depth -e OMP_NUM_THREADS=2 python3 "${project_dir}/scripts/create_model_data.py" "${each_individual}" &
    sleep 1
done

wait

status=$?
echo "Exit status of training run is: $status"

# aprun -n ${NTOTRANKS} -N ${NRANKS_PER_NODE} -d ${NDEPTH} -j ${NTHREADS_PER_CORE} -cc depth -e OMP_NUM_THREADS=2 ~/miniconda3/envs/r-env/bin/Rscript /projects/covid-ct/imlab/users/temi/projects/running-knl/scripts/xgboost-minimal-2.R

# status=$?

# echo "Exit status of training run is: $status"
# exit $status