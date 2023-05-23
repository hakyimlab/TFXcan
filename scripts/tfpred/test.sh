#!/bin/bash

function nodefinder (){
    printf "\n Argument ${1} is running on host `hostname`"
    printf "\nDone"
}

export -f nodefinder

"${@}"



#!/bin/bash

# echo Working directory is $PBS_O_WORKDIR
# cd $PBS_O_WORKDIR

# echo Jobid: $PBS_JOBID
# echo Running on host `hostname`
# echo Running on nodes `cat $PBS_NODEFILE`

# NNODES=`wc -l < $PBS_NODEFILE`
# NRANKS=1
# NDEPTH=48
# NTHREADS=1

# #NTOTRANKS=$(( NNODES * NRANKS ))
# NTOTRANKS=$(( NNODES * NRANKS_PER_NODE ))

# echo "NUM_OF_NODES=${NNODES}  TOTAL_NUM_RANKS=${NTOTRANKS}  RANKS_PER_NODE=${NRANKS}  THREADS_PER_RANK=${NTHREADS}"

# echo "PBS_JOBID = " $PBS_JOBID

# echo "NUM_OF_NODES=${NNODES}  TOTAL_NUM_RANKS=${NTOTRANKS}  RANKS_PER_NODE=${NRANKS}  THREADS_PER_RANK=${NTHREADS}"

# mpiexec="/opt/cray/pe/pals/1.1.7/bin/mpiexec"

# args=( `seq 1 5` )
# printf "Arguments are: ${args[@]} \n"

# module load gnu-parallel


# # ${mpiexec} --hostfile ${PBS_NODEFILE} -n ${NRANKS} --ppn ${NRANKS} --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS="${NTHREADS}" '/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/scripts/tfpred/test.sh' nodefinder 

# #parallel="parallel -j 1 --joblog runtask.log"

# # $parallel "${mpiexec} --hostfile ${PBS_NODEFILE} -np `cat ${PBS_NODEFILE} | wc -l` --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS="${NTHREADS}" '/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/scripts/tfpred/test.sh' nodefinder " ::: ${args[@]}


# #echo ${args[@]} | xargs -I {}  #-I {} ' echo {} '

# #printf '%s\n' "${args[@]}" | xargs

# NN=2 #`cat ${PBS_NODEFILE} | wc -l`

# printf '%s\n' "${args[@]}" | xargs -I {} ${mpiexec} --hostfile ${PBS_NODEFILE} -n ${NN} --ppn 1 --bynode --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS=${NTHREADS} /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/scripts/tfpred/test.sh nodefinder {}

#echo ${args[@]} | xargs -I {} ${mpiexec} --hostfile ${PBS_NODEFILE} -np `cat ${PBS_NODEFILE} | wc -l` --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS="${NTHREADS}" '/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/scripts/tfpred/test.sh' nodefinder {}"