#!/bin/bash
#PBS -l select=1:system=polaris
#PBS -l walltime=01:00:00,filesystems=grand
#PBS -A TFXcan
#PBS -q debug    
#PBS -N liftover
#PBS -k doe
#PBS -o /lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/baca_cwas/log/liftover.out
#PBS -e /lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/baca_cwas/log/liftover.err

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

echo "PBS_JOBID = " $PBS_JOBID

mpiexec="/opt/cray/pe/pals/1.1.7/bin/mpiexec"

${mpiexec} -n ${NRANKS} --ppn ${NRANKS} --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS="${NTHREADS}" "${liftover_cmd}" "${input_bed}" "${chain_file}" "${output_bed}" "${unmapped_bed}"

status=$?
echo "Exit status of lifting over is: $status"

# qsub -v 'data_file=/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/data/train-test-val/kawakami/data_2022-12-12/kawakami_aggByCenter_FOXA1_old.csv.gz,metainfo=old' train_enet_model_pbs.sh