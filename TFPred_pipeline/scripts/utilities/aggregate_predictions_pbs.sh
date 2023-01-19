#!/bin/bash
#PBS -l select=10:system=polaris
#PBS -l walltime=06:00:00,filesystems=grand
#PBS -A covid-ct
#PBS -q preemptable  
#PBS -N aggregate_predictions
#PBS -k doe
#PBS -o /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/logs/aggregate_predictions.out
#PBS -e /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/logs/aggregate_predictions.err

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

source ~/.bashrc
conda activate dl-tools
#printf "Aggregation type: ${agg_type}\n"

python3 "${aggregate_py}" --metadata_file "${aggregation_config}" --agg_types "aggByCenter aggByPreCenter aggByPostCenter aggByUpstreamDownstream aggByDownstream aggByUpstream aggByMean"

# python3 "${collect_py}" --metadata_file "${aggregation_config}" --agg_types "${agg_types}"
# #mpiexec -n 1 --ppn 1 --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS="${NTHREADS}" python3 "${collect_py}" --metadata_file "${aggregation_config}" --agg_types "${agg_types}"

status=$?
echo "Exit status of aggregating step is: $status"


# == A way to run this script ===
#qsub -v prediction_path=/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_pipeline/enformer_predictions/kawakami/predictions_2023-01-05/kawakami,log_dir=/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_pipeline/predictions_log/kawakami/predictions_log_2023-01-05,id=kawakami,agg_type=aggByCenter /lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/scripts/collect_model_data/collect_kawakami_pbs.sh 

#qsub -v prediction_path=/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_pipeline/enformer_predictions/kawakami/predictions_2022-12-16/kawakami,log_dir=/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_pipeline/predictions_log/kawakami/predictions_log_2022-12-16,id=kawakami,agg_type=aggByMean /lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/scripts/collect_model_data/collect_kawakami_pbs.sh 

# aprun -n ${NTOTRANKS} -N ${NRANKS_PER_NODE} -d ${NDEPTH} -j ${NTHREADS_PER_CORE} -cc depth -e OMP_NUM_THREADS=2 ~/miniconda3/envs/r-env/bin/Rscript /projects/covid-ct/imlab/users/temi/projects/running-knl/scripts/xgboost-minimal-2.R

# status=$?

# echo "Exit status of training run is: $status"
# exit $status


# qsub -v collect_py=/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/scripts/utilities/aggregate_predictions.py,aggregation_config=/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/metadata/aggregation_config_freedman_FOXA1.json,agg_types="aggByCenter aggByPreCenter aggByPostCenter aggByUpstreamDownstream aggByDownstream aggByUpstream aggByMean" /lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/scripts/utilities/aggregate_predictions_pbs.sh

# python3 ../scripts/utilities/aggregate_predictions.py --metadata_file /lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/metadata/aggregation_config_freedman_FOXA1.json --agg_types aggByCenter aggByPreCenter aggByPostCenter