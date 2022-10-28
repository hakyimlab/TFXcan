#!/bin/bash
#COBALT -n 1
#COBALT -t 04:00:00
#COBALT -A covid-ct
#COBALT --jobname=xgboost-cv
#COBALT -o /projects/covid-ct/imlab/users/temi/projects/TFXcan/log/xgboost-cv_log.out
#COBALT -e /projects/covid-ct/imlab/users/temi/projects/TFXcan/log/xgboost-cv_log.err
#COBALT --debuglog /projects/covid-ct/imlab/users/temi/projects/TFXcan/log/xgboost-cv_log.cobalt
#COBALT --mode script
#COBALT --attrs filesystems=theta-fs0
#COBALT -q full-node

echo 'Starting to run\n'

#bash /projects/covid-ct/imlab/users/temi/projects/TFXcan/scripts/call-xgboost-model.sh
# load conda
source ~/.bashrc

conda activate r-env
#module load cobalt/cobalt-gpu

~/miniconda3/envs/r-env/bin/Rscript /projects/covid-ct/imlab/users/temi/projects/TFXcan/scripts/xgboost-model.R

status=$?

echo "Exit status of training run is: $status"
exit $status