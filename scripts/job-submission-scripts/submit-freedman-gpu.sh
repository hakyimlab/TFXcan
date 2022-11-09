#!/bin/bash
#COBALT -n 1
#COBALT -t 02:00:00
#COBALT -A covid-ct
#COBALT --jobname=PREDICT-FREEDMAN
#COBALT -o /projects/covid-ct/imlab/users/temi/projects/TFXcan/log/personalized_predict.out
#COBALT -e /projects/covid-ct/imlab/users/temi/projects/TFXcan/log/personalized_predict.err
#COBALT --debuglog /projects/covid-ct/imlab/users/temi/projects/TFXcan/log/personalized_predict.cobalt
#COBALT --mode script
#COBALT -q full-node

echo 'Starting to run\n'
#module load cobalt/cobalt-gpu

source ~/.bashrc
conda activate dl-tools
module load cobalt/cobalt-gpu

python3 /projects/covid-ct/imlab/users/temi/projects/TFXcan/scripts/enformer-predict-personalized.py 

# qsub   -i /home/temi/imlab_folder/users/temi/projects/TFXcan/scripts/submit-reference-predictions.sh -O 
# qsub-gpu ./projects/TFXcan/scripts/submit-reference-predictions.sh
# qsub-gpu ./projects/TFXcan/scripts/submit-reference-predictions.sh