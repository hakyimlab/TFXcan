#!/bin/bash
#COBALT -n 1
#COBALT -t 02:00:00
#COBALT -A covid-ct
#COBALT --jobname=PREDICT-FREEDMAN
#COBALT -O /projects/covid-ct/imlab/users/temi/projects/TFXcan/log 
#COBALT --mode script
#COBALT --attrs filesystems=theta-fs0
#COBALT -q full-node

echo 'Starting to run\n'
module load cobalt/cobalt-gpu

source ~/.bashrc
conda activate dl-tools

python3 /projects/covid-ct/imlab/users/temi/projects/TFXcan/scripts/enformer-predict-reference.py 

# qsub   -i /home/temi/imlab_folder/users/temi/projects/TFXcan/scripts/submit-reference-predictions.sh -O 
# qsub-gpu ./projects/TFXcan/scripts/submit-reference-predictions.sh
# qsub-gpu ./projects/TFXcan/scripts/submit-reference-predictions.sh