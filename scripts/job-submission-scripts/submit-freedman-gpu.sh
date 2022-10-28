#!/bin/bash
#COBALT -n 1
#COBALT -t 12:00:00
#COBALT -A covid-ct
#COBALT --jobname=PREDICT-FREEDMAN
#COBALT -O /home/temi/imlab_folder/users/temi/projects/TFXcan/log/predict_freedman 
#COBALT --mode script
#COBALT --attrs filesystems=theta-fs0
#COBALT -q full-node

echo 'Starting to run\n'
module load cobalt/cobalt-gpu

source ~/.bashrc
conda activate dl-tools

python3 /home/temi/imlab_folder/users/temi/projects/TFXcan/scripts/enformer_predict.py 

# qsub   -i /home/temi/imlab_folder/users/temi/projects/TFXcan/scripts/submit-reference-predictions.sh -O 
# qsub-gpu ./projects/TFXcan/scripts/submit-reference-predictions.sh
# qsub-gpu ./projects/TFXcan/scripts/submit-reference-predictions.sh