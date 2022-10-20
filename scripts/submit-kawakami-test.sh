#!/bin/bash
#COBALT -n 1
#COBALT -t 1:00:00
#COBALT -A covid-ct
#COBALT --jobname=PREDICT-REFERENCE-KAWAKAMI-TEST
#COBALT -O /home/temi/imlab_folder/users/temi/projects/TFXcan/log/predict_kawakami_test_reference 
#COBALT --mode script
#COBALT --attrs filesystems=theta-fs0
#COBALT -q single-gpu


echo 'Starting to run\n'
module load cobalt/cobalt-gpu

~/miniconda3/bin/python3 /home/temi/imlab_folder/users/temi/projects/TFXcan/scripts/enformer-predict-reference.py '/projects/covid-ct/imlab/users/temi/projects/TFXcan/data/train-test-val/kawakami_test_motif_regions.csv' 'FOXA1' 'LuCaP' 'kawakami-test'

# qsub   -i /home/temi/imlab_folder/users/temi/projects/TFXcan/scripts/submit-reference-predictions.sh -O 
# qsub-gpu ./projects/TFXcan/scripts/submit-reference-predictions.sh
# qsub-gpu ./projects/TFXcan/scripts/submit-reference-predictions.sh