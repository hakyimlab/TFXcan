#!/bin/bash
#COBALT -n 1
#COBALT -t 1:00:00
#COBALT -A covid-ct
#COBALT --jobname=PREDICT-REFERENCE-KAWAKAMI
#COBALT -O /home/temi/imlab_folder/users/temi/projects/TFXcan/log/predict_kawakami_reference 
#COBALT --mode script
#COBALT --attrs filesystems=theta-fs0
#COBALT -q single-gpu


echo 'Starting to run\n'

module load cobalt/cobalt-gpu

# arguments are 
# motif_regions = sys.argv[1]: a file for the regions I want to predict on 
# TF = sys.argv[2] # a supplied transcription factor
# cell_line = sys.argv[3]
# experiment_details = sys.argv[4]

# first enter gpu node : run `ssh thetagpusn1`
# run: source ~/.bashrc
# second, request an hour allocation: run `qsub -I -n 1 -t 1:00:00 -q single-gpu -A covid-ct` and `~/.bashrc`

~/miniconda3/bin/python3 /home/temi/imlab_folder/users/temi/projects/TFXcan/scripts/enformer-predict-reference.py '/projects/covid-ct/imlab/users/temi/projects/TFXcan/data/train-test-val/kawakami_training_motif_regions.csv' 'FOXA1' 'LuCaP' 'kawakami-training'

# qsub   -i /home/temi/imlab_folder/users/temi/projects/TFXcan/scripts/submit-reference-predictions.sh -O 
# qsub-gpu ./projects/TFXcan/scripts/submit-reference-predictions.sh
# qsub-gpu ./projects/TFXcan/scripts/submit-reference-predictions.sh