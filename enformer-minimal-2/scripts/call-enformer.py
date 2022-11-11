import sys, os, json

# get the path of the script as well as parameters         
whereis_script = os.path.dirname(sys.argv[0])  
script_path = os.path.abspath(whereis_script) 

sys.path.append(f'{script_path}/utilities')

import kipoiseq # same as above, really
import subprocess
import pandas as pd
import kipoiseq 
import warnings

from enformerUsageCodes import *
import personalEnformerUtilities

personal_enformer = f'{script_path}/utilities/personalEnformerUtilities.py'
exec(open(personal_enformer).read(), globals(), globals())

with open(f'{script_path}/../metadata/enformer_parameters.json') as f:

    parameters = json.load(f)

    model_path = parameters['model_path']
    intervals_dir = parameters['interval_list_dir']
    fasta_file = parameters['hg38_fasta_file']
    output_dir = parameters['output_dir']
    individuals = parameters['individuals']
    vcf_file = parameters['vcf_file']
    path_to_bcftools = parameters['path_to_bcftools']
    path_to_tabix = parameters['path_to_tabix']
    temporary_vcf_dir = parameters['temporary_vcf_dir']
    TF = parameters['TF']
    sequence_folder = parameters['sequence_folder']

# # load the model
# enformer_model = Enformer(model_path)

parsl_config = f'{script_path}/utilities/parslConfiguration.py'
exec(open(parsl_config).read(), globals(), globals())

for each_individual in individuals:

    print(sys.path)

    # create the directories for this individual 
    if not os.path.exists(f'{output_dir}/{each_individual}'):
        print(f'\n[CREATING OUTPUT DIRECTORY] at {output_dir}/{each_individual}')
        os.makedirs(f'{output_dir}/{each_individual}')

    #open the sequences files
    with open(f'{sequence_folder}/{each_individual}_regions_sequences.txt', 'r') as f:

        sequence_list = [sequence.rstrip() for sequence in f]

    sequence_list_predictions_run = []
    for i, each_sequence in enumerate(sequence_list):
        sequence_list_predictions_run.append(run_predictions(sequence=each_sequence, region=str(i), sample=each_individual, seq_type='var', model_path=model_path, output_dir=output_dir))
    
    prediction_output = [s.result() for s in sequence_list_predictions_run]

    print(f'[DONE] prediction_output')