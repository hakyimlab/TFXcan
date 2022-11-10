

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

import createSequencesList
from enformerUsageCodes import *

# read the parameters file
with open(f'{script_path}/../metadata/enformer_parameters.json') as f:

    parameters = json.load(f)

    intervals_dir = parameters['interval_list_dir']
    fasta_file = parameters['hg38_fasta_file']
    individuals = parameters['individuals']
    vcf_file = parameters['vcf_file']
    path_to_bcftools = parameters['path_to_bcftools']
    path_to_tabix = parameters['path_to_tabix']
    temporary_vcf_dir = parameters['temporary_vcf_dir']
    TF = parameters['TF']
    sequence_folder = parameters['sequence_folder']

#enf_model = tf.saved_model.load(model_path).model
sequence_extractor = FastaStringExtractor(fasta_file)

if not os.path.exists(sequence_folder):
    os.makedirs(sequence_folder)

for each_individual in individuals:
    # read the intervals file
    # create the zip object for all the intervals
    a = pd.read_table(f'{intervals_dir}/{each_individual}_{TF}.txt', sep=' ', header=None)
    each_individual_intervals = a[0].tolist() # a list of queries


    with open(f'{sequence_folder}/{each_individual}_regions_sequences.txt', 'w') as f:
        
        for each_region in each_individual_intervals:

            qsplit = each_region.split('_')
            region = [qsplit[0], int(qsplit[1]), int(qsplit[2]), each_region]

            output = createSequencesList.create_sequences(sample=each_individual, region=region, fasta_file_path=fasta_file, fasta_extractor=sequence_extractor, open_vcf_file=vcf_file, temporary_vcf_dir=temporary_vcf_dir, software_paths=[path_to_bcftools, path_to_tabix], script_path=script_path)

        f.write(output['sequence'][each_individual] + '\n')