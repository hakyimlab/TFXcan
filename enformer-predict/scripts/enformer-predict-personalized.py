# === This script contains all codes to test if tensorflow and ENFORMER are working accurately  
# Created by Temi
# DATE: Sunday Nov 13 2022

from __future__ import absolute_import, division, print_function, unicode_literals
import os, re, sys, json, h5py, csv, warnings
import parsl
from parsl.app.app import python_app
import pandas as pd
import numpy as np
import tensorflow as tf
import tensorflow_hub as hub # for interacting with saved models and tensorflow hub
import joblib
import gzip # for manipulating compressed files
from kipoiseq import Interval # same as above, really
import pyfaidx # to index our reference genome file
import pandas as pd # for manipulating dataframes
import numpy as np # for numerical computations
import os, sys, re # functions for interacting with the operating system
from functools import lru_cache
import math

from tqdm import tqdm

# get the path of the script as well as parameters         
whereis_script = os.path.dirname(sys.argv[0])  
script_path = os.path.abspath(whereis_script)

sys.path.append(f'{script_path}/utilities')

import enformerUsageCodes
import runPredictionUtilities
import parslConfiguration


def generate_batch(lst, batch_size):
    """  Yields bacth of specified size """
    if batch_size <= 0:
        return
    for i in range(0, len(lst), batch_size):
        yield lst[i:(i + batch_size)]

def main():

    # personal_enformer = f'{script_path}/personal-enformer.py'
    # exec(open(personal_enformer).read(), globals(), globals())

    # read the parameters file
    with open(f'{script_path}/../metadata/enformer_parameters.json') as f:

        parameters = json.load(f)

        intervals_dir = parameters['interval_list_dir']
        model_path = parameters['model_path']
        fasta_file = parameters['hg38_fasta_file']
        output_dir = parameters['output_dir']
        individuals = parameters['individuals']
        vcf_file = parameters['vcf_file']
        path_to_bcftools = parameters['path_to_bcftools']
        path_to_tabix = parameters['path_to_tabix']
        temporary_vcf_dir = parameters['temporary_vcf_dir']
        TF = parameters['TF']
        logfile_path = parameters['logfile_path']
        batch_size = int(parameters['batch_size'])
    
    # individuals can be a list or a txt file of individuals per row

    if isinstance(individuals, list):
        pass
    else:
        individuals = pd.read_table(individuals, header=None)[0].tolist()
        print(individuals)

    # load the parsl config file here since you want to distribute across individuals
    # parsl_config = f'{script_path}/utilities/parslConfiguration.py'
    # exec(open(parsl_config).read(), globals(), globals()) 

    config_directives = parslConfiguration.create_parsl_configuration()
    parsl.load(config_directives)

    for each_individual in individuals:

        @lru_cache(1)
        def get_fastaExtractor(fasta_file_path=fasta_file):
            fasta_extractor = enformerUsageCodes.FastaStringExtractor(fasta_file_path)
            return fasta_extractor
        
        @lru_cache(1)
        def get_model(model_class, model_path):
            return model_class(model_path)

        #fasta_extractor = get_fastaExtractor(fasta_file_path)
        # create the directories for this individual 
        if not os.path.exists(f'{output_dir}/{each_individual}'):
            print(f'[INFO] Creating output directory at {output_dir}/{each_individual}')
            os.makedirs(f'{output_dir}/{each_individual}')

        print(f'\n[INFO] Loading intervals for {each_individual}')
        a = pd.read_table(f'{intervals_dir}/{each_individual}_{TF}_400000.txt', sep=' ', header=None)
        list_of_regions = a[0].tolist()[1:100] # a list of queries

        # I need a log file
        # read in the log file for this individual ; doing this so that the log file is not opened everytime
        logfile_csv = f'{logfile_path}/{each_individual}_predictions_log.csv'
        if os.path.isfile(logfile_csv):
            logfile = pd.read_csv(logfile_csv)
        else:
            logfile = None

        predict_batches = generate_batch(list_of_regions, batch_size=batch_size)

        run_predictions_tools = f'{script_path}/utilities/runPredictionUtilities.py'
        exec(open(run_predictions_tools).read(), globals(), globals())

        # # load the parsl config file here since you want to distribute across individuals
        # parsl_config = f'{script_path}/utilities/parslConfiguration.py'
        # exec(open(parsl_config).read(), globals(), globals()) 

        count = 0
        for query_batch in predict_batches:

            # run_predictions_tools = f'{script_path}/utilities/runPredictionUtilities.py'
            # exec(open(run_predictions_tools).read(), globals(), globals())

            print(f'[INFO] Starting on batch {count + 1} of {math.ceil(len(list_of_regions)/batch_size)}')

            sequence_list_appfutures = (create_input(region=query, individual=each_individual, vcf_file=vcf_file, subset_vcf_dir=temporary_vcf_dir, fasta_file_path=fasta_file, script_path=script_path, fasta_func=get_fastaExtractor, output_dir=output_dir, logfile=logfile, software_paths=[path_to_bcftools, path_to_tabix]) for query in tqdm(query_batch, desc='[INFO] Creating futures for inputs'))

            sequence_list = [s.result() for s in sequence_list_appfutures]
            
            #filter sequence_list to remove nones
            sequence_list = [s for s in sequence_list if s is not None]

            if not sequence_list:
                print(f'[INFO] Nothing to do. All predictions are available for batch {count + 1} for {each_individual}')
            else:
                predictions_appfutures = (enformer_predict(sequence['sequence'], region=sequence['region'], sample=each_individual, seq_type=sequence['sequence_source'], model_path=model_path, model_func=get_model, output_dir=output_dir, script_path=script_path, logfile_path=logfile_path) for sequence in tqdm(sequence_list, desc='[INFO] Creating futures of predictions'))

                print(list(predictions_appfutures))
                #predictions_output = [p.result() for p in predictions_appfutures]
            count = count + 1
        print(f'[INFO] Finished predictions for {each_individual}')

        #parsl.clear()

    print(f'[INFO] Finished all predictions')                  
                    
if __name__ == '__main__':
    main()


 # # for each interval check if prediction has been done; will only return these ==> not a parsl app
        # query_status = (runPredictionUtilities.check_query(sample=each_individual, query=query, output_dir=output_dir, logfile=logfile) for query in tqdm(list_of_regions, desc='[INFO] Checking query'))

        # filter this query_status to remove Nones
        #query_status = (l for l in query_status if not isinstance(l, type(None)))
        #query_status_peeked = runPredictionUtilities.peek(query_status)

        # if query_status_peeked is None:
        #     print('[INFO] Nothing to do. All predictions are available.')
        # else:
            
        #query_status = query_status_peeked
        #print('[INFO] Something to do. Making predictions.')