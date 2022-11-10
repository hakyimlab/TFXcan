# === This script contains all codes to test if tensorflow and ENFORMER are working accurately  
# Created by Temi
# DATE: Thursday Oct 20 2022

from __future__ import absolute_import, division, print_function, unicode_literals
import os, re, sys, json, h5py, csv, warnings
import parsl
from parsl.app.app import python_app
import pandas as pd
import numpy as np

def main():

    # get the path of the script as well as parameters         
    whereis_script = os.path.dirname(sys.argv[0])  
    script_path = os.path.abspath(whereis_script) 

    usage_codes = f'{script_path}/enformer-usage-codes.py'
    parsl_config = f'{script_path}/parsl-configuration.py'
    personal_enformer = f'{script_path}/personal-enformer.py'

    # import the enformer-usage_codes.py file
    exec(open(usage_codes).read(), globals(), globals())
    exec(open(parsl_config).read(), globals(), globals())
    exec(open(personal_enformer).read(), globals(), globals())

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
        sequence_folder = parameters['sequence_folder']

    for sam in individuals:

        # create the directories for this individual 
        if not os.path.exists(f'{output_dir}/{sam}'):
            print(f'\n[CREATING OUTPUT DIRECTORY] at {output_dir}/{sam}')
            os.makedirs(f'{output_dir}/{sam}')

        if not os.path.exists(f'{sequence_folder}/{sam}'):
            print(f'\n[CREATING SEQUENCE FOLDER] at {sequence_folder}/{sam}')
            os.makedirs(f'{sequence_folder}/{sam}')

        print(f'\n[LOADING INTERVALS] {sam}')

        # create the zip object for all the intervals
        a = pd.read_table(f'{intervals_dir}/{sam}_{TF}.txt', sep=' ', header=None)
        b = a[0].tolist() # a list of queries

        # I need a log file
        # read in the log file for this individual ; doing this so that the log file is not opened everytime
        logfile_csv = f'{logfile_path}/{sam}_predictions_log.csv'
        if os.path.isfile(logfile_csv):
            logfile = pd.read_csv(logfile_csv)
        else:
            logfile = None

        # for each interval check if prediction has been done
        q_status = [check_query(sample=sam, query=query, output_dir=output_dir, logfile=logfile) for q in b]
        for i, query in enumerate(b):
            query_status.append(check_query(sample=sam, query=query, output_dir=output_dir, logfile=logfile))
        
        q_status = [q.result() for q in query_status] # evaluate the results >> this should return a list of those that don't have predictions
        #q_status = [q[0] for q in q_status]
        print(f'Query results are: {q_status}n')

        # make predictions
        call_script = f"{script_path}/single-enformer-bashapp.sh"

        predict_status = [call_single_enformer_run(call_script=call_script, sequence_region=sreg, sam=sam) for sreg in q_status]

        prediction_finished = [p.result() for p in predict_status]

        print(f'Predict results are: {prediction_finished}')

        print(f'[FINISHED] {sam}\n')
    print(f'[FINISHED BOTH INDIVIDUALS]')                  
                    
if __name__ == '__main__':
    main()