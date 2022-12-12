# Usage: This script is used to predict on batches using ENFORMER on individuals' regions
# Author: Temi
# Date: Sunday Nov 13 2022

# "/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_predict/metadata/individuals.txt"

from __future__ import absolute_import, division, print_function, unicode_literals
import os, sys, json
import pandas as pd # for manipulating dataframes
import math, time, tqdm
import parsl
from parsl.app.app import python_app
from functools import lru_cache
import glob

whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)
script_path = os.path.abspath(whereis_script)
fpath = os.path.join(script_path, 'batch_utils')
sys.path.append(fpath)

def main():
    # get the path of the script as well as parameters         
    #print(sys.path)

    global use_parsl

    with open(f'{script_path}/../metadata/enformer_parameters.json') as f:

        parameters = json.load(f)
        intervals_dir = script_path + '/../' + parameters['interval_list_dir']
        output_dir = script_path + '/../' + parameters['output_dir'] + '/' + parameters['prediction_data_name'] + '/predictions_' + parameters['date']
        individuals = script_path + '/../' + parameters['individuals']
        vcf_file = script_path + '/../' + parameters['vcf_file']
        TF = parameters['TF']
        predictions_log_dir = script_path + '/../' + parameters['predictions_log_dir']
        log_dir = script_path + '/../' + parameters['log_dir']
        batch_size = int(parameters['batch_size'])
        use_parsl = True if parameters['use_parsl'] == 'true' else False
        n_regions = parameters["predict_on_n_regions"]
        parsl_parameters = parameters['parsl_parameters']

        if int(n_regions) == -1:
            predict_on_n_regions = None
        elif int(n_regions) > 0:
            predict_on_n_regions = (n_regions + 1) if isinstance(n_regions, int) else None

    # modify parsl parameters to add the working directory
    parsl_parameters['working_dir'] = f'{script_path}/../'

    if not os.path.isdir(predictions_log_dir):
        os.makedirs(predictions_log_dir)

    if not os.path.isdir(log_dir):
        os.makedirs(log_dir)

    if use_parsl == True:
        import parslConfiguration
        parsl.load(parslConfiguration.theta_htParslConfig(params=parsl_parameters))

    #individuals can be a given list or a txt file of individuals per row or a single string
    if isinstance(individuals, list):
        pass
    elif isinstance(individuals, type('str')):
        if os.path.isfile(individuals):
            individuals = pd.read_table(individuals, header=None)[0].tolist()
        else:
            individuals = [individuals]
        
    print(f'[INFO] Predicting for these individuals: {individuals}')
    
    predict_utils_one = f'{script_path}/batch_utils/predictUtils_one.py'
    exec(open(predict_utils_one).read(), globals(), globals())

    for each_individual in individuals:
        
        # this is specific to an individual but is cached per individual
        # I want to cache this but it is a bit tricky to do for now
        def make_cyvcf_object(vcf_file=vcf_file, sample=each_individual):
            import cyvcf2
            return(cyvcf2.cyvcf2.VCF(vcf_file, samples=sample))

        # create the directories for this individual 
        if not os.path.exists(f'{output_dir}/{each_individual}'):
            print(f'[INFO] Creating output directory at {output_dir}/{each_individual}')
            os.makedirs(f'{output_dir}/{each_individual}')

        interval_list_file = glob.glob(f'{intervals_dir}/{each_individual}_{TF}_*.txt')[0]
        a = pd.read_table(interval_list_file, sep=' ', header=None)
        list_of_regions = a[0].tolist()[0:(predict_on_n_regions)] # a list of queries

        # I need a log file
        # read in the log file for this individual ; doing this so that the log file is not opened everytime
        logfile_csv = f'{predictions_log_dir}/{each_individual}_predictions_log.csv'
        logfile = pd.read_csv(logfile_csv) if os.path.isfile(logfile_csv) else None

        tic_prediction = time.perf_counter() # as opposed to process_time

        batches = generate_batch(list_of_regions, batch_n=batch_size)
        count = 0
        app_futures = []
        for batch_query in tqdm.tqdm(batches, desc=f"[INFO] Creating futures for batch {count+1} of {batch_size}"):
            app_futures.append(run_batch_predictions(batch_regions=batch_query, batch_num = count+1, id=each_individual, vcf_func=make_cyvcf_object, script_path=script_path, output_dir=output_dir, logfile=logfile, predictions_log_dir=predictions_log_dir))

            count = count + 1

        if use_parsl == True:
            exec_futures = [q.result() for q in tqdm.tqdm(app_futures, desc=f'[INFO] Executing futures for all {len(app_futures)} batches')] #for q in tqdm.tqdm(query_futures, desc=f'[INFO] Executing futures for {len(query_futures)} input regions')]

        toc_prediction = time.perf_counter()

        print(f'[INFO] (time) to create inputs and predict on {len(list_of_regions)} queries is {toc_prediction - tic_prediction}')

        if use_parsl == True:
            print(f'[INFO] Finished predictions for {each_individual}: {exec_futures} ...\n')
        elif use_parsl == False:
            print(f'[INFO] Finished predictions for {each_individual}: {app_futures} ...\n')

    print(f'[INFO] Finished all predictions')                  
                    
if __name__ == '__main__':
    main()
