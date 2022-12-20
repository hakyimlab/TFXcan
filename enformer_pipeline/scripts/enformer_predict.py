# Usage: This script is used to predict on batches using ENFORMER on individuals' regions
# Author: Temi
# Date: Sunday Nov 13 2022

from __future__ import absolute_import, division, print_function, unicode_literals
import os, sys, json
import pandas as pd # for manipulating dataframes
import time, tqdm
import parsl
from parsl.app.app import python_app
from functools import lru_cache
import glob
from datetime import date

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
        TF = parameters['TF']
        predictions_log_dir = script_path + '/../' + parameters['predictions_log_dir']
        log_dir = script_path + '/../' + parameters['log_dir']
        batch_size = int(parameters['batch_size'])
        use_parsl = parameters['use_parsl']
        n_regions = parameters["predict_on_n_regions"]
        parsl_parameters = parameters['parsl_parameters']
        dataset_type = parameters['dataset_type']
        prediction_data_name = parameters['prediction_data_name']
        run_date = parameters['date'] if parameters['date'] is not None else date.today().strftime("%Y-%m-%d")


        output_dir = script_path + '/../' + parameters['output_dir'] + '/' + parameters['prediction_data_name'] + '/predictions_' + run_date

        if int(n_regions) == -1:
            predict_on_n_regions = None
        elif int(n_regions) > 0:
            predict_on_n_regions = (n_regions + 1) if isinstance(n_regions, int) else None


        individuals = script_path + '/../' + parameters['individuals'] if dataset_type == 'personalized' else None
        vcf_file = script_path + '/../' + parameters['vcf_file'] if dataset_type == 'personalized' else None

    # modify parsl parameters to add the working directory
    parsl_parameters['working_dir'] = f'{script_path}/../'

    predictions_log_dir = f'{predictions_log_dir}/{prediction_data_name}/predictions_log_{run_date}'
    if not os.path.isdir(predictions_log_dir):
        os.makedirs(predictions_log_dir)

    if not os.path.isdir(log_dir):
        os.makedirs(log_dir)

    if use_parsl == True:
        print(f'[INFO] Using parsl configuration: {use_parsl}')
        import parslConfiguration
        parsl.load(parslConfiguration.polaris_htParslConfig(params=parsl_parameters))

    predict_utils_one = f'{script_path}/batch_utils/predictUtils_one.py'
    exec(open(predict_utils_one).read(), globals(), globals())

    prediction_fxn = return_prediction_function(use_parsl)

    if dataset_type == 'personalized':
        if isinstance(individuals, list):
            id_list = individuals
            pass
        elif isinstance(individuals, type('str')):
            if os.path.isfile(individuals):
                id_list = pd.read_table(individuals, header=None)[0].tolist()
            else:
                id_list = [individuals]
        print(f'[INFO] Predicting for these individuals: {id_list}')
    elif dataset_type == 'reference':
        id_list = [prediction_data_name]
        print(f'[INFO] Predicting on a reference set')
        
    for each_id in id_list:
        
        # this is specific to an individual but is cached per individual
        # I want to cache this but it is a bit tricky to do for now
        #global make_cyvcf_object
        if dataset_type == 'personalized':
            def make_cyvcf_object(vcf_file=vcf_file, sample=each_id):
                import cyvcf2
                return(cyvcf2.cyvcf2.VCF(vcf_file, samples=sample))
        elif dataset_type == 'reference':
            make_cyvcf_object = None
            pass

        # create the directories for this individual 
        if not os.path.exists(f'{output_dir}/{each_id}'):
            print(f'[INFO] Creating output directory at {output_dir}/{each_id}')
            os.makedirs(f'{output_dir}/{each_id}')

        #print(f'{intervals_dir}/{each_id}_{TF}_*.txt')
        interval_list_file = glob.glob(f'{intervals_dir}/{each_id}_{TF}_*.txt')[0]
        a = pd.read_table(interval_list_file, sep=' ', header=None)
        list_of_regions = a[0].tolist()[0:(predict_on_n_regions)] # a list of queries

        # I need a log file
        # read in the log file for this individual ; doing this so that the log file is not opened everytime
        logfile_csv = f'{predictions_log_dir}/{each_id}_predictions_log.csv'
        logfile = pd.read_csv(logfile_csv) if os.path.isfile(logfile_csv) else None

        tic_prediction = time.perf_counter() # as opposed to process_time

        batches = generate_batch(list_of_regions, batch_n=batch_size)
        count = 0
        app_futures = []
        for batch_query in tqdm.tqdm(batches, desc=f"[INFO] Creating futures for batch {count+1} of {batch_size}"):
            app_futures.append(prediction_fxn(batch_regions=batch_query, batch_num = count+1, id=each_id, vcf_func=make_cyvcf_object, script_path=script_path, output_dir=output_dir, logfile=logfile, predictions_log_dir=predictions_log_dir, dataset_type=dataset_type))

            count = count + 1

        if use_parsl == True:
            exec_futures = [q.result() for q in tqdm.tqdm(app_futures, desc=f'[INFO] Executing futures for all {len(app_futures)} batches')] #for q in tqdm.tqdm(query_futures, desc=f'[INFO] Executing futures for {len(query_futures)} input regions')]

        toc_prediction = time.perf_counter()

        print(f'[INFO] (time) to create inputs and predict on {len(list_of_regions)} queries is {toc_prediction - tic_prediction}')

        if use_parsl == True:
            print(f'[INFO] Finished predictions for {each_id}: {exec_futures} ...\n')
        elif use_parsl == False:
            print(f'[INFO] Finished predictions for {each_id}: {app_futures} ...\n')

    print(f'[INFO] Finished all predictions')                  
                    
if __name__ == '__main__':
    main()
