# This script is used to predict sequentially using ENFORMER on individuals' regions
# AUTHOR: Temi
# DATE: Sunday Nov 13 2022

from __future__ import absolute_import, division, print_function, unicode_literals
import os, sys, json
import pandas as pd # for manipulating dataframes
from functools import lru_cache
import math, time, tqdm
import parsl

def main():

    # get the path of the script as well as parameters         
    whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)
    script_path = os.path.abspath(whereis_script)
    fpath = os.path.join(script_path, 'sequential_utils')
    sys.path.append(fpath)
    #print(sys.path)

    # load enformer parameters
    with open(f'{script_path}/../metadata/enformer_parameters.json') as f:

        parameters = json.load(f)
        intervals_dir = parameters['interval_list_dir']
        output_dir = parameters['output_dir']
        individuals = parameters['individuals']
        vcf_file = parameters['vcf_file']
        TF = parameters['TF']
        predictions_log_dir = parameters['predictions_log_dir']
        batch_size = int(parameters['batch_size'])
        use_parsl = True if parameters['use_parsl'] == 'true' else False
        n_regions = parameters["predict_on_n_regions"]

        if int(n_regions) == -1:
            predict_on_n_regions = None # i.e predict on all
        elif int(n_regions) > 0:
            predict_on_n_regions = (n_regions + 1) #if isinstance(n_regions, int) else None

    if use_parsl == True:
        import parslConfiguration
        parsl.load(parslConfiguration.htParslConfig())

    #individuals can be a given list or a txt file of individuals per row or a single string
    if isinstance(individuals, list):
        pass
    elif isinstance(individuals, type('str')):
        if os.path.isfile(individuals):
            individuals = pd.read_table(individuals, header=None)[0].tolist()[0:2]
        else:
            individuals = [individuals]
        
    print(f'[INFO] Predicting for these individuals: {individuals}')

    predict_utils_one = f'{script_path}/sequential_utils/predictUtils_one.py'
    exec(open(predict_utils_one).read(), globals(), globals())

    for each_individual in individuals:
        # this is specific to an individual but is cached per individual
        # I want to cache this but it is a bit tricky to do for now
        @lru_cache(1)
        def make_cyvcf_object(vcf_file=vcf_file, sample=each_individual):
            import cyvcf2
            return(cyvcf2.cyvcf2.VCF(vcf_file, samples=sample))

        # create the directories for this individual 
        if not os.path.exists(f'{output_dir}/{each_individual}'):
            print(f'[INFO] Creating output directory at {output_dir}/{each_individual}')
            os.makedirs(f'{output_dir}/{each_individual}')

        a = pd.read_table(f'{intervals_dir}/{each_individual}_{TF}_400000.txt', sep=' ', header=None)
        list_of_regions = a[0].tolist()[0:(predict_on_n_regions)] # a list of queries

        # I need a log file
        # read in the log file for this individual ; doing this so that the log file is not opened everytime
        logfile_csv = f'{predictions_log_dir}/{each_individual}_predictions_log.csv'
        logfile = pd.read_csv(logfile_csv) if os.path.isfile(logfile_csv) else None

        tic_prediction = time.process_time() # as opposed to perf_counter

        batches = generate_batch(list_of_regions, batch_size=batch_size)
        count = 0
        for batch_query in batches: # although I use batches here, predictions are not actually done in batches

            #enformer_model = get_model() # >> results in serialization errors

            query_futures = [run_single_predictions(region=query, individual=each_individual, vcf_func=make_cyvcf_object, script_path=script_path, output_dir=output_dir, logfile=logfile, predictions_log_dir=predictions_log_dir) for query in tqdm.tqdm(batch_query, desc=f'[INFO] Creating futures for batch {count + 1} of {math.ceil(len(list_of_regions)/batch_size)}')]

            if use_parsl == True:
                query_exec = [q.result() for q in tqdm.tqdm(query_futures, desc=f'[INFO] Executing futures for {len(query_futures)} input regions')]

            count = count + 1

        toc_prediction = time.process_time()

        print(f'[INFO] (time) to create inputs and predict on {predict_on_n_regions} queries is {toc_prediction - tic_prediction}')

        if use_parsl == True:
            print(f'[INFO] Finished predictions for {each_individual}: {query_exec} ...\n')
        else:
            print(f'[INFO] Finished predictions for {each_individual}: {query_futures[0:11]} ...\n')

    print(f'[INFO] Finished all predictions')                  
                    
if __name__ == '__main__':
    main()