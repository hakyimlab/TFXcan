# This script is used to predict on batches using ENFORMER on individuals' regions
# AUTHOR: Temi
# DATE: Sunday Nov 13 2022

# "/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_predict/metadata/individuals.txt"

from __future__ import absolute_import, division, print_function, unicode_literals
import os, sys, json
import pandas as pd # for manipulating dataframes
import math, time, tqdm
import parsl
from parsl.app.app import python_app

whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)
script_path = os.path.abspath(whereis_script)
fpath = os.path.join(script_path, 'utilities')
sys.path.append(fpath)

@python_app
def run_batch_predictions(batch_regions, batch_num, individual, vcf_func, script_path, output_dir, logfile, predictions_log_dir): #
  
    import sys, os, tqdm, faulthandler
    mpath = os.path.join(script_path, 'utilities') #os.path.dirname(__file__) #
    sys.path.append(mpath)
    
    faulthandler.enable() # to figure out where segmentation error is coming from

    try:
        import predictUtils_two
    except ModuleNotFoundError as merr:
        print(f'[ERROR] of type {type(merr)} at run_batch_predictions')

    # first check the queries
    check_result = (predictUtils_two.check_query(sample = individual, query = each_region, output_dir=output_dir, logfile=logfile) for each_region in tqdm.tqdm(batch_regions, desc=f'[INFO] Checking query for batch'))
    # filter for nones
    filtered_check_result = [r for r in check_result if r is not None]

    if not filtered_check_result: # i.e. if the list is empty
        return(1)
    else:
        reg_prediction = predictUtils_two.enformer_predict(batch_region=filtered_check_result, sample=individual, model=None, output_dir=output_dir, vcf_func=vcf_func, predictions_log_dir=predictions_log_dir, batch_num=batch_num)
        return(reg_prediction) # returns 0 returned by enformer_predict

# else:
#     print(f"[WARNING] {check_result}: Either length of input sequence is invalid (NoneType) or too long or too short")
            

def generate_batch(lst, batch_size):
    """  
    Given a list, this function yields batches of a specified size
    
    Parameters:
        lst: list
        batch_size: int
            Number of items in each batch.

    Yields
        Batches of the list containing `batch_size` elements.
    """
    if batch_size <= 0:
        return None
    for i in range(0, len(lst), batch_size):
        yield lst[i:(i + batch_size)]

def main():
    # get the path of the script as well as parameters         
    #print(sys.path)

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
            predict_on_n_regions = None
        elif int(n_regions) > 0:
            predict_on_n_regions = (n_regions + 1) if isinstance(n_regions, int) else None

    if use_parsl == True:
        import parslConfiguration
        parsl.load(parslConfiguration.htParslConfig())

    #individuals can be a given list or a txt file of individuals per row or a single string
    if isinstance(individuals, list):
        pass
    elif isinstance(individuals, type('str')):
        if os.path.isfile(individuals):
            individuals = pd.read_table(individuals, header=None)[0].tolist()
        else:
            individuals = [individuals]
        
    print(f'[INFO] Predicting for these individuals: {individuals}')

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

        a = pd.read_table(f'{intervals_dir}/{each_individual}_{TF}_400000.txt', sep=' ', header=None)
        list_of_regions = a[0].tolist()[0:(predict_on_n_regions)] # a list of queries

        # I need a log file
        # read in the log file for this individual ; doing this so that the log file is not opened everytime
        logfile_csv = f'{predictions_log_dir}/{each_individual}_predictions_log.csv'
        logfile = pd.read_csv(logfile_csv) if os.path.isfile(logfile_csv) else None

        tic_prediction = time.perf_counter() # as opposed to process_time

        batches = generate_batch(list_of_regions, batch_size=batch_size)
        count = 0
        app_futures = []
        for batch_query in tqdm.tqdm(batches, desc=f"[INFO] Creating futures for batch {count+1} of {math.ceil(len(list_of_regions)/batch_size)}"):
            app_futures.append(run_batch_predictions(batch_regions=batch_query, batch_num = count+1, individual=each_individual, vcf_func=make_cyvcf_object, script_path=script_path, output_dir=output_dir, logfile=logfile, predictions_log_dir=predictions_log_dir))

            count = count + 1

        if use_parsl == True:
            exec_futures = [q.result() for q in tqdm.tqdm(app_futures, desc=f'[INFO] Executing futures for all {len(app_futures)} batches')] #for q in tqdm.tqdm(query_futures, desc=f'[INFO] Executing futures for {len(query_futures)} input regions')]

        toc_prediction = time.perf_counter()

        print(f'[INFO] (time) to create inputs and predict on {predict_on_n_regions} queries is {toc_prediction - tic_prediction}')

        if use_parsl == True:
            print(f'[INFO] Finished predictions for {each_individual}: {exec_futures} ...\n')
        elif use_parsl == False:
            print(f'[INFO] Finished predictions for {each_individual}: {app_futures} ...\n')

    print(f'[INFO] Finished all predictions')                  
                    
if __name__ == '__main__':
    main()
