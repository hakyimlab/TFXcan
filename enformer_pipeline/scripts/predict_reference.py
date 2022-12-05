# USAGE: This script is used to predict on a reference sequence in batches using ENFORMER 
# AUTHOR: Temi
# DATE: Sunday Dec 4 2022

from __future__ import absolute_import, division, print_function, unicode_literals
import os, sys, json
import pandas as pd # for manipulating dataframes
import time, tqdm
import parsl

whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)
script_path = os.path.abspath(whereis_script)
fpath = os.path.join(script_path, 'batch_utils')
sys.path.append(fpath)

def main():
    # get the path of the script as well as parameters         
    #print(sys.path)

    global use_parsl

    with open(f'{script_path}/../metadata/reference_parameters.json') as f:

        parameters = json.load(f)
        intervals_dir = script_path + '/../' + parameters['interval_list_dir']
        output_dir = script_path + '/../' + parameters['output_dir']
        id = parameters['id']
        TF = parameters['TF']
        predictions_log_dir = script_path + '/../' + parameters['predictions_log_dir']
        log_dir = script_path + '/../' + parameters['log_dir']
        batch_size = int(parameters['batch_size'])
        use_parsl = True if parameters['use_parsl'] == 'true' else False
        n_regions = parameters["predict_on_n_regions"]

        if int(n_regions) == -1:
            predict_on_n_regions = None
        elif int(n_regions) > 0:
            predict_on_n_regions = (n_regions) if isinstance(n_regions, int) else None
    
    predict_id = f'{id}_{TF}'

    if not os.path.isdir(predictions_log_dir):
        os.makedirs(predictions_log_dir)

    if not os.path.isdir(log_dir):
        os.makedirs(log_dir)

    if use_parsl == True:
        import parslConfiguration
        parsl.load(parslConfiguration.theta_htParslConfig(job_name='enformer-predict-reference', workingdir=None, num_full_nodes=2))

    predict_utils_one = f'{script_path}/reference_utils/predictUtils_one.py'
    exec(open(predict_utils_one).read(), globals(), globals())

    # create the directories for this individual 
    if not os.path.exists(f'{output_dir}/{predict_id}'):
        print(f'[INFO] Creating output directory at {output_dir}/{predict_id}')
        os.makedirs(f'{output_dir}/{predict_id}')

    a = pd.read_table(f'{intervals_dir}/{predict_id}_40000.txt', sep=' ', header=None)
    list_of_regions = a[0].tolist()[0:(predict_on_n_regions)] # a list of queries

    print(f'[INFO] Predicting on {len(list_of_regions)} regions.')

    # I need a log file
    # read in the log file for this individual ; doing this so that the log file is not opened everytime
    logfile_csv = f'{predictions_log_dir}/{predict_id}_predictions_log.csv'
    logfile = pd.read_csv(logfile_csv) if os.path.isfile(logfile_csv) else None

    tic_prediction = time.perf_counter() # as opposed to process_time

    batches = generate_batch(list_of_regions, batch_n=batch_size)
    count = 1
    app_futures = []
    for batch_query in tqdm.tqdm(batches, desc=f"[INFO] Creating futures for batch {count} of {batch_size}"):
        l_batch_query = list(batch_query)
        app_futures.append(run_batch_predictions(batch_regions=l_batch_query, batch_num = count+1, id=predict_id, script_path=script_path, output_dir=output_dir, logfile=logfile, predictions_log_dir=predictions_log_dir))

        print(f'[INFO] length of batch {count} is {len(l_batch_query)}')

        count = count + 1

    if use_parsl == True:
        exec_futures = [q.result() for q in tqdm.tqdm(app_futures, desc=f'[INFO] Executing futures for all {len(app_futures)} batches')] #for q in tqdm.tqdm(query_futures, desc=f'[INFO] Executing futures for {len(query_futures)} input regions')]

    toc_prediction = time.perf_counter()

    print(f'[INFO] (time) to create inputs and predict on {len(list_of_regions)} queries is {toc_prediction - tic_prediction}')

    if use_parsl == True:
        print(f'[INFO] Finished predictions for {predict_id}: {exec_futures} ...\n')
    elif use_parsl == False:
        print(f'[INFO] Finished predictions for {predict_id}: {app_futures} ...\n')

    print(f'[INFO] Finished all predictions')                  
                    
if __name__ == '__main__':
    main()