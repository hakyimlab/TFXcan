
import pandas as pd
from functools import lru_cache
import os, time, tqdm, math, json
import sys
import parsl

# if use_parsl == True:
#     import parslConfiguration
#     parsl.load(parslConfiguration.htParslConfig())

print(script_path)
print(region_range)

global output_dir, TF, individuals, intervals_dir, log_dir, batch_size, is_ref, logfile

with open(f'{script_path}/../metadata/enformer_parameters.json') as f:

    parameters = json.load(f)

    intervals_dir = parameters['interval_list_dir']
    output_dir = parameters['output_dir']
    individuals = parameters['individuals']
    vcf_file = parameters['vcf_file']
    TF = parameters['TF']
    log_dir = parameters['log_dir']
    batch_size = int(parameters['batch_size'])
    is_ref = parameters['is_ref']

if isinstance(individuals, list):
    pass
elif isinstance(individuals, type('str')):
    if os.path.isfile(individuals):
        individuals = pd.read_table(individuals, header=None)[0].tolist()[0:2]
    else:
        individuals = [individuals]
    
print(f'[INFO] Predicting for these individuals: {individuals}')

global each_individual, make_cyvcf_object
for each_individual in individuals:
    #print(each_individual + ' ===')

    run_predictions_tools = f'{script_path}/utilities/runPredictionUtilities.py'
    exec(open(run_predictions_tools).read(), globals(), globals())
    
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
    list_of_regions = a[0].tolist()[0:(region_range + 1)] # a list of queries

    # I need a log file
    # read in the log file for this individual ; doing this so that the log file is not opened everytime
    logfile_csv = f'{log_dir}/{each_individual}_predictions_log.csv'
    logfile = pd.read_csv(logfile_csv) if os.path.isfile(logfile_csv) else None

    tic_prediction = time.process_time() # as opposed to perf_counter

    batches = generate_batch(list_of_regions, batch_size=batch_size)
    count = 0
    for batch_query in batches:

        #enformer_model = get_model() # >> results in serialization errors

        query_futures = [run_single_predictions(region=query, individual=each_individual, vcf_func=make_cyvcf_object, script_path=script_path, output_dir=output_dir, logfile=logfile, logfile_path=log_dir) for query in tqdm.tqdm(batch_query, desc=f'[INFO] Creating futures for batch {count + 1} of {math.ceil(len(list_of_regions)/batch_size)}')]

        if use_parsl == True:
            query_exec = [q.result() for q in tqdm.tqdm(query_futures, desc=f'[INFO] Executing futures for {len(query_futures)} input regions')]

        count = count + 1

    toc_prediction = time.process_time()

    print(f'[INFO] (time) to create inputs and predict on {region_range} queries is {toc_prediction - tic_prediction}')

    if use_parsl == True:
        print(f'[INFO] Finished predictions for {each_individual}: {query_exec[0:11]} ...\n')
    else:
        print(f'[INFO] Finished predictions for {each_individual}: {query_futures[0:11]} ...\n')

    print(f'[INFO] Finished all predictions')   