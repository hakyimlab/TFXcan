# This script is used to predict using ENFORMER on individuals' regions
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
    fpath = os.path.join(script_path, 'utilities')
    sys.path.append(fpath)
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
        predict_on_n_regions = (n_regions + 1) if isinstance(n_regions, int) else None

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

    for each_individual in individuals:

        predict_utils_one = f'{script_path}/utilities/predictUtils_one.py'
        exec(open(predict_utils_one).read(), globals(), globals())
        
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
        for batch_query in batches:

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


   #for query in tqdm(list_of_regions, desc='[INFO] Creating futures for input regions'):

            # run_predictions_tools = f'{script_path}/utilities/predictUtils_two.py'
            # exec(open(run_predictions_tools).read(), globals(), globals())

            #print(f'[INFO] Starting on batch {count + 1} of {math.ceil(len(list_of_regions)/batch_size)}')

            #tic_prediction = time.process_time()

            #query_parsl_appfutures.append(run_single_predictions(region=query, individual=each_individual, vcf_file=vcf_file, subset_vcf_dir=temporary_vcf_dir, fasta_file_path=fasta_file, script_path=script_path, fasta_func=get_fastaExtractor, output_dir=output_dir, logfile=logfile, model_path=model_path, model_func=get_model, logfile_path=logfile_path, software_paths=[path_to_bcftools, path_to_tabix]))

            #toc_prediction = time.process_time()

            #print(f'[INFO] (time) to predict on a single query is {toc_prediction - tic_prediction}')
# def main():

#     # personal_enformer = f'{script_path}/personal-enformer.py'
#     # exec(open(personal_enformer).read(), globals(), globals())

#     # read the parameters file
#     with open(f'{script_path}/../metadata/enformer_parameters.json') as f:

#         parameters = json.load(f)

#         intervals_dir = parameters['interval_list_dir']
#         model_path = parameters['model_path']
#         fasta_file = parameters['hg38_fasta_file']
#         output_dir = parameters['output_dir']
#         individuals = parameters['individuals']
#         vcf_file = parameters['vcf_file']
#         path_to_bcftools = parameters['path_to_bcftools']
#         path_to_tabix = parameters['path_to_tabix']
#         temporary_vcf_dir = parameters['temporary_vcf_dir']
#         TF = parameters['TF']
#         logfile_path = parameters['logfile_path']
#         batch_size = int(parameters['batch_size'])
    
#     # individuals can be a list or a txt file of individuals per row

#     if isinstance(individuals, list):
#         pass
#     else:
#         individuals = pd.read_table(individuals, header=None)[0].tolist()
#         print(individuals)

#     # load the parsl config file here since you want to distribute across individuals
#     # parsl_config = f'{script_path}/utilities/parslConfiguration.py'
#     # exec(open(parsl_config).read(), globals(), globals()) 

#     config_directives = parslConfiguration.create_parsl_configuration()
#     parsl.load(config_directives)

#     for each_individual in individuals:

#         @lru_cache(1)
#         def get_fastaExtractor(fasta_file_path=fasta_file):
#             fasta_extractor = enformerUsageCodes.FastaStringExtractor(fasta_file_path)
#             return fasta_extractor
        
#         @lru_cache(1)
#         def get_model(model_class, model_path):
#             return model_class(model_path)

#         #fasta_extractor = get_fastaExtractor(fasta_file_path)
#         # create the directories for this individual 
#         if not os.path.exists(f'{output_dir}/{each_individual}'):
#             print(f'[INFO] Creating output directory at {output_dir}/{each_individual}')
#             os.makedirs(f'{output_dir}/{each_individual}')

#         print(f'\n[INFO] Loading intervals for {each_individual}')
#         a = pd.read_table(f'{intervals_dir}/{each_individual}_{TF}_400000.txt', sep=' ', header=None)
#         list_of_regions = a[0].tolist()[1:100] # a list of queries

#         # I need a log file
#         # read in the log file for this individual ; doing this so that the log file is not opened everytime
#         logfile_csv = f'{logfile_path}/{each_individual}_predictions_log.csv'
#         if os.path.isfile(logfile_csv):
#             logfile = pd.read_csv(logfile_csv)
#         else:
#             logfile = None

#         predict_batches = generate_batch(list_of_regions, batch_size=batch_size)

#         run_predictions_tools = f'{script_path}/utilities/predictUtils_two.py'
#         exec(open(run_predictions_tools).read(), globals(), globals())

#         # # load the parsl config file here since you want to distribute across individuals
#         # parsl_config = f'{script_path}/utilities/parslConfiguration.py'
#         # exec(open(parsl_config).read(), globals(), globals()) 

#         count = 0
#         for query_batch in predict_batches:

#             # run_predictions_tools = f'{script_path}/utilities/predictUtils_two.py'
#             # exec(open(run_predictions_tools).read(), globals(), globals())

#             print(f'[INFO] Starting on batch {count + 1} of {math.ceil(len(list_of_regions)/batch_size)}')

#             sequence_list_appfutures = (create_input(region=query, individual=each_individual, vcf_file=vcf_file, subset_vcf_dir=temporary_vcf_dir, fasta_file_path=fasta_file, script_path=script_path, fasta_func=get_fastaExtractor, output_dir=output_dir, logfile=logfile, software_paths=[path_to_bcftools, path_to_tabix]) for query in tqdm(query_batch, desc='[INFO] Creating futures for inputs'))

#             sequence_list = [s.result() for s in sequence_list_appfutures]
            
#             #filter sequence_list to remove nones
#             sequence_list = [s for s in sequence_list if s is not None]

#             if not sequence_list:
#                 print(f'[INFO] Nothing to do. All predictions are available for batch {count + 1} for {each_individual}')
#             else:
#                 predictions_appfutures = (enformer_predict(sequence['sequence'], region=sequence['region'], sample=each_individual, seq_type=sequence['sequence_source'], model_path=model_path, model_func=get_model, output_dir=output_dir, script_path=script_path, logfile_path=logfile_path) for sequence in tqdm(sequence_list, desc='[INFO] Creating futures of predictions'))

#                 #print(list(predictions_appfutures))
#                 predictions_output = [p.result() for p in predictions_appfutures]
#             count = count + 1
#         print(f'[INFO] Finished predictions for {each_individual}')

#         #parsl.clear()

#     print(f'[INFO] Finished all predictions')                  
                    
# if __name__ == '__main__':
#     main()


 # # for each interval check if prediction has been done; will only return these ==> not a parsl app
        # query_status = (predictUtils_two.check_query(sample=each_individual, query=query, output_dir=output_dir, logfile=logfile) for query in tqdm(list_of_regions, desc='[INFO] Checking query'))

        # filter this query_status to remove Nones
        #query_status = (l for l in query_status if not isinstance(l, type(None)))
        #query_status_peeked = predictUtils_two.peek(query_status)

        # if query_status_peeked is None:
        #     print('[INFO] Nothing to do. All predictions are available.')
        # else:
            
        #query_status = query_status_peeked
        #print('[INFO] Something to do. Making predictions.')