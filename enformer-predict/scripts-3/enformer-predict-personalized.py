# === This script contains all codes to test if tensorflow and ENFORMER are working accurately  
# Created by Temi
# DATE: Sunday Nov 13 2022

#from __future__ import absolute_import, division, print_function, unicode_literals

import os, sys, json
import pandas as pd # for manipulating dataframes
from functools import lru_cache
import math, time, tqdm
import parsl
from parsl.app.app import python_app
import tensorflow as tf
import memory_profiler
from memory_profiler import LogFile

# how many regions should I predict on 
region_range = 5
use_parsl = False

# get the path of the script as well as parameters         
whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)
script_path = os.path.abspath(whereis_script)

fpath = os.path.join(script_path, 'utilities')
sys.path.append(fpath)
print(sys.path)

# #import enformerUsageCodes
# import runPredictionUtilities


@memory_profiler.profile(stream=open('/lus/theta-fs0/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer-predict/cobalt-log/mprof_enformer_predict.log', 'w+'))
def enformer_predict(sequence, region, sample, seq_type, model_path, model_func, output_dir, script_path, logfile_path):

    import numpy as np # for numerical computations
    #import sys # functions for interacting with the operating system
    import tensorflow as tf
    import kipoiseq, gc

    #tf.debugging.set_log_device_placement(True)
    if tf.config.list_physical_devices('GPU'):
        print(f"[GPU MEMORY] (at start of prediction) {sample} | {region}: {tf.config.experimental.get_memory_info('GPU:0')['current']}")

    #sys.path.append(f'{script_path}/utilities')
    #import runPredictionUtilities
    #try:
    enformer_model = tf.saved_model.load(model_path).model #model_func() # this load the model
    #print(f'[INFO] Model loaded successfully.')

    sequence_encoded = kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32) #one_hot_encode(sequence)
    sequence_tensor = tf.convert_to_tensor(sequence_encoded)[np.newaxis]
    predictions = enformer_model.predict_on_batch(sequence_tensor)
    del sequence_encoded
    #del enformer_model 
    prediction_dict = {k: v.numpy() for k, v in predictions.items()}
    target_prediction = prediction_dict['human'][0]
    obj_to_save = target_prediction[range(448 - 8, (448 + 8 + 1)), : ].squeeze()
    #print(f'[INFO] Shape of target prediction: {target_prediction.shape}')
    del target_prediction
    h5result = save_h5_prediction(obj_to_save, sample, region, seq_type, output_dir)
    write_logfile(logfile_path=logfile_path, each_individual=sample, what_to_write=h5result)

    if tf.config.list_physical_devices('GPU'):
        print(f"[GPU MEMORY] (at end of prediction, before clearing session) {sample} | {region}: {tf.config.experimental.get_memory_info('GPU:0')['current']}")

    tf.keras.backend.clear_session()
    gc.collect()

    if tf.config.list_physical_devices('GPU'):
        print(f"[GPU MEMORY] (at end of prediction, after clearing session) {sample} | {region}: {tf.config.experimental.get_memory_info('GPU:0')['current']}")

    return(0)

    # try:
    #     enformer_model = tf.saved_model.load(model_path).model #model_func() # this load the model
    #     print(f'[INFO] Model loaded successfully.')

    #     sequence_encoded = kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32) #one_hot_encode(sequence)
    #     sequence_tensor = tf.convert_to_tensor(sequence_encoded)[np.newaxis]
    #     predictions = enformer_model.predict_on_batch(sequence_tensor)
    #     del sequence_encoded
    #     #del enformer_model 
    #     prediction_dict = {k: v.numpy() for k, v in predictions.items()}
    #     target_prediction = prediction_dict['human'][0]
    #     obj_to_save = target_prediction[range(448 - 8, (448 + 8 + 1)), : ].squeeze()
    #     del target_prediction

    #     print(f'[INFO] Shape of target prediction: {target_prediction.shape}')

    #     if tf.config.list_physical_devices('GPU'):
    #         print(f"[MEMORY] (usage normal) {sample} | {region}: {tf.config.experimental.get_memory_info('GPU:0')['current']}")

    #     h5result = save_h5_prediction(obj_to_save, sample, region, seq_type, output_dir)
    #     write_logfile(logfile_path=logfile_path, each_individual=sample, what_to_write=h5result)

    #     tf.keras.backend.clear_session()
    #     gc.collect()

    #     return(0)
        
    # except Exception as tfe:

    #     if tf.config.list_physical_devices('GPU'):
    #         print(f"[MEMORY] (usage error) {sample} | {region}: {tf.config.experimental.get_memory_info('GPU:0')['current']}")

    #     print(f'[ERROR] (of some type) for {sample} at {region} of {type(tfe)}')
    #     #print(f'[CACHE ERROR INFO] {model_func.cache_info()}')
    #     #pass
    #     return(1)

# print(runPredictionUtilities.x)
# fp_single_predictions=open(f'{script_path}/../cobalt-log/memory_profile_run_single_predictions.log', 'w+')
@memory_profiler.profile(stream=open('/lus/theta-fs0/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer-predict/cobalt-log/mprof_run_single_predictions.log', 'w+'))
def run_single_predictions(region, individual, vcf_func, fasta_func, script_path, output_dir, logfile, model_path, model_func, logfile_path): #

    import sys, os
    import tensorflow as tf

    #log the device placement
    #tf.debugging.set_log_device_placement(True)

    gpus = tf.config.experimental.list_physical_devices('GPU')
    if gpus:
        try:
            # Currently, memory growth needs to be the same across GPUs
            for gpu in gpus:
                #if tf.config.list_physical_devices('GPU'):
                print(f"[GPU MEMORY] (at very start) {sample} | {region}: {tf.config.experimental.get_memory_info('GPU:0')['current']}")
                tf.config.experimental.set_memory_growth(gpu, True)
            print(f'[RUNTIME INFO] Memory growth set for {gpu} ({individual} at {region})')
            #logical_gpus = tf.config.experimental.list_logical_devices('GPU')
            #print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
        except RuntimeError as e:
            # Memory growth must be set before GPUs have been initialized
            print(f'[RUNTIME ERROR] Memory growth failed to set for {gpu} | Error of type {type(e)} ({individual} at {region})')


    #sys.path.append(f'{script_path}/utilities')
    # mpath = os.path.dirname(__file__) #os.path.join(script_path, 'utilities')
    # sys.path.append(mpath)
    # print(sys.path)

    # import runPredictionUtilities

    #enformer_model = model_func()

    run_predictions_tools = f'{script_path}/utility.py'
    exec(open(run_predictions_tools).read(), globals(), globals())

    b = create_individual_input_for_enformer(region=region, individual=individual, vcf_func=vcf_func, fasta_func=fasta_func, hap_type = 'hap1', resize_for_enformer=True, resize_length=None)

    if (b is not None) and (len(b['sequence']) == 393216): #(b['sequence'] is not None) and (len(b['sequence']) == 393216):
        reg_prediction = enformer_predict(b['sequence'], region=b['region'], sample=individual, seq_type=b['sequence_source'], model_path = model_path, model_func=model_func, output_dir=output_dir, script_path=script_path, logfile_path=logfile_path)

        return(reg_prediction)
    else:
        print(f"[WARNING] {region}: Either length of input sequence is invalid (NoneType) or too long or too short")


run_predictions_tools = f'{script_path}/utility.py'
exec(open(run_predictions_tools).read(), globals(), globals())

fp_main=open(f'{script_path}/../cobalt-log/mprof_main.log','w+')
@memory_profiler.profile(stream=fp_main)
def main():

    if use_parsl == True:
        import parslConfiguration
        parsl.load(parslConfiguration.localParslConfig())

    #import enformerUsageCodes
    #import runPredictionUtilities
    #read the parameters file
    with open(f'{script_path}/../metadata/enformer_parameters.json') as f:

        parameters = json.load(f)

        intervals_dir = parameters['interval_list_dir']
        model_path = parameters['model_path']
        fasta_file = parameters['hg38_fasta_file']
        output_dir = parameters['output_dir']
        individuals = parameters['individuals']
        vcf_file = parameters['vcf_file']
        TF = parameters['TF']
        logfile_path = parameters['logfile_path']
        batch_size = int(parameters['batch_size'])

    @lru_cache(1)
    def get_fastaExtractor(fasta_file_path=fasta_file, script_path=script_path):

        import sys
        sys.path.append(f'{script_path}/utilities')
        import enformerUsageCodes

        fasta_extractor = enformerUsageCodes.FastaStringExtractor(fasta_file_path)
        return fasta_extractor

    # @lru_cache(12)
    # def get_model(model_path=model_path):
    #     import tensorflow as tf
    #     #print(f'[INFO] Tensorflow found {tf.config.list_physical_devices()}')
    #     return tf.saved_model.load(model_path).model

    # individuals can be a given list or a txt file of individuals per row or a single string

    if isinstance(individuals, list):
        pass
    elif isinstance(individuals, type('str')):
        if os.path.isfile(individuals):
            individuals = pd.read_table(individuals, header=None)[0].tolist()[0:1]
        else:
            individuals = [individuals]
        
    print(f'[INFO] Predicting for these individuals: {individuals}')

    for each_individual in individuals:
        
        # this is specific to an individual but is cached per individual
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
        logfile_csv = f'{logfile_path}/{each_individual}_predictions_log.csv'
        logfile = pd.read_csv(logfile_csv) if os.path.isfile(logfile_csv) else None

        tic_prediction = time.process_time()

        batches = generate_batch(list_of_regions, batch_size=batch_size)
        #print(len(list(batches)))

        count = 0
        #enformer_model = get_model()

        for batch_query in batches:

            #enformer_model = get_model() # >> results in serialization errors

            # check if the ... exists
            #first check the query
            check_result = [check_query(sample = each_individual, query = query, output_dir=output_dir, logfile=logfile) for query in batch_query]
            check_result = [c for c in check_result if c is not None]
            print(check_result)

            if tf.config.list_physical_devices('GPU'):
                print(f"[MEMORY] (GPU usage pre parsl) {each_individual} | batch {count + 1}: {tf.config.experimental.get_memory_info('GPU:0')['current']}")
                #print(f"[MEMORY] (CPU usage pre parsl) {each_individual} | batch {count + 1}: {tf.config.experimental.get_memory_info('CPU:0')['current']}")
            else:
                #print(f"[MEMORY] (usage pre parsl) {each_individual} | batch {count + 1}: {tf.config.experimental.get_memory_info('GPU:0')['current']}")
                print(f'[MEMORY] Not using the GPU yet')

            if not check_result:
                pass
            else:
                query_futures = [run_single_predictions(region=query, individual=each_individual, vcf_func=make_cyvcf_object, script_path=script_path, fasta_func=get_fastaExtractor, output_dir=output_dir, logfile=logfile, model_path=model_path, model_func=None, logfile_path=logfile_path) for query in tqdm.tqdm(check_result, desc=f'[INFO] Creating futures for batch {count + 1} of {math.ceil(len(list_of_regions)/batch_size)}')]

                #query_futures = runPredictionUtilities.run_batch_predictions(batch_regions=batch_query, individual=each_individual, vcf_func=make_cyvcf_object, script_path=script_path, fasta_func=get_fastaExtractor, output_dir=output_dir, logfile=logfile, model_path=model_path, logfile_path=logfile_path).result()

                print(query_futures)

                if use_parsl == True:
                    query_exec = [q.result() for q in tqdm.tqdm(query_futures, desc=f'[INFO] Executing inputs and predicting on {len(query_futures)} input regions')]

                #print(f'[CACHE NORMAL INFO] {get_model.cache_info()}')
                print(f'[CACHE] (fasta) {get_fastaExtractor.cache_info()}')
                print(f'[CACHE] (vcf) {make_cyvcf_object.cache_info()}')


                count = count + 1
                #time.sleep(3)

        toc_prediction = time.process_time()

        print(f'[INFO] (time) to create inputs and predict on {region_range} queries is {toc_prediction - tic_prediction}')

        if use_parsl == True:
            print(f'[INFO] Finished predictions for {each_individual}: {query_exec[0:11]} ...\n')
        else:
            print(f'[INFO] Finished predictions for {each_individual}: {query_futures[0:11]} ...\n')

        #parsl.clear()

    print(f'[INFO] Finished all predictions')                  
                    
if __name__ == '__main__':
    main()


   #for query in tqdm(list_of_regions, desc='[INFO] Creating futures for input regions'):

            # run_predictions_tools = f'{script_path}/utilities/runPredictionUtilities.py'
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

#         run_predictions_tools = f'{script_path}/utilities/runPredictionUtilities.py'
#         exec(open(run_predictions_tools).read(), globals(), globals())

#         # # load the parsl config file here since you want to distribute across individuals
#         # parsl_config = f'{script_path}/utilities/parslConfiguration.py'
#         # exec(open(parsl_config).read(), globals(), globals()) 

#         count = 0
#         for query_batch in predict_batches:

#             # run_predictions_tools = f'{script_path}/utilities/runPredictionUtilities.py'
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
        # query_status = (runPredictionUtilities.check_query(sample=each_individual, query=query, output_dir=output_dir, logfile=logfile) for query in tqdm(list_of_regions, desc='[INFO] Checking query'))

        # filter this query_status to remove Nones
        #query_status = (l for l in query_status if not isinstance(l, type(None)))
        #query_status_peeked = runPredictionUtilities.peek(query_status)

        # if query_status_peeked is None:
        #     print('[INFO] Nothing to do. All predictions are available.')
        # else:
            
        #query_status = query_status_peeked
        #print('[INFO] Something to do. Making predictions.')