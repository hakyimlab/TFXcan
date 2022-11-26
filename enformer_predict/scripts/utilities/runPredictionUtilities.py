import parsl
from parsl.app.app import python_app
from functools import lru_cache
import logging

import os, sys
whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)
script_path = os.path.abspath(whereis_script)
fpath = os.path.join(script_path, 'utilities')
sys.path.append(fpath)

@lru_cache(5)
def get_fastaExtractor(script_path=script_path):

    """
    Create a fasta extractor object.

    Parameters:
        script_path: str (path), default is the path to where this file is.
            The path to the script directory
    """

    import sys, json, os
    fpath = os.path.join(script_path, 'utilities')
    sys.path.append(fpath)

    with open(f'{script_path}/../../metadata/enformer_parameters.json') as f:
        parameters = json.load(f)
        fasta_file = parameters['hg38_fasta_file']

    import enformerUsageCodes

    fasta_extractor = enformerUsageCodes.FastaStringExtractor(fasta_file)
    return fasta_extractor

@lru_cache(5)
def get_model(script_path=script_path):

    import sys, json, os
    fpath = os.path.join(script_path, 'utilities')
    sys.path.append(fpath)

    with open(f'{script_path}/../../metadata/enformer_parameters.json') as f:
        parameters = json.load(f)
        model_path = parameters['model_path']

    import tensorflow as tf
    
    return tf.saved_model.load(model_path).model

def get_gpu_name():
    import subprocess
    cmd = "cat $COBALT_NODEFILE"
    a = str(subprocess.run(cmd, shell=True, capture_output=True).stdout, encoding='utf-8').strip('\n')

    # cmd = "nvidia-smi -L" #"nvidia-smi --query-gpu=gpu_bus_id --format=csv"
    # b = str(subprocess.run(cmd, shell=True, capture_output=True).stdout, encoding='utf-8').strip('\n')

    # return(f"{a}-{b}")
    return(a)

def get_gpu_memory():
    import subprocess
    command = "nvidia-smi --query-gpu=memory.free,memory.used --format=csv"
    memory_info = subprocess.check_output(command.split()).decode('ascii').split('\n')[1].split(',')
    memory_values = [int(x.strip().split(' ')[0]) for i, x in enumerate(memory_info)]
    return memory_values

def setup_logger(logger_name, log_file, level=logging.INFO):

    log_setup = logging.getLogger(logger_name)
    formatter = logging.Formatter('[%(levelname)s: %(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    fileHandler = logging.FileHandler(log_file, mode='a')
    fileHandler.setFormatter(formatter)
    streamHandler = logging.StreamHandler()
    streamHandler.setFormatter(formatter)
    log_setup.setLevel(level)
    log_setup.addHandler(fileHandler)
    log_setup.addHandler(streamHandler)

    return None

def logger(msg, level, logfile):
 
    if logfile == 'memory'   : log = logging.getLogger('memory_log')
    if logfile == 'cache'   : log = logging.getLogger('cache_log') 
    if logfile == 'run_error'    : log = logging.getLogger('error_log')

    if level == 'info'    : log.info(msg) 
    if level == 'warning' : log.warning(msg)
    if level == 'error'   : log.error(msg)

    return None

def check_query(sample, query, output_dir, logfile):
    """
    Check whether a given region, for an individual has been predicted and logged.

    Parameters:
        sample: str 
            The name/id of an individual
        query: str
            A region in the genome in the form `chr_start_end`.
        output_dir: str (path)
            The folder where the predictions should have been logged. 
        logfile: pd.DataFrame 
            A dataframe of a log file or `None` if the log file does not exist. 
    
    Returns:
        dict of (1) the query region if it has not been logged or predictions don't exist and (2) whether it should be logged if it has not been logged. 

    If predictions exist and the query has been logged, this function returns None.
    """
    import pandas as pd
    import os

    if isinstance(logfile, type(None)):
        return({'query':query, 'logtype':'y'})

    elif isinstance(logfile, pd.DataFrame):# should have read it>
        motif_in_logfile = logfile.motif.values
        individual_in_logfile = logfile.individual.values
        query_saved = str(f'{output_dir}/{sample}/{query}_predictions.h5')

        # four conditions 
        a = (query in motif_in_logfile) # is the query already written
        b = (sample in individual_in_logfile)
        if (logfile.motif.empty) | (logfile.individual.empty):
            c = False
        elif logfile.loc[(logfile.motif==query) & (logfile.individual==sample)].empty == True:
            c = False
        else:
            c = (logfile.loc[(logfile.motif==query) & (logfile.individual==sample), 'status'].values[0] == 'completed')
        d = os.path.isfile(query_saved)

        if not (a and b and c and d):
            #return({'query':query, 'logtype':'y'})
            if (a and b) and (d == False):
                #print(f'Query exists in logfile')
                return({'query':query, 'logtype':'n'}) # log type is 'y' i.e to log or 'n' i.e not to log
            else:
                return({'query':query, 'logtype':'y'})

def extract_reference_sequence(region, fasta_func=None, resize_for_enformer=True, resize_length=None, write_log=True, script_path=script_path):

    """
    Given a region, extract the reference sequence from a fasta file through a fastaextractor object

    Parameters:
        region: str
            A region in the genome in the form `chr_start_end`.
        fasta_func:
            A function that returns a fastaExtractor object. This is currently not needed and depends on the `get_fastaExtractor` function in this module.
        resize_for_enformer: bool
            Should a given region be resized to an appropriate input of 393216 bp for ENFORMER?
        resize_length: int or None
            If `resize_for_enformer` is false, resize the sequence up to the supplied argument. If None, resizing is not done. 
        write_log: bool
            Should info/error/warnings be written to a log file? Depends on the `setup_logger` and `logger` functions in this module as well as the `logging` module. 
        script_path: str (path), default is the path to where this file is.
            The path to this module.

    Returns: dict
        sequences: str
            The sequences extracted from the fastaExtractor object
        interval_object: kipoiseq.Interval object
            The interval object giving the chr, start, end, and other information.
    """

    import kipoiseq

    SEQUENCE_LENGTH = 393216

    region_split = region.split('_')
    region_chr = region_split[0]
    region_start = int(region_split[1])
    region_end = int(region_split[2])

    if fasta_func is None:
        fasta_object = get_fastaExtractor()

    if resize_for_enformer == True:
        reg_interval = kipoiseq.Interval(region_chr, region_start, region_end).resize(SEQUENCE_LENGTH)
    else:
        reg_interval = kipoiseq.Interval(region_chr, region_start, region_end).resize(resize_length)

    # extract the sequence 
    ref_sequences = fasta_object.extract(interval=reg_interval, anchor=[])

    if write_log == True:
        msg_cac_log = f'[CACHE INFO] (fasta) [{get_fastaExtractor.cache_info()}, {get_gpu_name()}, {region}]'
        CACHE_LOG_FILE = f"{script_path}/../../cobalt-log/cache_usage.log"
        setup_logger('cache_log', CACHE_LOG_FILE)
        logger(msg_cac_log, 'info', 'cache')

    #print(f'[CACHE INFO] (fasta) {get_fastaExtractor.cache_info()}')
    return({'sequences': ref_sequences, 'interval_object': reg_interval})
    
def find_variants_in_vcf_file(cyvcf2_object, interval_object):

    '''
    Find variants for a region in a vcf file
    Returns an empty string if no variant is present
    
    '''
    query = f'{interval_object.chrom}:{interval_object.start}-{interval_object.end}'
    return([[variant.CHROM, variant.POS, variant.genotypes[0][0:2], variant.gt_bases[0].split('|')] for variant in cyvcf2_object(query)])

def create_mapping_dictionary(variants_array, interval_start, haplotype='hap1'):
    if haplotype == 'hap1':
        haplotype_map = {(variants_array[i][1] - interval_start): variants_array[i][3][0] for i in range(0, len(variants_array)) if variants_array[i][2][0] == 1}
        return(haplotype_map)
        #return({(k - interval_start): v for k, v in haplotype_map.items() if v is not None})
    elif haplotype == 'hap2':
        haplotype_map = {(variants_array[i][1] - interval_start): variants_array[i][3][1] for i in range(0, len(variants_array)) if variants_array[i][2][1] == 1}
        return(haplotype_map)
    elif haplotype == 'both':
        haplotype_map = {variants_array[i][1]: variants_array[i][3][0] if variants_array[i][2][0] == 1 else variants_array[i][3][1] if variants_array[i][2][1] == 1 else None for i in range(0, len(variants_array))}
        haplotype_map = {(k - interval_start): v for k, v in haplotype_map.items() if v is not None}
        return(haplotype_map)

def replace_variants_in_reference_sequence(query_sequences, mapping_dict):

    sequence_list = list(query_sequences)
    a = map(lambda i: mapping_dict.get(i, sequence_list[i]), range(len(sequence_list)))
    return(''.join(list(a)))

def create_individual_input_for_enformer(region, individual, fasta_func, hap_type = 'hap1', vcf_func=None, resize_for_enformer=True, resize_length=None, script_path=script_path, write_log=True):

    vcf_object = vcf_func()

    a = extract_reference_sequence(region, fasta_func, resize_for_enformer)

    if write_log:
        msg_cac_log = f'[CACHE INFO] (vcf) [{vcf_func.cache_info()}, {get_gpu_name()}, {individual}, {region}]'
        CACHE_LOG_FILE = f"{script_path}/../../cobalt-log/cache_usage.log"
        setup_logger('cache_log', CACHE_LOG_FILE)
        logger(msg_cac_log, 'info', 'cache')

    # check that all the sequences in a are valid
    if all(i == 'N' for i in a['sequences']):
        #print(f'[INPUT ERROR] {region} is invalid; all nucleotides are N.')
        if write_log:
            err_msg = f'[INPUT ERROR] {region} is invalid; all nucleotides are N.'
            MEMORY_ERROR_FILE = f"{script_path}/../../cobalt-log/error_details.log"
            setup_logger('error_log', MEMORY_ERROR_FILE)
            logger(err_msg, 'error', 'run_error')

        vcf_object.close()
        return(None)
    else:
        b = find_variants_in_vcf_file(vcf_object, a['interval_object'])
        vcf_object.close()
        
        if b: # go on and change the variants by position
            c = create_mapping_dictionary(b, a['interval_object'].start, haplotype=hap_type)
            variant_sequence = replace_variants_in_reference_sequence(a['sequences'], c)
            return({'sequence':variant_sequence, 'sequence_source':'var', 'region':region})
        else: # return the reference
            return({'sequence': a['sequences'], 'sequence_source':'ref', 'region':region})

def save_h5_prediction(prediction, sample, region, seq_type, output_dir):
    import h5py
    h5save = str(f'{output_dir}/{sample}/{region}_predictions.h5')
    with h5py.File(h5save, 'w') as hf:
        hf.create_dataset(region, data=prediction)
    return([region, sample, 'completed', seq_type])

def write_logfile(log_dir, each_individual, what_to_write):

    import os, csv
    logfile_csv = f'{log_dir}/{each_individual}_predictions_log.csv'
    open_mode = 'a' if os.path.isfile(logfile_csv) else 'w'

    #if not query_status: # i.e. if the list is not empty
    with open(logfile_csv, open_mode, encoding='UTF8') as running_log_file:
        logwriter = csv.writer(running_log_file)
        if open_mode == 'w':
            logwriter.writerow(['motif', 'individual', 'status', 'sequence_type']) # write the headers
        logwriter.writerow(what_to_write)

        running_log_file.flush()
        os.fsync(running_log_file)

def enformer_predict(sequence, region, sample, seq_type, model_func, output_dir, log_dir, logtype, script_path = script_path, grow_memory=True, write_log=True):

    import numpy as np # for numerical computations
    import sys # functions for interacting with the operating system
    import tensorflow as tf
    import kipoiseq, gc

    #sys.path.append(f'{script_path}/utilities')
    
    # try:
    #     import runPredictionUtilities
    # except ModuleNotFoundError as merr:
    #     err_msg = f"f'[MODULE ERROR] At enformer_predict, the module runPredictionUtilties was not found."
    #     MEMORY_ERROR_FILE = f"{script_path}/../cobalt-log/error_details.log"
    #     setup_logger('error_log', MEMORY_ERROR_FILE)
    #     logger(err_msg, 'error', 'run_error')

    if grow_memory == True:
        gpus = tf.config.experimental.list_physical_devices('GPU')
        if gpus:
            try:
                # Currently, memory growth needs to be the same across GPUs
                for gpu in gpus:
                    tf.config.experimental.set_memory_growth(gpu, True)
                #logical_gpus = tf.config.experimental.list_logical_devices('GPU')
                #print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
            except RuntimeError as e:
                # Memory growth must be set before GPUs have been initialized
                print(f'[RUNTIME ERROR] {sample} at {region} of {type(e)}')

    try:
        enformer_model = get_model() #tf.saved_model.load(model_path).model #model_func() # this load the model
        sequence_encoded = kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32) #one_hot_encode(sequence)
        sequence_tensor = tf.convert_to_tensor(sequence_encoded)[np.newaxis]
        predictions = enformer_model.predict_on_batch(sequence_tensor)
        del sequence_encoded
        del enformer_model 
        prediction_dict = {k: v.numpy() for k, v in predictions.items()}
        target_prediction = prediction_dict['human'][0]
        obj_to_save = target_prediction[range(448 - 8, (448 + 8 + 1)), : ].squeeze()
        #print(f'[INFO] Shape of target prediction: {target_prediction.shape}')
        del target_prediction
        h5result = save_h5_prediction(obj_to_save, sample, region, seq_type, output_dir)
        if logtype == 'n':
            pass
        elif logtype == 'y':
            write_logfile(log_dir=log_dir, each_individual=sample, what_to_write=h5result)

        # if tf.config.list_physical_devices('GPU'):
        #     print(f"[GPU MEMORY] (at end of prediction, before clearing session) {sample} | {region}: {tf.config.experimental.get_memory_info('GPU:0')['current']}")

        tf.keras.backend.clear_session()
        gc.collect()

        if write_log:
            if tf.config.list_physical_devices('GPU'):
                mem_use = get_gpu_memory()
                #msg_mem_log = f"[GPU MEMORY] (at end of prediction, after clearing session) {sample} | {region}: {tf.config.experimental.get_memory_info('GPU:0')['current']/1e+9}"
                msg_mem_log = f"[GPU MEMORY] (at end of prediction on {sample}, after clearing session for region {region}): free {mem_use[0]} mb, used {mem_use[1]} mb on {get_gpu_name()}"
                MEMORY_LOG_FILE = f"{script_path}/../../cobalt-log/memory_usage.log"
                temp = setup_logger('memory_log', MEMORY_LOG_FILE)
                temp = logger(msg_mem_log, 'info', 'memory')

            msg_cac_log = f'[CACHE INFO] (model) [{get_model.cache_info()}, {get_gpu_name()}, {sample}, {region}]'
            CACHE_LOG_FILE = f"{script_path}/../../cobalt-log/cache_usage.log"
            temp = setup_logger('cache_log', CACHE_LOG_FILE)
            temp = logger(msg_cac_log, 'info', 'cache')
            
        return(0)
        
    except (TypeError, AttributeError) as tfe:

        if write_log:
            if tf.config.list_physical_devices('GPU'):
                mem_use = get_gpu_memory()
                err_mem_log = f"[GPU MEMORY] (error type {type(tfe) } for individual {sample} for region {region}): free {mem_use[0]} mb, used {mem_use[1]} mb on {get_gpu_name()}"
                MEMORY_ERROR_FILE = f"{script_path}/../cobalt-log/error_details.log"
                setup_logger('error_log', MEMORY_ERROR_FILE)
                logger(err_mem_log, 'error', 'run_error')

            mem_use = get_gpu_memory()
            err_msg = f"[PREDICTION ERROR] (error type {type(tfe)} for individual {sample} for region {region}): free {mem_use[0]} mb, used {mem_use[1]} mb on {get_gpu_name()}"
            MEMORY_ERROR_FILE = f"{script_path}/../cobalt-log/error_details.log"
            setup_logger('error_log', MEMORY_ERROR_FILE)
            logger(err_msg, 'error', 'run_error')
        else:
            print(f'[ERROR] of {type(tfe)} at predictions')

        return(1)

@python_app
def run_single_predictions(region, individual, vcf_func, script_path, output_dir, logfile, log_dir): #

    import sys, os
    #sys.path.append(f'{script_path}/utilities')
    mpath = os.path.join(script_path, 'utilities') #os.path.dirname(__file__) #
    sys.path.append(mpath)
    # print(sys.path)
    
    try:
        import runPredictionUtilities
    except ModuleNotFoundError as merr:
        print(f'[MODULE NOT FOUND ERROR] at run_single_predictions')

    #first check the query
    check_result = runPredictionUtilities.check_query(sample = individual, query = region, output_dir=output_dir, logfile=logfile)

    if check_result is not None:
        #print(check_result)
        b = runPredictionUtilities.create_individual_input_for_enformer(region=check_result['query'], individual=individual, vcf_func=vcf_func, fasta_func=None, hap_type = 'hap1', resize_for_enformer=True, resize_length=None)

        #print(f'[CACHE INFO] (vcf) {vcf_func.cache_info()}')

        if (b is not None) and (len(b['sequence']) == 393216): #(b['sequence'] is not None) and (len(b['sequence']) == 393216):
            #print(type(b))
            reg_prediction = runPredictionUtilities.enformer_predict(b['sequence'], region=b['region'], sample=individual, seq_type=b['sequence_source'], model_func=None, output_dir=output_dir, log_dir=log_dir, logtype=check_result['logtype'])
            
            return(reg_prediction)
        else:
            print(f"[WARNING] {check_result}: Either length of input sequence is invalid (NoneType) or too long or too short")
            


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

@python_app
def run_batch_predictions(batch_regions, individual, vcf_func, fasta_func, script_path, output_dir, logfile, model_path, log_dir): #

    import sys, os
    #sys.path.append(f'{script_path}/utilities')
    mpath = os.path.join(script_path, 'utilities')
    sys.path.append(mpath)
    #print(sys.path)

    try:
        import runPredictionUtilities
    except ModuleNotFoundError as merr:
        print(f'[MODULE NOT FOUND ERROR] at run_batch_predictions')

    #first check the query
    check_result = [runPredictionUtilities.check_query(sample = individual, query = region, output_dir=output_dir, logfile=logfile) for region in batch_regions]
    check_result = [c for c in check_result if c is not None]
    #print(check_result)

    output = []
    for reg in check_result:
        if reg is not None:
            #print(f'{chk} is not None')
            b = runPredictionUtilities.create_individual_input_for_enformer(region=reg, individual=individual, vcf_func=vcf_func, fasta_func=fasta_func, hap_type = 'hap1', resize_for_enformer=True, resize_length=None)

            if (b is not None) and (len(b['sequence']) == 393216): #(b['sequence'] is not None) and (len(b['sequence']) == 393216):
                reg_prediction = runPredictionUtilities.enformer_predict(b['sequence'], region=b['region'], sample=individual, seq_type=b['sequence_source'], model_path = model_path, model_func=get_model(), output_dir=output_dir, script_path=script_path, log_dir=log_dir)

                print(f'[CACHE NORMAL INFO] {get_model.cache_info()}')

                output.append(reg_prediction)
            else:
                print(f"[WARNING] {reg}: Either length of input sequence is invalid (NoneType) or too long or too short")
    return(output)




 #import tensorflow as tf

    #log the device placement
    # tf.debugging.set_log_device_placement(True)

    # gpus = tf.config.experimental.list_physical_devices('GPU')
    # if gpus:
    #     try:
    #         # Currently, memory growth needs to be the same across GPUs
    #         for gpu in gpus:
    #             tf.config.experimental.set_memory_growth(gpu, True)
    #         #logical_gpus = tf.config.experimental.list_logical_devices('GPU')
    #         #print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
    #     except RuntimeError as e:
    #         # Memory growth must be set before GPUs have been initialized
    #         print(f'[RUNTIME ERROR] {individual} at {region} of {type(e)}')


# def one_hot_encode(sequence):
#     import kipoiseq
#     import numpy as np
#     return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)

# def convert_to_tensor(one_hot_encoding):
#     import tensorflow as tf
#     return(tf.convert_to_tensor(one_hot_encoding))

# def model_predict(input, model):
#     predictions = model.predict_on_batch(input)
#     prediction_dict = {k: v.numpy() for k, v in predictions.items()}

#     return(prediction_dict['human'][0])