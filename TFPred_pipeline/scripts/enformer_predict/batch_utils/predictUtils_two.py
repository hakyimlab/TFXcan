# Usage: This module is used to predict with ENFORMER on individual genomes
# Author: Temi
# Date:

from functools import lru_cache
import logging, json
import os, sys
import argparse

global module_path, log_dir, write_log, sequence_source

# read in the config_file
whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)
module_path = os.path.abspath(whereis_script)

#if __name__ == '__main__':
if __name__ == 'predictUtils_two':

    with open(f'{whereis_script}/config.json') as f:
        parameters = json.load(f)
        parameters_file = parameters['params_path']
    
    with open(f'{parameters_file}') as f:
        parameters = json.load(f)

        if parameters['sub_dir'] == True:
            project_dir = os.path.join(parameters['project_dir'], 'prediction_folder')
        else:
            project_dir = parameters['project_dir']

        sequence_source = parameters['sequence_source']
        # if sequence_source == 'personalized':
        #     vcf_file = parameters['vcf_file']

        log_dir = os.path.join(project_dir, parameters['log_dir'])
        write_log = parameters["write_log"]
        fasta_file = parameters['hg38_fasta_file']
        model_path = parameters['model_path']
        #vcf_file = parameters['vcf_file']


# === Function and class definitions start here =========================

#=== Fasta extractor ==================
class FastaStringExtractor:
    def __init__(self, fasta_file):
        import pyfaidx

        self.fasta = pyfaidx.Fasta(fasta_file)
        self._chromosome_sizes = {k: len(v) for k, v in self.fasta.items()}

    def extract(self, interval, **kwargs) -> str:
        # Truncate interval if it extends beyond the chromosome lengths.

        import kipoiseq
        chromosome_length = self._chromosome_sizes[interval.chrom]
        trimmed_interval = kipoiseq.Interval(interval.chrom,
                                    max(interval.start, 0),
                                    min(interval.end, chromosome_length),
                                    )
        # pyfaidx wants a 1-based interval
        sequence = str(self.fasta.get_seq(trimmed_interval.chrom,
                                            trimmed_interval.start + 1,
                                            trimmed_interval.stop).seq).upper()
        # Fill truncated values with N's.
        pad_upstream = 'N' * max(-interval.start, 0)
        pad_downstream = 'N' * max(interval.end - chromosome_length, 0)
        return pad_upstream + sequence + pad_downstream

    def close(self):
        return self.fasta.close()

# @lru_cache(5)
# def make_cyvcf_object(vcf_file=vcf_file, sample=vcf_id):
#     import cyvcf2
#     return(cyvcf2.cyvcf2.VCF(vcf_file, samples=sample))

@lru_cache(5)
def get_fastaExtractor(fasta_file=fasta_file):
    """
    Create a fasta extractor object.

    Parameters:
        script_path: str (path), default is the path to where this file is.
            The path to the script directory
    """
    fasta_extractor = FastaStringExtractor(fasta_file)
    return fasta_extractor

@lru_cache(5)
def get_model(model_path=model_path):
    import tensorflow as tf
    return tf.saved_model.load(model_path).model

def get_gpu_name():
    import subprocess
    cmd = "cat $COBALT_NODEFILE"
    #cmd = "cat $PBS_NODEFILE"
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

def generate_random_sequence_inputs(size=393216):
    import numpy as np
    r_seq_list = np.random.choice(['A', 'G', 'T', 'C'], size)
    return(''.join(r_seq_list))
  
def check_queries(sample, queries, output_dir, logfile):
    """
    Check whether a given region, for an individual has been predicted and logged.

    Parameters:
        sample: str 
            The name/id of an individual
        queries: A batch
            A region in the genome in the form `chr_start_end`.
        output_dir: str (path)
            The folder where the predictions should have been logged. 
        logfile: pd.DataFrame 
            A dataframe of a log file or `None` if the log file does not exist. 
    
    Returns: dict
        'query': the query region if it has not been logged or predictions don't exist
        'logtype': whether it should be logged if it has not been logged

    If predictions exist and the query has been logged, this function returns None.
    """
    import pandas as pd
    import os
    import numpy as np

    if isinstance(logfile, type(None)):
        output = [{'query':query, 'logtype':'y'} for query in queries]
        return(output)

    elif isinstance(logfile, pd.DataFrame):# should have read it>

        id_logfile = logfile.loc[logfile['individual'] == sample, : ]
        # check if the file is saved
        if sequence_source == 'personalized': # prediction must be present in two folders
            queries_saved_h1 = [str(f'{output_dir}/{sample}/haplotype1/{query}_predictions.h5') for query in queries]
            queries_saved_h2 = [str(f'{output_dir}/{sample}/haplotype2/{query}_predictions.h5') for query in queries]
            queries_saved = np.array([os.path.isfile(q) for q in queries_saved_h1]) * np.array([os.path.isfile(q) for q in queries_saved_h2])
        elif sequence_source == 'reference':
            queries_saved = [str(f'{output_dir}/{sample}/haplotype0/{query}_predictions.h5') for query in queries]
            queries_saved = np.array([os.path.isfile(q) for q in queries_saved])
        # check if the file is in the logfile
        query_logged = np.array([(query in id_logfile.motif.values) for query in queries])
        # 
        if not query_logged.shape[0] == queries_saved.shape[0]:
            raise Exception("Lengths of queries logged and saved conditions are not the same")
        queries_condition = queries_saved * query_logged
        queries_condition = queries_condition.tolist()
        # 
        output = [{'query':queries[i], 'logtype':'n'} if qc is True else {'query':queries[i], 'logtype':'y'} for i, qc in enumerate(queries_condition)]

        return(output)


def one_hot_encode(sequence, shap=None):
    import kipoiseq
    import tensorflow as tf
    import numpy as np

    try:
        sequence_encoded = kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)[np.newaxis]
        return(sequence_encoded)

    except ValueError as ve:
        print(f'[INPUT ERROR] {sequence} {len(shap)}')
        print(f'[INPUT ERROR] {sequence} {type(shap)}')


def extract_reference_sequence(region, fasta_func=None, resize_for_enformer=True, resize_length=None, write_log=write_log):

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

    try:
        region_split = region.split('_')
        region_chr = region_split[0]
        region_start = int(region_split[1])
        region_end = int(region_split[2])
    except ValueError:
        if write_log['error']:
            err_msg = f'[REGION ERROR] {region} input start or end is invalid.'
            MEMORY_ERROR_FILE = f"{log_dir}/error_details.log"
            setup_logger('error_log', MEMORY_ERROR_FILE)
            logger(err_msg, 'error', 'run_error')
        return(None)

    if fasta_func is None:
        fasta_object = get_fastaExtractor()
    else:
        fasta_object=fasta_func

    if resize_for_enformer == True:
        reg_interval = kipoiseq.Interval(region_chr, region_start, region_end).resize(SEQUENCE_LENGTH)
    else:
        reg_interval = kipoiseq.Interval(region_chr, region_start, region_end).resize(resize_length)

    # extract the sequence 
    ref_sequences = fasta_object.extract(interval=reg_interval, anchor=[])

    assert len(ref_sequences) == 393216, f'[ERROR] Length of reference sequence for region {region} is not 393216'

    if write_log['cache'] == True:
        msg_cac_log = f'[CACHE INFO] (fasta) [{get_fastaExtractor.cache_info()}, {get_gpu_name()}, {region}]'
        CACHE_LOG_FILE = f"{log_dir}/cache_usage.log"
        setup_logger('cache_log', CACHE_LOG_FILE)
        logger(msg_cac_log, 'info', 'cache')

    return({'sequence': one_hot_encode(ref_sequences), 'interval_object': reg_interval})
    
def find_variants_in_vcf_file(cyvcf2_object, interval_object):

    '''
    Find variants for a region in a vcf file

    Parameters:
        cyvcf2_object: cyvcf2_object
            A cyvcf2_object object
        interval_object: Interval object
            An interval object returned by kipoiseq.Interval

    Returns:
        A list of variants for that region in the form [chr, pos, [genotypes], [genotype bases]]
        An empty string if no variant is present
    
    '''

    import os
    query = f'{interval_object.chrom}:{interval_object.start}-{interval_object.end}'

    try: 
        return([[variant.CHROM, variant.POS, variant.genotypes[0][0:2], variant.gt_bases[0].split('|')] for variant in cyvcf2_object(query)])

    except OSError as OSE:
        #print(f'[ERROR] Cannot open {query} in vcf object.')
        if write_log['error']:
            err_mem_log = f'[ERROR] Cannot open {query} in vcf object.'
            MEMORY_ERROR_FILE = f"{log_dir}/error_details.log"
            setup_logger('error_log', MEMORY_ERROR_FILE)
            logger(err_mem_log, 'error', 'run_error')
        raise OSError(f'{query} cannot be extracted from vcf file.')

def create_mapping_dictionary(variants_array, interval_start, haplotype='both'):
    '''
    Create a map of which variants should be switched 

    Parameters:
        variants_array: a list/numpy array
            A cyvcf2_object object
        interval_start: Interval object start position
            The start position of the interval object returned by kipoiseq.Interval
        haplotype: str; default: 'hap1'; options: hap2, both
            Should haplotype 1 or 2 or both be used?

    Returns: dict
        {position: new genotype base}
    '''
    import numpy as np

    A = np.array([1,0,0,0], dtype=np.float32)
    C = np.array([0,1,0,0], dtype=np.float32)
    G = np.array([0,0,1,0], dtype=np.float32)
    T = np.array([0,0,0,1], dtype=np.float32)
    seq_dict = {'A': A, 'C': C, 'G': G, 'T': T}

    if haplotype == 'hap1':
        haplotype_map = {(variants_array[i][1] - interval_start): variants_array[i][3][0] for i in range(0, len(variants_array))}
        return(haplotype_map)
        #return({(k - interval_start): v for k, v in haplotype_map.items() if v is not None})
    elif haplotype == 'hap2':
        haplotype_map = {(variants_array[i][1] - interval_start): variants_array[i][3][1] for i in range(0, len(variants_array))}
        return(haplotype_map)
    elif haplotype == 'both':
        hap1 = {(variants_array[i][1] - interval_start): seq_dict[variants_array[i][3][0]] for i in range(0, len(variants_array))}
        hap2 = {(variants_array[i][1] - interval_start): seq_dict[variants_array[i][3][1]] for i in range(0, len(variants_array))}
        haplotype_map = {'haplotype1': hap1, 'haplotype2': hap2}
        return(haplotype_map)

def replace_variants_in_reference_sequence(query_sequences_encoded, mapping_dict):
    '''
    Using the mapping dictionary, mutate many variants in a given sequence 

    Parameters:
        query_sequences: str
            The query, perhaps, reference sequence
        mapping_dict: a numpy array or list
            List of the regions that should be modified

    Returns:
        A new sequence with all the changes applied.
    '''

    import copy
    import numpy as np

    if len(mapping_dict) != 2:
        raise Exception('Length of mapping dictionary is not 2')
    if query_sequences_encoded.shape != (1, 393216, 4):
        raise Exception(f'Shape of input numpy array is {query_sequences_encoded.shape}')

    try:
        output = {}
        haps = ['haplotype1', 'haplotype2']
        for j in range(len(mapping_dict)):
            ref_i = copy.deepcopy(query_sequences_encoded) # very important
            indices = np.array(list(mapping_dict[haps[j]].keys()))
            bases = list(mapping_dict[haps[j]].values())

            # print(indices)
            # print(bases)
            ref_i[:, indices, : ] = bases
            output[haps[j]] = ref_i

        return(output)

    except IndexError as ie:
        # if write_log['error']:
        #     err_mem_log = f'[ERROR] Cannot find variants in {ref_i.shape}.'
        #     MEMORY_ERROR_FILE = f"{log_dir}/error_details.log"
        #     setup_logger('error_log', MEMORY_ERROR_FILE)
        #     logger(err_mem_log, 'error', 'run_error')
        #raise IndexError(f'Cannot change variants in {ref_i.shape}.')
        return(None)


def create_input_for_enformer(region_details, fasta_func, hap_type = 'both', vcf_func=None, resize_for_enformer=True, resize_length=None, write_log=None, sequence_source=None):
    '''
    Given a region in the genome, a reference genome (fasta) and a VCF file, create an individual's sequence for that region

    Parameters:
        query_sequences: str
            The query, perhaps, reference sequence
        mapping_dict: a numpy array or list
            List of the regions that should be modified

    Returns:
        A new sequence with all the changes applied.
    '''

    import numpy as np

    region = region_details['query']
    logtype = region_details['logtype']

    if sequence_source == 'random':
        result = {'sequence': {'haplotype0': one_hot_encode(generate_random_sequence_inputs())}, 'sequence_source':'random', 'region':region, 'logtype': logtype}
    else:
        reference_sequence = extract_reference_sequence(region, fasta_func, resize_for_enformer)
        if np.all(reference_sequence['sequence'] == 0.25): # check if all the sequence are "NNNNNNNNNNN..."
            if write_log['error']:
                err_msg = f'[INPUT ERROR] {region} is invalid; all nucleotides are N.'
                MEMORY_ERROR_FILE = f"{log_dir}/error_details.log"
                setup_logger('error_log', MEMORY_ERROR_FILE)
                logger(err_msg, 'error', 'run_error')

            result = None
            return(result)

        else:
            if sequence_source == 'reference':
                result = {'sequence': {'haplotype0': reference_sequence['sequence']}, 'sequence_source':'ref', 'region':region, 'logtype': logtype}
                # check that result is intact and return it
                if result['sequence']['haplotype0'].shape != (1, 393216, 4):
                    if write_log['error']:
                        err_mem_log = f"[ERROR] Reference sequences for {region} is wrong. Shape is {result['sequence']['haplotype0'].shape}" # {result['sequence']['haplotype1'].shape}
                        MEMORY_ERROR_FILE = f"{log_dir}/error_details.log"
                        setup_logger('error_log', MEMORY_ERROR_FILE)
                        logger(err_mem_log, 'error', 'run_error')
                    else:
                        raise Exception('Shapes of input reference sequence is wrong.')
                else:
                    return(result)

            elif sequence_source == 'personalized':
                vcf_object = vcf_func()
                variants_coordinates = find_variants_in_vcf_file(vcf_object, reference_sequence['interval_object'])
                
                if variants_coordinates: # go on and change the variants by position
                    mapping_dictionary = create_mapping_dictionary(variants_coordinates, reference_sequence['interval_object'].start, haplotype=hap_type)

                    variant_sequence = replace_variants_in_reference_sequence(reference_sequence['sequence'], mapping_dictionary)

                    if variant_sequence is not None:
                        result = {'sequence':{'haplotype1': variant_sequence['haplotype1'], 'haplotype2': variant_sequence['haplotype2']}, 'sequence_source':'var', 'region':region, 'logtype': logtype}
                    else:
                        if write_log['error']:
                            err_mem_log = f'[ERROR] Cannot find variants in {region}. Using the reference.'
                            MEMORY_ERROR_FILE = f"{log_dir}/error_details.log"
                            setup_logger('error_log', MEMORY_ERROR_FILE)
                            logger(err_mem_log, 'error', 'run_error')

                        result = {'sequence':{'haplotype1': reference_sequence['sequence'], 'haplotype2': reference_sequence['sequence']}, 'sequence_source':'ref', 'region':region, 'logtype': logtype}

                    # check that result is intact and return it
                    if result['sequence']['haplotype1'].shape != (1, 393216, 4) or result['sequence']['haplotype2'].shape != (1, 393216, 4):
                        if write_log['error']:
                            err_mem_log = f"[ERROR] Variant sequences for {region} is wrong. Shape is {result['sequence']['haplotype1'].shape}" # {result['sequence']['haplotype1'].shape}
                            MEMORY_ERROR_FILE = f"{log_dir}/error_details.log"
                            setup_logger('error_log', MEMORY_ERROR_FILE)
                            logger(err_mem_log, 'error', 'run_error')
                        else:
                            raise Exception('Shapes of input variant sequence is wrong.')
                    else:
                        return(result)

        

def write_logfile(predictions_log_file, what_to_write):
    '''
    Write the prediction status to a log file

    Parameters:
        predictions_log_dir: str (path)
            A folder wherein to create the file in which to write the log
        each_individual: str
            The unique id of an individual
        what_to_write: list
            A list, comma separated, of what to write to the log file
        
    Returns:
        None
    '''

    import os, csv
    #logfile_csv = f'{predictions_log_dir}/{id}_predictions_log.csv'
    open_mode = 'a' if os.path.isfile(predictions_log_file) else 'w'

    #if not query_status: # i.e. if the list is not empty
    with open(predictions_log_file, open_mode, encoding='UTF8') as running_log_file:
        logwriter = csv.writer(running_log_file)
        if open_mode == 'w':
            logwriter.writerow(['motif', 'sample', 'status', 'sequence_type']) # write the headers
        logwriter.writerows(what_to_write)

        running_log_file.flush()
        os.fsync(running_log_file)

    return(0)

def enformer_predict_on_sequence(model, prepared_region):
    
    prediction_output = {}
    for haplotype, sequence_encoding in prepared_region['sequence'].items():
        prediction = model.predict_on_batch(sequence_encoding)['human'].numpy()[: , range(448 - 8, (448 + 8 + 1)), : ]

        prediction_output[haplotype] = prediction

    # delete `sequence` from the dictionary and return it
    del prepared_region['sequence']

    return( prediction_output, prepared_region )


def save_haplotypes_h5_prediction(haplotype_predictions, metadata, output_dir, sample):

    import h5py
    import os

    #output = []
    for key, values in batch_predictions.items():
        houtput = os.path.join(output_dir, sample, key)
        if not os.path.exists(houtput): os.makedirs(houtput)
        for i in range(0, values.shape[0]):
            region = metadata['region']
            h5save = str(f'{houtput}/{region}_predictions.h5')
            with h5py.File(h5save, 'w') as hf:
                hf.create_dataset(region, data=values[i, :])
    
    if metadata['logtype'] == 'y':
        output = [metadata['region'], sample, 'completed', metadata['sequence_source']]
    else:
        output = None

    return(output)

def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]


def enformer_predict_on_batch(batch_regions, sample, model, output_dir, predictions_log_file, vcf_func, batch_num, grow_memory=True, write_log=write_log):
    '''
    Use ENFORMER to predict on a batch of regions

    '''
    import numpy as np # for numerical computations
    import os, sys # functions for interacting with the operating system
    import tensorflow as tf
    import kipoiseq, gc

    if grow_memory == True:
        gpus = tf.config.experimental.list_physical_devices('GPU')
        if gpus:
            try:
                # Currently, memory growth needs to be the same across GPUs
                for gpu in gpus:
                    tf.config.experimental.set_memory_growth(gpu, True)
            except RuntimeError as e:
                # Memory growth must be set before GPUs have been initialized
                print(f'[RUNTIME ERROR] {sample} of {type(e)}')

    try:
        enformer_model = get_model(model_path) #tf.saved_model.load(model_path).model #model_func() # this load the model
        #print(f'[CACHE INFO] (model) [{get_model.cache_info()}, {get_gpu_name()}, {sample}]')

        for query_region in batch_regions:
            query_region_prepared = create_input_for_enformer(region_details=query_region, fasta_func=None, hap_type = 'both', vcf_func=vcf_func, resize_for_enformer=True, resize_length=None, write_log=write_log, sequence_source=sequence_source)

            if query_region_prepared is not None:

                pred, metadata = enformer_predict_on_sequence(enformer_model, query_region_prepared)
                saving_output = save_batch_h5_prediction(pred, metadata, output_dir, sample)

                if saving_output is not None:
                    writing_output = write_logfile(predictions_log_file, what_to_write=saving_output)

            if write_log['memory']:
                if tf.config.list_physical_devices('GPU'):
                    mem_use = get_gpu_memory()
                    #msg_mem_log = f"[GPU MEMORY] (at end of prediction, after clearing session) {sample} | {region}: {tf.config.experimental.get_memory_info('GPU:0')['current']/1e+9}"
                    msg_mem_log = f"[GPU MEMORY] (at end of prediction on {sample}\'s batch {batch_num}): free {mem_use[0]} mb, used {mem_use[1]} mb on {get_gpu_name()}"
                    MEMORY_LOG_FILE = f"{log_dir}/memory_usage.log"
                    setup_logger('memory_log', MEMORY_LOG_FILE)
                    logger(msg_mem_log, 'info', 'memory')

            if write_log['cache']:
                msg_cac_log = f'[CACHE INFO] (model) on batch {batch_num}: [{get_model.cache_info()}, {get_gpu_name()}, {sample}]'
                CACHE_LOG_FILE = f"{log_dir}/cache_usage.log"
                setup_logger('cache_log', CACHE_LOG_FILE)
                logger(msg_cac_log, 'info', 'cache')
            
        return(0) # for that batch
        
    except (TypeError, AttributeError) as tfe:
        if write_log['error']:
            if tf.config.list_physical_devices('GPU'):
                mem_use = get_gpu_memory()
                err_mem_log = f"[GPU MEMORY] (error type {type(tfe)} for individual {sample}): free {mem_use[0]} mb, used {mem_use[1]} mb on {get_gpu_name()}"
                MEMORY_ERROR_FILE = f"{log_dir}/error_details.log"
                setup_logger('error_log', MEMORY_ERROR_FILE)
                logger(err_mem_log, 'error', 'run_error')

        else:
            print(f'[ERROR] of {type(tfe)} at predictions')

        return(1)



  # for each_region in batch_region:
        #     each_region_seq_info = create_individual_input_for_enformer(region=each_region['query'], vcf_func=vcf_func, fasta_func=None, hap_type = 'both', resize_for_enformer=True, resize_length=None, sequence_source=sequence_source)

        #     if each_region_seq_info is not None:
        #         for haplotype, v_sequence in each_region_seq_info['sequence'].items():
        #             # ensure sequence validity
        #             if len(v_sequence) == 393216:
        #                 obj_to_save = predict_on_haplotype(model=enformer_model, sequence=v_sequence)
        #                 h5result = save_h5_prediction(obj_to_save, sample, each_region_seq_info['region'], each_region_seq_info['sequence_source'], os.path.join(output_dir, sample, haplotype))

        #         if each_region['logtype'] == 'n':
        #             pass
        #         elif each_region['logtype'] == 'y':
        #             write_logfile(predictions_log_file=predictions_log_file, id=sample, what_to_write=h5result)
        #     elif each_region_seq_info is None:
        #         pass



f#or bbatch in batch(batch_region, n=10):

        #     if sequence_source == 'personalized':
        #         t1, t2, variant_details = collect_inputs_into_dict(bbatch, vcf_func=vcf_func)

        #         if (t1 is not None) and (t2 is not None):
        #             predict_s1 = enformer_model.predict_on_batch(t1['haplotype1'])['human'].numpy()[: , range(448 - 8, (448 + 8 + 1)), : ]
        #             predict_s2 = enformer_model.predict_on_batch(t2['haplotype2'])['human'].numpy()[: , range(448 - 8, (448 + 8 + 1)), : ]

        #             print(f'[INFO] Final shape is {predict_s1.shape}')

        #             assert predict_s1.shape[0] == len(variant_details), '[ERROR] The shape of the predictions don\'t match the variant details'

        #             bpredictions = {'haplotype1': predict_s1, 'haplotype2': predict_s2}

        #             del t1
        #             del t2
        #             del predict_s1
        #             del predict_s2
        #         else:
        #             return(0)

        #     elif sequence_source == 'reference':
        #         t0, _, variant_details = collect_inputs_into_dict(bbatch,vcf_func=vcf_func)
        #         predict_s0 = enformer_model.predict_on_batch(t0['haplotype0'])['human'].numpy()[: , range(448 - 8, (448 + 8 + 1)), : ]

        #         assert predict_s0.shape[0] == len(variant_details), '[ERROR] The shape of the predictions don\'t match the variant details'

        #         bpredictions = {'haplotype0': predict_s0}

        #         del t0
        #         del predict_s0


# def collect_inputs_into_dict(batch_region, vcf_func):

#     import tensorflow as tf

#     output_sequences = [create_individual_input_for_enformer(region_details=each_region, vcf_func=vcf_func, fasta_func=None, hap_type = 'both', resize_for_enformer=True, resize_length=None, sequence_source=sequence_source) for each_region in batch_region] #>> [{sequence:{}, region: '', sequence_source: ''}, ..., {...}]
#     #output_metadata = [{'logtype': each_region['logtype']} for each_region in batch_region]
#     output = [out for out in output_sequences if out is not None]

#     if not output:
#         return(None, None, None)

#     else:
#         if sequence_source == 'personalized':
#             hap1 = [out_dict['sequence']['haplotype1'] for out_dict in output]
#             hap2 = [out_dict['sequence']['haplotype2'] for out_dict in output]

#             variant_details = [{'sequence_source': out_dict['sequence_source'],
#                                 'region': out_dict['region'],
#                                 'logtype': out_dict['logtype']} for out_dict in output]

#             shap1 = [create_enformer_input(s) for s in hap1]
#             shap2 = [create_enformer_input(s) for s in hap2]

#             # print(f'{len(shap1)} ========= {len(shap2)}')
#             # print(f'{type(shap1)} ========= {type(shap2)}')


#             if (len(hap1) == 1) and (len(hap2) == 1):
#                 tensor1 = {'haplotype1': shap1[0]}
#                 tensor2 = {'haplotype2': shap2[0]}
            
#             elif (len(hap1) > 1) and (len(hap2) > 1):
#                 tensor1 = {'haplotype1': tf.concat(shap1, axis=0)}
#                 tensor2 = {'haplotype2': tf.concat(shap2, axis=0)}

#             #query_names = [out_dict['query'] for out_dict in output]
#         elif sequence_source == 'reference':
#             hap0 = [out_dict['sequence']['haplotype0'] for out_dict in output]
#             variant_details = [{'sequence_source': out_dict['sequence_source'],
#                                 'region': out_dict['region'],
#                                 'logtype': out_dict['logtype']} for out_dict in output]

#             tensor1 = {'haplotype0': tf.concat([create_enformer_input(s) for s in hap0], axis=0)}
#             tensor2 = None 

#         return( tensor1, tensor2, variant_details )

