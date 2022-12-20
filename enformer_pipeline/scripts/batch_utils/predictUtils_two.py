# Usage: This module is used to predict with ENFORMER on individual genomes
# Author: Temi
# Date:

from functools import lru_cache
import logging, json
import os, sys
whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)

global module_path # this variable should be a global one and will be used by the functions defined below
module_path = os.path.abspath(whereis_script)

global log_dir, write_log, dataset_type

with open(f'{module_path}/../../metadata/enformer_parameters.json') as f:
    parameters = json.load(f)
    project_dir = parameters['project_dir']
    log_dir = module_path + '/../../' + parameters['log_dir']
    vcf_file = module_path + '/../../' + parameters['vcf_file']
    #write_log = True if parameters["write_log"] == 'true' else False
    write_log = parameters["write_log"]
    dataset_type = parameters['dataset_type']

#print(module_path, log_dir, vcf_file)

# === Fasta extractor ==================
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
def get_fastaExtractor(script_path=module_path):

    """
    Create a fasta extractor object.

    Parameters:
        script_path: str (path), default is the path to where this file is.
            The path to the script directory
    """

    with open(f'{script_path}/../../metadata/enformer_parameters.json') as f:
        parameters = json.load(f)
        fasta_file = parameters['hg38_fasta_file']

    fasta_extractor = FastaStringExtractor(fasta_file)
    return fasta_extractor

@lru_cache(5)
def get_model(script_path=module_path):

    with open(f'{script_path}/../../metadata/enformer_parameters.json') as f:
        parameters = json.load(f)
        model_path = parameters['model_path']

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
    
    Returns: dict
        'query': the query region if it has not been logged or predictions don't exist
        'logtye': whether it should be logged if it has not been logged

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
            if (a and b) and (d == False):
                return({'query':query, 'logtype':'n'}) # log type is 'y' i.e to log or 'n' i.e not to log
            else:
                return({'query':query, 'logtype':'y'})

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

    if resize_for_enformer == True:
        reg_interval = kipoiseq.Interval(region_chr, region_start, region_end).resize(SEQUENCE_LENGTH)
    else:
        reg_interval = kipoiseq.Interval(region_chr, region_start, region_end).resize(resize_length)

    # extract the sequence 
    ref_sequences = fasta_object.extract(interval=reg_interval, anchor=[])

    if write_log['cache'] == True:
        msg_cac_log = f'[CACHE INFO] (fasta) [{get_fastaExtractor.cache_info()}, {get_gpu_name()}, {region}]'
        CACHE_LOG_FILE = f"{log_dir}/cache_usage.log"
        setup_logger('cache_log', CACHE_LOG_FILE)
        logger(msg_cac_log, 'info', 'cache')

    return({'sequences': ref_sequences, 'interval_object': reg_interval})
    
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
    query = f'{interval_object.chrom}:{interval_object.start}-{interval_object.end}'
    return([[variant.CHROM, variant.POS, variant.genotypes[0][0:2], variant.gt_bases[0].split('|')] for variant in cyvcf2_object(query)])

def create_mapping_dictionary(variants_array, interval_start, haplotype='hap1'):
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

    sequence_list = list(query_sequences)
    a = map(lambda i: mapping_dict.get(i, sequence_list[i]), range(len(sequence_list)))
    return(''.join(list(a)))

def create_individual_input_for_enformer(region, fasta_func, hap_type = 'hap1', vcf_func=None, resize_for_enformer=True, resize_length=None, write_log=write_log, dataset_type=dataset_type):
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

    if dataset_type == 'reference':
        a = extract_reference_sequence(region, fasta_func, resize_for_enformer)
    elif dataset_type == 'personalized':
        #vcf_object = vcf_func() # make_cyvcf_object() #
        a = extract_reference_sequence(region, fasta_func, resize_for_enformer)

    # check that all the sequences in a are valid
    if all(i == 'N' for i in a['sequences']) or (a is None):
        #print(f'[INPUT ERROR] {region} is invalid; all nucleotides are N.')
        if write_log['error']:
            err_msg = f'[INPUT ERROR] {region} is invalid; all nucleotides are N.'
            MEMORY_ERROR_FILE = f"{log_dir}/error_details.log"
            setup_logger('error_log', MEMORY_ERROR_FILE)
            logger(err_msg, 'error', 'run_error')

        #vcf_object.close()
        return(None)
    else:
        if dataset_type == 'reference':
            return({'sequence': a['sequences'], 'sequence_source':'ref', 'region':region})
        elif dataset_type == 'personalized':
            vcf_object = vcf_func()
            b = find_variants_in_vcf_file(vcf_object, a['interval_object'])
            vcf_object.close()
        
            if b: # go on and change the variants by position
                c = create_mapping_dictionary(b, a['interval_object'].start, haplotype=hap_type)
                variant_sequence = replace_variants_in_reference_sequence(a['sequences'], c)
                return({'sequence':variant_sequence, 'sequence_source':'var', 'region':region})
            else: # return the reference
                return({'sequence': a['sequences'], 'sequence_source':'ref', 'region':region})

def save_h5_prediction(prediction, sample, region, seq_type, output_dir):
    """
    Save ENFORMER predictions as an hdf5 file.

    Parameters:
        prediction: numpy array 
            A numpy array of predictions
        sample: str
            unique id or name of an individual
        region: str
            a region in the genome (used to provide a name for the file to be saved)
        seq_type: str
            'ref': if the sequence was from the reference genome i.e. no variants were found
            'var': if the sequence was modified based on the VCF file

    Returns:
        A list of what to log [region, sample, 'completed', sequence type] if the prediction has been saved 
    """
    import h5py
    h5save = str(f'{output_dir}/{sample}/{region}_predictions.h5')
    with h5py.File(h5save, 'w') as hf:
        hf.create_dataset(region, data=prediction)
    return([region, sample, 'completed', seq_type])

def write_logfile(predictions_log_dir, id, what_to_write):
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
    logfile_csv = f'{predictions_log_dir}/{id}_predictions_log.csv'
    open_mode = 'a' if os.path.isfile(logfile_csv) else 'w'

    #if not query_status: # i.e. if the list is not empty
    with open(logfile_csv, open_mode, encoding='UTF8') as running_log_file:
        logwriter = csv.writer(running_log_file)
        if open_mode == 'w':
            logwriter.writerow(['motif', 'individual', 'status', 'sequence_type']) # write the headers
        logwriter.writerow(what_to_write)

        running_log_file.flush()
        os.fsync(running_log_file)

def enformer_predict_on_batch(batch_region, sample, model, output_dir, predictions_log_dir, vcf_func, batch_num, grow_memory=True, write_log=write_log):
    '''
    Use ENFORMER to predict on a batch of regions

    '''

    import numpy as np # for numerical computations
    import sys # functions for interacting with the operating system
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
        enformer_model = get_model(module_path) #tf.saved_model.load(model_path).model #model_func() # this load the model
        #print(f'[CACHE INFO] (model) [{get_model.cache_info()}, {get_gpu_name()}, {sample}]')

        for each_region in batch_region:
            each_region_seq_info = create_individual_input_for_enformer(region=each_region['query'], vcf_func=vcf_func, fasta_func=None, hap_type = 'hap1', resize_for_enformer=True, resize_length=None, dataset_type=dataset_type)

            if (each_region_seq_info is not None) and (len(each_region_seq_info['sequence']) == 393216): #(b['sequence'] is not None) and (len(b['sequence']) == 393216):
                sequence_encoded = kipoiseq.transforms.functional.one_hot_dna(each_region_seq_info['sequence']).astype(np.float32) #one_hot_encode(sequence)
                sequence_tensor = tf.convert_to_tensor(sequence_encoded)[np.newaxis] # (1000, 393216, 4)
                predictions = enformer_model.predict_on_batch(sequence_tensor)
                del sequence_encoded
                prediction_dict = {k: v.numpy() for k, v in predictions.items()}
                target_prediction = prediction_dict['human'][0]
                obj_to_save = target_prediction[range(448 - 8, (448 + 8 + 1)), : ].squeeze()
                #print(f'[INFO] Shape of target prediction: {target_prediction.shape}')
                del target_prediction
                h5result = save_h5_prediction(obj_to_save, sample, each_region_seq_info['region'], each_region_seq_info['sequence_source'], output_dir)
                if each_region['logtype'] == 'n':
                    pass
                elif each_region['logtype'] == 'y':
                    write_logfile(predictions_log_dir=predictions_log_dir, id=sample, what_to_write=h5result)
                    
        tf.keras.backend.clear_session()
        gc.collect()

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
                err_mem_log = f"[GPU MEMORY] (error type {type(tfe) } for individual {sample}): free {mem_use[0]} mb, used {mem_use[1]} mb on {get_gpu_name()}"
                MEMORY_ERROR_FILE = f"{log_dir}/error_details.log"
                setup_logger('error_log', MEMORY_ERROR_FILE)
                logger(err_mem_log, 'error', 'run_error')

        else:
            print(f'[ERROR] of {type(tfe)} at predictions')

        return(1)