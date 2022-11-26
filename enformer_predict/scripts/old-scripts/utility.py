
import parsl
from parsl.app.app import python_app

def check_query(sample, query, output_dir, logfile):
    '''
    Checks of predictions are available for an individual or sample
    '''
    import pandas as pd
    import os

    if isinstance(logfile, type(None)):
        return(query)

    elif isinstance(logfile, pd.DataFrame):# should have read it>
        motif_in_logfile = logfile.motif.values
        individual_in_logfile = logfile.individual.values
        query_saved = str(f'{output_dir}/{sample}/{query}_predictions.h5')

        # four conditions 
        a = (query in motif_in_logfile)
        b = (sample in individual_in_logfile)
        if (logfile.motif.empty) | (logfile.individual.empty):
            c = False
        elif logfile.loc[(logfile.motif==query) & (logfile.individual==sample)].empty == True:
            c = False
        else:
            c = (logfile.loc[(logfile.motif==query) & (logfile.individual==sample), 'status'].values[0] == 'completed')
        d = os.path.isfile(query_saved)

        if not (a and b and c and d):
            return(query)

def extract_reference_sequence(region, fasta_func, resize_for_enformer=True, resize_length=None):
    import kipoiseq

    SEQUENCE_LENGTH = 393216

    region_split = region.split('_')
    region_chr = region_split[0]
    region_start = int(region_split[1])
    region_end = int(region_split[2])

    #fasta_object = fasta_func()

    if resize_for_enformer == True:
        reg_interval = kipoiseq.Interval(region_chr, region_start, region_end).resize(SEQUENCE_LENGTH)
    else:
        reg_interval = kipoiseq.Interval(region_chr, region_start, region_end).resize(resize_length)

    # extract the sequence 
    ref_sequences = fasta_func().extract(interval=reg_interval, anchor=[])
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

def create_individual_input_for_enformer(region, individual, vcf_func, fasta_func, hap_type = 'hap1', resize_for_enformer=True, resize_length=None):

    #vcf_object = vcf_func()

    a = extract_reference_sequence(region, fasta_func, resize_for_enformer)

    # check that all the sequences in a are valid
    if all(i == 'N' for i in a['sequences']):
        print(f'[ERROR] {region} is invalid; all nucleotides are N.')
        return(None)
    else:
        b = find_variants_in_vcf_file(vcf_func(), a['interval_object'])

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

def one_hot_encode(sequence):
    import kipoiseq
    import numpy as np
    return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)

def convert_to_tensor(one_hot_encoding):
    import tensorflow as tf
    return(tf.convert_to_tensor(one_hot_encoding))


def write_logfile(logfile_path, each_individual, what_to_write):

    import os, csv
    import pandas as pd
    logfile_csv = f'{logfile_path}/{each_individual}_predictions_log.csv'
    open_mode = 'a' if os.path.isfile(logfile_csv) else 'w'

    #if not query_status: # i.e. if the list is not empty
    with open(logfile_csv, open_mode, encoding='UTF8') as running_log_file:
        logwriter = csv.writer(running_log_file)
        if open_mode == 'w':
            logwriter.writerow(['motif', 'individual', 'status', 'sequence_type']) # write the headers
        logwriter.writerow(what_to_write)

        running_log_file.flush()
        os.fsync(running_log_file)

def enformer_predict(sequence, region, sample, seq_type, model_path, model_func, output_dir, script_path, logfile_path):

    import numpy as np # for numerical computations
    import sys # functions for interacting with the operating system
    import tensorflow as tf
    import kipoiseq
    from functools import lru_cache

    #tf.debugging.set_log_device_placement(True)

    sys.path.append(f'{script_path}/utilities')
    #import runPredictionUtilities

    try:
        enformer_model = tf.saved_model.load(model_path).model #model_func() # this load the model
        #print(f'[INFO] Model loaded successfully.')

        sequence_encoded = kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32) #one_hot_encode(sequence)
        sequence_tensor = tf.convert_to_tensor(sequence_encoded)[np.newaxis]
        predictions = enformer_model.predict_on_batch(sequence_tensor)
        del sequence_encoded
        prediction_dict = {k: v.numpy() for k, v in predictions.items()}
        target_prediction = prediction_dict['human'][0]
        obj_to_save = target_prediction[range(448 - 8, (448 + 8 + 1)), : ].squeeze()
        del target_prediction

        if tf.config.list_physical_devices('GPU'):
            print(f"[MEMORY NORMAL USAGE] {sample} | {region}: {tf.config.experimental.get_memory_info('GPU:0')['current']}")

        h5result = save_h5_prediction(obj_to_save, sample, region, seq_type, output_dir)
        write_logfile(logfile_path=logfile_path, each_individual=sample, what_to_write=h5result)

        return(0)
        
    except Exception as tfe:

        if tf.config.list_physical_devices('GPU'):
            print(f"[MEMORY ERROR USAGE] {sample} | {region}: {tf.config.experimental.get_memory_info('GPU:0')['current']}")

        print(f'[MEMORY ERROR USAGE] (device placement) for {sample} at {region} of {type(tfe)}')
        #print(f'[CACHE ERROR INFO] {model_func.cache_info()}')
        #pass
        return(1)

def generate_batch(lst, batch_size):
    """  Yields batc of specified size """
    if batch_size <= 0:
        return
    for i in range(0, len(lst), batch_size):
        yield lst[i:(i + batch_size)]
