import parsl
from parsl.app.app import python_app
from functools import lru_cache

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

    fasta_object = fasta_func()

    if resize_for_enformer == True:
        reg_interval = kipoiseq.Interval(region_chr, region_start, region_end).resize(SEQUENCE_LENGTH)
    else:
        reg_interval = kipoiseq.Interval(region_chr, region_start, region_end).resize(resize_length)

    # extract the sequence 
    ref_sequences = fasta_object.extract(interval=reg_interval, anchor=[])
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

    #import types
    #import cyvcf2

    # if (isinstance(vcf_object, types.FunctionType)) and (type(vcf_object) != cyvcf2.cyvcf2.VCF): # not a string
    #     print('VCF object is a function to call load the vcf.')
    #     vcf_object = vcf_object()
    # # elif isinstance(vcf_object, type('a')): # or a path to a vcf file
    # #     import cyvcf2
    # #     vcf_object = cyvcf2.cyvcf2.VCF(vcf_object, samples=individual)
    # elif callable(vcf_object) and (type(vcf_object) == cyvcf2.cyvcf2.VCF):
    #     pass

    vcf_object = vcf_func()

    a = extract_reference_sequence(region, fasta_func, resize_for_enformer)

    # check that all the sequences in a are valid
    if all(i == 'N' for i in a['sequences']):
        print(f'[ERROR] {region} is invalid; all nucleotides are N.')
        return(None)
    else:
        b = find_variants_in_vcf_file(vcf_object, a['interval_object'])

        if b: # go on and change the variants by position
            c = create_mapping_dictionary(b, a['interval_object'].start, haplotype=hap_type)
            variant_sequence = replace_variants_in_reference_sequence(a['sequences'], c)
            return({'sequence':variant_sequence, 'sequence_source':'var', 'region':region})
        else: # return the reference
            return({'sequence': a['sequences'], 'sequence_source':'ref', 'region':region})

#define the class
# class Enformer:
    
#     def __init__(self, tfhub_url):
#         #self._model = hub.load(tfhub_url).model
#         import tensorflow as tf
#         self._model = tf.saved_model.load(tfhub_url).model

#     def predict_on_batch(self, inputs):
#         predictions = self._model.predict_on_batch(inputs)
#         return {k: v.numpy() for k, v in predictions.items()}
# class Enformer:

#     from functools import lru_cache

#     @classmethod
#     @lru_cache(maxsize=1)
#     def __init__(self, tfhub_url):
#         #self._model = hub.load(tfhub_url).model
#         import tensorflow as tf
#         self._model = tf.saved_model.load(tfhub_url).model

#     @classmethod
#     @lru_cache(maxsize=1)
#     def predict_on_batch(self, inputs):
#         predictions = self._model.predict_on_batch(inputs)
#         return {k: v.numpy() for k, v in predictions.items()}

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

def model_predict(input, model):
    predictions = model.predict_on_batch(input)
    prediction_dict = {k: v.numpy() for k, v in predictions.items()}

    return(prediction_dict['human'][0])

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

    sys.path.append(f'{script_path}/utilities')
    import runPredictionUtilities

    try:
        #enformer_model = model_func()
        #print(f'[INFO] Model loaded successfully.')
        #print('Model loaded')

        #with tf.device('/physical_device:CPU:0'):
        sequence_encoded = runPredictionUtilities.one_hot_encode(sequence)
        sequence_tensor = runPredictionUtilities.convert_to_tensor(sequence_encoded)[np.newaxis]
        del sequence_encoded
        #print(f'[INFO] {region}: input matrix shape is {sequence_encoded.shape}')
    
        target_prediction = runPredictionUtilities.model_predict(input=sequence_tensor, model=model_func) #['human'][0]
        #with tf.device('/physical_device:CPU:0'):
        obj_to_save = target_prediction[range(448 - 8, (448 + 8 + 1)), : ].squeeze()

        del target_prediction

        if tf.config.list_physical_devices('GPU'):
            print(f"[MEMORY NORMAL USAGE] {sample} | {region}: {tf.config.experimental.get_memory_info('GPU:0')['current']}")
        
        #print(f'[CACHE NORMAL INFO] {model_func.cache_info()}')

        h5result = runPredictionUtilities.save_h5_prediction(obj_to_save, sample, region, seq_type, output_dir)
        runPredictionUtilities.write_logfile(logfile_path=logfile_path, each_individual=sample, what_to_write=h5result)

        return(0)
        
    except (TypeError, AttributeError) as tfe:

        if tf.config.list_physical_devices('GPU'):
            print(f"[MEMORY ERROR USAGE] {sample} | {region}: {tf.config.experimental.get_memory_info('GPU:0')['current']}")

        print(f'[MEMORY ERROR USAGE] (device placement) for {sample} at {region} of {type(tfe)}')
        #print(f'[CACHE ERROR INFO] {model_func.cache_info()}')
        #pass
        return(1)

    
@python_app
def run_single_predictions(region, individual, vcf_func, fasta_func, script_path, output_dir, logfile, model_path, model_func, logfile_path): #

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

    import sys, os
    #sys.path.append(f'{script_path}/utilities')
    mpath = os.path.join(script_path, 'utilities')
    sys.path.append(mpath)
    #print(sys.path)
    import runPredictionUtilities

    #first check the query
    check_result = runPredictionUtilities.check_query(sample = individual, query = region, output_dir=output_dir, logfile=logfile) #for region in region_batch
    #print(check_result)

    if check_result is not None:
        #print(f'{chk} is not None')
        b = runPredictionUtilities.create_individual_input_for_enformer(region=check_result, individual=individual, vcf_func=vcf_func, fasta_func=fasta_func, hap_type = 'hap1', resize_for_enformer=True, resize_length=None)

        if (b is not None) and (len(b['sequence']) == 393216): #(b['sequence'] is not None) and (len(b['sequence']) == 393216):
            reg_prediction = runPredictionUtilities.enformer_predict(b['sequence'], region=b['region'], sample=individual, seq_type=b['sequence_source'], model_path = model_path, model_func=model_func, output_dir=output_dir, script_path=script_path, logfile_path=logfile_path)
            return(reg_prediction)
        else:
            print(f"[WARNING] {check_result}: Either length of input sequence is invalid (NoneType) or too long or too short")


def generate_batch(lst, batch_size):
    """  Yields batc of specified size """
    if batch_size <= 0:
        return
    for i in range(0, len(lst), batch_size):
        yield lst[i:(i + batch_size)]





















# ======================================================================
# import parsl
# from parsl.app.app import python_app, join_app

# def check_query(sample, query, output_dir, logfile):
#     '''
#     Checks of predictions are available for an individual or sample
#     '''
#     import pandas as pd
#     import os

#     if isinstance(logfile, type(None)):
#         return(query)

#     elif isinstance(logfile, pd.DataFrame):# should have read it>
#         motif_in_logfile = logfile.motif.values
#         individual_in_logfile = logfile.individual.values
#         query_saved = str(f'{output_dir}/{sample}/{query}_predictions.h5')

#         # four conditions 
#         a = (query in motif_in_logfile)
  
#         b = (sample in individual_in_logfile)
       
#         if (logfile.motif.empty) | (logfile.individual.empty):
#             c = False
#         elif logfile.loc[(logfile.motif==query) & (logfile.individual==sample)].empty == True:
#             c = False
#         else:
#             c = (logfile.loc[(logfile.motif==query) & (logfile.individual==sample), 'status'].values[0] == 'completed')
        
#         d = os.path.isfile(query_saved)

#         if not (a and b and c and d):
#             return(query)

# def create_region_file(vcf_file, region, subset_vcf_dir, individual, software_paths=[]):

#     '''
#     Creates a subsetted vcf file per region
#     Arguments:
#         vcf_file: A vcf that is already opened
#         region: region to query, a list [chr, start, end]
#         subset_vcf_dir: the directory to save the subsetted file
#         individual: the individual to subset the file for
#         software_paths: a list of paths to software to use [bcftools, tabix]
#     '''

#     import kipoiseq
#     import subprocess

#     path_to_bcftools = software_paths[0]
#     path_to_tabix = software_paths[1]

#     SEQUENCE_LENGTH = 393216
    
#     reg_split = region.split('_')
#     region = [reg_split[0], int(reg_split[1]), int(reg_split[2]), region]

#     # Center the interval at the region
#     interval = kipoiseq.Interval(region[0], region[1], region[2]).resize(SEQUENCE_LENGTH) # resizing will change the regions
#     path = f'{subset_vcf_dir}/{individual}_{region[3]}_subset_genotypes.vcf.gz'
#     region_interval = f'{interval.chr}:{interval.start}-{interval.end}'
#     view_cmd = f"{path_to_bcftools} view {vcf_file} -r {region_interval} -s {individual} --output-type z --output-file {path} && {path_to_tabix} -p vcf --force {path}"
#     out = subprocess.run(view_cmd, shell=True)

#     return {'subset_path':path, 'interval':interval, 'individual':individual, 'region':region}


# #@python_app
# def extract_individual_sequence(subset_dict, fasta_file_path, script_path, fasta_func, delete_region=False):

#     '''
#     Extracts a sequence from a reference for a region with the variants of an individual applied
#     Arguments:
#         subset_dict: dict - the result of `create_region_file()`
#         fasta_file_path: string/path - the path to the fasta file
#         fasta_extractor: extractor object - in case there are not variants to apply to that region, this helps extract the reference sequence instead
#         delete_region: bool - after the sequence has been applied, should the temporary vcf file be deleted?
#     '''
#     import kipoiseq, pyfaidx
#     import warnings
#     import os
#     from functools import lru_cache

#     import sys
#     sys.path.append(f'{script_path}/utilities')
#     import runPredictionUtilities
#     import enformerUsageCodes

#     SEQUENCE_LENGTH = 393216

#     warnings.filterwarnings('error')

#     try:
#         kseq_extractor = kipoiseq.extractors.SingleSeqVCFSeqExtractor(fasta_file=fasta_file_path, vcf_file=subset_dict['subset_path'])
#         center = subset_dict['interval'].center() - subset_dict['interval'].start
#         individual = subset_dict['individual']
#         individual_sequences = kseq_extractor.extract(interval=subset_dict['interval'], anchor=center, sample_id=individual)
#         seq_source = 'var'
        
#     except Warning:
#         warnings.simplefilter("always", category=UserWarning)
#         print(f"[WARNING] (from pyfaidx) No variants for {subset_dict['region'][3]}. Extracting from the reference genome.")
#         # load a fasta extractor by calling the lru_cached fasta extractor
#         fasta_extractor = fasta_func()
#         individual_sequences = fasta_extractor.extract(interval=subset_dict['interval'], anchor=[])
#         fasta_extractor.close()
#         seq_source = 'ref'

#     except (ValueError, pyfaidx.FetchError) as err: #Exception as e:
#         individual_sequences = None
#         seq_source = None
#         print(f"[ERROR] (from pyfaidx) {type(err)} for {subset_dict['region'][3]}")

#     if delete_region == True:
#         os.remove(subset_dict['subset_path'])
#         os.remove(f"{subset_dict['subset_path']}.tbi")

#     # ensure that the length of the sequence is 
#     #if len(individual_sequences) == SEQUENCE_LENGTH:
#     return {'sequence':individual_sequences, 'sequence_source':seq_source, 'region':subset_dict['region'][3]}
#     # else:
#     #     print(f"[SEQUENCE LENGTH ERROR] Length of {subset_dict['region'][3]} ({len(individual_sequences)}) after extracting is not equal to {SEQUENCE_LENGTH}.")

# def one_hot_encode(sequence):
#     import kipoiseq
#     import numpy as np
#     return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)

# def model_predict(input, model):
#     predictions = model.predict_on_batch(input)
#     prediction_dict = {k: v.numpy() for k, v in predictions.items()}

#     return(prediction_dict['human'][0])

# def save_h5_prediction(prediction, sample, region, seq_type, output_dir):
#     h5save = str(f'{output_dir}/{sample}/{region}_predictions.h5')
#     with h5py.File(h5save, 'w') as hf:
#         hf.create_dataset(region, data=prediction)

#     return([region, sample, 'completed', seq_type])


# #@python_app
# def enformer_predict(sequence, region, sample, seq_type, model_path, model_func, output_dir, script_path, logfile_path):

#     import tensorflow as tf
#     import tensorflow_hub as hub # for interacting with saved models and tensorflow hub
#     import joblib
#     import gzip # for manipulating compressed files
#     import numpy as np # for numerical computations
#     import os, sys, re # functions for interacting with the operating system
#     from functools import lru_cache

#     #define the class
#     class Enformer:
#         def __init__(self, tfhub_url):
#             #self._model = hub.load(tfhub_url).model
#             self._model = tf.saved_model.load(tfhub_url).model

#         def predict_on_batch(self, inputs):
#             predictions = self._model.predict_on_batch(inputs)
#             return {k: v.numpy() for k, v in predictions.items()}

#     def save_h5_prediction(prediction, sample, region, seq_type, output_dir):
#         import h5py
#         h5save = str(f'{output_dir}/{sample}/{region}_predictions.h5')
#         with h5py.File(h5save, 'w') as hf:
#             hf.create_dataset(region, data=prediction)
#         return([region, sample, 'completed', seq_type])

#     def one_hot_encode(sequence):
#         import kipoiseq
#         import numpy as np
#         return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)

#     def model_predict(input, model):
#         predictions = model.predict_on_batch(input)
#         prediction_dict = {k: v.numpy() for k, v in predictions.items()}

#         return(prediction_dict['human'][0])
    
#     def write_logfile(logfile_path, each_individual, what_to_write):

#         import os, csv
#         import pandas as pd
#         logfile_csv = f'{logfile_path}/{each_individual}_predictions_log.csv'
#         if os.path.isfile(logfile_csv):
#             open_mode = 'a'
#         else:
#             open_mode = 'w'

#         #if not query_status: # i.e. if the list is not empty
#         with open(logfile_csv, open_mode, encoding='UTF8') as running_log_file:
#             logwriter = csv.writer(running_log_file)
#             if open_mode == 'w':
#                 logwriter.writerow(['motif', 'individual', 'status', 'sequence_type']) # write the headers
#             logwriter.writerow(what_to_write)

#             running_log_file.flush()
#             os.fsync(running_log_file)

#     # @lru_cache(1)
#     # def get_model(model_class, model_path):
#     #     return model_class(model_path)

#     enformer_model = model_func(Enformer, model_path)
#     #print('Model loaded')

#     sequence_encoded = one_hot_encode(sequence)[np.newaxis]
#     #print(f'[INFO] {region}: input matrix shape is {sequence_encoded.shape}')
#     target_prediction = enformer_model.predict_on_batch(sequence_encoded)['human'][0]
#     obj_to_save = target_prediction[range(448 - 8, (448 + 8 + 1)), : ].squeeze()
#     h5result = save_h5_prediction(obj_to_save, sample, region, seq_type, output_dir)
#     write_logfile(logfile_path=logfile_path, each_individual=sample, what_to_write=h5result)

#     return(0)


# @python_app
# def run_single_predictions(region, individual, vcf_file, subset_vcf_dir, fasta_file_path, script_path, fasta_func, output_dir, logfile, model_path, model_func, logfile_path, software_paths=[]): #

#     # load the parsl config file here since you want to distribute across individuals
#     # run_predictions_tools = f'{script_path}/utilities/runPredictionUtilities.py'
#     # exec(open(run_predictions_tools).read(), globals(), globals()) 

#     import sys
#     sys.path.append(f'{script_path}/utilities')
#     import runPredictionUtilities

#     # first check the query
#     check_result = runPredictionUtilities.check_query(sample = individual, query = region, output_dir=output_dir, logfile=logfile) #for region in region_batch

#     if check_result is not None:
#         #print(f'{chk} is not None')
#         a = runPredictionUtilities.create_region_file(vcf_file=vcf_file, region=check_result, subset_vcf_dir=subset_vcf_dir, individual=individual, software_paths=software_paths)
#         b = runPredictionUtilities.extract_individual_sequence(subset_dict=a, fasta_file_path=fasta_file_path, script_path=script_path, fasta_func=fasta_func, delete_region=True)

#         if (b['sequence'] is not None) and (len(b['sequence']) == 393216):
#             return runPredictionUtilities.enformer_predict(b['sequence'], region=b['region'], sample=individual, seq_type=b['sequence_source'], model_path=model_path, model_func=model_func, output_dir=output_dir, script_path=script_path, logfile_path=logfile_path)
#         else:
#             print(f"[WARNING] {check_result}: Either length of input sequence is invalid (NoneType) or too long or too short")


# # @python_app
# # def create_input(region, individual, vcf_file, subset_vcf_dir, fasta_file_path, script_path, fasta_func, output_dir, logfile, model_path, model_func, logfile_path, software_paths=[]): #(model_path, model_func, )

# #     # load the parsl config file here since you want to distribute across individuals
# #     # run_predictions_tools = f'{script_path}/utilities/runPredictionUtilities.py'
# #     # exec(open(run_predictions_tools).read(), globals(), globals()) 

# #     import sys
# #     sys.path.append(f'{script_path}/utilities')
# #     import runPredictionUtilities

# #     # first check the query
# #     chk_result = runPredictionUtilities.check_query(sample = individual, query = region, output_dir=output_dir, logfile=logfile) #for region in region_batch
# #     #chk_result = chk_appfuture.result()
# #     # print(f'[INFO] {region}: length of input sequence is {len(chk['sequence'])}')

# #     if chk_result is not None:
# #         #print(f'{chk} is not None')
# #         a = runPredictionUtilities.create_region_file(vcf_file=vcf_file, region=chk_result, subset_vcf_dir=subset_vcf_dir, individual=individual, software_paths=software_paths)
# #         b = runPredictionUtilities.extract_individual_sequence(subset_dict=a, fasta_file_path=fasta_file_path, script_path=script_path, fasta_func=fasta_func, delete_region=True)

# #         #check that length of sequence is compatible
# #         l_sequence = len(b['sequence'])
# #         if l_sequence == 393216:
# #             #return(b) 
# #             runPredictionUtilities.enformer_predict(b['sequence'], region=b['region'], sample=individual, seq_type=b['sequence_source'], model_path=model_path, model_func=model_func, output_dir=output_dir, script_path=script_path, logfile_path=logfile_path)
# #         else:
# #             print(f"[INFO] {chk_result}: length of input sequence is {l_sequence}")

# #         runPredictionUtilities.enformer_predict()