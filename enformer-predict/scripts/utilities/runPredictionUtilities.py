import parsl
from parsl.app.app import python_app

# def peek(iterable):
#     import itertools
#     try:
#         first = next(iterable)
#     except StopIteration:
#         return None

#     return itertools.chain([first], iterable)

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

def create_region_file(vcf_file, region, subset_vcf_dir, individual, software_paths=[]):

    '''
    Creates a subsetted vcf file per region
    Arguments:
        vcf_file: A vcf that is already opened
        region: region to query, a list [chr, start, end]
        subset_vcf_dir: the directory to save the subsetted file
        individual: the individual to subset the file for
        software_paths: a list of paths to software to use [bcftools, tabix]
    '''

    import kipoiseq
    import subprocess

    path_to_bcftools = software_paths[0]
    path_to_tabix = software_paths[1]

    SEQUENCE_LENGTH = 393216
    
    reg_split = region.split('_')
    region = [reg_split[0], int(reg_split[1]), int(reg_split[2]), region]

    # Center the interval at the region
    interval = kipoiseq.Interval(region[0], region[1], region[2]).resize(SEQUENCE_LENGTH) # resizing will change the regions
    path = f'{subset_vcf_dir}/{individual}_{region[3]}_subset_genotypes.vcf.gz'
    region_interval = f'{interval.chr}:{interval.start}-{interval.end}'
    view_cmd = f"{path_to_bcftools} view {vcf_file} -r {region_interval} -s {individual} --output-type z --output-file {path} && {path_to_tabix} -p vcf {path}"
    out = subprocess.run(view_cmd, shell=True)

    return {'subset_path':path, 'interval':interval, 'individual':individual, 'region':region}


def extract_individual_sequence(subset_dict, fasta_file_path, script_path, fasta_func, delete_region=False):

    '''
    Extracts a sequence from a reference for a region with the variants of an individual applied
    Arguments:
        subset_dict: dict - the result of `create_region_file()`
        fasta_file_path: string/path - the path to the fasta file
        fasta_extractor: extractor object - in case there are not variants to apply to that region, this helps extract the reference sequence instead
        delete_region: bool - after the sequence has been applied, should the temporary vcf file be deleted?
    '''
    import kipoiseq
    import warnings
    import os
    from functools import lru_cache

    import sys
    sys.path.append(f'{script_path}/utilities')
    import runPredictionUtilities
    import enformerUsageCodes

    warnings.filterwarnings('error')

    try:
        kseq_extractor = kipoiseq.extractors.SingleSeqVCFSeqExtractor(fasta_file=fasta_file_path, vcf_file=subset_dict['subset_path'])
        center = subset_dict['interval'].center() - subset_dict['interval'].start
        individual = subset_dict['individual']
        individual_sequences = kseq_extractor.extract(interval=subset_dict['interval'], anchor=center, sample_id=individual)
        seq_source = 'var'
    except Warning:
        try:
            warnings.simplefilter("always", category=UserWarning)
            #print('No variants for this region. Using reference genome.\n')
            # load a fasta extractor by calling the lru_cached fasta extractor

            fasta_extractor = fasta_func()

            individual_sequences = fasta_extractor.extract(interval=subset_dict['interval'], anchor=[])
            fasta_extractor.close()
            seq_source = 'ref'

        except Exception as e:
            individual_sequences = 'NA'
            seq_source = 'NA'
            print(f"[ERROR] pyfaidx FetchError for {subset_dict['region'][3]}")

    if delete_region == True:
        os.remove(subset_dict['subset_path'])
        os.remove(f"{subset_dict['subset_path']}.tbi")

    return {'sequence':individual_sequences, 'sequence_source':seq_source, 'region':subset_dict['region'][3]}

def one_hot_encode(sequence):
    import kipoiseq
    import numpy as np
    return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)

def model_predict(input, model):
    predictions = model.predict_on_batch(input)
    prediction_dict = {k: v.numpy() for k, v in predictions.items()}

    return(prediction_dict['human'][0])

def save_h5_prediction(prediction, sample, region, seq_type, output_dir):
    h5save = str(f'{output_dir}/{sample}/{region}_predictions.h5')
    with h5py.File(h5save, 'w') as hf:
        hf.create_dataset(region, data=prediction)

    return([region, sample, 'completed', seq_type])


@python_app
def enformer_predict(sequence, region, sample, seq_type, model_path, output_dir, script_path, logfile_path):

    import tensorflow as tf
    import tensorflow_hub as hub # for interacting with saved models and tensorflow hub
    import joblib
    import gzip # for manipulating compressed files
    import numpy as np # for numerical computations
    import os, sys, re # functions for interacting with the operating system
    from functools import lru_cache

    #define the class
    class Enformer:
        def __init__(self, tfhub_url):
            #self._model = hub.load(tfhub_url).model
            self._model = tf.saved_model.load(tfhub_url).model

        def predict_on_batch(self, inputs):
            predictions = self._model.predict_on_batch(inputs)
            return {k: v.numpy() for k, v in predictions.items()}

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

    def model_predict(input, model):
        predictions = model.predict_on_batch(input)
        prediction_dict = {k: v.numpy() for k, v in predictions.items()}

        return(prediction_dict['human'][0])
    
    def write_logfile(logfile_path, each_individual, what_to_write):

        import os, csv
        import pandas as pd
        logfile_csv = f'{logfile_path}/{each_individual}_predictions_log.csv'
        if os.path.isfile(logfile_csv):
            open_mode = 'a'
        else:
            open_mode = 'w'

        #if not query_status: # i.e. if the list is not empty
        with open(logfile_csv, open_mode, encoding='UTF8') as running_log_file:
            logwriter = csv.writer(running_log_file)
            if open_mode == 'w':
                logwriter.writerow(['motif', 'individual', 'status', 'sequence_type']) # write the headers
            logwriter.writerow(what_to_write)

            running_log_file.flush()
            os.fsync(running_log_file)

    @lru_cache(1)
    def get_model(model_class, model_path):
        return model_class(model_path)

    enformer_model = get_model(Enformer, model_path)
    #print('Model loaded')

    sequence_encoded = one_hot_encode(sequence)[np.newaxis]
    print(f'[INFO] {region}: input matrix shape is {sequence_encoded.shape}')
    target_prediction = enformer_model.predict_on_batch(sequence_encoded)['human'][0]
    obj_to_save = target_prediction[range(448 - 8, (448 + 8 + 1)), : ].squeeze()
    h5result = save_h5_prediction(obj_to_save, sample, region, seq_type, output_dir)
    write_logfile(logfile_path=logfile_path, each_individual=sample, what_to_write=h5result)

    return(0)

@python_app
def create_input(region, individual, vcf_file, subset_vcf_dir, fasta_file_path, script_path, fasta_func, output_dir, logfile, software_paths=[]):

    # load the parsl config file here since you want to distribute across individuals
    # run_predictions_tools = f'{script_path}/utilities/runPredictionUtilities.py'
    # exec(open(run_predictions_tools).read(), globals(), globals()) 

    import sys
    sys.path.append(f'{script_path}/utilities')
    import runPredictionUtilities

    # first check the query
    chk = runPredictionUtilities.check_query(sample = individual, query = region, output_dir=output_dir, logfile=logfile)
    # #chk_result = chk.result()
    # print(f'[INFO] {region}: length of input sequence is {len(chk['sequence'])}')

    if chk is not None:
        #print(f'{chk} is not None')
        a = runPredictionUtilities.create_region_file(vcf_file=vcf_file, region=chk, subset_vcf_dir=subset_vcf_dir, individual=individual, software_paths=software_paths)
        b = runPredictionUtilities.extract_individual_sequence(subset_dict=a, fasta_file_path=fasta_file_path, script_path=script_path, fasta_func=fasta_func, delete_region=True)
        print(f"[INFO] {region}: length of input sequence is {len(b['sequence'])}")
        return(b)