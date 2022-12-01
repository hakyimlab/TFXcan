import parsl
from parsl.app.app import python_app

# def check_query(sample, query, output_dir, logfile):
#     '''
#     Checks of predictions are available for an individual or sample
#     '''
#     import pandas as pd
#     import os

#     if isinstance(logfile, pd.DataFrame):# should have read it>
#         motif_in_logfile = logfile.motif.values
#         individual_in_logfile = logfile.individual.values
#         query_saved = str(f'{output_dir}/{sample}/{query}_predictions.h5')

#         # four conditions 
#         a = (query in motif_in_logfile)
#         b = (sample in individual_in_logfile)

#         if (logfile.motif.empty) | (logfile.individual.empty):
#             c = False
#         else:
#             c = (logfile.loc[(logfile.motif==query) & (logfile.individual==sample), 'status'].values[0] == 'completed')

#         d = os.path.isfile(query_saved)

#         if a and b and c and d:
#             pass
#         else:
#             out = query
#     elif isinstance(logfile, type(None)):
#         out = query

#     return(out)

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
    path = f'{subset_vcf_dir}/{individual}_{interval.chr}_{interval.start}_{interval.end}_subset_genotypes.vcf.gz'
    region_interval = f'{interval.chr}:{interval.start}-{interval.end}'
    view_cmd = f"{path_to_bcftools} view {vcf_file} -r {region_interval} -s {individual} --output-type z --output-file {path} && {path_to_tabix} -p vcf {path}"
    out = subprocess.run(view_cmd, shell=True)

    return {'subset_path':path, 'interval':interval, 'individual':individual, 'region':region}

def extract_individual_sequence(subset_dict, fasta_file_path, fasta_extractor, delete_region=False):

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

    kseq_extractor = kipoiseq.extractors.SingleSeqVCFSeqExtractor(fasta_file=fasta_file_path, vcf_file=subset_dict['subset_path'])
    center = subset_dict['interval'].center() - subset_dict['interval'].start

    individual = subset_dict['individual']
    individual_sequences = {}
    for ind in [individual]:

        warnings.filterwarnings('error')
        
        try:
            individual_sequences[ind] = kseq_extractor.extract(interval=subset_dict['interval'], anchor=center, sample_id=ind)
            seq_source = 'var'
        except Warning:
            warnings.simplefilter("always", category=UserWarning)
            print('No variants for this region. Using reference genome.\n')
            individual_sequences[ind] = fasta_extractor.extract(interval=subset_dict['interval'], anchor=[])
            seq_source = 'ref'

    if delete_region == True:
        os.remove(subset_dict['subset_path'])
        os.remove(f"{subset_dict['subset_path']}.tbi")

    return {'sequence':individual_sequences, 'sequence_source':seq_source, 'region':subset_dict['region'][3]}

def one_hot_encode(sequence):
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

def enformer_predict(sequence, region, sample, seq_type, model_path, output_dir):

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
        return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)

    def model_predict(input, model):
        predictions = model.predict_on_batch(input)
        prediction_dict = {k: v.numpy() for k, v in predictions.items()}

        return(prediction_dict['human'][0])

    @lru_cache(1)
    def get_model(model_class, model_path):
        return model_class(model_path)

    enformer_model = get_model(Enformer, model_path)
    print('Model loaded')

    #sequence_encoded = one_hot_encode(sequence)[np.newaxis]
    target_prediction = enformer_model.predict_on_batch(one_hot_encode(sequence)[np.newaxis])['human'][0]
    obj_to_save = target_prediction[range(448 - 8, (448 + 8 + 1)), : ].squeeze()
    h5result = save_h5_prediction(obj_to_save, sample, region, seq_type, output_dir)

    return(h5result)

# @python_app
# def run_predictions(region, individual, vcf_file, subset_vcf_dir, fasta_file_path, fasta_extractor, model_path, output_dir, software_paths=[]):

#     a = create_region_file(vcf_file=vcf_file, region=region, subset_vcf_dir=subset_vcf_dir, individual=individual, software_paths=software_paths)
#     b = extract_individual_sequence(subset_dict=a, fasta_file_path=fasta_file_path, fasta_extractor=fasta_extractor, delete_region=False)
#     c = enformer_predict(b['sequence'][individual], region=region, sample=individual, seq_type=b['sequence_source'], model_path=model_path, output_dir=output_dir)


#     return(c)




#@python_app
# def create_sequences(sample, region, fasta_file_path, fasta_extractor, open_vcf_file, temporary_vcf_dir, software_paths, script_path):
#     import pandas as pd
#     import os
#     import h5py
#     import numpy as np
#     import kipoiseq

#     usage_codes = f'{script_path}/enformer-usage-codes.py'
#     exec(open(usage_codes).read(), globals(), globals())

#     # create the region file
#     a = create_region_file(open_vcf=open_vcf_file, region=region, subset_vcf_dir=temporary_vcf_dir, individual=sample, software_paths=software_paths)
#     print(f'Region file is created')
#     try:
#         print('Extracting individual sequences')
#         b = extract_individual_sequence(subset_dict=a, fasta_file_path=fasta_file_path, fasta_extractor=fasta_extractor, delete_region=True)
#         #b['sequence'][sample] = one_hot_encode(b['sequence'][sample])[np.newaxis]
#         return(b)
#     except ValueError:
#         #status = [0, 'NA'] #evaluate_check(is_query_done, b={'sequence_source':'NA'})
#         return({'sequence':None, 'sequence_source':'NA', 'region':'NA'})

# def create_input(sample, region, fasta_file_path, fasta_extractor, open_vcf_file, temporary_vcf_dir, software_paths):

#     # create the region file
#     a = create_region_file(open_vcf=open_vcf_file, region=region, subset_vcf_dir=temporary_vcf_dir, individual=sample, software_paths=software_paths)
#     print(f'Region file is created')
#     try:
#         print('Extracting individual sequences')
#         b = extract_individual_sequence(subset_dict=a, fasta_file_path=fasta_file_path, fasta_extractor=fasta_extractor, delete_region=True)
#         b['sequence'][sample] = one_hot_encode(b['sequence'][sample])[np.newaxis]
#         return(b)
#     except ValueError:
#         #status = [0, 'NA'] #evaluate_check(is_query_done, b={'sequence_source':'NA'})
#         return({'sequence':None, 'sequence_source':'NA', 'region':'NA'})
    
# @bash_app
# def call_single_enformer_run(call_script, sequence_region, sam, stderr='/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer-minimal/log/bashapp_error.err'):
#     return(' '.join(['bash', call_script, sequence_region, sam]))


# @bash_app
# def call_single_enformer_run(call_script, model_path, output_dir, sequence_folder, sequence_info, sam, logfile_path, stderr='/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer-minimal/log/bashapp_error.err'):

#     region_id = sequence_info['region']
#     sequence_type = sequence_info['sequence_source']

#     sequence_file = f'{sequence_folder}/{sam}/{region_id}_{sequence_type}.txt'

#     with open(sequence_file, 'w') as f:
#         f.write(sequence_info['sequence'][sam])

#     return(' '.join(['bash', call_script, model_path, output_dir, sequence_file, region_id, sequence_type, sam, logfile_path]))