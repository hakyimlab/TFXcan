
import parsl
from parsl.app.app import python_app, join_app, bash_app
import kipoiseq

#@python_app
def create_sequences(sample, region, fasta_file_path, fasta_extractor, open_vcf_file, temporary_vcf_dir, software_paths, script_path):
    import pandas as pd
    import os
    import h5py
    import numpy as np
    import kipoiseq

    usage_codes = f'{script_path}/enformer-usage-codes.py'
    exec(open(usage_codes).read(), globals(), globals())

    # create the region file
    a = create_region_file(open_vcf=open_vcf_file, region=region, subset_vcf_dir=temporary_vcf_dir, individual=sample, software_paths=software_paths)
    print(f'Region file is created')
    try:
        print('Extracting individual sequences')
        b = extract_individual_sequence(subset_dict=a, fasta_file_path=fasta_file_path, fasta_extractor=fasta_extractor, delete_region=True)
        #b['sequence'][sample] = one_hot_encode(b['sequence'][sample])[np.newaxis]
        return(b)
    except ValueError:
        #status = [0, 'NA'] #evaluate_check(is_query_done, b={'sequence_source':'NA'})
        return({'sequence':None, 'sequence_source':'NA', 'region':'NA'})

def create_input(sample, region, fasta_file_path, fasta_extractor, open_vcf_file, temporary_vcf_dir, software_paths):

    # create the region file
    a = create_region_file(open_vcf=open_vcf_file, region=region, subset_vcf_dir=temporary_vcf_dir, individual=sample, software_paths=software_paths)
    print(f'Region file is created')
    try:
        print('Extracting individual sequences')
        b = extract_individual_sequence(subset_dict=a, fasta_file_path=fasta_file_path, fasta_extractor=fasta_extractor, delete_region=True)
        b['sequence'][sample] = one_hot_encode(b['sequence'][sample])[np.newaxis]
        return(b)
    except ValueError:
        #status = [0, 'NA'] #evaluate_check(is_query_done, b={'sequence_source':'NA'})
        return({'sequence':None, 'sequence_source':'NA', 'region':'NA'})

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

@python_app
def run_predictions(sequence, region, sample, seq_type, model_path, output_dir):

    import tensorflow as tf
    import tensorflow_hub as hub # for interacting with saved models and tensorflow hub
    import joblib
    import gzip # for manipulating compressed files
    import numpy as np # for numerical computations
    import os, sys, re # functions for interacting with the operating system
    
    # define the class
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

    enformer_model = tf.saved_model.load(model_path).model

    sequence_encoded = one_hot_encode(sequence)[np.newaxis]
    target_prediction = model_predict(sequence_encoded, enformer_model)
    obj_to_save = target_prediction[range(448 - 8, (448 + 8 + 1)), : ].squeeze()
    h5result = save_h5_prediction(obj_to_save, sample, region, seq_type, output_dir)

    return(h5result)


@bash_app
def call_single_enformer_run(call_script, sequence_region, sam, stderr='/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer-minimal/log/bashapp_error.err'):
    return(' '.join(['bash', call_script, sequence_region, sam]))


# @bash_app
# def call_single_enformer_run(call_script, model_path, output_dir, sequence_folder, sequence_info, sam, logfile_path, stderr='/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer-minimal/log/bashapp_error.err'):

#     region_id = sequence_info['region']
#     sequence_type = sequence_info['sequence_source']

#     sequence_file = f'{sequence_folder}/{sam}/{region_id}_{sequence_type}.txt'

#     with open(sequence_file, 'w') as f:
#         f.write(sequence_info['sequence'][sam])

#     return(' '.join(['bash', call_script, model_path, output_dir, sequence_file, region_id, sequence_type, sam, logfile_path]))