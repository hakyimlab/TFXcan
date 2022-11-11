

import tensorflow as tf
import tensorflow_hub as hub # for interacting with saved models and tensorflow hub
import joblib
import gzip # for manipulating compressed files
from kipoiseq import Interval # same as above, really
import pyfaidx # to index our reference genome file
import pandas as pd # for manipulating dataframes
import numpy as np # for numerical computations
import os, sys, re

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

def run_predictions(sequence, region, sample, seq_type, model_path, output_dir):

    import tensorflow as tf
    import tensorflow_hub as hub # for interacting with saved models and tensorflow hub
    import joblib
    import gzip # for manipulating compressed files
    import numpy as np # for numerical computations
    import os, sys, re # functions for interacting with the operating system
    
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

    gpus = tf.config.list_physical_devices('GPU')
    print(f'GPU devices found: {gpus}')

    # gpus = tf.config.experimental.list_physical_devices('GPU')
    # if gpus:
    #     try:
    #         # Currently, memory growth needs to be the same across GPUs
    #         for gpu in gpus:
    #             tf.config.experimental.set_memory_growth(gpu, True)
    #         logical_gpus = tf.config.experimental.list_logical_devices('GPU')
    #         print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
    #     except RuntimeError as e:
    #         # Memory growth must be set before GPUs have been initialized
    #         print(e)

    # #enformer_model = tf.saved_model.load(model_path).model
    # from functools import lru_cache

    # @lru_cache(1)
    def get_model(model_class, model_path):
        return model_class(model_path)

    enformer_model = get_model(Enformer, model_path)
    print('Model loaded')

    #sequence_encoded = one_hot_encode(sequence)[np.newaxis]
    target_prediction = enformer_model.predict_on_batch(one_hot_encode(sequence)[np.newaxis])['human'][0]
    obj_to_save = target_prediction[range(448 - 8, (448 + 8 + 1)), : ].squeeze()
    h5result = save_h5_prediction(obj_to_save, sample, region, seq_type, output_dir)

    return(h5result)

def main():
    sequence_folder = "/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer-minimal-2/sequence_folder"
    each_individual = 'LuCaP_145'

    model_path = "/projects/covid-ct/imlab/data/enformer/raw"
    output_dir = "/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer-minimal-2/outputs"
    enformer_model = tf.saved_model.load(model_path).model

    with open(f'{sequence_folder}/{each_individual}_regions_sequences.txt', 'r') as f:
        sequence_list = [sequence.rstrip() for sequence in f]

    for i, sequence in enumerate(sequence_list):
        out = run_predictions(sequence, str(i), each_individual, 'var', model_path, output_dir)
        print(out)


                    
if __name__ == '__main__':
    main()


# @bash_app
# def call_single_enformer_run(call_script, model_path, output_dir, sequence_folder, sequence_info, sam, logfile_path, stderr='/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer-minimal/log/bashapp_error.err'):

#     region_id = sequence_info['region']
#     sequence_type = sequence_info['sequence_source']

#     sequence_file = f'{sequence_folder}/{sam}/{region_id}_{sequence_type}.txt'

#     with open(sequence_file, 'w') as f:
#         f.write(sequence_info['sequence'][sam])

#     return(' '.join(['bash', call_script, model_path, output_dir, sequence_file, region_id, sequence_type, sam, logfile_path]))