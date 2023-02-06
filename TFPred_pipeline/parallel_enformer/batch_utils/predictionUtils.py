# Utilities to predict using enformer
# Author: Temi
# Date: Thurs Feb 2 2023

import functools

@functools.lru_cache(5)
def get_model(model_path):
    import tensorflow as tf
    return tf.saved_model.load(model_path).model

def enformer_predict_on_sequence(model, sample_input):
    
    prediction_output = {}
    for haplotype, sequence_encoding in sample_input.items():
        if not sequence_encoding.shape == (1, 393216, 4):
            raise Exception(f'[ERROR] Fatal. Input sequence shape is not appropriate')
        prediction = model.predict_on_batch(sequence_encoding)['human'].numpy()[: , range(448 - 8, (448 + 8 + 1)), : ]
        prediction_output[haplotype] = prediction
        
    return(prediction_output)


def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]