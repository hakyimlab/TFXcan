
# === This script contains all codes used to apply ENFORMER on genome sequences to predict ENCODE tracks
# === Modified by Temi from Deepmind people

import tensorflow as tf
import tensorflow_hub as hub # for interacting with saved models and tensorflow hub
from kipoiseq import Interval # same as above, really
import pyfaidx # to index our reference genome file
import pandas as pd # for manipulating dataframes
import numpy as np # for numerical computations

# class Enformer:

#     from functools import lru_cache

#     @lru_cache(maxsize=1)
#     def __init__(self, tfhub_url):
#         #self._model = hub.load(tfhub_url).model
#         import tensorflow as tf
#         self._model = tf.saved_model.load(tfhub_url).model

#     @lru_cache(maxsize=1)
#     def predict_on_batch(self, inputs):
#         predictions = self._model.predict_on_batch(inputs)
#         return {k: v.numpy() for k, v in predictions.items()}

class FastaStringExtractor:
    
    def __init__(self, fasta_file):
        self.fasta = pyfaidx.Fasta(fasta_file)
        self._chromosome_sizes = {k: len(v) for k, v in self.fasta.items()}

    def extract(self, interval: Interval, **kwargs) -> str:
        # Truncate interval if it extends beyond the chromosome lengths.
        chromosome_length = self._chromosome_sizes[interval.chrom]
        trimmed_interval = Interval(interval.chrom,
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