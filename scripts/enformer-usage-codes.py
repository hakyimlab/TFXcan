

# === This script contains all codes used to apply ENFORMER on genome sequences to predict ENCODE tracks
# === Modified by Temi from Deepmind people

import subprocess
import kipoiseq # for manipulating fasta files # for some reason, kipoiseq has to be loaded first
import time
import tensorflow as tf
import tensorflow_hub as hub # for interacting with saved models and tensorflow hub
import joblib
import gzip # for manipulating compressed files
from kipoiseq import Interval # same as above, really
import pyfaidx # to index our reference genome file
import pandas as pd # for manipulating dataframes
import numpy as np # for numerical computations
import pickle # for saving large objects
import os, sys, re # functions for interacting with the operating system

# @title `Enformer`, `EnformerScoreVariantsNormalized`, `EnformerScoreVariantsPCANormalized`,
SEQUENCE_LENGTH = 393216

class Enformer:

    def __init__(self, tfhub_url):
        #self._model = hub.load(tfhub_url).model
        self._model = tf.saved_model.load(tfhub_url).model

    def predict_on_batch(self, inputs):
        predictions = self._model.predict_on_batch(inputs)
        return {k: v.numpy() for k, v in predictions.items()}

    @tf.function
    def contribution_input_grad(self, input_sequence,
                              target_mask, output_head='human'):
        input_sequence = input_sequence[tf.newaxis]

        target_mask_mass = tf.reduce_sum(target_mask)
        with tf.GradientTape() as tape:
            tape.watch(input_sequence)
            prediction = tf.reduce_sum(target_mask[tf.newaxis] * self._model.predict_on_batch(input_sequence)[output_head]) / target_mask_mass

        input_grad = tape.gradient(prediction, input_sequence) * input_sequence
        input_grad = tf.squeeze(input_grad, axis=0)
        return tf.reduce_sum(input_grad, axis=-1)
    
class EnformerScoreVariantsRaw:

    def __init__(self, tfhub_url, organism='human'):
        self._model = Enformer(tfhub_url)
        self._organism = organism

    def predict_on_batch(self, inputs):
        ref_prediction = self._model.predict_on_batch(inputs['ref'])[self._organism]
        alt_prediction = self._model.predict_on_batch(inputs['alt'])[self._organism]

        return alt_prediction.mean(axis=1) - ref_prediction.mean(axis=1)


class EnformerScoreVariantsNormalized:

    def __init__(self, tfhub_url, transform_pkl_path,
               organism='human'):
        assert organism == 'human', 'Transforms only compatible with organism=human'
        self._model = EnformerScoreVariantsRaw(tfhub_url, organism)
        with tf.io.gfile.GFile(transform_pkl_path, 'rb') as f:
            transform_pipeline = joblib.load(f)
        self._transform = transform_pipeline.steps[0][1]  # StandardScaler.
    
    def predict_on_batch(self, inputs):
        scores = self._model.predict_on_batch(inputs)
        return self._transform.transform(scores)


class EnformerScoreVariantsPCANormalized:

    def __init__(self, tfhub_url, transform_pkl_path,
               organism='human', num_top_features=500):
        self._model = EnformerScoreVariantsRaw(tfhub_url, organism)
        with tf.io.gfile.GFile(transform_pkl_path, 'rb') as f:
            self._transform = joblib.load(f)
        self._num_top_features = num_top_features
    
    def predict_on_batch(self, inputs):
        scores = self._model.predict_on_batch(inputs)
        return self._transform.transform(scores)[:, :self._num_top_features]


# @title `variant_centered_sequences`

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

def one_hot_encode(sequence):
    return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)

# @title `plot_tracks`

def plot_tracks(tracks, interval, height=1.5, vline=True):
    fig, axes = plt.subplots(len(tracks), 1, figsize=(20, height * len(tracks)), sharex=True)
    for ax, (title, y) in zip(axes, tracks.items()):
        ax.fill_between(np.linspace(interval.start, interval.end, num=len(y)), y)
        if not isinstance(vline, type(None)):
            a = np.linspace(interval.start, interval.end, num=len(y))
            center = a[0] + (a[-1] - a[0])/2
            ax.axvline(x=center, color='red')
        ax.set_title(title)
        sns.despine(top=True, right=True, bottom=True)
    ax.set_xlabel(str(interval))
    plt.tight_layout()


# Temi's defined functions ===
def create_region_file(vcf_file, region, output_dir, individual, software_paths=[]):
    '''
    Creates a subsetted vcf file per region
    '''

    import subprocess
    import kipoiseq # for manipulating fasta files # for some reason, kipoiseq has to be loaded first
    import time
    import tensorflow as tf
    import tensorflow_hub as hub # for interacting with saved models and tensorflow hub
    import joblib
    import gzip # for manipulating compressed files
    from kipoiseq import Interval # same as above, really
    import pyfaidx # to index our reference genome file
    import pandas as pd # for manipulating dataframes
    import numpy as np # for numerical computations
    import pickle # for saving large objects
    import os, sys, re # functions for interacting with the operating system

    path_to_bcftools = software_paths[0]
    path_to_tabix = software_paths[1]

    # Center the interval at the region
    interval = kipoiseq.Interval(region[0], region[1], region[2]).resize(SEQUENCE_LENGTH) # resizing will change the regions

    path = f'{output_dir}/{individual}_{interval.chr}_{interval.start}_{interval.end}_subset_genotypes.vcf.gz'

    region = f'{interval.chr}:{interval.start}-{interval.end}'

    view_cmd = f"{path_to_bcftools} filter {vcf_file} -r {region} --output-type z --output {path} && {path_to_tabix} -p vcf {path}"

    out = subprocess.run(view_cmd, shell=True)

    return {'subset_path':path, 'interval':interval}

def extract_individual_sequence(region_details, individual, fasta_file_path, fasta_extractor, delete_region=False):

    '''
    Extracts a sequence from a reference for a region with the variants of an individual applied
    '''

    import subprocess
    import kipoiseq # for manipulating fasta files # for some reason, kipoiseq has to be loaded first
    import time
    import tensorflow as tf
    import tensorflow_hub as hub # for interacting with saved models and tensorflow hub
    import joblib
    import gzip # for manipulating compressed files
    from kipoiseq import Interval # same as above, really
    import pyfaidx # to index our reference genome file
    import pandas as pd # for manipulating dataframes
    import numpy as np # for numerical computations
    import pickle # for saving large objects
    import os, sys, re # functions for interacting with the operating system

    kseq_extractor = kipoiseq.extractors.SingleSeqVCFSeqExtractor(fasta_file=fasta_file_path, vcf_file=region_details['subset_path'])

    center = region_details['interval'].center() - region_details['interval'].start

    individual_sequences = {}
    for ind in [individual]:

        warnings.filterwarnings('error')
        
        try:
            individual_sequences[ind] = kseq_extractor.extract(interval=region_details['interval'], anchor=center, sample_id=ind)
            seq_source = 'var'
        except Warning:
            warnings.simplefilter("always", category=UserWarning)
            
            print('No variants for this region. Using reference genome.\n')
            individual_sequences[ind] = fasta_extractor.extract(interval=region_details['interval'], anchor=[])
            seq_source = 'ref'

    if delete_region == True:
        os.remove(region_details['subset_path'])
        os.remove(f"{region_details['subset_path']}.tbi")

    return {'sequence':individual_sequences, 'sequence_source':seq_source}