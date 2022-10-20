

from __future__ import absolute_import, division, print_function, unicode_literals

import tensorflow as tf
print(f'\nTensorflow version is {tf.__version__}.\n')

print(tf.config.list_physical_devices('GPU'))

import tensorflow_hub as hub # for interacting with saved models and tensorflow hub
import joblib
import gzip # for manipulating compressed files
import kipoiseq # for manipulating fasta files
from kipoiseq import Interval # same as above, really
import pyfaidx # to index our reference genome file
import pandas as pd # for manipulating dataframes
import numpy as np # for numerical computations
import pickle # for saving large objects
import os, sys # functions for interacting with the operating system

transform_path = 'gs://dm-enformer/models/enformer.finetuned.SAD.robustscaler-PCA500-robustscaler.transform.pkl'
model_path = 'https://tfhub.dev/deepmind/enformer/1'
fasta_file_directory = '/projects/covid-ct/imlab/data/hg19_genome'

#bed_files_directory = '/projects/covid-ct/imlab/users/temi/projects/TFXcan/processed-data/motif-bed'

#if not os.path.isfile(fasta_file_directory + '/genome.fa') & os.path.isfile(fasta_file_directory + '/genome.fa.fai'):
#    !wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz -P {fasta_file_directory}
#    !gunzip -c {fasta_file_directory + '/hg19.fa.gz'} > {fasta_file_directory + '/genome.fa'}
#    pyfaidx.Faidx(fasta_file_directory + '/genome.fa')
#else:
#    print('genome.fa and genome.fa.fai files exist.')
    
    
 # Enformer architecture
# @title `Enformer`, `EnformerScoreVariantsNormalized`, `EnformerScoreVariantsPCANormalized`,
SEQUENCE_LENGTH = 393216

class Enformer:

    def __init__(self, tfhub_url):
        self._model = hub.load(tfhub_url).model

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

def create_intervals(motif_region_path):
    #collect_intervals()
    motif_regions = pd.read_table(motif_region_path, sep=' ', names=['chr', 'start', 'end', 'id', 'score', 'strand', 'set', 'num'])
    motif_regions['motif_name'] = motif_regions['set'] + motif_regions['num'].astype(str)

    temp_chr = [re.sub(pattern='chr', repl='', string=cc) for cc in motif_regions.chr.tolist()]
    temp_motif_list = motif_regions.motif_name.tolist()
    
    return (temp_chr, temp_motif_list)

def collect_intervals(chromosomes = ["1"], motif_list=['TP1']):

    '''
    Parameters : 
      chromosomes : a list of chromosome numbers; each element should be a string format
      gene_list : a list of genes; the genes should be located on those chromosomes

    Returns :
      A dictionary of genes (from gene_list) and their intervals within their respective chromosomes
    '''

    motif_intervals = {} # Collect intervals for our genes of interest
    motif_regions = pd.read_table(Gata3_motif_regions, sep=' ', names=['chr', 'start', 'end', 'id', 'score', 'strand', 'set', 'num'])
    motif_regions['motif_name'] = motif_regions['set'] + motif_regions['num'].astype(str)
    for i, motif in enumerate(motif_list):
        temp = motif_regions.loc[motif_regions['motif_name'] == motif]
        motif_intervals[motif] = [chromosomes[i], temp['start'].values[0], temp['end'].values[0]]
    return(motif_intervals)


def run_predictions(motif_intervals, fasta_extractor, model):
    '''
    Parameters :
    gene_intervals : the results from calling `collect_intervals`
    tss_dataframe : a list of the TSSs dataframes i.e. the TSS for the genes in the chromosomes
    individuals_list : a list of individuals on which we want to make predictions; defaults to None

    Returns :
    A list of predictions; the first element is the predictions around the TSS for each gene. The second is the prediction across CAGE tracks
    '''
    
    output = {}

    for motif in motif_intervals.keys():
        motif_interval = motif_intervals[motif]
        target_interval = kipoiseq.Interval("chr" + motif_interval[0],
                                        motif_interval[1],
                                        motif_interval[2]) # creates an interval to select the right sequences
                               
        target_fa = fasta_extractor.extract(target_interval.resize(SEQUENCE_LENGTH))  # extracts the fasta sequences, and resizes such that it is compatible with the sequence_length
        
        #window_coords = target_interval.resize(SEQUENCE_LENGTH) # we also need information about the start and end locations after resizing
        
        # predict on the individual's two haplotypes
        output[motif] = model.predict_on_batch(one_hot_encode(target_fa)[np.newaxis])['human'][0]

    return(output)

def Save_predictions_as_pickle(obj_to_save, path_to_save='/projects/covid-ct/imlab/users/temi/projects/TFXcan/output/reference_predictions'):
    ## Save the results as a pickle object

    #print(str(f'Saving predictions to {path_to_save}/GATA3_four_motifs_enformer_predictions_2022-07-14.pkl ===\n'))

    with open(str(f'{path_to_save}/GATA3_reference_motifs_predictions_2022-07-19.pkl'), 'wb') as obj:
        pickle.dump(obj_to_save, obj)

    print('Done saving the predictions')
    
def main():
    
#     TF_motif_regions = sys.argv[1]
#     print(TF_motif_regions)
    
#     individuals = sys.argv[2:]
#     print(individuals)

    # for Gata3
    Gata3_motif_regions = '/projects/covid-ct/imlab/users/temi/projects/TFXcan/processed-data/Gata3_motif_regions_TEMP.txt'
    
    # load the model
    model = Enformer(model_path)
    print('Model loaded.\n')

    # instantiate fasta string extractor
    fasta_extractor = FastaStringExtractor(fasta_file_directory + '/genome.fa')
    print(f'Fasta strings extracted from "{fasta_file_directory}/genome.fa".\n')
    
    # create the intervals
    int_ = create_intervals(Gata3_motif_regions)
    print('Intervals created\n')

    # collect the intervals
    motif_int = collect_intervals(int_[0], int_[1])
    print('Intervals collected\n')
    
    # make predictions
    print('\nMaking predictions ===\n')
    motif_preds = run_predictions(motif_intervals=motif_int, fasta_extractor=fasta_extractor, model=model)

    # save the predictions
    Save_predictions_as_pickle(motif_preds)

    #return(motif_preds)
    
    print('Everything done!')
    
    



def main2():
    print(__name__)

if __name__ == '__main__':
    main()
    #Save_predictions_as_pickle(res)
    
    
