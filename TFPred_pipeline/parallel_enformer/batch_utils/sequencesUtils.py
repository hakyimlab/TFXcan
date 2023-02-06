
# Usage: Utilities to extract and manipulate DNA sequences
# Author: Temi
# Date: Thurs Feb 2 2023

import functools
import loggerUtils

class FastaStringExtractor:
    def __init__(self, fasta_file):
        import pyfaidx

        self.fasta = pyfaidx.Fasta(fasta_file)
        self._chromosome_sizes = {k: len(v) for k, v in self.fasta.items()}

    def extract(self, interval, **kwargs) -> str:
        # Truncate interval if it extends beyond the chromosome lengths.

        import kipoiseq
        chromosome_length = self._chromosome_sizes[interval.chrom]
        trimmed_interval = kipoiseq.Interval(interval.chrom,
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


@functools.lru_cache(5)
def get_fastaExtractor(fasta_file):
    """
    Create a fasta extractor object.

    Parameters:
        script_path: str (path), default is the path to where this file is.
            The path to the script directory
    """
    fasta_extractor = FastaStringExtractor(fasta_file)
    return fasta_extractor

def generate_random_sequence_inputs(size=393216):
    import numpy as np
    r_seq_list = np.random.choice(['A', 'G', 'T', 'C'], size)
    return(''.join(r_seq_list))
  
def one_hot_encode(sequence):

    import kipoiseq
    import numpy as np

    if not isinstance(sequence, str):
        raise Exception(f'[ERROR] Input to be one-hot encoded must be a str type. You provided a {type(sequence_encoded).__name__} type.')

    #try:
    sequence_encoded = kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)[np.newaxis]
    return(sequence_encoded)


def extract_reference_sequence(region, fasta_func=None, resize_for_enformer=True, resize_length=None, print_sequence = False, write_log=None):

    """
    Given a region, extract the reference sequence from a fasta file through a fastaextractor object

    Parameters:
        region: str
            A region in the genome in the form `chr_start_end`.
        fasta_func:
            A function that returns a fastaExtractor object. This is currently not needed and depends on the `get_fastaExtractor` function in this module.
        resize_for_enformer: bool
            Should a given region be resized to an appropriate input of 393216 bp for ENFORMER?
        resize_length: int or None
            If `resize_for_enformer` is false, resize the sequence up to the supplied argument. If None, resizing is not done. 
        write_log: bool
            Should info/error/warnings be written to a log file? Depends on the `setup_logger` and `logger` functions in this module as well as the `logging` module. 
        script_path: str (path), default is the path to where this file is.
            The path to this module.

    Returns: dict
        sequences: str
            The sequences extracted from the fastaExtractor object
        interval_object: kipoiseq.Interval object
            The interval object giving the chr, start, end, and other information.
    """

    import kipoiseq

    SEQUENCE_LENGTH = 393216

    if not region.startswith('chr'):
        raise Exception(f'[ERROR] Input region does not start with chromosome name. e.g. `chr1_start_end`, but you supplied {region}.')
    if fasta_func is None:
        raise Exception(f'[ERROR] Fatal for `fasta_func`. Please, pass in a fasta extractor object.')

    try:
        region_split = region.split('_')
        region_chr = region_split[0]
        region_start = int(region_split[1])
        region_end = int(region_split[2])
    except ValueError:
        if (write_log is not None) and write_log['logtypes']['error']:
            err_msg = f'[REGION ERROR] {region} input start or end is invalid.'
            MEMORY_ERROR_FILE = f"{write_log['logdir']}/error_details.log"
            loggerUtils.write_logger(log_msg_type='error', logfile=MEMORY_ERROR_FILE, message=err_msg)
        return(None)
    
    if resize_for_enformer == True:
        reg_interval = kipoiseq.Interval(region_chr, region_start, region_end).resize(SEQUENCE_LENGTH)
    else:
        reg_interval = kipoiseq.Interval(region_chr, region_start, region_end).resize(resize_length)

    # extract the sequence 
    #print(f'extracting reference sequence')
    ref_sequences = fasta_func.extract(interval=reg_interval, anchor=[])
    #print(f'finished extracting reference sequence of length {len(ref_sequences)}')

    if print_sequence == True:
        print(ref_sequences)

    if (write_log is not None) and (write_log['logtypes']['cache'] == True):
        msg_cac_log = f'[CACHE INFO] (fasta) [{fasta_func.cache_info()}, {loggerUtils.get_gpu_name()}, {region}]'
        CACHE_LOG_FILE = f"{write_log['logdir']}/cache_usage.log"
        loggerUtils.write_logger(log_msg_type='cache', logfile=CACHE_LOG_FILE, message=msg_cac_log)

    return({'sequence': one_hot_encode(ref_sequences), 'interval_object': reg_interval})


def find_variants_in_vcf_file(cyvcf2_object, interval_object, samples):

    #n_samples = len(samples)

    # check that samples are in the vcf file
    # if not set(samples).issubset(cyvcf2_object.samples):
    #     raise Exception(f'[ERROR] Fatal. Some samples are not in the VCF file.')
        
    query = f'{interval_object.chrom}:{interval_object.start}-{interval_object.end}'

    variants_dictionary = {}
    variants_dictionary['chr'] = interval_object.chrom
    variants_dictionary['positions'] = tuple(variant.POS for variant in cyvcf2_object(query))
    if not variants_dictionary['positions']:
        return(None)
        
    for i, sample in enumerate(samples):
        try:
            if sample in cyvcf2_object.samples:
                variants_dictionary[sample] = tuple([variant.genotypes[i][0:2], variant.gt_bases[i].split('|')] for variant in cyvcf2_object(query))
        except UserWarning:
            print(f'[WARNING] {sample} is not in the VCF file.')
            continue

    return(variants_dictionary)


def create_mapping_dictionary(variants_array, interval_start, haplotype='both', sequence_source = 'reference', samples=None):

    # sequence_source is needed to do the right check

    import numpy as np

    # never change this whatsoever
    A = np.array([1,0,0,0], dtype=np.float32)
    C = np.array([0,1,0,0], dtype=np.float32)
    G = np.array([0,0,1,0], dtype=np.float32)
    T = np.array([0,0,0,1], dtype=np.float32)
    seq_dict = {'A': A, 'C': C, 'G': G, 'T': T}

    if haplotype == 'hap1':
        haplotype_map = {(variants_array[i][1] - interval_start): variants_array[i][3][0] for i in range(0, len(variants_array))}
        return(haplotype_map)
        #return({(k - interval_start): v for k, v in haplotype_map.items() if v is not None})
    elif haplotype == 'hap2':
        haplotype_map = {(variants_array[i][1] - interval_start): variants_array[i][3][1] for i in range(0, len(variants_array))}
        return(haplotype_map)

    elif haplotype == 'both':
        samples_haplotype_map = {}
        samples_haplotype_map['positions'] = tuple(variants_array['positions'][i] - interval_start for i in range(len(variants_array['positions'])))

        if samples is None:
            samples = list(variants_array.keys())
            samples.remove('chr')
            samples.remove('positions')
        else:
            pass
        
        for i, sample in enumerate(samples):
            samples_haplotype_map[sample] = {}
            samples_haplotype_map[sample]['haplotype1'] = tuple(seq_dict[variants_array[sample][i][1][0]] for i in range(0, len(variants_array[sample])))
            samples_haplotype_map[sample]['haplotype2'] = tuple(seq_dict[variants_array[sample][i][1][1]] for i in range(0, len(variants_array[sample])))

        # check 
        for sample in samples:
            # reference 
            if sequence_source == 'reference':
                condition = len(samples_haplotype_map['positions']) == len(samples_haplotype_map[sample]['haplotype0'])
            elif sequence_source == 'personalized':
                condition = len(samples_haplotype_map['positions']) == len(samples_haplotype_map[sample]['haplotype1']) == len(samples_haplotype_map[sample]['haplotype2'])

            if not condition:
                raise Exception('[ERROR] Fatal. {sample} positions and haplotypes do not match.')
                
        return(samples_haplotype_map)

    
def replace_variants_in_reference_sequence(query_sequences_encoded, mapping_dict, samples=None):

    import copy
    import numpy as np

    variant_encoded = {}

    if samples is None:
        samples = list(mapping_dict.keys())
        samples.remove('positions')
    else:
        pass

    indices = mapping_dict['positions']
    for i, sample in enumerate(samples):
        sample_output = {} # will store haplotype1 and haplotype2
        if len(mapping_dict[sample]) != 2:
            raise Exception(f'[ERROR] Number of haplotypes sequences for {sample} is not equal to 2.')
        else:
            haps = ['haplotype1', 'haplotype2']
            for j in range(len(mapping_dict[sample])):
                ref_i = copy.deepcopy(query_sequences_encoded) # very important

                # check that ref_i is of the right shape
                if ref_i.shape != query_sequences_encoded.shape:
                    raise Exception(f'[ERROR] Fatal. The shape of copied encoded sequence is {ref_i.shape}')

                bases = mapping_dict[sample][haps[j]]

                try:
                    ref_i[:, np.array(indices), : ] = bases
                except IndexError:
                    # this is a hack - I need to find a way around this later
                    try:
                        if 393216 in indices:
                            
                            indexerror_index = indices.index(393216)
                            popped_index = indices.pop(indexerror_index)
                            popped_base = bases.pop(indexerror_index)
                            print(f'[INFO] Encountered an index error for {sample}. Correcting... Removed {popped_index} & {popped_base}')
                    except:
                        raise Exception(f'[ERROR] Fatal. The shape of copied encoded sequence is {ref_i.shape} for {sample}: max {np.max(indices)} min {np.min(indices)}')

                sample_output[haps[j]] = ref_i

                if ref_i.shape != query_sequences_encoded.shape:
                    raise Exception(f'[ERROR] Fatal. {sample} encoded variant sequence shape is {ref_i.shape}, and not {query_sequences_encoded.shape}')
    
        variant_encoded[sample] = sample_output
    return(variant_encoded)



def create_input_for_enformer(region_details, samples, path_to_vcf, fasta_func, hap_type = 'both', resize_for_enformer=True, resize_length=None, write_log=None, sequence_source=None):
    '''
    Given a region in the genome, a reference genome (fasta) and a VCF file, create an individual's sequence for that region

    Parameters:
        query_sequences: str
            The query, perhaps, reference sequence
        mapping_dict: a numpy array or list
            List of the regions that should be modified

    Returns:
        A new sequence with all the changes applied.
    '''

    import numpy as np
    import cyvcf2
    import loggerUtils

    # print('Within create input function')
    # print(f'Regions are: {region_details}')
    # print(f'samples are: {samples}')
    # print(f'path_to_vcf are: {path_to_vcf}')
    # print(f'sequence_source are: {sequence_source}')

    region = region_details['query']
    logtype = region_details['logtype']

    if sequence_source is None:
        raise Exception(f'[ERROR] Fatal. Please, pass in a genome sequence source e.g. personalized, reference, or random.')

    if sequence_source == 'random':
        return({'sequence': {'haplotype0': one_hot_encode(generate_random_sequence_inputs())}, 'metadata': {'sequence_source':'random', 'region':region, 'logtype': logtype}})
    else:
        reference_sequence = extract_reference_sequence(region=region, fasta_func=fasta_func, resize_for_enformer=resize_for_enformer, write_log=write_log, resize_length=resize_length)
        #print(f'Region {region} sequences successfully created within create input function')
        if np.all(reference_sequence['sequence'] == 0.25): # check if all the sequence are "NNNNNNNNNNN..."
            if (write_log is not None) and (write_log['logtypes']['error']):
                err_msg = f'[INPUT ERROR] {region} is invalid; all nucleotides are N.'
                MEMORY_ERROR_FILE = f"{write_log['logdir']}/error_details.log"
                loggerUtils.write_logger(log_msg_type = 'error', logfile = MEMORY_ERROR_FILE, message = err_msg)
            return(None)

        else:
            if sequence_source == 'reference':
                return({'sequence': {'haplotype0': reference_sequence['sequence']}, 'metadata': {'sequence_source':'ref', 'region':region, 'logtype': logtype}})
            elif sequence_source == 'personalized':
                # which vcf file
                vcf_chr = cyvcf2.cyvcf2.VCF(path_to_vcf, samples=samples)
                samples_variants = find_variants_in_vcf_file(cyvcf2_object=vcf_chr, interval_object=reference_sequence['interval_object'], samples=samples)

                if samples_variants: # go on and change the variants by position
                    samples_mapping_dictionary = create_mapping_dictionary(variants_array=samples_variants, interval_start = reference_sequence['interval_object'].start, haplotype=hap_type, sequence_source=sequence_source, samples=None)
                    samples_variants_encoded = replace_variants_in_reference_sequence(query_sequences_encoded = reference_sequence['sequence'], mapping_dict = samples_mapping_dictionary, samples=None)
                    return({'sequence': samples_variants_encoded, 'metadata': {'sequence_source':'var', 'region':region, 'logtype': logtype}})

                else: # return the reference and no need for having the same reference genome for all samples since one is enough 
                    return({'sequence':{'haplotype1': reference_sequence['sequence'], 'haplotype2': reference_sequence['sequence']}, 'metadata': {'sequence_source':'ref', 'region':region, 'logtype': logtype, 'samples': samples}})
