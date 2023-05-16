

write_log = {'memory':False, 'error':False, 'time':False, 'cache':False}

def generate_random_sequence_inputs(size=393216):
    import numpy as np
    r_seq_list = np.random.choice(['A', 'G', 'T', 'C'], size)
    return(''.join(r_seq_list))


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


def one_hot_encode(sequence, shap=None):
    import kipoiseq
    import tensorflow as tf
    import numpy as np

    try:
        sequence_encoded = kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)[np.newaxis]
        return(sequence_encoded)

    except ValueError as ve:
        print(f'[INPUT ERROR] {sequence} {len(shap)}')
        print(f'[INPUT ERROR] {sequence} {type(shap)}')

def extract_reference_sequence(region, fasta_func=None, resize_for_enformer=True, resize_length=None, print_sequence=False, write_log=write_log):

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

    try:
        region_split = region.split('_')
        region_chr = region_split[0]
        region_start = int(region_split[1])
        region_end = int(region_split[2])
    except ValueError:
        if write_log['error']:
            err_msg = f'[REGION ERROR] {region} input start or end is invalid.'
            MEMORY_ERROR_FILE = f"{log_dir}/error_details.log"
            setup_logger('error_log', MEMORY_ERROR_FILE)
            logger(err_msg, 'error', 'run_error')
        return(None)


    if fasta_func is None:
        fasta_object = get_fastaExtractor()
    else:
        fasta_object=fasta_func

    if resize_for_enformer == True:
        reg_interval = kipoiseq.Interval(region_chr, region_start, region_end).resize(SEQUENCE_LENGTH)
    else:
        reg_interval = kipoiseq.Interval(region_chr, region_start, region_end).resize(resize_length)

    # extract the sequence 
    ref_sequences = fasta_object.extract(interval=reg_interval, anchor=[])

    if print_sequence == True:
        print(ref_sequences)

    if write_log['cache'] == True:
        msg_cac_log = f'[CACHE INFO] (fasta) [{get_fastaExtractor.cache_info()}, {get_gpu_name()}, {region}]'
        CACHE_LOG_FILE = f"{log_dir}/cache_usage.log"
        setup_logger('cache_log', CACHE_LOG_FILE)
        logger(msg_cac_log, 'info', 'cache')

    return({'sequence': one_hot_encode(ref_sequences), 'interval_object': reg_interval})
    
def find_variants_in_vcf_file(cyvcf2_object, interval_object):

    '''
    Find variants for a region in a vcf file

    Parameters:
        cyvcf2_object: cyvcf2_object
            A cyvcf2_object object
        interval_object: Interval object
            An interval object returned by kipoiseq.Interval

    Returns:
        A list of variants for that region in the form [chr, pos, [genotypes], [genotype bases]]
        An empty string if no variant is present
    
    '''
    query = f'{interval_object.chrom}:{interval_object.start}-{interval_object.end}'
    return([[variant.CHROM, variant.POS, variant.genotypes[0][0:2], variant.gt_bases[0].split('|')] for variant in cyvcf2_object(query)])

def create_mapping_dictionary(variants_array, interval_start, haplotype='both'):
    '''
    Create a map of which variants should be switched 

    Parameters:
        variants_array: a list/numpy array
            A cyvcf2_object object
        interval_start: Interval object start position
            The start position of the interval object returned by kipoiseq.Interval
        haplotype: str; default: 'hap1'; options: hap2, both
            Should haplotype 1 or 2 or both be used?

    Returns: dict
        {position: new genotype base}
    '''
    import numpy as np

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
        hap1 = {(variants_array[i][1] - interval_start): seq_dict[variants_array[i][3][0]] for i in range(0, len(variants_array))}
        hap2 = {(variants_array[i][1] - interval_start): seq_dict[variants_array[i][3][1]] for i in range(0, len(variants_array))}
        haplotype_map = {'haplotype1': hap1, 'haplotype2': hap2}
        return(haplotype_map)

def replace_variants_in_reference_sequence(query_sequences_encoded, mapping_dict):
    '''
    Using the mapping dictionary, mutate many variants in a given sequence 

    Parameters:
        query_sequences: str
            The query, perhaps, reference sequence
        mapping_dict: a numpy array or list
            List of the regions that should be modified

    Returns:
        A new sequence with all the changes applied.
    '''

    import copy
    import numpy as np

    output = {}

    if len(mapping_dict) == 2:
        haps = ['haplotype1', 'haplotype2']
        for j in range(len(mapping_dict)):
            ref_i = copy.deepcopy(query_sequences_encoded) # very important

            #(ref_i.shape)

            indices = list(mapping_dict[haps[j]].keys())
            bases = list(mapping_dict[haps[j]].values())

            ref_i[:, np.array(indices), : ] = bases
            output[haps[j]] = ref_i

    # if len(mapping_dict) == 2:
    #     haps = ['haplotype1', 'haplotype2']
    #     for j in range(len(mapping_dict)):
    #         sequence_list = list(query_sequences)
    #         a = map(lambda i: mapping_dict[haps[j]].get(i, sequence_list[i]), range(len(sequence_list)))
    #         output[haps[j]] = ''.join(list(a))

    
    return(output)



def create_input_for_enformer(region_details, fasta_func, hap_type = 'both', vcf_func=None, resize_for_enformer=True, resize_length=None, write_log=None, sequence_source=None):
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

    #print(region_details)
    region = region_details['query']
    logtype = region_details['logtype']

    if sequence_source == 'random':
        return({'sequence': {'haplotype0': one_hot_encode(generate_random_sequence_inputs())}, 'sequence_source':'random', 'region':region, 'logtype': logtype})
    else:
        reference_sequence = extract_reference_sequence(region, fasta_func, resize_for_enformer)
        if np.all(reference_sequence['sequence'] == 0.25): # check if all the sequence are "NNNNNNNNNNN..."
            if write_log['error']:
                err_msg = f'[INPUT ERROR] {region} is invalid; all nucleotides are N.'
                MEMORY_ERROR_FILE = f"{log_dir}/error_details.log"
                setup_logger('error_log', MEMORY_ERROR_FILE)
                logger(err_msg, 'error', 'run_error')

            return(None)

        else:
            if sequence_source == 'reference':
                return({'sequence': {'haplotype0': reference_sequence['sequence']}, 'sequence_source':'ref', 'region':region, 'logtype': logtype})
            elif sequence_source == 'personalized':
                vcf_object = vcf_func
                variants_coordinates = find_variants_in_vcf_file(vcf_object, reference_sequence['interval_object'])
                
                if variants_coordinates: # go on and change the variants by position
                    mapping_dictionary = create_mapping_dictionary(variants_coordinates, reference_sequence['interval_object'].start, haplotype=hap_type)

                    variant_sequence = replace_variants_in_reference_sequence(reference_sequence['sequence'], mapping_dictionary)

                    return({'sequence':{'haplotype1': variant_sequence['haplotype1'], 
                                        'haplotype2': variant_sequence['haplotype2']}, 
                        'sequence_source':'var', 'region':region, 'logtype': logtype})
                else: # return the reference
                    return({'sequence':{'haplotype1': reference_sequence['sequence'], 
                                        'haplotype2': reference_sequence['sequence']}, 
                        'sequence_source':'ref', 'region':region, 'logtype': logtype})