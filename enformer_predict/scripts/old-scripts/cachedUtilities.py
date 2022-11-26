from functools import lru_cache

@lru_cache(1)
def get_fastaExtractor():

    import sys, json, os
    #sys.path.append(f'{script_path}/utilities')
    whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)
    script_path = os.path.abspath(whereis_script)
    fpath = os.path.join(script_path, 'utilities')
    sys.path.append(fpath)

    with open(f'{script_path}/../../metadata/enformer_parameters.json') as f:
        parameters = json.load(f)
        fasta_file = parameters['hg38_fasta_file']

    import enformerUsageCodes

    fasta_extractor = enformerUsageCodes.FastaStringExtractor(fasta_file)
    return fasta_extractor

@lru_cache(1)
def get_model():

    import sys, json, os
    #sys.path.append(f'{script_path}/utilities')
    whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)
    script_path = os.path.abspath(whereis_script)
    fpath = os.path.join(script_path, 'utilities')
    sys.path.append(fpath)

    with open(f'{script_path}/../../metadata/enformer_parameters.json') as f:
        parameters = json.load(f)
        model_path = parameters['model_path']

    import tensorflow as tf
    
    return tf.saved_model.load(model_path).model

def print_cache_status(m):
    return(f'[CACHE NORMAL INFO] (get_model) {m.cache_info()}')