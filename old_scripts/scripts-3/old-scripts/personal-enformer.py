
import parsl
from parsl.app.app import python_app, join_app
import kipoiseq

@python_app
def check_query(sample, query, output_dir, logfile):
    '''
    Checks of predictions are available for an individual or sample
    '''
    import pandas as pd
    import os

    if isinstance(logfile, pd.DataFrame):# should have read it>
        print('Logfile is a dataframe')
        motif_in_logfile = logfile.motif.values
        individual_in_logfile = logfile.individual.values
        query_saved = str(f'{output_dir}/{sample}/{query}_predictions.h5')

        # four conditions 
        a = (query in motif_in_logfile)
        b = (sample in individual_in_logfile)

        if (logfile.motif.empty) | (logfile.individual.empty):
            pass
        else:
            c = (logfile.loc[(logfile.motif==query) & (logfile.individual==sample), 'status'].values[0] == 'completed')
        d = os.path.isfile(query_saved)

        if a and b and c and d:
            pass
        else:
            out = query
    elif isinstance(logfile, type(None)):
        out = query

    return(out)


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

def run_predictions(sequence_encoded, region, sample, seq_type, model, output_dir):  
    #target_sequence = one_hot_encode(sequence)[np.newaxis]
    target_prediction = model_predict(sequence_encoded, model)
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