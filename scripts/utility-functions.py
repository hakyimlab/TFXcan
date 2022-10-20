

import zipfile, gzip
import pandas as pd

def get_bed_files(dcids=[1], data_type='TF', cistrome_dir='/projects/covid-ct/imlab/data/cistrome/compressed', data_info=None):

    '''
    Extract the bedfiles from each ...
    
    Returns:
        A dictionary
    '''
    
    bedfiles = {}

    # check if the file exists in the data_info : do this later

    if data_type not in ['TF', 'chromatin', 'histone']:
        print('Data type must be of `TF`, `chromatin`, or `histone`.')
        exit()
    
    if data_type == 'TF':
        dt_zipfile = zipfile.ZipFile(f'{cistrome_dir}/human_fa.zip', mode = 'r')
        for dcid in dcids:
            try:
                bedfiles[dcid] = pd.read_table(dt_zipfile.open(f'human_factor/{dcid}_sort_peaks.narrowPeak.bed.gz'), compression='gzip', header=None)
            except:
                bedfiles[dcid] = 'File not found.'
            
    elif data_type == 'chromatin':
        dt_zipfile = zipfile.ZipFile(f'{cistrome_dir}/human_ca.zip', mode = 'r')
        for dcid in dcids:
            try:
                bedfiles[dcid] = pd.read_table(dt_zipfile.open(f'human_ca/{dcid}_sort_peaks.narrowPeak.bed.gz'), compression='gzip', header=None)
            except:
                bedfiles[dcid] = 'File not found.'
            
    elif data_type == 'histone':
        dt_zipfile = zipfile.ZipFile(f'{cistrome_dir}/human_hm.zip', mode = 'r')
        for dcid in dcids:
            try:
                bedfiles[dcid] = pd.read_table(dt_zipfile.open(f'human_hm/{dcid}_sort_peaks.narrowPeak.bed.gz'), compression='gzip', header=None)
            except:
                bedfiles[dcid] = 'File not found.'
    
    # close the zipfile
    dt_zipfile.close()
            
    return bedfiles
