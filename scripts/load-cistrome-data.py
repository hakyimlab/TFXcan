

<<<<<<< HEAD
import zipfile, os, re
import pandas as pd
from itertools import compress

# %%
=======
import zipfile
import pandas as pd

a = 2+2
print(a)

>>>>>>> e15d2c58b1ffd0f4fd0508b4b72f37aad9dbebf7
def get_bed_file(dcids=[1], data_type='TF', cistrome_dir='/projects/covid-ct/imlab/data/cistrome/compressed', data_info=None):

    '''
    Extract the bedfiles from each ...
    
    Returns:
        A dictionary
    '''
    
    bedfiles = {}

    #bedcolumns = ['chrom']

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
<<<<<<< HEAD

def read_kawakami_data(TF, kawakami_dir='../data/kawakami-human/detailed_info/', read_data=True):

    # check if the data exists
    mask = [bool(re.match(f'^({TF}).+', each)) for each in os.listdir(kawakami_dir)]

    if not any(mask):
        return(f'{TF} information not available.')
    else:
        tf_data = list(compress(os.listdir(kawakami_dir), mask))
        if read_data == True:
            #mask = [each.startswith(TF) for each in os.listdir(kawakami_dir)]
            tf_dataframe = [pd.read_table(kawakami_dir + which_bed) for which_bed in tf_data]
            tf_dt = ['_'.join(info.split('_')[0:2]) for info in tf_data]
            out = dict(zip(tf_dt, tf_dataframe))

            return(out)
        else:
            return(tf_data)
=======
# %%
>>>>>>> e15d2c58b1ffd0f4fd0508b4b72f37aad9dbebf7
