

import zipfile, gzip
import pandas as pd
import numpy as np
import re

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


# these are the bins
upstream = list(range(0, 8))
center = [8]
downstream = list(range(9, 17))

# can aggregate by the mean of all bins, mean of the upstream and/or downstream alone, or just select the center
def agg_by_mean(pred_tracks, use_bins=None):

    y = []
    X = []

    for k, v in pred_tracks.items():
        y.append(1) if k.startswith('pos') else y.append(0)

        if isinstance(use_bins, type(None)):
            v = v.mean(axis=0)
        elif isinstance(use_bins, type([])):
            v = v[use_bins, :].mean(axis=0)
        v = np.expand_dims(v, axis=1).T
        X.append(v)

    y = np.expand_dims(np.array(y), axis=1)
    dt = np.hstack((y, np.vstack(X)))

    return dt

def agg_by_center(pred_tracks, center=8):

    y = []
    X = []

    for k, v in pred_tracks.items():
        y.append(1) if k.startswith('pos') else y.append(0)
        v = v[center, :]
        v = np.expand_dims(v, axis=1).T
        X.append(v)

    y = np.expand_dims(np.array(y), axis=1)
    dt = np.hstack((y, np.vstack(X)))

    return dt


def agg_byall(pred_tracks, center=8):
    
    return(agg_by_mean(pred_tracks), agg_by_center(pred_tracks), agg_by_mean(pred_tracks, use_bins=upstream), agg_by_mean(pred_tracks, use_bins=downstream), agg_by_mean(pred_tracks, use_bins=upstream + downstream))

def select_tracks(list_or_np_array, select_tracks=None, select_bins=None, motif_bin=None):

    # make sure that only one is evaluated in select_tracks and remove_tracks
    # I think it is many times faster if it removes 2 tracks vs if it selects 5311 tracks out of 5313
    # or if it 

    if isinstance(list_or_np_array, type([])):
        for i, np_array in enumerate(list_or_np_array):
            if isinstance(select_bins, type(8)):
                temp = np_array[:, select_tracks]
                list_or_np_array[i] = temp[range(motif_bin - select_bins, (motif_bin + select_bins + 1)), : ]
            else:
                list_or_np_array[i] = np_array[:, select_tracks]
            
    if isinstance(list_or_np_array, type(np.empty((2, 2)))):
        if isinstance(select_bins, type(8)):
            temp = list_or_np_array[:, select_tracks]
            list_or_np_array = temp[range(motif_bin - select_bins, (motif_bin + select_bins + 1)), : ].squeeze()
        else:
            list_or_np_array = list_or_np_array[:, select_tracks].squeeze()

    return list_or_np_array







# ===== 
def faster_parse_vcf():
    sample_names = sys.argv[1]
    chr_names = list(range(1, 23)) + ['Y']
    #print(chr_names)

    vcf_folder = '/projects/covid-ct/imlab/data/GEUVADIS/vcf_files'
    output_folder = '/projects/covid-ct/imlab/users/temi/projects/TFXcan/output'

    # header: #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT

    with open(sample_names, 'r') as sn:
        for line in sn:
            for chrom in chr_names:
                # each line is an individual
                query = line.split()
                
                # the regions for this individual
                regions = query[1:]
                
                # filter for those regions in that chromosome >> a boolean list is return
                filt_bool = [True if str(chrom) == re.split(':', reg)[0] else False for reg in regions]
                
                if not any(filt_bool):
                    print(f'\nchr{chrom} not present. Moving on === \n')
                    continue
                
                print(f'\nchr{chrom} present. Extracting variant information === \n')
                
                chrom_regions = [x for x, y in zip(regions, filt_bool) if y == True]
                chrom_regions = ','.join(chrom_regions)
            
            # create the sh call
                call_ = str(f'~/miniconda3/bin/bcftools view -v snps -H -r {chrom_regions} -s {query[0]} {vcf_folder}/ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz > {output_folder}/{query[0]}_chr{chrom}_subsetted_vcf_file.vcf')
                
                #call_ = str(f'~/miniconda3/bin/bcftools view -H -r {chr_name}:{",".join(query[1:])} -s {query[0]} {vcf_folder}/ALL.chr{chr_name}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz > {output_folder}/{query[0]}_subsetted_vcf_file.vcf')
            
                #print(call_)
                
                # call the shell script
                subprocess.run(call_, shell=True)
                subprocess.run(str(f'gzip -f {output_folder}/{query[0]}_chr{chrom}_subsetted_vcf_file.vcf'), shell=True) # removed the txt file
        

