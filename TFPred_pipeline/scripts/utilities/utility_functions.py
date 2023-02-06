
import pandas as pd
import numpy as np

# def agg_byall(pred_tracks):
    
#     return(agg_by_mean(pred_tracks), agg_by_center(pred_tracks), agg_by_mean(pred_tracks, use_bins=upstream), agg_by_mean(pred_tracks, use_bins=downstream), agg_by_mean(pred_tracks, use_bins=upstream + downstream))

def aggregate_enformer_predictions(each_id, log_data, predictions_path, TF, save_dir, agg_types, batch_num=None):

    import h5py
    import numpy as np
    import os, sys, subprocess, tqdm
    import pandas as pd
    import multiprocessing
    import itertools

    # these are the bins/ positions
    upstream = list(range(0, 7)) # or 0 to 7
    center = 8
    pre_center = 7
    post_center = 9
    mean_center=[7,8,9]
    downstream = list(range(10, 17)) # or 10 to 17

    global read_file

    # can aggregate by the mean of all bins, mean of the upstream and/or downstream alone, or just select the center
    def agg_by_mean(pred_tracks, use_bins=None):

        y = []
        X = []

        for k, v in pred_tracks.items():
            y.append(k) #if k.startswith('pos') else y.append(0)

            if isinstance(use_bins, type(None)):
                v = v.mean(axis=0)
            elif isinstance(use_bins, type([])):
                v = v[use_bins, :].mean(axis=0)
            v = np.expand_dims(v, axis=1).T
            X.append(v)

        y = np.expand_dims(np.array(y), axis=1)
        dt = np.hstack((y, np.vstack(X)))
        #dt = np.vstack(X)

        return dt

    def agg_by_center(pred_tracks, center):

        y = []
        X = []

        for k, v in pred_tracks.items():
            y.append(k)
            v = v[center, :]
            v = np.expand_dims(v, axis=1).T
            X.append(v)

        y = np.expand_dims(np.array(y), axis=1)
        dt = np.hstack((y, np.vstack(X)))

        return dt

    if each_id in ['kawakami', 'cistrome']:
        pred_type = 'ref'
        haplotypes = ['haplotype0']
    elif each_id in ['random']:
        pred_type = 'random'
        haplotypes = ['haplotype0']
    else:
        pred_type = 'var'
        haplotypes  = ['haplotype1', 'haplotype2']

    print(f'[INFO] Seeing {multiprocessing.cpu_count()} CPUs')
    print(f'[INFO] Starting to collect predictions for {each_id}')

    def read_file(motif, dir, haplotype):
        import os
        import h5py
        import numpy as np

        output = {}
        fle = os.path.join(dir, haplotype, f'{motif}_predictions.h5')  #f'{dir}//{motif}_predictions.h5'
        if os.path.isfile(fle):
            with h5py.File(fle, 'r') as f:
                filekey = list(f.keys())[0]
                output[motif] = np.vstack(list(f[filekey]))
        else:
            print(f'[ERROR] {motif} predictions file does not exist.')

        return(output)

    if __name__ == '__main__':

        motifs_list = log_data.loc[log_data['sequence_type'] == pred_type, ].motif.values.tolist()
        print(len(motifs_list))

        pooled_dictionary = {}
        for haplotype in haplotypes:
            pool = multiprocessing.Pool(16)
            outputs_list = pool.starmap(read_file, itertools.product(motifs_list, [predictions_path], [haplotype])) #pool.map(read_file, motifs_list) # 'haplotype1

            predictions =  {k: v for d in outputs_list for k, v in d.items()}
            pooled_dictionary[haplotype] = predictions

        if len(haplotypes) == 2:
            # check that the keys match
            match_condition = sorted(list(pooled_dictionary[haplotypes[0]].keys())) == sorted(list(pooled_dictionary[haplotypes[1]].keys()))
            if not match_condition:
                raise Exception(f'[ERROR] Fatal: Haplotypes 1 and 2 regions are different.Haplotype1 length is {len(list(pooled_dictionary[haplotypes[0]].keys()))} abd Haplotype2 length is {len(list(pooled_dictionary[haplotypes[1]].keys()))}')
            else:
                # one one of them
                ids = list(pooled_dictionary[haplotypes[0]].keys())
                summed_pooled_dictionary = {id: np.add(pooled_dictionary[haplotypes[0]][id], pooled_dictionary[haplotypes[1]][id]) for id in ids}

        elif len(haplotypes) == 1:
            summed_pooled_dictionary = pooled_dictionary[haplotypes[0]]
        print(f'[INFO] Successfully read all files into pooled dictionary.')


        for agg_type in agg_types:
            try:
                data_dict = {}
                if agg_type == 'aggByMean': data_dict[agg_type] = agg_by_mean(summed_pooled_dictionary)
                if agg_type == 'aggByMeanCenter': data_dict[agg_type] = agg_by_mean(summed_pooled_dictionary, use_bins=mean_center)
                if agg_type == 'aggByCenter': data_dict[agg_type] = agg_by_center(summed_pooled_dictionary, center=center)
                if agg_type == 'aggByPreCenter': data_dict[agg_type] = agg_by_center(summed_pooled_dictionary, center=pre_center)
                if agg_type == 'aggByPostCenter': data_dict[agg_type] = agg_by_center(summed_pooled_dictionary, center=post_center)
                if agg_type == 'aggByUpstream': data_dict[agg_type] = agg_by_mean(summed_pooled_dictionary, use_bins=upstream)
                if agg_type == 'aggByDownstream': data_dict[agg_type] = agg_by_mean(summed_pooled_dictionary, use_bins=downstream)
                if agg_type == 'aggByUpstreamDownstream': data_dict[agg_type] = agg_by_mean(summed_pooled_dictionary, use_bins=upstream + downstream)
            
            except ValueError as ve:
                raise ValueError(f'[ERROR] Problem with arrays for {each_id}, {agg_type}')

            #ty = pd.concat([pd.Series(list(predictions.keys())), pd.DataFrame(dt)], axis=1)
            ty = pd.DataFrame(data_dict[agg_type])
            print(f'[INFO] Dimension of collected data is {ty.shape[0]} by {ty.shape[1]}')

            column_names = ['id']
            column_names.extend([f'f_{i}' for i in range(1, ty.shape[1])])

            #print(len(column_names))
            ty = ty.set_axis(column_names, axis=1, inplace=False)
            print(ty.iloc[0:5, 0:5])

            print(f'[INFO] Saving file to {save_dir}/{each_id}_{agg_type}_{TF}.csv.gz')
            if batch_num is None:
                ty.to_csv(path_or_buf=f'{save_dir}/{each_id}_{agg_type}_{TF}.csv.gz', index=False, compression='gzip')
                print(f'[INFO] Finished saving data for {each_id}')
            else:
                ty.to_csv(path_or_buf=f'{save_dir}/{each_id}_{agg_type}_{TF}_batch_{batch_num}.csv.gz', index=False, compression='gzip')
                print(f'[INFO] Finished saving data for {each_id} for batch {batch_num}')
            
    return(0)

def return_prediction_function(use_parsl, fxn=aggregate_enformer_predictions):
    '''
    Decorate or not the `run_batch_predictions` function based on whether `use_parsl` is true or false
    Returns: 
        function object
        The function if parsl is not used
        The parsl decorated function if parsl is meant to be used
    
    '''
    from parsl.app.app import python_app
    if use_parsl == True:
        return python_app(fxn)
    elif use_parsl == False:
        return fxn

def generate_batch(lst, batch_n, len_lst = None):
    """
    Given a list, this function yields batches of an unspecified size but the number of batches is equal to `batch_n`
    E.g. generate_batch([0, 1, 2, 3, 4, 5, 6], batch_n=2) -> (0, 1, 2, 3), (4, 5, 6)
    
    Parameters:
        lst: list
        batch_n: int
            Number of batches to return
        len_lst: None or num (length of the input list)
    Yields
        `batch_n` batches of the list
    """
    import math
    # how many per batch
    if len_lst is not None:
        n_elems = math.ceil(len_lst/batch_n)
    else:
        n_elems = math.ceil(len(lst)/batch_n)
    for i in range(0, len(lst), n_elems):
        yield lst[i:(i + n_elems)]



# def pool_read_files(log_dt = log_data):
    #     import multiprocessing
    #     import pandas as pd

    #     pool = multiprocessing.Pool(64)
    #     motifs_list = log_dt.loc[log_dt['sequence_type'] == pred_type, ].motif.values.tolist()
    #     outputs_list = pool.map(read_file, motifs_list)
    #     return {k: v for d in outputs_list for k, v in d.items()}

    # for dt in tqdm.tqdm(log_data.loc[log_data['sequence_type'] == pred_type, ].motif.values.tolist()):
    #     fle = f'{predictions_path}/{dt}_predictions.h5'
    #     if os.path.isfile(fle):
    #         with h5py.File(fle, 'r') as f:
    #             filekey = list(f.keys())[0]
    #             predictions[dt] = np.vstack(list(f[filekey]))
    #     else:
    #         print(f'[ERROR] {dt} predictions file does not exist.')

    # print(f'[INFO] Finished collecting {len(predictions)} predictions for {each_id}')

    #dt_aggbycenter = agg_by_center(predictions, center=8)
    #data_list = [dt_aggbycenter]


# def get_bed_files(dcids=[1], data_type='TF', cistrome_dir='/projects/covid-ct/imlab/data/cistrome/compressed', data_info=None):

#     '''
#     Extract the bedfiles from each ...
    
#     Returns:
#         A dictionary
#     '''
    
#     bedfiles = {}

#     # check if the file exists in the data_info : do this later

#     if data_type not in ['TF', 'chromatin', 'histone']:
#         print('Data type must be of `TF`, `chromatin`, or `histone`.')
#         exit()
    
#     if data_type == 'TF':
#         dt_zipfile = zipfile.ZipFile(f'{cistrome_dir}/human_fa.zip', mode = 'r')
#         for dcid in dcids:
#             try:
#                 bedfiles[dcid] = pd.read_table(dt_zipfile.open(f'human_factor/{dcid}_sort_peaks.narrowPeak.bed.gz'), compression='gzip', header=None)
#             except:
#                 bedfiles[dcid] = 'File not found.'
            
#     elif data_type == 'chromatin':
#         dt_zipfile = zipfile.ZipFile(f'{cistrome_dir}/human_ca.zip', mode = 'r')
#         for dcid in dcids:
#             try:
#                 bedfiles[dcid] = pd.read_table(dt_zipfile.open(f'human_ca/{dcid}_sort_peaks.narrowPeak.bed.gz'), compression='gzip', header=None)
#             except:
#                 bedfiles[dcid] = 'File not found.'
            
#     elif data_type == 'histone':
#         dt_zipfile = zipfile.ZipFile(f'{cistrome_dir}/human_hm.zip', mode = 'r')
#         for dcid in dcids:
#             try:
#                 bedfiles[dcid] = pd.read_table(dt_zipfile.open(f'human_hm/{dcid}_sort_peaks.narrowPeak.bed.gz'), compression='gzip', header=None)
#             except:
#                 bedfiles[dcid] = 'File not found.'
    
#     # close the zipfile
#     dt_zipfile.close()
            
#     return bedfiles
# def select_tracks(list_or_np_array, select_tracks=None, select_bins=None, motif_bin=None):

#     # make sure that only one is evaluated in select_tracks and remove_tracks
#     # I think it is many times faster if it removes 2 tracks vs if it selects 5311 tracks out of 5313
#     # or if it 

#     if isinstance(list_or_np_array, type([])):
#         for i, np_array in enumerate(list_or_np_array):
#             if isinstance(select_bins, type(8)):
#                 temp = np_array[:, select_tracks]
#                 list_or_np_array[i] = temp[range(motif_bin - select_bins, (motif_bin + select_bins + 1)), : ]
#             else:
#                 list_or_np_array[i] = np_array[:, select_tracks]
            
#     if isinstance(list_or_np_array, type(np.empty((2, 2)))):
#         if isinstance(select_bins, type(8)):
#             temp = list_or_np_array[:, select_tracks]
#             list_or_np_array = temp[range(motif_bin - select_bins, (motif_bin + select_bins + 1)), : ].squeeze()
#         else:
#             list_or_np_array = list_or_np_array[:, select_tracks].squeeze()

#     return list_or_np_array




   #test_aggbymean, test_aggbycenter, test_aggbymean_upstream, test_aggbymean_downstream, test_aggbymean_upstream_downstream = agg_byall(predictions)
    #data_list = [test_aggbymean, test_aggbycenter, test_aggbymean_upstream, test_aggbymean_downstream, test_aggbymean_upstream_downstream]

            #for i, agg_type in enumerate(data_dict):


# # ===== 
# def faster_parse_vcf():
#     sample_names = sys.argv[1]
#     chr_names = list(range(1, 23)) + ['Y']
#     #print(chr_names)

#     vcf_folder = '/projects/covid-ct/imlab/data/GEUVADIS/vcf_files'
#     output_folder = '/projects/covid-ct/imlab/users/temi/projects/TFXcan/output'

#     # header: #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT

#     with open(sample_names, 'r') as sn:
#         for line in sn:
#             for chrom in chr_names:
#                 # each line is an individual
#                 query = line.split()
                
#                 # the regions for this individual
#                 regions = query[1:]
                
#                 # filter for those regions in that chromosome >> a boolean list is return
#                 filt_bool = [True if str(chrom) == re.split(':', reg)[0] else False for reg in regions]
                
#                 if not any(filt_bool):
#                     print(f'\nchr{chrom} not present. Moving on === \n')
#                     continue
                
#                 print(f'\nchr{chrom} present. Extracting variant information === \n')
                
#                 chrom_regions = [x for x, y in zip(regions, filt_bool) if y == True]
#                 chrom_regions = ','.join(chrom_regions)
            
#             # create the sh call
#                 call_ = str(f'~/miniconda3/bin/bcftools view -v snps -H -r {chrom_regions} -s {query[0]} {vcf_folder}/ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz > {output_folder}/{query[0]}_chr{chrom}_subsetted_vcf_file.vcf')
                
#                 #call_ = str(f'~/miniconda3/bin/bcftools view -H -r {chr_name}:{",".join(query[1:])} -s {query[0]} {vcf_folder}/ALL.chr{chr_name}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz > {output_folder}/{query[0]}_subsetted_vcf_file.vcf')
            
#                 #print(call_)
                
#                 # call the shell script
#                 subprocess.run(call_, shell=True)
#                 subprocess.run(str(f'gzip -f {output_folder}/{query[0]}_chr{chrom}_subsetted_vcf_file.vcf'), shell=True) # removed the txt file
        

