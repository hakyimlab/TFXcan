
import pandas as pd
import numpy as np

def localParslConfig(params):
    from tqdm import tqdm
    from parsl import python_app
    import parsl
    import json

    # Make a config that runs on two nodes
    from parsl.executors import HighThroughputExecutor
    from parsl.providers import LocalProvider
    from parsl.config import Config
    from parsl.channels import LocalChannel
    from parsl.launchers import MpiExecLauncher

    import os

    workingdir = params['working_dir'] #'/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_predict'
    rundir = os.path.join(workingdir, 'runinfo')
    #rundir = '/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_predict/runinfo'

    #parsl.clear()

    local_htex = Config(
        executors=[
            HighThroughputExecutor(
                label="htex_Local",
                max_workers=8, # vs max_workers
                available_accelerators=8,
                worker_debug=True,
                cores_per_worker=1, # how many cores per worker #nodes_per_block, 2 is usually enough or 1.
                working_dir=workingdir,
                provider=LocalProvider(
                    channel=LocalChannel(),
                    init_blocks=1,
                    max_blocks=1,
                    launcher=MpiExecLauncher(),
                    worker_init='source ~/.bashrc; conda activate dl-tools; which python; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/temi/miniconda3/envs/dl-tools/lib'
                ),
            )
        ],
        strategy=None,
        run_dir=rundir
    )

    return(local_htex)

def localParslConfig_threadpool(params):
 
    import parsl
    # Make a config that runs on two nodes
    from parsl.executors import ThreadPoolExecutor
    from parsl.config import Config

    import os
    workingdir = params['working_dir']
    rundir = os.path.join(workingdir, 'runinfo')

    parsl.clear()

    local_tpex = Config(
        executors=[
            ThreadPoolExecutor(
                label="tpex_Local",
                max_threads=16,
                thread_name_prefix='tpex_run',
                working_dir=workingdir,
            )
        ],
        strategy=None,
        run_dir=rundir
    )

    return(local_tpex)

# these are the bins
upstream = list(range(0, 8))
center = [8]
downstream = list(range(9, 17))

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

def agg_by_center(pred_tracks, center=8):

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

def agg_byall(pred_tracks, center=8):
    
    return(agg_by_mean(pred_tracks), agg_by_center(pred_tracks), agg_by_mean(pred_tracks, use_bins=upstream), agg_by_mean(pred_tracks, use_bins=downstream), agg_by_mean(pred_tracks, use_bins=upstream + downstream))

def collect_modeling_data_for_kawakami(each_id, log_data, predictions_path, TF, base_path, save_dir, agg_types):

    import h5py
    import numpy as np
    import os
    import pandas as pd
    # read in one of the files

    #exec(open(f'{base_path}/modeling_pipeline/scripts/collect_model_data/utility-functions.py').read(), globals(), globals())
    # bpath = os.path.join(base_path, 'modeling_pipeline')
    # localParslConfig_threadpool({'working_dir': bpath})

    kawakami_predictions = {}

    for dt in tqdm.tqdm(log_data.loc[log_data['sequence_type'] == 'ref', ].motif.values.tolist()):
        fle = f'{predictions_path}/{dt}_predictions.h5'
        if os.path.isfile(fle):
            with h5py.File(fle, 'r') as f:
                filekey = list(f.keys())[0]
                kawakami_predictions[dt] = np.vstack(list(f[filekey]))
        else:
            print(f'[ERROR] {dt} predictions file does not exist.')

    print(f'[INFO] Finished collecting {len(kawakami_predictions)} predictions for {each_id}')

    #dt_aggbycenter = agg_by_center(kawakami_predictions, center=8)
    #data_list = [dt_aggbycenter]

    data_dict = {}
    for agg_type in agg_types:
        if agg_type == 'aggByMean': data_dict[agg_type] = agg_by_mean(kawakami_predictions)
        if agg_type == 'aggByCenter': data_dict[agg_type] = agg_by_center(kawakami_predictions)
        if agg_type == 'aggByUpstream': data_dict[agg_type] = agg_by_mean(kawakami_predictions, use_bins=upstream)
        if agg_type == 'aggByDownstream': data_dict[agg_type] = agg_by_mean(kawakami_predictions, use_bins=downstream)
        if agg_type == 'aggByUpstreamDownstream': data_dict[agg_type] = agg_by_mean(kawakami_predictions, use_bins=upstream + downstream)

    #test_aggbymean, test_aggbycenter, test_aggbymean_upstream, test_aggbymean_downstream, test_aggbymean_upstream_downstream = agg_byall(kawakami_predictions)
    #data_list = [test_aggbymean, test_aggbycenter, test_aggbymean_upstream, test_aggbymean_downstream, test_aggbymean_upstream_downstream]

    for i, agg_type in enumerate(data_dict):

        #ty = pd.concat([pd.Series(list(kawakami_predictions.keys())), pd.DataFrame(dt)], axis=1)
        ty = pd.DataFrame(data_dict[agg_type])
        print(f'[INFO] Dimension of collected data is {ty.shape[0]} by {ty.shape[1]}')

        column_names = ['id']
        column_names.extend([f'f_{i}' for i in range(1, ty.shape[1])])

        #print(len(column_names))
        ty = ty.set_axis(column_names, axis=1, inplace=False)
        print(ty.iloc[0:5, 0:5])

        ty.to_csv(path_or_buf=f'{save_dir}/{each_id}_{agg_type}_{TF}.csv.gz', index=False, compression='gzip')
    print(f'[INFO] Finished saving data for {each_id}')

    return(0)





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
        

