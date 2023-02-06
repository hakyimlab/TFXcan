# Usage: This script is used to predict on batches using ENFORMER on individuals' regions
# Author: Temi
# Date: Wed 25 Jan 2023

from __future__ import absolute_import, division, print_function, unicode_literals
import os, sys, json
import pandas as pd # for manipulating dataframes
import time, tqdm
# import parsl
# from parsl.app.app import python_app
from functools import lru_cache
import glob
from datetime import date
import argparse

# needed arguments 
parser = argparse.ArgumentParser()
parser.add_argument("--param_config", help="Path to JSON file of parameters and directives to be used by ENFORMER", type=str)
args = parser.parse_args()

# some locations and folders
whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)
script_path = os.path.abspath(whereis_script)
batch_utils_path = os.path.join(script_path, 'batch_utils')
sys.path.append(batch_utils_path)

import check_input_parameters



# main 
def main():

    global use_parsl
    params_path = args.param_config

    with open(f'{params_path}') as f:

        parameters = json.load(f)

        if parameters['sub_dir'] == True:
            project_dir = os.path.join(parameters['project_dir'], 'prediction_folder')
        else:
            project_dir = parameters['project_dir']

        interval_list_file = parameters['interval_list_file']
        TF = parameters['TF']
        predictions_log_dir = os.path.join(project_dir, parameters['predictions_log_dir'])
        log_dir = os.path.join(project_dir, parameters['log_dir'])
        batch_size = int(parameters['batch_size'])
        n_individuals = int(parameters['n_individuals'])
        use_parsl = parameters['use_parsl']
        n_regions = parameters["predict_on_n_regions"]
        parsl_parameters = parameters['parsl_parameters']
        sequence_source = parameters['sequence_source']
        prediction_data_name = parameters['prediction_data_name']
        create_hdf5_file = parameters["create_hdf5_file"]
        run_date = parameters['date'] if parameters['date'] is not None else date.today().strftime("%Y-%m-%d")

        metadata_dir = parameters['metadata_dir']
        if not os.path.isdir(metadata_dir):
            os.makedirs(metadata_dir)

        output_dir = os.path.join(project_dir, parameters['output_dir'], parameters['prediction_data_name'] + '_' + parameters['TF'], 'predictions_' + run_date)

        if int(n_regions) == -1:
            predict_on_n_regions = None
        elif int(n_regions) > 0:
            predict_on_n_regions = (n_regions) if isinstance(n_regions, int) else None

        # personalized parameters 
        individuals = parameters['individuals'] if sequence_source == 'personalized' else None
        #vcf_file = parameters['vcf_file']if sequence_source == 'personalized' else None
        vcf_files_dict = parameters['vcf_files'] if sequence_source == 'personalized' else None

        # personalized parameters
        if sequence_source == 'personalized':
             # use only the chromosomes that have been made available in the config file vcf params
            chromosomes = list(vcf_files_dict['files'].keys())
        # list of chromosomes (if the sequence source is reference)
        elif sequence_source == 'reference':
            chromosomes = [f'chr{i}' for i in range(1, 23)]
            chromosomes.extend(['chrX'])

    # write the params_path to a config.json file in a predefined folder
    config_data = {'params_path': params_path}
    with open(f'{batch_utils_path}/config.json', mode='w') as cj:
        json.dump(config_data, cj)

    # modify parsl parameters to add the working directory
    parsl_parameters['working_dir'] = project_dir

    # create necessary folders if not available
    predictions_log_dir = f'{predictions_log_dir}/{prediction_data_name}_{TF}'
    if not os.path.isdir(predictions_log_dir):
        os.makedirs(predictions_log_dir)

    if not os.path.isdir(log_dir):
        os.makedirs(log_dir)

    if use_parsl == True:
        print(f'[INFO] Using parsl configuration: {use_parsl}')
        import parslConfiguration
        parsl.load(parslConfiguration.polaris_htParslConfig(params=parsl_parameters))

    predict_utils_one = f'{script_path}/batch_utils/predictUtils_one.py'
    exec(open(predict_utils_one).read(), globals(), globals())

    # decorate the prediction function with or without parsl
    prediction_fxn = return_prediction_function(use_parsl)

    # determine what individuals to predict on and all that
    if sequence_source == 'personalized':
        if isinstance(individuals, list):
            id_list = individuals
            pass
        elif isinstance(individuals, type('str')):
            if os.path.isfile(individuals):
                if n_individuals == -1:
                    id_list = pd.read_table(individuals, header=None)[0].tolist()[0:]
                elif n_individuals > 0:
                    id_list = pd.read_table(individuals, header=None)[0].tolist()[0:(n_individuals)]
            else:
                id_list = [individuals]
        #print(f'[INFO] Predicting for these individuals: {id_list}')
    elif sequence_source == 'reference':
        id_list = [prediction_data_name]
        print(f'[INFO] Predicting on a reference set')
    elif sequence_source == 'random':
        id_list = [prediction_data_name]
        print(f'[INFO] Predicting on a randomly generated set')

    # I need a log file ==> for personalized predictions, all logs should be in the same file
    # read in the log file for this individual ; doing this so that the log file is not opened everytime
    logfile_csv = f'{predictions_log_dir}/predictions_log_{run_date}.csv'
    # better to read the file in and pass it around than to just pass the file path adn read it everytime
    logfile = pd.read_csv(logfile_csv) if os.path.isfile(logfile_csv) else None
        
    # list of intervals to be predicted on
    a = pd.read_table(interval_list_file, sep=' ', header=None).dropna(axis=0)
    #print(a.head())
    list_of_regions = a[0].tolist()[0:(predict_on_n_regions)] # a list of queries
    #print(list_of_regions[0:10])

    # filter the list of chromosomes to be compatible with the available regions
    chromosomes = list(set([r.split('_')[0] for r in list_of_regions]))
    print(f'[INFO] Available chromosomes are: {chromosomes}')

    #chr_list_of_regions = [r for r in list_of_regions if r.startswith(f"{chromosome}_")]
    for each_id in id_list:
        app_futures = [] # collect futures here
        
        # filter where each_id is in the log file
        if not logfile is None:
            id_logfile = logfile.loc[logfile['individual'] == each_id, : ]
        elif logfile is None:
            id_logfile = logfile
        
        if sequence_source == 'personalized':
            if not os.path.exists(f'{output_dir}/{each_id}'):
                print(f'[INFO] Creating output directory at {output_dir}/{each_id}')
                os.makedirs(f'{output_dir}/{each_id}/haplotype1')
                os.makedirs(f'{output_dir}/{each_id}/haplotype2')
        elif sequence_source == 'reference':
            if not os.path.exists(f'{output_dir}/{each_id}'):
                print(f'[INFO] Creating output directory at {output_dir}/{each_id}')
                os.makedirs(f'{output_dir}/{each_id}/haplotype0')
        elif sequence_source == 'random':
            if not os.path.exists(f'{output_dir}/{each_id}'):
                print(f'[INFO] Creating output directory at {output_dir}/{each_id}')
                os.makedirs(f'{output_dir}/{each_id}/haplotype0')

        print(f'[INFO] Collecting appfutures for {each_id}')
        for chromosome in chromosomes:
            chr_list_of_regions = [r for r in list_of_regions if r.startswith(f"{chromosome}_")]

            if not chr_list_of_regions:
                continue

            chr_vcf_file = os.path.join(vcf_files_dict['folder'], vcf_files_dict['files'][chromosome])
            if sequence_source == 'personalized':
                def make_cyvcf_object(vcf_file=chr_vcf_file, sample=each_id):
                    import cyvcf2
                    return(cyvcf2.cyvcf2.VCF(vcf_file, samples=sample))
            elif sequence_source == 'reference':
                make_cyvcf_object = None
            elif sequence_source == 'random':
                make_cyvcf_object = None

            #print(f'[INFO] Collecting parsl appfuture batches for {chromosome}: {len(chr_list_of_regions)}')
            #for batch_query in batches:

            batches = generate_batch(chr_list_of_regions, batch_n=batch_size)
            count = 0
            #app_futures = []
            for batch_query in batches:
                count = count + 1
                app_futures.append(prediction_fxn(batch_regions=batch_query, batch_num = count, id=each_id, vcf_func=make_cyvcf_object, script_path=script_path, output_dir=output_dir, logfile=id_logfile, predictions_log_file=logfile_csv))

        if use_parsl == True:
            print(f'[INFO] Executing parsl futures for {len(app_futures)} batches')
            exec_futures = [q.result() for q in app_futures] 
            print(f'[INFO] Finished predictions for {each_id}: {exec_futures[0:10]} ...\n')
        elif use_parsl == False:
            print(f'[INFO] Finished predictions for {each_id}: {app_futures} ...\n')

    print(f'[INFO] Finished predicting for all inputs.')
    
    print(f'[INFO] Success for all.')
        
    # == After predictions are complete, a json file will be written out to help with aggregation
    if sequence_source == 'reference':
        print(f'[INFO] Writing `aggregation_config_{prediction_data_name}_{TF}.json` file to {metadata_dir}')

        agg_dt = {'predictions_folder': project_dir, 'enformer_prediction_path': f'{output_dir}', 'log_file':logfile_csv, 'prediction_data_name':prediction_data_name, 'sequence_source': sequence_source, 'run_date':run_date, 'transcription_factor':TF, 'individuals':None, 'n_individuals':n_individuals}

        with(open(f'{metadata_dir}/aggregation_config_{prediction_data_name}_{TF}.json', mode='w')) as wj:
            json.dump(agg_dt, wj) 

    elif sequence_source == 'personalized':
        print(f'[INFO] Writing `aggregation_config_{prediction_data_name}_{TF}.json` file to {metadata_dir}')

        agg_dt = {'predictions_folder': project_dir, 'enformer_prediction_path': f'{output_dir}', 'log_file':logfile_csv, 'prediction_data_name':prediction_data_name, 'sequence_source': sequence_source, 'run_date':run_date, 'transcription_factor':TF, 'individuals':individuals, 'n_individuals':n_individuals}

        print(agg_dt)

        with(open(f'{metadata_dir}/aggregation_config_{prediction_data_name}_{TF}.json', mode='w')) as wj:
            json.dump(agg_dt, wj)     

    elif sequence_source == 'random':
        print(f'[INFO] Writing `aggregation_config_{prediction_data_name}_{TF}.json` file to {metadata_dir}')

        agg_dt = {'predictions_folder': project_dir, 'enformer_prediction_path': f'{output_dir}', 'log_file':logfile_csv, 'prediction_data_name':prediction_data_name, 'sequence_source': sequence_source, 'run_date':run_date, 'transcription_factor':TF, 'individuals':None, 'n_individuals':n_individuals}

        with(open(f'{metadata_dir}/aggregation_config_{prediction_data_name}_{TF}.json', mode='w')) as wj:
            json.dump(agg_dt, wj)         
                    
if __name__ == '__main__':
    check_input_parameters.check_inputs(args.param_config)
    #main()



# if create_hdf5_file == True:
#         print(f'[INFO] Creating HDF5 database(s)')
#         finished_predictions = pd.read_csv(logfile_csv)
#         make_db = make_h5_db_parsl(use_parsl = use_parsl)

#         db_parsl = []
#         for each_id in id_list:
#             motifs_list = finished_predictions.loc[finished_predictions['individual'] == each_id, ].motif.values.tolist()
#             motifs_list = list(set(motifs_list))

#             print(f'[INFO] Creating HDF5 database for {each_id} for {len(motifs_list)} predictions.')

#             motifs_list_paths = [f'{output_dir}/{each_id}/{i}_predictions.h5' for i in motifs_list]
#             csv_file = f'{output_dir}/{each_id}_{TF}_predictions.csv'
#             h5_file = f'{output_dir}/{each_id}_{TF}_predictions.hdf5'
#             db_parsl.append(make_db(h5_file = h5_file, csv_file = csv_file, files_list = motifs_list, files_path = motifs_list_paths, dataset = each_id))

        
#         print(db_parsl)
#         if use_parsl == True:
#             exec_parsl = [q.result() for q in tqdm.tqdm(db_parsl, desc=f'[INFO] Executing database futures.')] 
#             print(exec_parsl)