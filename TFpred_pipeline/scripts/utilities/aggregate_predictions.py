
# import os
# os.chdir('/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/scripts/collect_model_data/')

import os
import pandas as pd
import os, sys, json
from datetime import date
import parsl
from parsl.app.app import python_app
import subprocess
import argparse

import multiprocessing



# needed arguments 
parser = argparse.ArgumentParser()
parser.add_argument("--metadata_file", help="Path to the metadata file", type=str)
parser.add_argument("--agg_types", nargs='+', help='<Required> aggregation type to be used', required=True)
args = parser.parse_args()

print(f'[INFO] Available CPUs are {multiprocessing.cpu_count()}')

# use parsl if the num of rows of log_data is more than 10000
use_parsl = True

#todays_date = date.today().strftime("%Y-%m-%d")

# some locations and folders
whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)
script_path = os.path.abspath(whereis_script)
utils_path = os.path.join(script_path, 'utilities')
sys.path.append(utils_path)

#print(sys.path)

try:
    import utility_functions
    import parslConfiguration
except ModuleNotFoundError:
    print(f'[ERROR] utility_functions module not found.')

with open(f'{args.metadata_file}') as f:
    parameters = json.load(f)

    enformer_predictions_path = parameters['enformer_prediction_path']
    log_file = parameters['log_file']
    sequence_source = parameters['sequence_source']# e.g. "kawakami" or "cistrome"
    todays_date = parameters['run_date']
    base_path = parameters['predictions_folder']
    TF = parameters['transcription_factor']
    individuals = parameters['individuals']

# determine what individuals to predict on and all that
if sequence_source in ['freedman']:
    if isinstance(individuals, list):
        pass
    elif isinstance(individuals, type('str')):
        if os.path.isfile(individuals):
            individuals = pd.read_table(individuals, header=None)[0].tolist()

            print(type(individuals))

agg_types = args.agg_types
agg_types = agg_types[0].split(' ')
print(f'[INFO] Aggregating these: {agg_types}')

print(f'[INFO] Currently on {sequence_source}')

save_dir = f'{base_path}/aggregated_predictions/{sequence_source}_{TF}_{todays_date}'
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)

if individuals is None:
    ids_names = [sequence_source]
elif isinstance(individuals, list):
    ids_names = individuals

print(ids_names)

log_data_all = pd.read_csv(log_file)
#log_data_all = log_data_all.drop_duplicates(subset=['motif']) #.iloc[1:1000, ]

#fpath = os.path.join(base_path, 'scripts', 'collect_model_data')
predict_utils_one = f'{script_path}/utility_functions.py'
exec(open(predict_utils_one).read(), globals(), globals())

if use_parsl == True:
    #bpath = os.path.join(base_path, 'modeling_pipeline')
    print(f'[INFO] Using parsl.')
    parsl_params = {'working_dir':base_path, 'job_name':'aggregate_predictions', 'queue':"preemptable", 'walltime':"06:00:00", 'num_of_full_nodes':4, 'min_num_blocks':0, 'max_num_blocks':4}
    #parsl.load(parslConfiguration.localParslConfig_htpool(parsl_params))
    #parsl.load(parslConfiguration.localParslConfig_htpool(parsl_params))
    parsl.load(parslConfiguration.polaris_htParslConfig(parsl_params))

collection_fxn = return_prediction_function(use_parsl)


app_futures = []
for each_id in ids_names:
    log_data = log_data_all.loc[log_data_all['individual'] == each_id, ]
    log_data = log_data.drop_duplicates(subset=['motif'])
    #print(log_data.iloc[0:5, ])
    print(f'[INFO] {each_id}\'s prediction log has {log_data.shape[0]} rows and {log_data.shape[1]} columns.')
    #log_data = log_data.iloc[0:3000, ]

    ind_path = os.path.join(enformer_predictions_path, each_id)

    app_futures.append(collection_fxn(each_id=each_id, log_data=log_data, predictions_path=ind_path, TF=TF, agg_types=agg_types, base_path=base_path, save_dir=save_dir, batch_num=None))

print(app_futures)
if use_parsl == True:
    app_execs = [r.result() for r in app_futures]

print(f'[INFO] Aggregation complete for all.')

# for each_id in ids_names:
#     for agg_type in agg_types:
#         cmd = f"for i in `ls {save_dir}/{each_id}_{agg_type}_{TF}_batch_*.csv.gz`; do zcat $i | sed '1d'; done | pigz -c >{save_dir}/{each_id}_{agg_type}_{TF}.csv.gz"    
#         p = subprocess.Popen(cmd, shell=True)
#         print(f'[INFO] Status: Aggregation complete for {log_data.shape[0]} predictions for {each_id}.')
# print(f'[INFO] Entire aggregation complete.')

    # cmd = f"for i in `ls {save_dir}/{each_id}_{agg_type}_{TF}_batch_*.csv.gz`; do zcat $i | sed '1d'; done | pigz -c >{save_dir}/{each_id}_{agg_type}_{TF}.csv.gz"    
    # p = subprocess.Popen(cmd, shell=True)
    # print(f'[INFO] Status: Aggregation complete for {log_data.shape[0]} predictions for {each_id}.')

# for i in `ls ./*.gz`; do zcat \"$i\" | sed '1d'; done > k_aggByMean_T_compiled.csv & gzip -f k_aggByMean_T_compiled.csv

# python3 ../scripts/utilities/aggregate_predictions.py --metadata_file /lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFpred_pipeline/metadata/aggregation_config_freedman_FOXA1.json --agg_types "aggByCenter aggByPreCenter aggByPostCenter aggByUpstreamDownstream aggByDownstream aggByUpstream aggByMean"