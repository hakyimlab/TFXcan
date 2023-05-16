
import os
os.chdir('/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/scripts/collect_model_data/')

import pandas as pd
import os, sys, tqdm, json
from datetime import date
import parsl
import subprocess
from parsl.app.app import python_app
import argparse

# needed arguments 
parser = argparse.ArgumentParser()
parser.add_argument("--metadata_file", help="Path to the metadata file", type=str)
parser.add_argument("--agg_type", help="aggregation type to be used", type=str)
args = parser.parse_args()

#todays_date = date.today().strftime("%Y-%m-%d")

# some locations and folders
whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)
script_path = os.path.abspath(whereis_script)
utils_path = os.path.join(script_path, 'utilities')
sys.path.append(utils_path)

try:
    import utility_functions
    import parslConfiguration
except ModuleNotFoundError:
    print(f'[ERROR] utility_functions module not found.')

#prediction_path=/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_pipeline/enformer_predictions/kawakami/predictions_2022-12-20/kawakami,log_dir=/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_pipeline/predictions_log/kawakami/predictions_log_2022-12-20,id=kawakami,agg_type=aggByMean

with open(f'{args.metadata_file}') as f:
    parameters = json.load(f)

    enformer_predictions_path = parameters['enformer_predictions_path']
    log_file = parameters['log_file']
    each_id = parameters['each_id']# e.g. "kawakami" or "cistrome"
    todays_date = parameters['run_date']
    base_path = parameters['predictions_folder']
    TF = parameters['transcription_factor']

agg_type = args.agg_type

print(f'[INFO] Currently on {each_id}')

#base_path = '/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline'

save_dir = f'{base_path}/aggregated_predictions/{each_id}_{TF}_{todays_date}'
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)

#logpath = f'{log_path}/{each_id}_{TF}_predictions_log.csv'
log_data = pd.read_csv(log_file)
log_data = log_data.drop_duplicates(subset=['motif']) #.iloc[1:1000, ]

# use parsl if the num of rows of log_data is more than 10000
use_parsl = True
print(f'[INFO] Using parsl. Aggregation will be split into multiple files, and batches will be collected into one later.')

# fpath = os.path.join(base_path, 'scripts', 'collect_model_data')
# predict_utils_one = f'{fpath}/utility_functions.py'
# exec(open(predict_utils_one).read(), globals(), globals())

#bpath = os.path.join(base_path, 'modeling_pipeline')
parsl_params = {'working_dir':base_path, 'job_name':'collect_kawakami_data', 'queue':"preemptable", 'walltime':"05:59:00", 'num_of_full_nodes':1, 'min_num_blocks':0, 'max_num_blocks':10}
parsl.load(parslConfiguration.localParslConfig_htpool(parsl_params))

collection_fxn = utility_functions.return_prediction_function(use_parsl)
range_batches = utility_functions.generate_batch(range(0, log_data.shape[0]), batch_n=10)

app_futures = []
for i, range_batch in enumerate(range_batches):

    print(i)
    app_futures.append(collection_fxn(each_id=each_id, log_data=log_data.iloc[range_batch, ], predictions_path=enformer_predictions_path, TF=TF, agg_types=[agg_type], base_path=base_path, save_dir=save_dir, batch_num=i))

print(app_futures)

app_execs = [r.result() for r in app_futures]

cmd = f"for i in `ls {save_dir}/*_{agg_type}_{TF}_batch_*.csv.gz`; do zcat $i | sed '1d'; done | pigz -c >{save_dir}/{each_id}_{agg_type}_{TF}.csv.gz"    

p = subprocess.Popen(cmd, shell=True)

# remove the files


print(f'[INFO] Status: Aggregation complete for {log_data.shape[0]} predictions for {each_id}.')

# for i in `ls ./*.gz`; do zcat \"$i\" | sed '1d'; done > k_aggByMean_T_compiled.csv & gzip -f k_aggByMean_T_compiled.csv