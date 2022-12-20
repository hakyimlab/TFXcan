import pandas as pd
import os, sys, tqdm
from datetime import date
import parsl
import subprocess
from parsl.app.app import python_app

todays_date = date.today().strftime("%Y-%m-%d")


import utility_functions

#prediction_path=/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_pipeline/enformer_predictions/kawakami/predictions_2022-12-20/kawakami,log_dir=/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_pipeline/predictions_log/kawakami/predictions_log_2022-12-20,id=kawakami,agg_type=aggByMean

# enformer_predictions_path = "/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_pipeline/enformer_predictions/kawakami/predictions_2022-12-20/kawakami" #sys.argv[1] 
# log_path = "/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_pipeline/predictions_log/kawakami/predictions_log_2022-12-20" #sys.argv[2]
# each_id = "kawakami" #sys.argv[3] # e.g. "kawakami" or "cistrome"
# agg_type = "aggByMean" #sys.argv[4]

enformer_predictions_path = sys.argv[1] 
log_path = sys.argv[2]
each_id = sys.argv[3] # e.g. "kawakami" or "cistrome"
agg_type = sys.argv[4]

print(f'[INFO] Currently on {each_id}')

upstream = list(range(0, 8))
center = [8]
downstream = list(range(9, 17))

base_path = '/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan'
# enformer_predictions_path = f'{base_path}/enformer_pipeline/enformer_predictions/{each_id}_reference/predictions_2022-12-11/{each_id}_FOXA1'
# log_path = f'{base_path}/enformer_pipeline/predictions_log/predictions_log_2022-12-11'

save_dir = f'{base_path}/modeling_pipeline/data/train-test-val/{each_id}/data_{todays_date}'
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)

TF = 'FOXA1'
#data_names = ['aggByCenter']
#data_names = ['aggByMean', 'aggByCenter', 'aggByUpstream', 'aggByDownstream', 'aggByUpstreamDownstream']

logpath = f'{log_path}/{each_id}_{TF}_predictions_log.csv'
log_data = pd.read_csv(logpath)
log_data = log_data.drop_duplicates(subset=['motif'])#.iloc[1:1000, ]

# use parsl if the num of rows of log_data is more than 10000
use_parsl = True
print(f'[INFO] Using parsl. Aggregation will be split into multiple files. Collect the batches into one later.')

bpath = os.path.join(base_path, 'modeling_pipeline')
parsl.load(utility_functions.localParslConfig_threadpool({'working_dir': bpath}))

collection_fxn = utility_functions.return_prediction_function(use_parsl)
range_batches = utility_functions.generate_batch(range(0, log_data.shape[0]), batch_n=5)

app_futures = []
for i, range_batch in enumerate(range_batches):

    app_futures.append(collection_fxn(each_id=each_id, log_data=log_data.iloc[range_batch, ], predictions_path=enformer_predictions_path, TF=TF, agg_types=[agg_type], base_path=base_path, save_dir=save_dir, batch_num=i))

app_execs = [r.result() for r in app_futures]

cmd = f"for i in `ls {save_dir}/*_batch_*.csv.gz`; do zcat $i | sed '1d'; done | pigz -c >{save_dir}/{each_id}_{agg_type}_{TF}.csv.gz"    

p = subprocess.Popen(cmd, shell=True)

# remove the files


print(f'[INFO] Status: Aggregation complete for {log_data.shape[0]} predictions for {each_id}.')

# for i in `ls ./*.gz`; do zcat \"$i\" | sed '1d'; done > k_aggByMean_T_compiled.csv & gzip -f k_aggByMean_T_compiled.csv