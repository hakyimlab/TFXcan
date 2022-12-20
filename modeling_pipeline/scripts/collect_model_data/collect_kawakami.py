import pandas as pd
import os, sys, tqdm
from datetime import date

todays_date = date.today().strftime("%Y-%m-%d")

#agg_type = sys.argv[1] #'aggByMean' # data_names = ['aggByMean', 'aggByCenter', 'aggByUpstream', 'aggByDownstream', 'aggByUpstreamDownstream']

import utility_functions


#prediction_path=/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_pipeline/enformer_predictions/kawakami/predictions_2022-12-20/kawakami,log_dir=/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_pipeline/predictions_log/kawakami/predictions_log_2022-12-20,id=kawakami,agg_type=aggByMean

enformer_predictions_path = "/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_pipeline/enformer_predictions/kawakami/predictions_2022-12-20/kawakami" #sys.argv[1] 
log_path = "/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_pipeline/predictions_log/kawakami/predictions_log_2022-12-20" #sys.argv[2]
each_id = "kawakami" #sys.argv[3] # e.g. "kawakami" or "cistrome"
agg_type = "aggByMean" #sys.argv[4]

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
log_data = log_data.drop_duplicates(subset=['motif'])

#exec(open(f'{base_path}/modeling_pipeline/scripts/collect_model_data/utility-functions.py').read(), globals(), globals())
collected = utility_functions.collect_modeling_data_for_kawakami(each_id=each_id, log_data=log_data, predictions_path=enformer_predictions_path, TF=TF, agg_types=[agg_type], base_path=base_path, save_dir=save_dir)
print(f'[INFO] Status: {collected} for {log_data.shape[0]} predictions for {each_id}.')