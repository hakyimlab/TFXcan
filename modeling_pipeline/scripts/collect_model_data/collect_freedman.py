import pandas as pd
import os, sys
from datetime import date

todays_date = date.today().strftime("%Y-%m-%d")
print(todays_date)

def collect_modeling_data_for_freedman(each_individual, log_data, individual_predictions_path, TF, data_names, base_path, save_dir):

    import h5py
    import numpy as np
    import os
    import pandas as pd
    # read in one of the files

    exec(open(f'{base_path}/modeling_pipeline/scripts/collect_model_data/utility-functions.py').read(), globals(), globals())

    freedman_predictions = {}

    for dt in log_data.loc[log_data['sequence_type'] == 'var', ].motif.values.tolist():
        fle = f'{individual_predictions_path}/{dt}_predictions.h5'
        if os.path.isfile(fle):
            with h5py.File(fle, 'r') as f:
                filekey = list(f.keys())[0]
                # should I select tracks? ; maybe not yet
                freedman_predictions[dt] = np.vstack(list(f[filekey]))
        else:
            print('File does not exist')

    print(f'[INFO] Finished collecting {len(freedman_predictions)} predictions for {each_individual}')

    dt_aggbycenter = agg_by_center(freedman_predictions, center=8)
    data_list = [dt_aggbycenter]

    # test_aggbymean, test_aggbycenter, test_aggbymean_upstream, test_aggbymean_downstream, test_aggbymean_upstream_downstream = agg_byall(freedman_predictions)
    # data_list = [test_aggbymean, test_aggbycenter, test_aggbymean_upstream, test_aggbymean_downstream, test_aggbymean_upstream_downstream]

    for i, dt in enumerate(data_list):

        ty = pd.concat([pd.Series(list(freedman_predictions.keys())), pd.DataFrame(dt)], axis=1)

        print(f'[INFO] Dimension of collected data is {ty.shape[0]} by {ty.shape[1]}')

        column_names = ['id']
        column_names.extend([f'f_{i}' for i in range(1, ty.shape[1])])

        ty = ty.set_axis(column_names, axis=1, inplace=False)

        print(ty.iloc[0:5, 0:5])

        ty.to_csv(path_or_buf=f'{save_dir}/{each_individual}_{data_names[i]}_{TF}.csv.gz', index=False, compression='gzip')
    print(f'[INFO] Finished saving data for {each_individual}')

    return(0)

each_individual = sys.argv[1]
print(f'[INFO] Currently on {each_individual}')

upstream = list(range(0, 8))
center = [8]
downstream = list(range(9, 17))

base_path = '/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan'
enformer_predictions_path = f'{base_path}/enformer_pipeline/enformer_predictions/freedman/predictions_2022-12-11'
individuals_file = f'{base_path}/enformer_pipeline/metadata/individuals.txt'
individuals_log_path = f'{base_path}/enformer_pipeline/predictions_log/predictions_log_2022-12-11'

save_dir = f'{base_path}/modeling_pipeline/data/train-test-val/freedman/data_{todays_date}'
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)

TF = 'FOXA1'
data_names = ['aggByCenter']
#data_names = ['aggByMean', 'aggByCenter', 'aggByUpstream', 'aggByDownstream', 'aggByUpstreamDownstream']

logpath = f'{individuals_log_path}/{each_individual}_predictions_log.csv'
individual_predictions_path = f'{enformer_predictions_path}/{each_individual}'
log_data = pd.read_csv(logpath)
log_data = log_data.drop_duplicates(subset=['motif'])

collected = collect_modeling_data_for_freedman(each_individual=each_individual, log_data=log_data, individual_predictions_path=individual_predictions_path, TF=TF, data_names=data_names, base_path=base_path, save_dir=save_dir)

#print(collected)

print(f'[INFO] Status: {collected} for {log_data.shape[0]} predictions for {each_individual}.')