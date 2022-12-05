import pandas as pd
import os, sys, h5py
import numpy as np
from parsl import python_app
import parsl
# these are the bins

print(sys.path)

import parslConfiguration
parsl.load(parslConfiguration.localParslConfig())

@python_app
def collect_modeling_data_for_freedman(each_individual, log_data, individual_predictions_path, TF, data_names, base_path, save_dir):

    import h5py
    import numpy as np
    import os
    import pandas as pd
    # read in one of the files

    exec(open(f'{base_path}/modeling_pipeline/scripts/utility-functions.py').read(), globals(), globals())

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

    test_aggbymean, test_aggbycenter, test_aggbymean_upstream, test_aggbymean_downstream, test_aggbymean_upstream_downstream = agg_byall(freedman_predictions)
    data_list = [test_aggbymean, test_aggbycenter, test_aggbymean_upstream, test_aggbymean_downstream, test_aggbymean_upstream_downstream]

    for i, dt in enumerate(data_list):

        ty = pd.concat([pd.Series(list(freedman_predictions.keys())), pd.DataFrame(dt)], axis=1)

        column_names = ['id', 'class']
        column_names.extend([f'f_{i}' for i in range(1, ty.shape[1] - 1)])

        ty = ty.set_axis(column_names, axis=1, copy=False)

        ty.to_csv(path_or_buf=f'{save_dir}/{each_individual}_{data_names[i]}_{TF}.csv.gz', index=False, compression='gzip')
    print(f'[INFO] Finished saving data for {each_individual}')

    return(0)


#each_individual = sys.argv[1]
upstream = list(range(0, 8))
center = [8]
downstream = list(range(9, 17))

base_path = '/projects/covid-ct/imlab/users/temi/projects/TFXcan'
enformer_predictions_path = f'{base_path}/enformer_pipeline/enformer_predictions/freedman'
individuals_file = f'{base_path}/enformer_pipeline/metadata/individuals.txt'
individuals_log_path = f'{base_path}/enformer_pipeline/predictions-log'

save_dir = f'{base_path}/modeling_pipeline/data/train-test-val/freedman'
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)

individuals = pd.read_table(individuals_file, header=None)[0].tolist()
#exec(open(f'{base_path}/modeling_pipeline/scripts/utility-functions.py').read(), globals(), globals())

TF = 'FOXA1'
data_names = ['aggByMean', 'aggByCenter', 'aggByUpstream', 'aggByDownstream', 'aggByUpstreamDownstream']

app_futures = []
for each_individual in individuals:
    logpath = f'{individuals_log_path}/{each_individual}_predictions_log.csv'
    individual_predictions_path = f'{enformer_predictions_path}/{each_individual}'
    log_data = pd.read_csv(logpath)

    app_futures.append(collect_modeling_data_for_freedman(each_individual=each_individual, log_data=log_data, individual_predictions_path=individual_predictions_path, TF=TF, data_names=data_names, base_path=base_path, save_dir=save_dir))

print(app_futures)

exec_futures = [a.result() for a in app_futures]

print('[INFO] All done.')