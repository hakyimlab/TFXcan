import h5py
import pandas as pd
import numpy as np
import os, sys, tqdm


data_dir = '/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_pipeline/enformer_predictions/kawakami-reference'

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

predictions_list = os.listdir(data_dir)
training_predictions = {}
for dt in tqdm.tqdm(predictions_list):
    fle = f'{data_dir}/{dt}'
    if os.path.isfile(fle):
        with h5py.File(fle, 'r') as f:
            filekey = list(f.keys())[0]
            #print(filekey)
            # should I select tracks? ; maybe not yet
            #kk = 0 if dt.startswith('neg') else 1
            training_predictions[dt] = np.vstack(list(f[filekey]))
    else:
        print('File does not exist')


full_dt = pd.DataFrame(agg_by_center(training_predictions))

