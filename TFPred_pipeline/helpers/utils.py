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
    import numpy as np

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
    import numpy as np

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

def read_file(motif, dir):
    import os
    import h5py
    import numpy as np

    output = {}
    fle = f'{dir}/{motif}_predictions.h5'
    if os.path.isfile(fle):
        with h5py.File(fle, 'r') as f:
            filekey = list(f.keys())[0]
            output[motif] = np.vstack(list(f[filekey]))
    else:
        print(f'[ERROR] {motif} predictions file does not exist.')

    return(output)