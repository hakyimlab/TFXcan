



import os, h5py, re
import numpy as np
import math
import warnings


def calculate_genomic_coordinates(locus, num_bins = 896):

    locus = locus.split('_')
    start = int(locus[1])
    end = int(locus[2])
    midn = math.ceil((start + end) / 2)

    # from position choose center bin
    center_ind = midn - 1
    center_bin = center_ind // 128
    
    half_bins = num_bins // 2
    start_bin = center_bin - half_bins
    end_bin = center_bin + half_bins
    if num_bins % 2 != 0: # if num_bins is odd
        end_bin += 1

    return([start_bin, end_bin])

def extract_enformer_predictions(enfref_dir, chr_num, locs):
    with h5py.File(os.path.join(enfref_dir, f"chr{chr_num}_cat.h5"), "r") as f:
        epigen = f[f'chr{chr_num}'][locs[0]:locs[1], :] 

    return epigen


def slice_bins(locus, bin_size=128, nbins=896):
    locus = locus.split('_')
    start = int(locus[1])
    end = int(locus[2])
    midn = math.ceil((start + end) / 2)
    nstart = midn - ((bin_size*nbins) / 2) # (128*896) / 2
    sstart = (start - nstart)
    send = sstart + ((end - start) - 1)
    bins = list(range(0, nbins*bin_size, bin_size))
    out = np.digitize([sstart, send], bins=bins).tolist()

    if((end - start) <= 128):
        pass
    elif(((end - start) > 128) or (end - start) % bin_size > 0):
        out[1] = out[1] + 1

    # if [448], make it [448, 449]
    # if [448, 448] make it [448, 449]
    
    if len(out) == 1:
        out.append(out[0] + 1)
    elif len(out) == 2:
        if out[0] == out[1]:
            out[1] = out[1] + 1
    #print(f"INFO - Bins to take: {out}")
    return(out)
    

def aggregate_and_collect_epigenome(locus, reference_epigenome_dir, pad_bins = 1):

    output = dict()

    # get the chrom
    chrom = locus.split('_')[0]
    chrom = int(re.sub('chr', '', chrom))

    # get bins coords
    bins_coords = calculate_genomic_coordinates(locus)
    enf_pred = extract_enformer_predictions(enfref_dir=reference_epigenome_dir, chr_num=chrom, locs=bins_coords)
    bins_to_take = slice_bins(locus)
    bins_to_aggregate = [bins_to_take[0] - pad_bins, bins_to_take[1] + pad_bins] # add one more bin upstream and downstrea
    print(f"INFO - After padding bins with {pad_bins}: {bins_to_aggregate}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        try:
            output['locus'] = locus
            output['values'] = enf_pred[bins_to_aggregate[0]:bins_to_aggregate[1], :].mean(axis = 0)
        except RuntimeWarning as rw: # this probably won't work since I am ignoring the warning anyway
            print(f'ERROR - {locus} has encountered an error')
            output['locus'] = None
            output['values'] = None
    return(output)