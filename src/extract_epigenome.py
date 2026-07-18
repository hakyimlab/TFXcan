# Author: Temi
# Date: Wednesday June 18 2025
# Description: helper functions for the bin/coordinate math that slices Enformer's 896-bin output around a genomic locus
# Usage: imported by aggregate_epigenomes.py, not run directly

import os, h5py, re
import numpy as np
import math
import warnings


def calculate_genomic_coordinates(locus, num_bins = 896):
    # locus is "chr_start_end"; find the whole-chromosome bin indices (128bp/bin,
    # matching Enformer's resolution) for a num_bins window centered on the locus midpoint
    locus = locus.split('_')
    start = int(locus[1])
    end = int(locus[2])
    midn = math.ceil((start + end) / 2)

    # from position choose center bin
    center_ind = midn - 1 # 0-based position (locus coords are 1-based)
    center_bin = center_ind // 128 # which 128bp bin this position falls in

    half_bins = num_bins // 2
    start_bin = center_bin - half_bins
    end_bin = center_bin + half_bins
    if num_bins % 2 != 0: # if num_bins is odd
        end_bin += 1 # put the extra bin on the end so the window still spans num_bins bins

    return([start_bin, end_bin])

def extract_enformer_predictions(enfref_dir, chr_num, locs):
    # chr{n}_cat.h5 holds that whole chromosome's Enformer predictions, one row per 128bp bin;
    # locs is the [start_bin, end_bin] window from calculate_genomic_coordinates
    with h5py.File(os.path.join(enfref_dir, f"chr{chr_num}_cat.h5"), "r") as f:
        epigen = f[f'chr{chr_num}'][locs[0]:locs[1], :]

    return epigen


def slice_bins(locus, bin_size=128, nbins=896):
    # find where the locus itself sits *within* the nbins-bin context window extracted above,
    # so aggregate_and_collect_epigenome can average over just the locus's own bins (plus padding)
    locus = locus.split('_')
    start = int(locus[1])
    end = int(locus[2])
    midn = math.ceil((start + end) / 2)
    nstart = midn - ((bin_size*nbins) / 2) # (128*896) / 2 -- bp coordinate of the window's left edge (mirrors calculate_genomic_coordinates, but in bp instead of bin index)
    sstart = (start - nstart) # locus start, as a bp offset from the window's left edge
    send = sstart + ((end - start) - 1) # locus end, as a bp offset from the window's left edge
    bins = list(range(0, nbins*bin_size, bin_size)) # bp boundary of each of the nbins bins within the window
    out = np.digitize([sstart, send], bins=bins).tolist() # convert those bp offsets back into local bin indices (0..nbins)

    if((end - start) <= 128):
        pass
    elif(((end - start) > 128) or (end - start) % bin_size > 0):
        out[1] = out[1] + 1 # locus doesn't fit cleanly in one bin, so extend the end bin by one

    # if [448], make it [448, 449]
    # if [448, 448] make it [448, 449]

    if len(out) == 1:
        out.append(out[0] + 1)
    elif len(out) == 2:
        if out[0] == out[1]:
            out[1] = out[1] + 1 # avoid a zero-width bin range
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