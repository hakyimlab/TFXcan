

# Description: Common snakefile
# Author: Temi
# Date: Wed Mar 29 2023


import pandas as pd
import os, glob, sys, re
from snakemake.io import glob_wildcards

import module

cfile = "config/pipeline.yaml"
configfile: cfile



rule all:
    input:
        os.path.join('data', 'genotypes'),
        expand(os.path.join(HOMERFILES_DIR, '{tf}', '{motif_file}.motif'), zip, tf = homer_motifs_wildcards.tfs, motif_file = homer_motifs_wildcards.motif_files),
        expand(os.path.join(HOMERFILES_DIR, '{tf}', 'scanMotifsGenomeWide_{motif_file}.txt'), zip, tf = homer_motifs_wildcards.tfs, motif_file = homer_motifs_wildcards.motif_files),
        expand(os.path.join(HOMERFILES_DIR, '{tf}', 'merged_motif_file.txt'), tf = set(homer_motifs_wildcards.tfs)),
        expand(os.path.join(PREDICTORS_DIR, '{tf}_{tissue}_predictors.txt'), zip, tf=TF_list, tissue=tissue_list),
        expand(os.path.join(PREDICTORS_DIR, '{tf}_{tissue}_ground_truth.txt'), zip, tf=TF_list, tissue=tissue_list),
        expand(os.path.join(PREDICTION_PARAMS_DIR, f'enformer_parameters_{config["dataset"]}_{{tf}}_{{tissue}}.json'), zip, tf=TF_list, tissue=tissue_list)