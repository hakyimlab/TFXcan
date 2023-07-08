


import pandas as pd
import cyvcf2
import numpy as np
import subprocess, re, os, sys
from functools import reduce
import utils
import parsl
import random
import parslConfiguration


parsl_parameters = {
    "job_name": "collect_l_tfpred_training_data",
    "num_of_full_nodes": 4,
    "walltime": "06:00:00",
    "init_blocks":1,
    "min_num_blocks": 0,
    "max_num_blocks": 10,
    "queue": "preemptable",
    "account": "covid-ct",
    "hpc": "polaris",
    "provider": "highthroughput",
    "worker_init":"source ~/.bashrc; conda activate /lus/grand/projects/TFXcan/imlab/users/temi/software/conda_envs/compbio-tools; which python; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lus/grand/projects/TFXcan/imlab/users/temi/software/conda_envs/compbio-tools/lib"
  }

parsl_parameters['working_dir'] = '/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/compare_predictors'

parsl.load(parslConfiguration.polaris_htParslConfig(params=parsl_parameters))
print(f"INFO - parsl loaded correctly")
print(f"INFO - Using {parsl_parameters['provider']} parsl configuration for polaris")
print(f"INFO - Using parsl with version {parsl.__version__}")

utils = os.path.join('.', 'utils.py')
exec(open(utils).read(), globals(), globals())

imlab_dir = '/lus/grand/projects/TFXcan/imlab'
data_dir = f'{imlab_dir}/data/freedman_data/peaks_liftover'
project_dir = f'{imlab_dir}/users/temi/projects/TFXcan'
hg38_snps_vcf = '/lus/grand/projects/TFXcan/imlab/data/variants_data/hg38_snps/GCF_000001405.40.gz'
save_dir = '/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/compare_predictors/prediction_matrices'
if not os.path.isdir(save_dir): os.makedirs(save_dir, exist_ok = True)

# read in the tfpred scores and collect
tfpred_dir = '/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/compare_predictors/tfpred_scores'
tfpred_scores = {i: pd.read_csv(f'{tfpred_dir}/aggByCollect_AR_Prostate_1KG.linear.0{i}.txt') for i in [1,2,3,4]}
# I only want loci that are present in all of them
common_loci = [a['locus'].tolist() for a in tfpred_scores.values()]
common_loci = list(set([a for b in common_loci for a in b]))
tfpred_scores = {k: v.loc[v['locus'].isin(common_loci)] for k, v in tfpred_scores.items()}
tfpred_scores_dt = reduce(lambda x, y: pd.merge(x, y, on = 'locus'), tfpred_scores.values())

# cwas loci ===
cwas_loci = pd.read_table(f'/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/compare_predictors/metadata/cwas_intervals.txt', header=None).iloc[:,0].tolist()
cwas_loci = [r for r in cwas_loci if isinstance(r, str)]

#random.shuffle(cwas_loci)
# individuals ===
individuals = pd.read_table(f'/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/compare_predictors/metadata/1000_genome_individuals.txt', header=None).iloc[:,0].tolist()

print(f'INFO - Collecting matrix for {len(individuals)} individuals and {len(cwas_loci)} loci.')

# create a dictionary of loci separated by chromosome
chromosomes = [f"chr{i}" for i in range(1, 23)]
chromosomes.extend(['chrX', 'chrY'])

available_loci_by_chr = {}
for chromosome in chromosomes:
    chr_list_of_regions = [r for r in cwas_loci if r.startswith(f"{chromosome}_")]
    if chr_list_of_regions:
        available_loci_by_chr[chromosome] = chr_list_of_regions

# step 2 - collect locus information (vcf file and location) by chromosome
ltfpred_obj = {k: create_locus_info(k, v) for k, v in available_loci_by_chr.items()}
print(f"INFO - created locus information per chromosome")

# step 3 - get prediction matrices for each locus by chromosome
# this is what wastes time the most and this is where I use parsl
alleles_matrices = {k: create_prediction_matrix(chr_=k, chr_info = ltfpred_obj[k], individuals=individuals, tfpred_matrix = tfpred_scores_dt, save_dir = save_dir).result() for k in list(ltfpred_obj.keys())}

# training_loci = []
# for chrK, chrV in prediction_matrices.items():
#     training_loci.extend(list(chrV.keys()))

# pd.DataFrame(training_loci).to_csv(os.path.join(chrF, '..', f'training_locus_1KG.csv'), sep = '\t', index=False, header=False)
print(f"INFO - Matrices of alleles, dosages and predictions have been saved.")

# for chrK, chrV in alleles_matrices.items():
#     chrF = os.path.join(save_dir, chrK)
#     if not os.path.isdir(chrF): os.makedirs(chrF, exist_ok = True)
#     for k, v in chrV.items():
#         v.to_csv(os.path.join(chrF, f'{k}_1KG_alleles.csv'), sep = '\t', index=False)


# dosages_matrices = {}
# for chrK, chrV in alleles_matrices.items():
#     dosages_matrices[chrK] = {k: v.apply(alleles_to_dosages, axis=1).set_index(v['locus']).T.reset_index() for k, v in chrV.items() if v is not None}

# for chrK, chrV in dosages_matrices.items():
#     chrF = os.path.join(save_dir, chrK)
#     if not os.path.isdir(chrF): os.makedirs(chrF, exist_ok = True)
#     for k, v in chrV.items():
#         v.to_csv(os.path.join(chrF, f'{k}_1KG_dosages.csv'), sep = '\t', index=False)


# prediction_matrices = {}
# for chrK, chrV in dosages_matrices.items():
#     prediction_matrices[chrK] = {k: create_per_locus_training_matrix(v, k, tfpred_scores_dt) for k, v in chrV.items() if v is not None}


# for chrK, chrV in prediction_matrices.items():
#     chrF = os.path.join(save_dir, chrK)
#     if not os.path.isdir(chrF): os.makedirs(chrF, exist_ok = True)
#     for k, v in chrV.items():
#         v.to_csv(os.path.join(chrF, f'{k}_1KG_prediction_matrix.csv'), sep = '\t', index=False)


# save the list of loci to be trained

