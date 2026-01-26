import os, sys
import pandas as pd, numpy as np
import argparse
import h5py

def list_of_ints(arg):
    return list(map(int, arg.split(',')))

# needed arguments 
parser = argparse.ArgumentParser()
parser.add_argument("--metadata", help="Path to file", type=str, default=None)
# parser.add_argument("--process_by_haplotype", help="Path to file", action=argparse.BooleanOptionalAction)
# parser.add_argument("--process_function", help="Path", type=str, default='sum')
parser.add_argument("--output_basename", help="", type=str)
args = parser.parse_args()

# read in the metadata
dt_metadata = pd.read_table(args.metadata)

print(dt_metadata.shape)

dt_temp = pd.read_csv(dt_metadata.path[0])
dt_temp.insert(loc=0, column='individual', value=dt_metadata.individual[0])
dt_temp.rename(columns={'id': 'locus'}, inplace=True)

tpp = list(dt_metadata.itertuples(index=False))

def process_file(tp):
    dt_temp = pd.read_csv(tp.path)
    dt_temp.insert(loc=0, column='individual', value=tp.individual)
    dt_temp.rename(columns={'id': 'locus'}, inplace=True)
    return dt_temp

processed_lsit = [process_file(tp) for tp in tpp]

eur_1kg_predictions = pd.concat(processed_lsit, axis = 0, ignore_index=True)

xmatrix = pd.DataFrame(eur_1kg_predictions.drop(columns = ['locus', 'individual']).to_numpy())
meta = eur_1kg_predictions[['locus', 'individual']]

# ensure that the rows are the same
assert xmatrix.shape[0] == meta.shape[0]

# write out the metadata
#bname = 'prostate_cancer_risk.2024-09-30.processed'
metadata_file = os.path.join(f'{args.output_basename}.metadata.tsv')
matrix_file = os.path.join(f'{args.output_basename}.matrix.h5.gz')
meta.to_csv(metadata_file, sep = '\t', index = False)
xmatrix.to_hdf(matrix_file, key = 'matrix', mode='w', complevel=9)








# # check if available
# if not os.path.exists(args.merged_h5_file):
#     print(f'ERROR - File {args.merged_h5_file} does not exist')
#     sys.exit(1)

# pfile = args.merged_h5_file
# with h5py.File(pfile, 'r') as hf:
#     mt = hf.get('metadata')[()] # metadata i.e. chr_pos:individual:haplotype:batch
#     dt = hf.get('predictions')[()]

# print(f"INFO - Read in the file {pfile} with shape {dt.shape}")

# # process the data
# if args.process_by_haplotype:
#     if args.process_function == 'sum':
#         dt = dt.reshape(-1, 2, dt.shape[-1]).sum(1)
#     elif args.process_function == 'mean':
#         dt = dt.reshape(-1, 2, dt.shape[-1]).mean(1)
#     else:
#         print(f'ERROR - Process function {args.process_function} not implemented')
#         sys.exit(1)

# elif not args.process_by_haplotype:
#     pass

# print(f"INFO - Processed the data to shape {dt.shape}")

# # process the metadata
# decoded_m = [m.decode("utf-8").split(":") for m in mt.tolist()]
# xd = pd.DataFrame(decoded_m)
# xd.columns = ['locus', 'individual', 'haplotype', 'batch']
# # remove duplicates
# df_metadata = xd.drop_duplicates(subset = ['locus', 'individual'])[['locus', 'individual']].reset_index(drop = True)

# # ensure that the rows are the same
# assert df_metadata.shape[0] == dt.shape[0]

# # write out the metadata
# metadata_file = f'{args.output_basename}.metadata.tsv'
# df_metadata.to_csv(metadata_file, sep = '\t', index = False)

# # write out the data
# xt = pd.DataFrame(dt)
# # data_file = f'{args.output_basename}.matrix.tsv.gz'
# data_file = f'{args.output_basename}.matrix.h5.gz'
# xt.to_hdf(data_file, key = 'matrix', mode='w', complevel=9)
# #xt.to_csv(data_file, sep = '\t', index = False, compression = 'gzip')
