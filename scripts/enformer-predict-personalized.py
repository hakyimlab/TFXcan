
# === This script contains all codes to test if tensorflow and ENFORMER are working accurately  
# Created by Temi
# DATE: Thursday Oct 20 2022

from __future__ import absolute_import, division, print_function, unicode_literals
import os, re, sys, json, h5py, csv, warnings
from cyvcf2 import VCF
import parsl
from parsl.app.app import python_app


@python_app
def check_and_predict_query(line, sample, query_name, output_dir, logfile, query, model, vcf_file, temporary_vcf_dir, fasta_extractor, fasta_file_path, software_paths):

    import subprocess
    import kipoiseq # for manipulating fasta files # for some reason, kipoiseq has to be loaded first
    import time
    import tensorflow as tf
    import tensorflow_hub as hub # for interacting with saved models and tensorflow hub
    import joblib
    import gzip # for manipulating compressed files
    from kipoiseq import Interval # same as above, really
    import pyfaidx # to index our reference genome file
    import pandas as pd # for manipulating dataframes
    import numpy as np # for numerical computations
    import seaborn as sns # for plotting
    import pickle # for saving large objects
    import os, sys, re # functions for interacting with the operating system

    def check_query(sample, query_name, output_dir, logfile):
        '''
        Checks of predictions are available for an individual or sample
        '''

        if isinstance(logfile, pd.DataFrame):# should have read it>
            motif_in_logfile = logfile.motif.values
            individual_in_logfile = logfile.individual.values

            if (query_name in motif_in_logfile) and (sample in individual_in_logfile): # both the individual and motifs are present

                # check if the motif has been saved and if the status is `completed`
                query_saved = str(f'{output_dir}/{sample}/{query_name}_predictions.h5')
                
                if (logfile.loc[(logfile.motif==query_name) & (logfile.individual==sample), 'status'].values[0] == 'completed') and os.path.isfile(query_saved):
                    return(0)
                elif (logfile.loc[(logfile.motif==query_name) & (logfile.individual==sample), 'status'].values[0] != 'completed') or (not os.path.isfile(query_saved)):
                    return(1)

        elif isinstance(logfile, type(None)):
            return(2)

    def predict_query(query, output_dir, model, fasta_extractor, vcf_file, temporary_vcf_dir, fasta_file_path, sample, software_paths=[], SEQUENCE_LENGTH = 393216):

        '''
        Parameters :
        query : a list [chr, start, end, unique_id]
        model : the path to the ENFORMER Model
        fasta_extractor : fasta extractor

        Returns :
        A list of predictions; the first element is the predictions around the TSS for each gene. The second is the prediction across CAGE tracks
        '''

        query_chr = query[0]
        query_start, query_end = int(query[1]), int(query[2])
        query_id = query[3] # i.e. the unique id for the motif

        region = [query_chr, query_start, query_end]
        vcf = VCF(vcf_file, samples=[sample]) 

        # create the region file
        a = create_region_file(vcf_file, region, output_dir=temporary_vcf_dir, software_paths=software_paths, individual=sample)

        try:
            b = extract_individual_sequence(region_details=a, individual=sample, fasta_file_path=fasta_file_path, delete_region=True, fasta_extractor=fasta_extractor)
            target_fa = b['sequence'][sample]
            obj_to_save = model.predict_on_batch(one_hot_encode(target_fa)[np.newaxis])['human'][0]
        except ValueError:
            return({'sequence_source':'NA'})

        # select 8 bins upstream and downstream the center
        obj_to_save = obj_to_save[range(448 - 8, (448 + 8 + 1)), : ].squeeze()

        # make sure to select the 17 bins needed

        if not os.path.isdir(f'{output_dir}/{sample}'):
            os.makedirs(f'{output_dir}/{sample}')

        h5_save = str(f'{output_dir}/{sample}/{query_id}_predictions.h5')
        with h5py.File(h5_save, 'w') as hf:
            hf.create_dataset(query_id, data=obj_to_save)

        # close the vcf file - be a responsible coder
        vcf.close()

        return(b['sequence_source'])

    query = line.strip().split(' ') 
    query_exists = check_query(sample=sam, query_name=query[3], output_dir=output_dir, logfile=logfile)
    print(f'Query status: {query_exists}')

    if query_exists == 0: # i.e. both 'completed' and file exists
        print(f'Query status: {query_exists}')
        print(f'[NOTHING TO DO]\n')
        return([0]) # to the next line
    elif query_exists != 0: # i.e. does not exist; and will create a log file if it does not exist, of course
        print(f'Query status: {query_exists}')
        print(f"[PREDICTING] {query[3]}\n")
        p = predict_query(query=query, output_dir=output_dir, model=model, vcf_file=vcf_file, temporary_vcf_dir=temporary_vcf_dir, fasta_extractor=fasta_extractor, fasta_file_path=fasta_file_path, sample=sample, software_paths=software_paths)
        return([1, p])

def main():

    # get the path of the script as well as parameters         
    whereis_script = os.path.dirname(sys.argv[0])  
    script_path = os.path.abspath(whereis_script) 

    usage_codes = f'{script_path}/enformer-usage-codes.py'
    parsl_config = f'{script_path}/parsl-configuration.py'

    # import the enformer-usage_codes.py file
    exec(open(usage_codes).read(), globals(), globals())
    exec(open(parsl_config).read(), globals(), globals())

    # read the parameters file
    with open(f'{script_path}/../code_arguments/enformer_parameters.json') as f:

        parameters = json.load(f)

        intervals_dir = parameters['interval_list_dir']
        model_path = parameters['model_path']
        fasta_file = parameters['hg38_fasta_file']
        output_dir = parameters['output_dir']
        individuals = parameters['individuals']
        vcf_file = parameters['vcf_file']
        path_to_bcftools = parameters['path_to_bcftools']
        path_to_tabix = parameters['path_to_tabix']
        temporary_vcf_dir = parameters['temporary_vcf_dir']
        TF = parameters['TF']
        logfile_path = parameters['logfile_path']


    # create directories as needed
    if not os.path.exists(output_dir):
        os.makedirs(output_dir) 

    enf_model = Enformer(model_path)
    sequence_extractor = FastaStringExtractor(fasta_file)

    for sam in individuals:
        print(f'\n[OPENING & CHECKING INTERVALS] {sam} === \n')
        with open(f'{intervals_dir}/{sam}_{TF}.txt', 'r') as ind_r:
            all_intervals = ind_r.readlines()

            # read in the log file for this individual
            # doind this so that the log file is not opened everytime
            logfile_csv = f'{logfile_path}/{sam}_predictions_log.csv'

            if os.path.isfile(logfile_csv):
                logfile = pd.read_csv(logfile_csv)
                open_mode = 'a'
            else:
                logfile = None
                open_mode = 'w'

            print('Opening log file\n')
            with open(logfile_csv, open_mode, encoding='UTF8') as running_log_file:
                logwriter = csv.writer(running_log_file)

                if open_mode == 'w':
                    logwriter.writerow(['motif', 'individual', 'status', 'sequence_type']) # write the headers
                
                print('[CREATING] appfutures list\n')
                appfuture_runs = []
                for line in all_intervals: # for each line (i.e. chr, start, end and query_id or motif)
                    query = line.strip().split(' ')
                    appfuture_runs.append(check_and_predict_query(line, sample=sam, query_name=query[3], output_dir=output_dir, logfile=logfile, query=query, model=enf_model, vcf_file=vcf_file, temporary_vcf_dir=temporary_vcf_dir, fasta_extractor=sequence_extractor, fasta_file_path=fasta_file, software_paths=[path_to_bcftools, path_to_tabix]))

                print(f'Length of appfuture list is {len(appfuture_runs)}\n')

                print(appfuture_runs)
                
                print('[DONE] predicting and writing to list\n')
                output = [logwriter.writerow([query[3], sam, 'completed', apf.result()[1]]) if apf.result()[0] == 1 else None for apf in appfuture_runs]
                    
                # if the file gets saved, write out the status
                #logwriter.writerow([query[3], sam, 'completed', ['sequence_source']])

                # always flush
                running_log_file.flush()
                os.fsync(running_log_file)
                        
                    
if __name__ == '__main__':
    main()







# query = line.strip().split(' ') 

                    # query_exists = check_query(sample=sam, query_name=query[3], output_dir=output_dir, logfile=logfile)

                    # if (query_exists.done() == True) & (query_exists.result() == 0): # i.e. both 'completed' and file exists
                    #     print(f'Query status: {query_exists.result()}')
                    #     print(f'[NOTHING TO DO]\n')
                    #     continue # to the next line
                    # elif (query_exists.done() == True) & (query_exists.result() != 0): # i.e. does not exist; and will create a log file if it does not exist, of course
                    #     print(f'Query status: {query_exists.result()}')
                    #     print(f"[PREDICTING] {query[3]}\n")
                    #     p = predict_query(query=query, output_dir=output_dir, model=enf_model, vcf_file=vcf_file, temporary_vcf_dir=temporary_vcf_dir, fasta_extractor=sequence_extractor, fasta_file_path=fasta_file, samples=[sam], software_paths=[path_to_bcftools, path_to_tabix])