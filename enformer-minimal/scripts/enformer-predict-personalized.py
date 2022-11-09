
# === This script contains all codes to test if tensorflow and ENFORMER are working accurately  
# Created by Temi
# DATE: Thursday Oct 20 2022

from __future__ import absolute_import, division, print_function, unicode_literals
import os, re, sys, json, h5py, csv, warnings
from cyvcf2 import VCF
import parsl
from parsl.app.app import python_app
import pandas as pd

def main():

    # get the path of the script as well as parameters         
    whereis_script = os.path.dirname(sys.argv[0])  
    script_path = os.path.abspath(whereis_script) 

    usage_codes = f'{script_path}/enformer-usage-codes.py'
    parsl_config = f'{script_path}/parsl-configuration.py'
    personal_enformer = f'{script_path}/personal-enformer.py'

    # import the enformer-usage_codes.py file
    exec(open(usage_codes).read(), globals(), globals())
    exec(open(parsl_config).read(), globals(), globals())
    exec(open(personal_enformer).read(), globals(), globals())

    # read the parameters file
    with open(f'{script_path}/../metadata/enformer_parameters.json') as f:

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

    enf_model = Enformer(model_path)
    sequence_extractor = FastaStringExtractor(fasta_file)
    
    for sam in individuals:
        # open the vcf file
        variants_vcf = VCF(vcf_file, samples=sam)
        # create the directory for this individual 
        if not os.path.exists(f'{output_dir}/{sam}'):
            print(f'\n[CREATING DIR] at {output_dir}/{sam}')
            os.makedirs(f'{output_dir}/{sam}')

        print(f'\n[LOADING INTERVALS] {sam}')

        # create the zip object for all the intervals
        a = pd.read_table(f'{intervals_dir}/{sam}_{TF}.txt', sep=' ', header=None)
        b = a[0].tolist() # a list of queries

        # I need a log file
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

            # for each interval check if prediction has been done
            query_status = []
            print(f'Type of query status before filling is: {type(query_status)}\n')
            for i, query in enumerate(b):
                query_status.append(check_query(sample=sam, query=query, output_dir=output_dir, logfile=logfile))
            # evaluate the results >> this should return a list of those that don't have predictions
            q_status = [q.result() for q in query_status]
            print(f'Query results are: {q_status} ======\n')

            # predict_status = []
            # for q in q_status:
            #     if q[0] == 'no':
            #         predict_status.append(predict_query(query = q[1:], sample=sam, logfile=logfile, output_dir=output_dir, model=enf_model, fasta_extractor=sequence_extractor, open_vcf_file=vcf_file, temporary_vcf_dir=temporary_vcf_dir, fasta_file_path=fasta_file, script_path=script_path, software_paths=[path_to_bcftools, path_to_tabix], SEQUENCE_LENGTH = 393216))

            # p_status = [p.result() for p in predict_status]
            # print(f'Predict results are: {p_status}')

                # if is_query_done.result() == 'yes':
                #     print(f'[NOTHING TO DO] {query[3]}\n')
                #     #return(evaluate_check(is_query_done))

                # elif is_query_done.result() == 'no':
                #     #print(f'[SOMETHING TO DO] {query[3]}\n')

                # return_value = predict_query(query, sample=sam, logfile=logfile, output_dir=output_dir, model=enf_model, fasta_extractor=sequence_extractor, open_vcf_file=vcf_file, temporary_vcf_dir=temporary_vcf_dir, fasta_file_path=fasta_file, software_paths=[path_to_bcftools, path_to_tabix], SEQUENCE_LENGTH = 393216)

                # predict_query_list.append(return_value.result())

                # print(f'Length of apps future list is {len(predict_query_list)}')
                #print(predict_query_list)

                # # if the file gets saved, write out the status
                # if return_value[0] == 1:
                #     logwriter.writerow([query[3], sam, 'completed', return_value[1]])

                # # always flush
                # running_log_file.flush()
                # os.fsync(running_log_file)

        print(f'[FINISHED] {sam}\n')

        # close the vcf file when done with both sams
        variants_vcf.close()
    print(f'[FINISHED BOTH INDIVIDUALS]')                  
                    
if __name__ == '__main__':
    main()