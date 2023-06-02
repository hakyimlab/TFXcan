
# === This script contains all codes to test if tensorflow and ENFORMER are working accurately  
# Modified by Temi from Deepmind people 
# DATE: Thursday Oct 20 2022

from __future__ import absolute_import, division, print_function, unicode_literals
import os, re, sys, json, h5py, csv, warnings
# import parsl
# from parsl.app.app import python_app, bash_app

def check_query_reference(sample, query_name, output_dir, logfile):
    '''
    Checks of predictions are available for an individual or sample
    '''

    if isinstance(logfile, pd.DataFrame):# should have read it>

        # check if the motif has been saved and if the status is `completed`
        query_saved = str(f'{output_dir}/{query_name}_predictions.h5')

        motif_in_logfile = logfile.motif.values
        individual_in_logfile = logfile.individual.values

        if (query_name in motif_in_logfile) and (sample in individual_in_logfile): # both the individual and motifs are present
            if (logfile.loc[(logfile.motif==query_name) & (logfile.individual==sample), 'status'].values[0] == 'completed') and os.path.isfile(query_saved):
                return(0)
            elif (logfile.loc[(logfile.motif==query_name) & (logfile.individual==sample), 'status'].values[0] != 'completed') or (not os.path.isfile(query_saved)):
                return(1)

    elif isinstance(logfile, type(None)):
        return(2)


def predict_query_reference(query, output_dir, model, fasta_extractor, SEQUENCE_LENGTH = 393216):

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
    query_id = query[3]


    try:
        target_interval = kipoiseq.Interval(query_chr, query_start, query_end) # creates an interval to select the right sequences
        target_fa = fasta_extractor.extract(target_interval.resize(SEQUENCE_LENGTH))  # extracts the fasta sequences, and resizes such that it is compatible with the sequence_length
        obj_to_save = model.predict_on_batch(one_hot_encode(target_fa)[np.newaxis])['human'][0]
    except ValueError:
        return('NA')
    
    # select 8 bins upstream and downstream the center
    obj_to_save = obj_to_save[range(448 - 8, (448 + 8 + 1)), : ].squeeze()

    # make sure to select the 17 bins needed

    h5_save = str(f'{output_dir}/{query_id}_predictions.h5')
    with h5py.File(h5_save, 'w') as hf:
        hf.create_dataset(query_chr, data=obj_to_save)

    return('ref')


def main():

    # get the path of the script as well as parameters         
    whereis_script = os.path.dirname(sys.argv[0])  
    script_path = os.path.abspath(whereis_script) 

    # read the parameters file
    with open(f'{script_path}/../data/freedman/enformer-reference-parameters.json') as f:

        parameters = json.load(f)

        model_path = parameters['model_path']
        fasta_file = parameters['hg38_fasta_file']
        output_dir = parameters['output_dir']
        TF = parameters['TF']
        logfile_path = parameters['logfile_path']
        interval_list_dir = parameters['interval_list_dir']

        # additionally
        base_name = parameters['base_name']

        # create directories as needed
        if not os.path.exists(output_dir):
            os.makedirs(output_dir) 

        if not os.path.exists(logfile_path):
            os.makedirs(logfile_path) 

        usage_codes = f'{script_path}/enformer-usage-codes.py'
        #parsl_config = f'{script_path}/parsl-configuration.py'

        # import the enformer-usage_codes.py file
        exec(open(usage_codes).read(), globals(), globals())
        #exec(open(parsl_config).read(), globals(), globals())

        enf_model = Enformer(model_path)
        sequence_extractor = FastaStringExtractor(fasta_file)

        with open(f'{interval_list_dir}/{base_name}_{TF}.txt', 'r') as ind_r:
            all_intervals = ind_r.readlines()

            # read in the log file for this individual
            # doind this so that the log file is not opened everytime
            logfile_csv = f'{logfile_path}/{base_name}_predictions_log.csv'

            if os.path.isfile(logfile_csv):
                logfile = pd.read_csv(logfile_csv)
                open_mode = 'a'
            else:
                logfile = None
                open_mode = 'w'

            with open(logfile_csv, open_mode, encoding='UTF8') as running_log_file:
                logwriter = csv.writer(running_log_file)

                if open_mode == 'w':
                    logwriter.writerow(['motif', 'individual', 'status', 'sequence_type']) # write the headers

                for line in all_intervals: # for each line (i.e. chr, start, end and query_id or motif)
                    query = line.strip().split(' ') 

                    print(f"[PREDICTING] {query[3]}")

                    # you should check that the file does not exist here
                    #print(logfile)
                    query_exists = check_query_reference(sample=base_name, query_name=query[3], output_dir=output_dir, logfile=logfile)

                    if query_exists == 0: # i.e. both 'completed' and file exists
                        print(f'[NOTHING TO DO]\n')
                        continue # to the next interval
                    else: # i.e. does not exist; and will create a log file if it does not exist, of course
                        p = predict_query_reference(query=query, output_dir=output_dir, model=enf_model, fasta_extractor=sequence_extractor)

                        print(f"[DONE]\n")

                        # if the file gets saved, write out the status
                        logwriter.writerow([query[3], base_name, 'completed', p])

                    # always flush
                    running_log_file.flush()
                    os.fsync(f)
                                
if __name__ == '__main__':
    main()