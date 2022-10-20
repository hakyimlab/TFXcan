
# import the enformer-usage_codes.py file
exec(open('/home/temi/imlab_folder/users/temi/projects/TFXcan/scripts/enformer-usage-codes.py').read())
<<<<<<< HEAD

transform_path = 'gs://dm-enformer/models/enformer.finetuned.SAD.robustscaler-PCA500-robustscaler.transform.pkl'
#model_path = 'https://tfhub.dev/deepmind/enformer/1'
model_path = '/projects/covid-ct/imlab/data/enformer/raw'
fasta_file_directory = '/projects/covid-ct/imlab/data/hg19_genome'
bed_files_directory = '/projects/covid-ct/imlab/users/temi/projects/TFXcan/processed-data/motif-bed'

# Helps to quickly generate variant bed files

import tabix, csv, h5py
=======
>>>>>>> e15d2c58b1ffd0f4fd0508b4b72f37aad9dbebf7

def create_intervals(motif_region_path):
    #collect_intervals()
    try:
        motif_regions = pd.read_table(motif_region_path, sep=' ')
    except:
        motif_regions = pd.read_csv(motif_region_path)
    print(motif_regions.head(5))
    #motif_regions['motif_name'] = motif_regions['set'] + motif_regions['num'].astype(str)

    temp_chr = [re.sub(pattern='chr', repl='', string=cc) for cc in motif_regions.chr.tolist()]
    temp_motif_list = motif_regions.motif_name.tolist()
    
    return (temp_chr, temp_motif_list)

def collect_intervals(reg, chromosomes = ["1"], motif_list=['TP1']):

    motif_intervals = {} # Collect intervals for our genes of interest
    try:
        motif_regions = pd.read_table(reg, sep=' ') # names=['chr', 'motif_center_start', 'motif_center_end', 'id', 'score', 'strand', 'start', 'end', 'motif_name']
    except:
        motif_regions = pd.read_csv(reg) # names=['chr', 'motif_center_start', 'motif_center_end', 'id', 'score', 'strand', 'start', 'end', 'motif_name']

    print(motif_regions.head(5))

    #motif_regions['motif_name'] = motif_regions['set'] + motif_regions['num'].astype(str)
    for i, motif in enumerate(motif_list):
        temp = motif_regions.loc[motif_regions['motif_name'] == motif]
        motif_intervals[motif] = [chromosomes[i], temp['motif_center_start'].values[0], temp['motif_center_end'].values[0]]
    return(motif_intervals)


def run_predictions(motif_intervals, model, fasta_extractor, experiment_details, TF='', cell_line='', path_to_save='/projects/covid-ct/imlab/users/temi/projects/TFXcan/enformer_predictions', SEQUENCE_LENGTH = 393216):
    '''
    Parameters :
    gene_intervals : the results from calling `collect_intervals`
    tss_dataframe : a list of the TSSs dataframes i.e. the TSS for the genes in the chromosomes
    individuals_list : a list of individuals on which we want to make predictions; defaults to None

    Returns :
    A list of predictions; the first element is the predictions around the TSS for each gene. The second is the prediction across CAGE tracks
    '''

    #current_motif = 0
    path_to_save = path_to_save + '/' + experiment_details
    # create the folder
    if not os.path.exists(path_to_save):
        os.makedirs(path_to_save)
    
    logfile_csv = f'/projects/covid-ct/imlab/users/temi/projects/TFXcan/log/log_ref_{experiment_details}_{TF}_{cell_line}.csv'
    
    if os.path.isfile(logfile_csv):
        # if a log file exists, read it use that to check the files that are completed
        print('Log file exists. Reading it in ===\n')
        logfile = pd.read_csv(logfile_csv)
        motif_in_logfile = logfile.motif.values
        print('Opening in append mode...\n')
        open_mode = 'a'
    else:
        logfile = None
        print('No log file. Moving on...')
        print('Opening in write mode...\n')
        open_mode = 'w'
        
    # if isinstance(logfile, type(None)):
    #     print('Opening in write mode...\n')
    #     open_mode = 'w'
    # elif isinstance(logfile, pd.DataFrame):
    #     print('Opening in append mode...\n')
    #     open_mode = 'a'
    
    with open(logfile_csv, open_mode, encoding='UTF8') as running_log_file:
        logwriter = csv.writer(running_log_file)
        
        if open_mode == 'w':
            logwriter.writerow(['motif', 'status']) # write the headers
        
        for motif, interval_info in motif_intervals.items():

            if not isinstance(interval_info, list): # if the value is 'Not found', just skip it
                continue

            # the file to be saved 
            h5_save = str(f'{path_to_save}/{TF}_{cell_line}_reference_{motif}_predictions.h5')
            
            # the file is done running
            #print(isinstance(logfile, pd.DataFrame))
            if isinstance(logfile, pd.DataFrame):# should have read it>
                if motif in motif_in_logfile:
                    print(f'{motif} in logfile already.\n')
                    h5_save = str(f'{path_to_save}/{TF}_{cell_line}_reference_{motif}_predictions.h5')
                    try:
                        if (logfile.loc[logfile.motif == motif, 'status'].values[0] == 'completed') and os.path.isfile(h5_save):
                            print(f'{motif} in logfile and predictions available. Moving on...\n')
                            continue
                    except:
                        continue

            #current_motif += 1
            print(f'Predicting on motif {motif} ===\n')

            motif_interval = motif_intervals[motif]
            target_interval = kipoiseq.Interval("chr" + interval_info[0],
                                            interval_info[1],
                                            interval_info[2]) # creates an interval to select the right sequences

            target_fa = fasta_extractor.extract(target_interval.resize(SEQUENCE_LENGTH))  # extracts the fasta sequences, and resizes such that it is compatible with the sequence_length

            #window_coords = target_interval.resize(SEQUENCE_LENGTH) # we also need information about the start and end locations after resizing

            # predict on the individual's two haplotypes
            #output[motif] = model.predict_on_batch(one_hot_encode(target_fa)[np.newaxis])['human'][0]
            obj_to_save = model.predict_on_batch(one_hot_encode(target_fa)[np.newaxis])['human'][0]

            # save the predictions
            # with open(str(f'{path_to_save}/{TF}_{cell_line}_reference_{motif}_predictions.pkl'), 'wb') as obj:
            #     pickle.dump(obj_to_save, obj)
            h5_save = str(f'{path_to_save}/{TF}_{cell_line}_reference_{motif}_predictions.h5')
            with h5py.File(h5_save, 'w') as hf:
                hf.create_dataset(experiment_details + ',' + motif, data=obj_to_save)

            # if the file gets saved
            logwriter.writerow([motif, 'completed'])
            running_log_file.flush()
            
        print('Done with prediction and/or saving\n')

    #return(output)


# def Save_predictions_as_pickle(obj_to_save, path_to_save='/projects/covid-ct/imlab/users/temi/projects/TFXcan/output/reference_predictions'):
#     ## Save the results as a pickle object

#     #print(str(f'Saving predictions to {path_to_save}/GATA3_four_motifs_enformer_predictions_2022-07-14.pkl ===\n'))

#     with open(str(f'{path_to_save}/GATA3_reference_motifs_predictions_2022-07-19.pkl'), 'wb') as obj:
#         pickle.dump(obj_to_save, obj)

#     print('Done saving the predictions')

# def save_predictions_as_hdf5(obj_to_save, path_to_save='/projects/covid-ct/imlab/users/temi/projects/TFXcan/output/reference_predictions'):
#     ## Save the results as a pickle object

#     #print(str(f'Saving predictions to {path_to_save}/GATA3_four_motifs_enformer_predictions_2022-07-14.pkl ===\n'))

#     with open(str(f'{path_to_save}/GATA3_reference_motifs_predictions_2022-07-19.pkl'), 'wb') as obj:
#         pickle.dump(obj_to_save, obj)

#     print('Done saving the predictions')
    
def main():
    
#     TF_motif_regions = sys.argv[1]
#     print(TF_motif_regions)
    
#     individuals = sys.argv[2:]
#     print(individuals)

    # for Gata3
    # '/projects/covid-ct/imlab/users/temi/projects/TFXcan/processed-data/Gata3_motif_regions_TEMP.txt'
    motif_regions = sys.argv[1] # a bed file containing the regions you want to predict for
    TF = sys.argv[2] # a supplied transcription factor
    cell_line = sys.argv[3] # cell line
    experiment_details = sys.argv[4] # another discriminating information
    
    # load the model
    model = Enformer(model_path)
    print('Model loaded.\n')

    # instantiate fasta string extractor
    fasta_extractor = FastaStringExtractor(fasta_file_directory + '/genome.fa')
    print(f'Fasta strings extracted from "{fasta_file_directory}/genome.fa".\n')
    
    # create the intervals
    motif_int = create_intervals(motif_regions)
    print('Intervals created\n')

    # collect the intervals
    motif_int = collect_intervals(motif_regions, motif_int[0], motif_int[1])
    print('Intervals collected\n')
    
    # make predictions
    print('\nMaking predictions ===\n')
    
    #path_to_save='/projects/covid-ct/imlab/users/temi/projects/TFXcan/output/reference_predictions'
    motif_preds = run_predictions(motif_intervals=motif_int, fasta_extractor=fasta_extractor, model=model, experiment_details=experiment_details, TF=TF, cell_line=cell_line)

    #return(motif_preds)
    
    print('Everything done!')


# def main2():
#     print(__name__)

if __name__ == '__main__':
    main()
    #Save_predictions_as_pickle(res)
    
    
