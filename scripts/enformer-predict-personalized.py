

# import the enformer-usage_codes.py file
exec(open('/home/temi/imlab_folder/users/temi/projects/TFXcan/scripts/enformer-usage-codes.py').read())

transform_path = 'gs://dm-enformer/models/enformer.finetuned.SAD.robustscaler-PCA500-robustscaler.transform.pkl'
model_path = 'https://tfhub.dev/deepmind/enformer/1'
fasta_file_directory = '/projects/covid-ct/imlab/data/hg19_genome'

bed_files_directory = '/projects/covid-ct/imlab/users/temi/projects/TFXcan/processed-data/motif-bed'

# Helps to quickly generate variant bed files

import tabix, csv

def extract_individual_snps_from_bed(bed_path, loci):
    bed_f = tabix.open(bed_path)
    cur_records = bed_f.query(loci[0],loci[1],loci[2])
    return(pd.DataFrame([x for x in cur_records], columns=['CHROM', 'POS_START', 'POS_END', 'REF', 'ALT', 'CHANGE']))

def collect_intervals(motif_list=['TP1', 'TN9464'], motif_regions_file=None):

    '''
    should take a dictionary of motifs/genes and chromosomes for them e.g {'TP1': ['1', '3'], 'TN5':['4', '6']} >> {'TP1':[['1', 184101981, 184495197], ['3', ..., ...], 'TN5':[[...], [...]]}
    Parameters : 
      gene_list : a list of genes; the genes should be located on those chromosomes

    Returns :
      A dictionary of genes (from gene_list) and their intervals within their respective chromosomes
    '''

    motif_intervals = {} # Collect intervals for our genes of interest

    # may supply a string or a file : should not matter
    if isinstance(motif_regions_file, type('a')):
      motif_regions = pd.read_table(motif_regions_file, sep=' ', names=['chr', 'start', 'end', 'id', 'score', 'strand', 'set', 'num'])
    elif isinstance(motif_regions_file, type(pd.DataFrame)):
      motif_regions = pd.DataFrame(motif_regions_file, columns = ['chr', 'start', 'end', 'id', 'score', 'strand', 'set', 'num'])

    # create motif name
    motif_regions['motif_name'] = motif_regions['set'] + motif_regions['num'].astype(str)
    # if `motif_list` is None, predict on all the motifs available
    if isinstance(motif_list, type(None)):
        print(f'Predicting on all the motifs because a list has not been explictly supplied\n')
        motif_list = motif_regions.motif_name.values

    for i, motif in enumerate(motif_list):
      try:
        temp = motif_regions.loc[motif_regions['motif_name'] == motif]
        motif_intervals[motif] = [re.sub('chr', '', temp['chr'].values[0]), temp['start'].values[0], temp['end'].values[0]]
      except:
        motif_intervals[motif] = 'Not found'

    return(motif_intervals)

def run_predictions(motif_intervals, model, fasta_extractor, individuals_list=['HG00096'], TF='GATA3', bfiles_path = '/projects/covid-ct/imlab/users/saideep/prepare_vcfs/hackathon/full_beds', SEQUENCE_LENGTH = 393216):
    '''
    Parameters :
    motif_intervals : the results from calling `collect_intervals`
    bfiles_path : a dataframe of the individual's variants, the path to where these file can be extracted quickly
    individuals_list : a list of individuals on which we want to make predictions; defaults to None

    Returns :
    A list of predictions; the first element is the predictions around the TSS for each gene. The second is the prediction across CAGE tracks
    '''

    motif_output = dict()
    motif_predictions = dict()

    logfile_csv = '/projects/covid-ct/imlab/users/temi/projects/TFXcan/log/logging_predictions_personalized.csv'
    path_to_save='/projects/covid-ct/imlab/users/temi/projects/TFXcan/output/personalized_predictions'

    if not os.path.exists(path_to_save):
        os.makedirs(path_to_save) 
    else:
        print('Saving folder exists. Moving on...\n')

    if os.path.isfile(logfile_csv):
        # if a log file exists, read it and use that to check the files that are completed
        print('Log file exists. Reading it in ===\n')
        logfile = pd.read_csv(logfile_csv)
        motif_in_logfile = logfile.motif.values
        individual_in_logfile = logfile.individual.values
        print(logfile.head())

        print('Opening in append mode...\n')
        open_mode = 'a'

        
    else:
        logfile = None
        print('No log file. Moving on...')

        print('Opening in write mode...\n')
        open_mode = 'w'

    # if isinstance(logfile, type(pd.DataFrame)):
    #     print('Opening in append mode...\n')
    #     open_mode = 'a'
    # elif isinstance(logfile, type(None)):
    #     print('Opening in write mode...\n')
    #     open_mode = 'w'

    with open(logfile_csv, open_mode, encoding='UTF8') as running_log_file:
        
        logwriter = csv.writer(running_log_file)
        
        if open_mode == 'w':
            logwriter.writerow(['motif', 'individual', 'status']) # write the headers

        for motif, interval_info in motif_intervals.items():
            if not isinstance(interval_info, list): # if the value is 'Not found', just skip it
                continue

            target_interval = kipoiseq.Interval("chr" + interval_info[0],
                                            interval_info[1],
                                            interval_info[2]) # creates an interval to select the right sequences
            #print(target_interval)
                                
            target_fa = fasta_extractor.extract(target_interval.resize(SEQUENCE_LENGTH))  # extracts the fasta sequences, and resizes such that it is compatible with the sequence_length
            #print(target_fa)

            window_coords = target_interval.resize(SEQUENCE_LENGTH) # we also need information about the start and end locations after resizing
            #print(f'Start: {window_coords.start} | End: {window_coords.end}')

            # this is where you read the bed files, per individual
        
            individual_predictions = dict() # where to save the individual's predictions

            # if an individual is supplied - actually, make sure individuals are supplied
            if isinstance(individuals_list, list) or isinstance(individuals_list, type(np.empty([1, 1]))):
                use_individuals = individuals_list
            # elif isinstance(individuals_list, type(None)):
            #     use_individuals = cur_gene_vars.columns[4:]

            for individual in use_individuals: # check if the predictions are available

                # the file is done running
                #print(isinstance(logfile, pd.DataFrame))
                if isinstance(logfile, pd.DataFrame):# should have read it>
                    if (motif in motif_in_logfile) and (individual in individual_in_logfile): # both the individual and motifs are present
                        print(f'{individual}\'s {motif} in logfile already.\n')

                        # check if the motif has been saved and if the status is `completed`
                        motif_saved = str(f'{path_to_save}/{individual}/{TF}_{motif}_personalized_predictions_2022-07-23.pkl')

                        #print(any(logfile.loc[logfile.motif == motif, 'status'].values[0] == 'completed'))
                        if (logfile.loc[(logfile.motif==motif) & (logfile.individual==individual), 'status'].values[0] == 'completed') and os.path.isfile(motif_saved):
                            print(f'{individual}\'s {motif} predictions already available. Moving on...\n')
                            continue
                
                # otherwise, the individual's predictions don't exist

                # create the folder
                if not os.path.exists(f'{path_to_save}/{individual}'):
                    os.makedirs(f'{path_to_save}/{individual}')

                print(f'Currently on individual {individual}, and predicting on motif {motif}...')
                try: # this can fail if there is not tabix file for the bed.gz file in the folder
                    curr_vars = extract_individual_snps_from_bed(bed_path=f'{bfiles_path}/{individual}.bed.gz', loci=interval_info)
                    print('Variant file extacted')
                except:
                    continue

                #print(curr_vars.head(5))

                #two haplotypes per individual
                haplo_1 = list(target_fa[:])
                haplo_2 = list(target_fa[:])

                for i, row in curr_vars.iterrows():
            
                    geno = row['CHANGE'].split("|")
          
                    if (int(row["POS_START"]) - window_coords.start-1) >= len(haplo_2):
                        continue
                    if (int(row["POS_START"]) - window_coords.start-1) < 0:
                        continue
                    if geno[0] == "1":
                        haplo_1[int(row["POS_START"]) - window_coords.start-1] = row["ALT"]
                    if geno[1] == "1":
                        haplo_2[int(row["POS_START"]) - window_coords.start-1] = row["ALT"]

                #predict on the individual's two haplotypes
                prediction_1 = model.predict_on_batch(one_hot_encode("".join(haplo_1))[np.newaxis])['human'][0]
                prediction_2 = model.predict_on_batch(one_hot_encode("".join(haplo_2))[np.newaxis])['human'][0]

                predictions = [prediction_1, prediction_2]

                # save the predictions
                with open(str(f'{path_to_save}/{individual}/{TF}_{motif}_personalized_predictions_2022-07-23.pkl'), 'wb') as obj:
                    pickle.dump(predictions, obj)

                # if the file gets saved, write out the status
                logwriter.writerow([motif, individual, 'completed'])
                running_log_file.flush()
                
                print('Done with prediction and/or saving\n')

def main():
    
    # '/projects/covid-ct/imlab/users/temi/projects/TFXcan/processed-data/Gata3_motif_regions_TEMP.txt'
    GATA3_motif_regions = sys.argv[1]

    TF = sys.argv[2] # a supplied transcription factor
    
    # load the model
    model = Enformer(model_path)
    print('Model loaded.\n')

    # instantiate fasta string extractor
    fasta_extractor = FastaStringExtractor(fasta_file_directory + '/genome.fa')
    print(f'Fasta strings extracted from "{fasta_file_directory}/genome.fa".\n')
    
    # create the intervals
    # int_ = create_intervals(Gata3_motif_regions)
    # print('Intervals created\n')

    # collect the intervals
    motif_int = collect_intervals(motif_list=None, motif_regions_file=GATA3_motif_regions)
    print('Intervals collected\n')
    
    # make predictions
    print('\nMaking predictions ===\n')
    
    #path_to_save='/projects/covid-ct/imlab/users/temi/projects/TFXcan/output/reference_predictions'
    motif_preds = run_predictions(motif_intervals=motif_int, model=model, fasta_extractor=fasta_extractor, individuals_list = ['HG00479', 'HG00590'])

    
    print('\nEverything done!')

if __name__ == '__main__':
    main()
    