
import parsl
from parsl.app.app import python_app, join_app

def check_query(sample, query_name, output_dir, logfile):
        '''
        Checks of predictions are available for an individual or sample
        '''

        import pandas as pd

        if isinstance(logfile, pd.DataFrame):# should have read it>
            print('Logfile is a dataframe')
            motif_in_logfile = logfile.motif.values
            individual_in_logfile = logfile.individual.values
            query_saved = str(f'{output_dir}/{sample}/{query_name}_predictions.h5')

            # four conditions 
            a = (query_name in motif_in_logfile)
            b = (sample in individual_in_logfile)
            c = (logfile.loc[(logfile.motif==query_name) & (logfile.individual==sample), 'status'].values[0] == 'completed')
            d = os.path.isfile(query_saved)

            if a and b and c and d: 
                return('yes')
            else:
                print('No predictions')
                return('no')
        elif isinstance(logfile, type(None)):
            print('Logfile does not exist')
            return('no')

def evaluate_check(is_query_done, b=None):

    '''
    b is either a dict or none
    '''
    if is_query_done == 'yes':
        return([0, 'NA'])
    elif is_query_done == 'no':
        return([1, b['sequence_source']])


@python_app
def predict_query(query, sample, logfile, output_dir, model, fasta_extractor, open_vcf_file, temporary_vcf_dir, fasta_file_path, software_paths=[], SEQUENCE_LENGTH = 393216):

    '''
    Parameters :
    query : a list [chr, start, end, unique_id]
    model : the path to the ENFORMER Model
    fasta_extractor : fasta extractor

    Returns :
    A numpy matrix of ENFORMER predictions for that region
    '''
    import pandas as pd

    # check if the query has been predicted first
    is_query_done = check_query(sample=sample, query_name=query[3], output_dir=output_dir, logfile=logfile)

    if is_query_done == 'yes':
        print(f'[NOTHING TO DO] {query[3]}\n')
        return(evaluate_check(is_query_done))

    elif is_query_done == 'no':
        print(f'[SOMETHING TO DO] {query[3]}\n')
        query_chr = query[0]
        query_start, query_end = int(query[1]), int(query[2])
        query_id = query[3] # i.e. the unique id for the motif
        region = [query_chr, query_start, query_end] 

        # create the region file
        a = create_region_file(open_vcf=open_vcf_file, region=region, subset_vcf_dir=temporary_vcf_dir, individual=sample, software_paths=software_paths)

        print(f'a\'s results are: {a}')

        try:
            b = extract_individual_sequence(subset_dict=a, fasta_file_path=fasta_file_path, fasta_extractor=fasta_extractor, delete_region=True)
            target_fa = b['sequence'][sample]
            obj_to_save = model.predict_on_batch(one_hot_encode(target_fa)[np.newaxis])['human'][0]
            # select 8 bins upstream and downstream the center
            obj_to_save = obj_to_save[range(448 - 8, (448 + 8 + 1)), : ].squeeze()
            h5_save = str(f'{output_dir}/{sample}/{query_id}_predictions.h5')
            with h5py.File(h5_save, 'w') as hf:
                hf.create_dataset(query_id, data=obj_to_save)

            return(evaluate_check(is_query_done, b))
        except ValueError:
            return(evaluate_check(is_query_done, b={'sequence_source':'NA'}))