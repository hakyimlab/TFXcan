import parsl
from parsl.app.app import python_app

@python_app
def run_single_predictions(region, individual, vcf_func, script_path, output_dir, logfile, predictions_log_dir): #
    """
    Given a region for an individual, check if predictions have been made and logged, and if not, retrieve the sequences, change the variants as directly by the vcf file, predict on that region, and save the predictions.

    Parameters:
        region: str
            A region in the genome in the form `chr_start_end`.
        individual: str
            The name/id of an individual
        vcf_func: cyvcf object
            A function that returns a cyvcf object
        script_path: str (path)
            The path to this script (This is needed to add the utilities directory to the path for the module to be imported)
        output_dir: str (path)
            The folder where predictions should be stored.
        logfile: pd.DataFrame or None (if no log file exists)
            The log file holding the logged predictions
        predictions_log_dir: str (path)
            The folder where
    
    Returns:

    """
    import sys, os
    #sys.path.append(f'{script_path}/utilities')
    mpath = os.path.join(script_path, 'utilities') #os.path.dirname(__file__) #
    sys.path.append(mpath)
    print(sys.path)
    
    try:
        import predictUtils_two
    except ModuleNotFoundError as merr:
        print(f'[MODULE NOT FOUND ERROR] at run_single_predictions')

    #first check the query
    check_result = predictUtils_two.check_query(sample = individual, query = region, output_dir=output_dir, logfile=logfile)

    if check_result is not None:
        #print(check_result)
        b = predictUtils_two.create_individual_input_for_enformer(region=check_result['query'], individual=individual, vcf_func=vcf_func, fasta_func=None, hap_type = 'hap1', resize_for_enformer=True, resize_length=None)

        #print(f'[CACHE INFO] (vcf) {vcf_func.cache_info()}')

        if (b is not None) and (len(b['sequence']) == 393216): #(b['sequence'] is not None) and (len(b['sequence']) == 393216):
            #print(type(b))
            reg_prediction = predictUtils_two.enformer_predict(b['sequence'], region=b['region'], sample=individual, seq_type=b['sequence_source'], model_func=None, output_dir=output_dir, predictions_log_dir=predictions_log_dir, logtype=check_result['logtype'])
            
            return(reg_prediction)
        else:
            print(f"[WARNING] {check_result}: Either length of input sequence is invalid (NoneType) or too long or too short")
            

def generate_batch(lst, batch_size):
    """  
    Given a list, this function yields batches of a specified size
    
    Parameters:
        lst: list
        batch_size: int
            Number of items in each batch.

    Yields
        Batches of the list containing `batch_size` elements.
    """
    if batch_size <= 0:
        return None
    for i in range(0, len(lst), batch_size):
        yield lst[i:(i + batch_size)]