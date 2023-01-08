# Usage: This module is used to predict with ENFORMER on individual genomes
# Author: Temi
# Date: 

import parsl
from parsl.app.app import python_app

#@python_app
def run_batch_predictions(batch_regions, batch_num, id, script_path, vcf_func, output_dir, logfile, predictions_log_file): #
    """
    Predict and save on a given batch of regions in the genome

    This function also filters the regions in the batch for those that have been predicted and logged

    Parameters:
        batch_regions: list
            A list of regions with each element (region) in the form `chr_start_end`.
        batch_num: num
            The number of the batch e.g. batch 1, batch 2, e.t.c.
        id: str
            The unique id for predictions or the name of the individual, used to save the and name the folders (e.g. 'kawakami', 'freedman', )
        output_dir: str (path)
            Where should the predictions be saved? Predictions are saved as `{id}/{region}_predictions.h5` 
        logfile: str (path)
            When predictions are made, they are logged in this file. Also useful for checking the log to prevent re-predicting an existing prediction. 
        script_path: str (path), default is the path to where this file is.
            The path to this module.
        predictions_log_dir: str (path)
            This is the folder wherein predictions will be logged.

    Returns: num
        A single value of either 0 (if predictions were successful) or 1 (if there was a error).
        Check the call logs or stacks for the source of the error. 
    """
  
    import sys, os, tqdm, faulthandler, time

    mpath = os.path.join(script_path, 'batch_utils') #os.path.dirname(__file__) #
    sys.path.append(mpath)
    
    faulthandler.enable() # to figure out where segmentation error is coming from

    global vcf_id
    vcf_id = id

    try:
        import predictUtils_two
    except ModuleNotFoundError as merr:
        print(f'[ERROR] of type {type(merr)} at run_batch_predictions')

    # first check the queries
    check_result = (predictUtils_two.check_query(sample = id, query = each_region, output_dir=output_dir, logfile=logfile) for each_region in tqdm.tqdm(batch_regions, desc=f'[INFO] Checking query for batch'))
    # filter for nones
    filtered_check_result = [r for r in check_result if r is not None]

    if not filtered_check_result: # i.e. if the list is empty
        return(1)
    else:

        tic = time.perf_counter()

        reg_prediction = predictUtils_two.enformer_predict_on_batch(batch_region=filtered_check_result, sample=id, model=None, vcf_func=vcf_func, output_dir=output_dir, predictions_log_file=predictions_log_file, batch_num=batch_num)
        
        toc = time.perf_counter()

        print(f'[INFO] (time) to predict on this batch is {toc - tic}')
        return(reg_prediction) # returns 0 returned by enformer_predict

def return_prediction_function(use_parsl, fxn=run_batch_predictions):
    '''
    Decorate or not the `run_batch_predictions` function based on whether `use_parsl` is true or false

    Returns: 
        function object
        The function if parsl is not used
        The parsl decorated function if parsl is meant to be used
    
    '''
    from parsl.app.app import python_app
    if use_parsl == True:
        return python_app(fxn)
    elif use_parsl == False:
        return fxn

def generate_batch(lst, batch_n, len_lst = None):
    """
    Given a list, this function yields batches of an unspecified size but the number of batches is equal to `batch_n`
    E.g. generate_batch([0, 1, 2, 3, 4, 5, 6], batch_n=2) -> (0, 1, 2, 3), (4, 5, 6)
    
    Parameters:
        lst: list
        batch_n: int
            Number of batches to return
        len_lst: None or num (length of the input list)

    Yields
        `batch_n` batches of the list
    """
    import math
    # how many per batch
    if len_lst is not None:
        n_elems = math.ceil(len_lst/batch_n)
    else:
        n_elems = math.ceil(len(lst)/batch_n)
        
    for i in range(0, len(lst), n_elems):
        yield lst[i:(i + n_elems)]