# Usage: This module is used to predict with ENFORMER on individual genomes
# Author: Temi
# Date: Thurs Feb 2 2023


#@python_app
def run_batch_predictions(batch_regions, samples, path_to_vcf, batch_num, script_path, output_dir, prediction_logfiles_folder, sequence_source): #
    """
    Predict and save on a given batch of regions in the genome

    This function also filters the regions in the batch for those that have been predicted and logged

    Parameters:
        batch_regions: list
            A list of regions with each element (region) in the form `chr_start_end`.
        batch_num: num
            The number of the batch e.g. batch 1, batch 2, e.t.c.
        samples: list
            a list of samples: [a, b, c]
        output_dir: str (path)
            Where should the predictions be saved? Predictions are saved as `{sample}/{haplotype0, haplotype1, haplotype2}/{region}_predictions.h5` 
        path_to_vcf: str (path)
            The path to the vcf file
        prediction_logfiles_folder: str (path)
            When predictions are made, they are logged in this folder and saved as {sample}_log.csv. Also useful for checking the log to prevent re-predicting an existing prediction. 
        script_path: str (path), default is the path to where this file is.
            The path to this module.
        sequence_source: one of 'personalized' or 'reference'

    Returns: num
        A single value of either 0 (if predictions were successful) or 1 (if there was a error).
        Check the call logs or stacks for the source of the error. 
    """

    import sys, os, faulthandler, time, itertools

    mpath = os.path.join(script_path, 'batch_utils') #os.path.dirname(__file__) #
    sys.path.append(mpath)
    
    faulthandler.enable() # to figure out where segmentation error is coming from

    try:
        import checksUtils
        import predictUtils_two
    except ModuleNotFoundError as merr:
        raise Exception(f'[ERROR] {type(merr).__name__} at run_batch_predictions. Cannot locate either of `checkUtils` or `predictUtils_two` modules.')

    # check_these = itertools.product(samples, [batch_regions])
    # check_results = [checksUtils.check_queries(sample=cq[0], queries=cq[1], output_dir=output_dir, prediction_logfiles_folder=prediction_logfiles_folder, sequence_source=sequence_source) for cq in check_these]

    check_results = {sample: checksUtils.check_queries(sample=sample, queries=batch_regions, output_dir=output_dir, prediction_logfiles_folder=prediction_logfiles_folder, sequence_source=sequence_source) for sample in samples}

    filtered_check_result = {k: v for k, v in check_results.items() if k is not None}
   # print(filtered_check_result)

    #print(filtered_check_result)

    # filter out nones
    # filtered_check_result = [r for r in check_results if r is not None]
    if not filtered_check_result: # i.e. if the list is empty
        return(1)
    else:
        # unlist the list of dictionaries
        # filtered_check_result = [f for f in filtered_check_result for f in f]
        # filtered_check_result = list({v['query']: v for v in filtered_check_result}.values())

        pqueries = [v for k, v in filtered_check_result.items()]
        pqueries = [l for l in pqueries for l in l]
        pqueries = list(set([d['query'] for d in pqueries]))
        # print(f'{len(pqueries)}')
        #print(pqueries)

        if pqueries:
            samples = list(filtered_check_result.keys())
            tic = time.perf_counter()

            reg_prediction = predictUtils_two.enformer_predict_on_batch(batch_regions=pqueries, samples=samples, logging_dictionary=filtered_check_result, path_to_vcf = path_to_vcf, output_dir=output_dir, prediction_logfiles_folder=prediction_logfiles_folder, batch_num=batch_num, sequence_source=sequence_source)
            
            toc = time.perf_counter()
            print(f'[INFO] (time) to predict on batch {batch_num} is {toc - tic}')
            return(reg_prediction) # returns 0 returned by enformer_predict
        else:
            return(1)

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

def generate_n_batches(lst, batch_n, len_lst = None):
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

def generate_batch_n_elems(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]


def make_h5_db(h5_file, csv_file, files_list, files_path, dataset):
    import h5py
    import pandas as pd

    with h5py.File(f"{h5_file}", "w") as f_dst:
        #h5files = [f for f in os.listdir(f'{project_dir}/') if f.endswith(".h5")]

        dset = f_dst.create_dataset(f'{dataset}_dataset', shape=(len(files_list), 17, 5313), dtype='f4')
        for i, filename in enumerate(files_list):
            with h5py.File(files_path[i]) as f_src:
                dset[i] = f_src[filename]

    pd.DataFrame(files_list, columns=['motif']).to_csv(f"{csv_file}")

    return(0)

def make_h5_db_parsl(use_parsl, fxn=make_h5_db):
    '''
    Decorate or not the `make_h5_db` function based on whether `use_parsl` is true or false

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