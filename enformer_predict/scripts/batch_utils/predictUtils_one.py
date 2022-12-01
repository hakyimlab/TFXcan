import parsl
from parsl.app.app import python_app

@python_app
def run_batch_predictions(batch_regions, batch_num, individual, vcf_func, script_path, output_dir, logfile, predictions_log_dir): #
  
    import sys, os, tqdm, faulthandler, time
    mpath = os.path.join(script_path, 'batch_utils') #os.path.dirname(__file__) #
    sys.path.append(mpath)
    
    faulthandler.enable() # to figure out where segmentation error is coming from

    try:
        import predictUtils_two
    except ModuleNotFoundError as merr:
        print(f'[ERROR] of type {type(merr)} at run_batch_predictions')

    # first check the queries
    check_result = (predictUtils_two.check_query(sample = individual, query = each_region, output_dir=output_dir, logfile=logfile) for each_region in tqdm.tqdm(batch_regions, desc=f'[INFO] Checking query for batch'))
    # filter for nones
    filtered_check_result = [r for r in check_result if r is not None]

    if not filtered_check_result: # i.e. if the list is empty
        return(1)
    else:

        tic = time.perf_counter()

        reg_prediction = predictUtils_two.enformer_predict(batch_region=filtered_check_result, sample=individual, model=None, output_dir=output_dir, vcf_func=vcf_func, predictions_log_dir=predictions_log_dir, batch_num=batch_num)
        
        toc = time.perf_counter()

        print(f'[INFO] (time) to predict on this batch is {toc - tic}')
        return(reg_prediction) # returns 0 returned by enformer_predict

# else:
#     print(f"[WARNING] {check_result}: Either length of input sequence is invalid (NoneType) or too long or too short")

def return_prediction_function(use_parsl, fxn=run_batch_predictions):
    from parsl.app.app import python_app
    if use_parsl == True:
        return python_app(fxn)
    elif use_parsl == False:
        return fxn

# def generate_batch(lst, batch_size):
#     """  
#     Given a list, this function yields batches of a specified size
    
#     Parameters:
#         lst: list
#         batch_size: int
#             Number of items in each batch.

#     Yields
#         Batches of the list containing `batch_size` elements.
#     """
#     if batch_size <= 0:
#         return None
#     for i in range(0, len(lst), batch_size):
#         yield lst[i:(i + batch_size)]

def generate_batch(lst, batch_n, len_lst = None):
    """
    Given a list, this function yields batches of a specified size
    
    Parameters:
        lst: list
        batch_size: int
            Number of items in each batch.

    Yields
        Batches of the list containing `batch_size` elements.
    """
    import math
    # how many per batch
    if len_lst is not None:
        n_elems = math.ceil(len_lst/batch_n)
    else:
        n_elems = math.ceil(len(lst)/batch_n)
        
    for i in range(0, len(lst), n_elems):
        yield lst[i:(i + n_elems)]

