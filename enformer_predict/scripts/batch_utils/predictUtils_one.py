# import parsl
# from parsl.app.app import python_app

from parsl.app.app import python_app

@python_app
def run_batch_predictions(batch_regions, batch_num, individual, vcf_func, script_path, output_dir, logfile, predictions_log_dir): #
  
    import sys, os, tqdm, faulthandler
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
        reg_prediction = predictUtils_two.enformer_predict(batch_region=filtered_check_result, sample=individual, model=None, output_dir=output_dir, vcf_func=vcf_func, predictions_log_dir=predictions_log_dir, batch_num=batch_num)
        return(reg_prediction) # returns 0 returned by enformer_predict

# else:
#     print(f"[WARNING] {check_result}: Either length of input sequence is invalid (NoneType) or too long or too short")
            

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

