


def check_queries(sample, queries, output_dir, prediction_logfiles_folder, sequence_source):
    """
    Check whether a given region, for an individual has been predicted and logged.

    Parameters:
        sample: str 
            The name/id of an individual
        queries: A batch
            A region in the genome in the form `chr_start_end`.
        output_dir: str (path)
            The folder where the predictions should have been logged. 
        logfile: pd.DataFrame 
            A dataframe of a log file or `None` if the log file does not exist. 
    
    Returns: dict
        'query': the query region if it has not been logged or predictions don't exist
        'logtype': whether it should be logged if it has not been logged

    If predictions exist and the query has been logged, this function returns None.
    """
    import pandas as pd
    import os
    import numpy as np

    if isinstance(prediction_logfiles_folder, type(None)):
        output = [{'query':query, 'logtype':'y'} for query in queries]
        return(output)

    else:
        id_logfile = os.path.join(prediction_logfiles_folder, f'{sample}_log.csv')
        if os.path.isfile(id_logfile):
            try:
                id_logfile = pd.read_csv(id_logfile) 
                id_logfile = id_logfile.loc[id_logfile['sample'] == sample, : ]
                # check if the file is saved
                if sequence_source == 'personalized': # prediction must be present in two folders
                    queries_saved_h1 = [str(f'{output_dir}/{sample}/haplotype1/{query}_predictions.h5') for query in queries]
                    queries_saved_h2 = [str(f'{output_dir}/{sample}/haplotype2/{query}_predictions.h5') for query in queries]
                    queries_saved = np.array([os.path.isfile(q) for q in queries_saved_h1]) * np.array([os.path.isfile(q) for q in queries_saved_h2])
                elif sequence_source == 'reference':
                    queries_saved = [str(f'{output_dir}/{sample}/haplotype0/{query}_predictions.h5') for query in queries]
                    queries_saved = np.array([os.path.isfile(q) for q in queries_saved])
                # check if the file is in the logfile
                query_logged = np.array([(query in id_logfile.motif.values) for query in queries])
                # 
                if not query_logged.shape[0] == queries_saved.shape[0]:
                    raise Exception("Lengths of queries logged and saved conditions are not the same")
                queries_condition = queries_saved * query_logged # true should not be predicted or logged; False should be
                queries_condition = queries_condition.tolist()
                # 
                output = [{'query':queries[i], 'logtype':'y'} for i, qc in enumerate(queries_condition) if qc is False]
            except pd.errors.EmptyDataError:
                id_logfile = None
                output = [{'query':query, 'logtype':'y'} for query in queries]
        else:
            output = [{'query':query, 'logtype':'y'} for query in queries]
              
        return(output)