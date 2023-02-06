


def check_inputs(params_path):
    import os
    import json
    import pandas as pd

    gatherer = {}

    with open(f'{params_path}') as f:
        parameters = json.load(f)

        project_dir = parameters['project_dir']
        if not os.path.exists(project_dir):
            gatherer['project_dir'] = False
        else:
            print(f'[INFO] Project directory ... good')

        interval_list_file = parameters['interval_list_file']
        if not os.path.isfile(interval_list_file):
            raise Exception(f'[ERROR] Interval list file cannot be found at {os.path.dirname(interval_list_file)}')
        else:
            tf = pd.read_table(interval_list_file, nrows=20, sep=' ', header=None)
            
            print(f'[INFO] Interval list file ... good')
        
        predictions_log_dir = os.path.join(project_dir, parameters['predictions_log_dir'])
        if not os.path.exists(predictions_log_dir):
            print(f'[INFO] Log folder for completed predictions will be created ... good')

        log_dir = os.path.join(project_dir, parameters['log_dir'])
        if not os.path.exists(predictions_log_dir):
            print(f'[INFO] Log folder for submitted job info will be created ... good')

        if int(parameters['batch_size']) is None:
            raise Exception(f'[ERROR] batch size for predictions is not provided. Please read the documentation.')
        
        if int(parameters['n_individuals']) is None:
            raise Exception(f'[ERROR] Number of individuals to predict on is not provided. Please read the documentation.')

        if int(parameters['predict_on_n_regions']) is None:
            raise Exception(f'[ERROR] Number of regions to predict on is not provided. Please read the documentation.')

        if parameters['sequence_source'] is None:
            raise Exception(f'[ERROR] Sequence source is not provided. Please read the documentation.')
        elif parameters['sequence_source'] == 'personalized':
            if not os.path.exists(parameters['hg38_fasta_file']):
                raise Exception(f'[ERROR] genome fasta file is not provided. Please read the documentation.')
            if parameters['vcf_files'] is None:
                raise Exception(f'[ERROR] vcf files paths are not provided. Please read the documentation.')
            if not os.path.exists(parameters['vcf_files']['folder']):
                raise Exception(f'[ERROR] vcf files folder does not exist. Please read the documentation.')
            if not parameters['vcf_files']['files']:
                raise Exception(f'[ERROR] vcf files have not been provided. Please read the documentation.')
        elif parameters['sequence_source'] == 'reference':
            if not os.path.exists(parameters['hg38_fasta_file']):
                raise Exception(f'[ERROR] genome fasta file is not provided. Please read the documentation.')
                
        if parameters['prediction_data_name'] is None:
            raise Exception(f'[ERROR] Prediction data name is not provided. Please read the documentation.')

    return(0)