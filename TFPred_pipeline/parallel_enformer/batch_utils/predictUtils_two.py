# Usage: This module is used to predict with ENFORMER on genomic sequences
# Author: Temi
# Date: Thurs Feb 2 2023

import json
import os
import tensorflow as tf

global module_path, write_log, sequence_source, grow_memory

# read in the config_file
whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)
#module_path = os.path.abspath(whereis_script)

#if __name__ == '__main__':
if __name__ == 'predictUtils_two':

    with open(f'{whereis_script}/config.json') as f:
        parameters = json.load(f)
        parameters_file = parameters['params_path']
    
    with open(f'{parameters_file}') as f:
        parameters = json.load(f)

        if parameters['sub_dir'] == True:
            project_dir = os.path.join(parameters['project_dir'], 'prediction_folder')
        else:
            project_dir = parameters['project_dir']

        sequence_source = parameters['sequence_source']
        write_logdir = parameters['write_log']['logdir']
        write_log = parameters["write_log"]
        fasta_file = parameters['hg38_fasta_file']
        model_path = parameters['model_path']
        #vcf_file = parameters['vcf_file']

# for ltype in ['memory', 'error', 'time', 'cache']:
#     if write_log['logtypes'][ltype]: os.makedirs

if any([write_log['logtypes'][ltype] for ltype in ['memory', 'error', 'time', 'cache']]):
    write_log['logdir'] = os.path.join(project_dir, write_logdir)


print(write_log)
print( write_log['logdir'])

#import checksUtils
import loggerUtils
import sequencesUtils
import predictionUtils

global enformer_model
global fasta_extractor

#enformer_model = predictionUtils.get_model(model_path)
grow_memory = True
print(f'GPU Memory before calling batch predict function is {loggerUtils.get_gpu_memory()}')

def enformer_predict_on_batch(batch_regions, samples, path_to_vcf, batch_num, output_dir, prediction_logfiles_folder, sequence_source):

    # this could mean
    # - check_queries returned nothing because
        # - nothing should be returned - good ; and it should return none
        # - something is wrong with check_queries ; and I should fix that
    if (not batch_regions) or (batch_regions is None):
        raise Exception(f'[INFO] There are no regions in this batch {batch_num}.')

    # print(f'batch_regions are: {batch_regions}')
    # print(f'samples are: {samples}')
    # print(f'path_to_vcf are: {path_to_vcf}')
    # print(f'output_dir are: {output_dir}')
    # print(f'prediction_logfiles_folder are: {prediction_logfiles_folder}')
    # print(f'sequence_source are: {sequence_source}')

    print(f'GPU Memory at start of batch {batch_num} predict function is {loggerUtils.get_gpu_memory()}')

    if grow_memory == True:
        gpus = tf.config.experimental.list_physical_devices('GPU')
        if gpus:
            try:
                # Currently, memory growth needs to be the same across GPUs
                for gpu in gpus:
                    tf.config.experimental.set_memory_growth(gpu, True)
            except RuntimeError as e:
                print(f'[RUNTIME ERROR] Batch number: {batch_num} of {type(e)} in module {__name__}')

    try:
        #model = predictionUtils.get_model(model_path)
        #model = enformer_model # check global definitions
        enformer_model = predictionUtils.get_model(model_path)
        fasta_extractor = sequencesUtils.get_fastaExtractor(fasta_file)
        #print('Fasta and model successfully loaded')
        logger_output = []
        for input_region in batch_regions:
            #print(f'Creating sequences for {input_region}')
            samples_enformer_inputs = sequencesUtils.create_input_for_enformer(region_details=input_region, samples=samples, path_to_vcf=path_to_vcf, fasta_func=fasta_extractor, hap_type = 'both', resize_for_enformer=True, resize_length=None, write_log=write_log, sequence_source=sequence_source)

            print(f'Region {input_region} sequences successfully created')

            # check that all the samples are accounted for
            #print(sorted(list(samples_enformer_inputs['sequence'].keys())))
            if samples_enformer_inputs is not None:
                if samples_enformer_inputs['metadata']['sequence_source'] == 'var':
                    if sorted(samples) != sorted(list(samples_enformer_inputs['sequence'].keys())):
                        raise Exception(f'[ERROR] Some samples cannot be found.')
            elif samples_enformer_inputs is None:
                continue

            #logging_info_list = [] # collect all logging information here
            for sample in samples:
                if samples_enformer_inputs['metadata']['sequence_source'] == 'var':
                    sample_predictions = predictionUtils.enformer_predict_on_sequence(model=enformer_model, sample_input=samples_enformer_inputs['sequence'][sample])
                elif samples_enformer_inputs['metadata']['sequence_source'] in ['ref', 'random']:
                    sample_predictions = predictionUtils.enformer_predict_on_sequence(model=enformer_model, sample_input=samples_enformer_inputs['sequence'])
                
                # check that the predictions have the appropriate shapes
                for hap in sample_predictions.keys():
                    print(f'Sample {sample} {hap} predictions of shape {sample_predictions[hap].shape} successfully created and saved.')
                    if sample_predictions[hap].shape != (1, 17, 5313):
                        raise Exception(f'[ERROR] {sample}\'s {hap} predictions shape is {sample_predictions[hap].shape} and is not the right shape.')
                    
                # otherwise, you can save the predictions ; prediction will be reshaped to (17, 5313) here
                sample_logging_info = loggerUtils.save_haplotypes_h5_prediction(haplotype_predictions=sample_predictions, metadata=samples_enformer_inputs['metadata'], output_dir=output_dir, sample=sample)

                if (sample_logging_info is not None) and (len(sample_logging_info) == 4):
                    predictions_log_file = os.path.join(prediction_logfiles_folder, f'{sample}_log.csv')
                    logger_output.append(loggerUtils.log_predictions(predictions_log_file=predictions_log_file, what_to_write=sample_logging_info))
                else:
                    logger_output.append(1)

        if write_log['logtypes']['memory']:
            if tf.config.list_physical_devices('GPU'):
                mem_use = loggerUtils.get_gpu_memory()
                msg_mem_log = f"[INFO] GPU memory at the end of prediction batch {batch_num}): free {mem_use[0]} mb, used {mem_use[1]} mb on {loggerUtils.get_gpu_name()}"
                MEMORY_LOG_FILE = os.path.join(write_log['logdir'], "memory_usage.log")
                loggerUtils.write_logger(log_msg_type = 'memory', logfile = MEMORY_LOG_FILE, message = msg_mem_log)

        if write_log['logtypes']['cache']:
            msg_cac_log = f'[CACHE] (model) at batch {batch_num}: [{predictionUtils.get_model.cache_info()}, {loggerUtils.get_gpu_name()}]'
            CACHE_LOG_FILE = os.path.join(write_log['logdir'], 'cache_usage.log')
            loggerUtils.write_logger(log_msg_type = 'cache', logfile = CACHE_LOG_FILE, message = msg_cac_log)

        return(logger_output)
    
    except (TypeError, AttributeError) as tfe:
        if write_log['logtypes']['error']:
            if tf.config.list_physical_devices('GPU'):
                mem_use = loggerUtils.get_gpu_memory()
                err_mem_log = f"[ERROR] GPU memory error of type {type(tfe).__name__} for batch {batch_num}): free {mem_use[0]} mb, used {mem_use[1]} mb on {loggerUtils.get_gpu_name()}"
                MEMORY_ERROR_FILE = os.path.join(write_log['logdir'], 'error_details.log')
                loggerUtils.write_logger(log_msg_type = 'error', logfile = MEMORY_ERROR_FILE, message = err_mem_log)
        else:
            raise Exception(f'[ERROR] of {type(tfe).__name__} at batch {batch_num}')



