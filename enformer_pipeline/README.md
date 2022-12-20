## Introduction
This folder contains scripts used to run ENFORMER to predict on a given a set of intervals/regions. 

## Usage
After editing the `enformer_parameters.json` file, simply call `python3 ./scripts/enformer_predict`

## Options
### General (and necessary)
- `project_dir`: (Absolute path) Provide a directory where files, logs, and outputs will be created.
- `model_path`: (Absolute path) Path to ENFORMER model.
- `hg38_fasta_file`: (Absolute path) Path to the hg38 fasta file.
- `interval_list_dir`: () Path to the folder where the intervals files are saved.
- `output_dir`: (Relative path) Path to where the outputs will be created.
- `dataset_type`: ("reference" or "personalized") A name for the type of dataset. 
- `TF`: (str) The transcription factor.
- `predict_on_n_regions`: (int: -1 or any number greater than 0): How many regions should be predicted on? If "-1" all regions present in the interval list are predicted on. 
- `prediction_data_name`: (str) A unique id for the predictions. This will be used to create folders mostly.
- `predictions_log_dir`: (path to folder): Where should predictions be logged? In the event of the job not completing and you intend to re-run, the file `{predictions_log_dir}/{id}_{TF}_predictions_log.csv` will be read and used such that if a prediction for a region exists and the `{region}_predicions.h5` file exists, predictions will not be made for that region.
- `batch_size`: 

### Personalized predictions
- `individuals`:(path, list or str) The unique ids of the individuals whose predictions are to be made. If providing a file, the ids should be written row-wise. A list of ids can be supplied and a single id can be supplied as a string. 
- `vcf_file`: (path) The path to the phased vcf file that contains the genotypes/variants of the individuals. 

### Optional
- `date`: (str in the form "YYYY-MM-DD") The date these predictions are made. If not supplied, today's date will be used.
- `use_parsl`: (bool: true or false) Should parsl be used to make predicting run faster? If `false` only one GPU is used. 
- `write_log`: (dict of bools) Controls what logs should be written. Any of `memory`, `cache`, `error` or `time`.
- `parsl_parameters`: (dict of bools) If `use_parsl` is true, these parameters are passed to parsl. Dictionary keys are `job_name`, `num_of_full_nodes`, `walltime`, `min_num_blocks`, and `max_num_blocks`, and possible values are `true` or `false`.

## To-do
