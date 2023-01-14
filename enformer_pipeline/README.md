
## Introduction
This folder contains scripts used to run ENFORMER to predict on a given a set of intervals/regions. There are to kinds of datasets to predict on: `reference` or `personalized`. 

## Usage
After editing the [`enformer_parameters.json`](./metadata/enformer_parameters.json) file , simply call [`python3 ./scripts/enformer_predict`](./scripts/enformer_predict.py)

## Options
### General (and necessary)
- `project_dir`: (Absolute path) Provide a directory where files, logs, and outputs will be created. **Currently not being used.** 
- `model_path`: (Absolute path) Path to ENFORMER model.
- `hg38_fasta_file`: (Absolute path) Path to the hg38 fasta file.
- `interval_list_dir`: (Relative path to the this directory) Path to the folder where the intervals files are saved.
- `output_dir`: (Relative path) Path to where the outputs will be created. The folder, as well as other necessary folder levels, will be created.
- `dataset_type`: (str: "reference" or "personalized") A name for the type of dataset. This is necessary so that the pipeline knows whether to use vcf files or not. 
- `TF`: (str) The transcription factor. Necessary for naming conventions.
- `predict_on_n_regions`: (int: -1 or any number greater than 0): How many regions should be predicted on? If `-1` all regions present in the interval list are predicted on. 
- `prediction_data_name`: (str) A unique id for the predictions that will be used to create folders. Can also use the name of the dataset e.g. "freedman", or "kawakami". If `dataset_type` is "reference", this will be the name of the folder within which predictions will be made. If `dataset_type` is personalized, the unique ids of the individuals will be used.
- `predictions_log_dir`: (path to folder): Where should predictions be logged? In the event of the job not completing and you intend to re-run, the file `{predictions_log_dir}/{id}_{TF}_predictions_log.csv` will be read and used such that if a prediction for a region exists and the `{region}_predictions.h5` file exists, predictions will not be made for that region anylonger.
- `batch_size`: (int) Predictions for intervals will be split into `batch_size`.  
- `log_data`: (relative path) All kinds of logs are placed here. 

### Personalized predictions
- `individuals`:(path, list or str) The unique ids of the individuals whose predictions are to be made. If providing a file, the ids should be written row-wise. A list of ids can be supplied and a single id can be supplied as a string. 
- `vcf_file`: (path) The path to the phased vcf file that contains the genotypes/variants of the individuals. 

### Optional
- `date`: (str in the form "YYYY-MM-DD") The date these predictions are made. If not supplied, today's date will be used.
- `use_parsl`: (bool: true or false) Should parsl be used to make predicting run faster? If `false` only one GPU is used. 
- `write_log`: (dict of bools) Controls what logs should be written. Any of `memory`, `cache`, `error` or `time`, and possible values are `true` or `false`.
- `parsl_parameters`: (dict of bools) If `use_parsl` is true, these parameters are passed to parsl. Dictionary keys are `job_name`, `num_of_full_nodes`, `walltime`, `min_num_blocks`, `queue`, and `max_num_blocks`
    - `job_name`: (str) e.g. "my_predictions"
    - `num_of_full_nodes`: (int) e.g. 10
    - `walltime`: (str) "HH:MM:SS" e.g. "01:10:59"
    - `min_num_blocks`: (int) e.g. 0
    - `max_num_blocks`: (int) e.g. 1
    - `queue`: (str) e.g. "preemptable", "full-node"
## To-do
- [ ] Check that `os.path.join` paths work correctly.
- [ ] User should supply the intervals list within the predictions folder.
- [ ] Provide checks to ensure that necessary files and folders are available.
- [ ] Change `batch_size` to `n_batches`. 

