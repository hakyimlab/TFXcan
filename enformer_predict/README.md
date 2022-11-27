This folder includes scripts used to make personalized, as well as, reference predictions using ENFORMER.

### Usage 
After updating the metadata json file, simply call: python3 [./scripts/enformer-predict-personalized.py](./scripts/enformer-predict-personalized.py)

### Paramaters
The following folders/files/changes are needed:

[ENFORMER parameters](./metadata/enformer_parameters.json) - This is a json file and is the only file that needs changing to suit your specific circumstance. All paths given to the json file should be absolute paths, please. Also, please don't add a forward slash, `/`, after the path names. Thank you!

- `output_dir`: Predictions are saved here as `{individual_id}/{region}_predictions.h5`, using an `hdf5` format. 
- `interval_list_dir`: This folder should contain a `.txt` file for each individual you want to predict for. The text file containing the intervals should be named this way: `{individual_id}_{transcription_factor}_{any other information you want}.txt`, and should contain just the regions of interest listed row-wise. Also each region should be named thus: `chr_start_end` e.g. "chr1_4000_4020". This region will be expanded or resized for ENFORMER. 
- `individiduals`: The `individuals.txt` file should be a similar file as above, but should contain the individuals' ids you want to make predictions for. Alternatively, you can pass a list of individuals directly to the json file. 
- Parsl is used to distribute jobs across multiple GPUs. If not needed, you can disable parsl by toggling the `use_parsl` parameter (either "true" or "false") in the json file. 
- `batch_size`: is used so that given batches are passed to the GPUs. Change this as needed. 
- `predictions_log_dir`: When predictions are made, they should be logged per individual. You will find the logs here. These prediction log files are important if there is an interruption during the running of the script. If the file is available, the script will not re-predict on regions whose `*.h5` files are available and whose predictions have been logged, to save time. If you need to force prediction, delete either the log file, or the `*.h5` prediction for that region.
- `log_dir`: Files detailing error messages, memory consumption, and cache usage are deposited here. 
- Besides predicting and saving in the user-defined folder, some log files are also returned and are found in the `cobalt-log` folder. 

### To-do
- I may need to find a way to toggle Parsl's `@python_app` decorator on/off, depending on if the `use_parsl` parameter is "true" or "false".
- I also want to find a way to toggle if logging should be done and at what level.
