This folder includes scripts used to make personalized, as well as, reference predictions. 

The following folders/files are needed:

[ENFORMER parameters](./metadata/enformer_parameters.json) - This is a json file and is the only file that needs changing to suit your specific circumstance. All paths given to the json file should be absolute paths, please.

    - The text file containing the intervals should be named this way: `{individual_id}_{transcription_factor}_{any other information you want}.txt`, and should contain just the regions of interest listed row-wise.

    - The `individuals.txt` file should be a similar file as above, but should contain the individuals you want to make predictions for. Alternatively, you can pass a list of individuals directly to the json file. 
    
    - Parsl is used to distribute jobs across multiple GPUs. If not needed, you can disable parsl by toggling the `use_parsl` parameter (either "true" or "false") in the json file. 
