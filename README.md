## TFXcan


This repo contains scripts and codes in the analysis of the TFXcan paper.

TFXcan (transcription factor binding prediction and correlation with complex traits) method

The notebook can be found [here](./notebooks/recreateEnpactResults.qmd)

In addtion, we provide access to the following pipelines to generate the data used in the paper:
1. To develop DL-based predictors of TF binding: use this pipeline: https://github.com/hakyimlab/TFpred-snakemake
2. To run TFXcan, use this pipeline: https://github.com/hakyimlab/TFXcan-snakemake


## Brief description of scripts and folders
- [src](./src/): contains standalone bash, r or python scripts. Used in the [notebook](./notebooks/recreateEnpactResults.qmd) to reproduce the results.
- [software](./software/): contains software used such as liftover, e.t.c.
- [metadata](./metadata/): contains some files and data used to analyse the results.
