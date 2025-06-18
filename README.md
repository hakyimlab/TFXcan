## TFXcan
TFXcan (transcription factor binding prediction and correlation with complex traits).

## Date
Wednesday, June 18 2025

## Author
Temi

## Description
In addtion, we provide access to the following pipelines to generate the data used in the paper:
1. To develop DL-based predictors of TF binding: use this pipeline: https://github.com/hakyimlab/TFpred-snakemake
2. To run TFXcan, use this pipeline: https://github.com/hakyimlab/TFXcan-snakemake


## Brief description of scripts and folders
- [recreateEnpactResults](./notebooks/recreateEnpactResults.qmd): contains step-by-step instructions to run scripts that reproduce the main results in the TFXcan paper. THIS IS NOT AN ANALYSIS SCRIPT. It is written in Quarto and can be rendered to HTML or PDF, or run interactively in RStudio.
- [src](./src/): contains standalone bash, r or python scripts. Used in the [recreateEnpactResults](./notebooks/recreateEnpactResults.qmd) notebook.
- [software](./software/): contains software used such as liftover, e.t.c.
- [metadata](./metadata/): contains some files and data used to analyse the results.
- [recreateAnalysis](./notebooks/recreateAnalysis.qmd): contains scripts to reproduce the analysis, figures, supplementary figures, data and results in the TFXcan paper.
- [misc](./misc/): contains some random scripts and files [to be deleted later].
