## TFXcan


This repo contains scripts and codes in the analysis of the TFXcan paper.

TFXcan (transcription factor binding prediction and correlation with complex traits) method



In addtion, we provide access to the following pipelines to generate the data used in the paper:
1. To develop DL-based predictors of TF binding: use this pipeline: https://github.com/hakyimlab/TFpred-snakemake
2. To run TFXcan, use this pipeline: https://github.com/hakyimlab/TFXcan-snakemake


## Brief description of scripts and folders
- [notebook](./notebooks/recreateEnpactResults.qmd): contains step by step instructions to reproduce the main results in the TFXcan paper. It is written in Quarto and can be rendered to HTML or PDF, or run interactively in RStudio.
- [src](./src/): contains standalone bash, r or python scripts. Used in the [notebook](./notebooks/recreateEnpactResults.qmd).
- [software](./software/): contains software used such as liftover, e.t.c.
- [metadata](./metadata/): contains some files and data used to analyse the results.
