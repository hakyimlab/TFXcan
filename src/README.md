# src/

Standalone R, Python and sbatch scripts used by the notebooks in [../notebooks](../notebooks/) (mainly
[reproduceData.qmd](../notebooks/reproduceData.qmd) and [recreateAnalysis.qmd](../notebooks/recreateAnalysis.qmd))
to reproduce the Enpact / TFXcan paper's data and results. Everything here is called with `sbatch <script>.sbatch <args...>`
on the RCC cluster; the `.sbatch` files are thin wrappers that `conda activate` the right environment and call the
matching `.R` / `.py` script with positional args (`${1}`, `${2}`, ...).

Every script now has a header comment (Author/Date/Description/Usage) and its `--flag` options have real
`help=` text. If you're unsure what a script does, read its header first, then check whether it's mentioned below.

Convention: an `.sbatch` file with the same base name as an `.R`/`.py` file is that script's cluster wrapper
(e.g. `train_enpact.R` + `train_enpact.sbatch`). Only wrappers with no same-named script, or scripts with more
than one wrapper, are called out explicitly below.

## Pipeline map

### 1. Train the DL-based (Enformer/Enpact) TF-binding models
The actual deep-learning training happens outside this repo, in
[TFPred-snakemake](https://github.com/hakyimlab/TFPred-snakemake). The scripts here train/evaluate the
lightweight elastic-net head on top of pre-computed Enformer features:

- `train_enpact.R` / `train_enpact.sbatch` — train the elastic-net model from a training matrix.
- `evaluateEnpactModel.R` — evaluate a trained model against a held-out test set.
- `trainAndEvaluateEnpactModel.sbatch` — runs `train_enpact.R` then `evaluateEnpactModel.R` back to back (script
  paths are passed in as args, so it can also run other train/eval script pairs).
- `create_training_data.R`, `merge_epigenome_ground_truth.R`, `train_enpact_different_nobs.sbatch` — ablation
  experiment: train models on increasing numbers of positive binding sites (50 to 18,000) to see how model
  quality scales with training-set size.
- `test_enpact_different_features.R` / `.sbatch`, `train_enpact_different_features.sbatch` — feature-ablation
  experiments (train/test after dropping specific feature columns). Note:
  `train_enpact_different_features.sbatch` currently calls `train_enpact.R` under the hood, not a script named
  `train_enpact_different_features.R` (that version lives in [archive/](archive/), superseded).

### 2. Extract & aggregate epigenomic features (Enformer predictions)
- `extract_epigenome.py` — helper module (bin/coordinate math); imported by `aggregate_epigenomes.py`, not run directly.
- `aggregate_epigenomes.py` / `aggregate_epigenomes.sbatch`, `gather_epigenomes.sbatch` — aggregate per-locus
  Enformer predictions from the reference epigenome HDF5s into a feature matrix/CSV per locus list. There are a
  few near-duplicate copies of this script across sibling repos (TFPred-snakemake, Enpact); the wrapper you use
  determines which copy actually runs — check the `python3 <path>` line in the `.sbatch` file.
- `create_h5database_from_csvs.py` / `.sbatch` — bundle many per-individual aggregated CSVs into a single HDF5 database.
- `collect_predictions_into_rds.R` / `.sbatch` — collect per-individual prediction files into one `.rds`.
- `enformer_predict.sbatch`, `enformer_merge.sbatch`, `enformer_process.sbatch` — run/merge/post-process raw
  Enformer predictions for custom regions (e.g. the gene-modules ACR experiments), as opposed to the
  reference-epigenome pipeline above.

### 3. Compute Enpact scores (apply trained weights to features)
- `collect_enpact_weights.R` / `.sbatch` — gather per-TF-tissue model weights into one compiled weights file.
- `calculateEnpactScores.R` / `.sbatch` — compute Enpact scores for a set of individuals from aggregated CSV features.
- `enpact_predict_hdf5.R` / `.sbatch` — same idea, reading features from an HDF5 database instead of CSVs.
- `predict_enpact_scores.sbatch`, `enpact_predict.sbatch` — other variants of "apply weights matrix to a feature
  matrix" used at different pipeline stages (see each script's Usage line for its exact flags).
- `cross_predict_enpact_scores.R` / `.sbatch` — apply one model's weights to another model's feature set (cross-prediction).

### 4. CWAS (Baca) weights and predictions
- `predictCWASscores.sbatch` — predict CWAS scores in Baca individuals from the CWAS-paper SQLite weight
  databases, via MetaXcan's `Predict.py`.
- `predictENPACTscores.sbatch` — same, for the linearized SNP-based Enpact ("lEnpact") models.

### 5. Train a SNP-based ("linearized") Enpact model + PredictDB
- `formatForPredictDB.R` / `.sbatch` — reformat Enpact scores into PredictDB-compatible weights/annotation files.
- `generateLinearModels.sbatch` — run the [PredictDb-nextflow](https://github.com/hakyimlab/PredictDb-nextflow)
  pipeline to fit the SNP-based elastic-net models (no matching `.R` script; it drives a Nextflow workflow).

### 6. TFXcan / S-PrediXcan association testing
- `sTFXcan.sbatch` — run S-PrediXcan ("TFXcan") for one model DB against GWAS summary stats.
- `run_summary_tfxcan.sbatch` — run S-PrediXcan across many TF-tissue models against one phenotype's GWAS.
- `gather_sTFXcan.sbatch` — collect/merge per-model S-PrediXcan outputs into a single results file.

### 7. Consensus matrix factorization ("tenerife")
Decomposes the loci-by-TF/tissue z-ratio matrix into latent factors, repeated across random subsets of loci for
stability, then clusters the repeated runs into consensus programs.

- `repeat_flash.R` / `repeat_flash.sbatch` — run Flashier factorization on one batch of random loci subsets (the
  active, parallelized worker — see the sbatch's `parallel` loop).
- `gather_repeats.R` — merge the per-batch Flashier outputs from `repeat_flash.R`.
- `cluster_programs.R` / `cluster_programs.sbatch` — consensus-cluster the merged factor loadings into "programs".
- `factorize_tfxcan.sbatch` — the end-to-end wrapper: `prepare_summary_matrices.R` → `repeat_flash.R` (parallel)
  → `gather_repeats.R` → `cluster_programs.R`. This is the generalized/current version of the manual steps
  walked through in `reproduceData.qmd`'s "Consensus matrix factorization" section.
- `prepare_summary_matrices.R` — build the z-ratio matrix and random-subset splits from GWAS + TFXcan summary stats.
- `runFlashier.R`, `consensus_flashier.R` — earlier single-run (non-batched) Flashier variants. Superseded by
  `repeat_flash.R` — see the commented-out block at the bottom of `repeat_flash.sbatch`. Kept for reference.
- `factorize_associations.R` — same option list as `prepare_summary_matrices.R` and not called from any
  `.sbatch` here; looks like an earlier draft that got forked into `prepare_summary_matrices.R`. Not part
  of the active pipeline — kept around as reference.

### 8. Liftover (hg19 → hg38)
- `liftover.sbatch`, `liftoverSingleBed.sbatch`, `liftoverMultipleBeds.sbatch` — wrappers around UCSC `liftOver`
  for converting BED coordinates between genome builds.

### 9. Motifs, accessibility, and binding/expression correlation (gene-module / ACR experiments)
- `calculateDistanceToMotifs.R` / `.sbatch` — distance from each locus to the nearest TF motif.
- `find_motifs_given_bed.sbatch`, `annotate_bed_with_motifs.sbatch` — motif scanning/annotation of a BED file.
- `accessibility_effect_on_binding.R` / `.sbatch`, `accessibility_effect_on_binding.ttest.R` — test whether
  chromatin accessibility affects predicted binding.
- `correlate_binding_expression.R` / `.sbatch`, `correlate_binding_expression.acrs.R` / `.sbatch` — correlate
  predicted TF binding with gene expression (the `.acrs` variant works on accessible chromatin regions).
- `find_snps_in_windows.R` / `.sbatch` — find SNPs falling inside a set of genomic windows.
- `sample_random_sites.R` — sample random genomic sites, e.g. as a null/background set for the above.

### 10. Heritability / GRM analysis
- `make_grm.sbatch`, `make_grm_single.sbatch`, `model_grms.sbatch` — build genetic relationship matrices (GRMs).
- `estimate_hsq.R` / `.sbatch`, `estimate_hsq_one_grm.R` / `.sbatch` — estimate SNP heritability from GRM(s).
- `estimate_null_ppi.R` / `.sbatch` — estimate a null distribution for protein-protein-interaction-based scores.

### 11. Model evaluation / stats
- `calculate_auprc.R` / `.sbatch` — AUPRC of Enpact model predictions vs. ground truth.
- `calculate_wilcoxon_test.R` / `.sbatch` — Wilcoxon test between two groups of scores.
- `count_training_test_overlaps.R` / `.sbatch` — check for train/test set leakage/overlap.

### 12. Promoter GRNs
- `promoter_grns.tree.R` / `.sbatch` — random-forest-based gene regulatory network inference from binding and
  expression matrices.

### Shared / misc
- `modules.R` — shared R helper functions (plot themes, small transforms) sourced by other scripts. Not runnable on its own.
- `test.sbatch` — placeholder/scratch sbatch script (just echoes "Running test"); not part of the pipeline.
- `infer.sbatch` — generic wrapper that runs an arbitrary epigenome-extraction Python script (path passed as
  arg 1); its header comment was copied from `trainAndEvaluateEnpactModel.sbatch` and described the wrong
  thing — fixed in the header now.
- `archive/` — older, superseded versions of scripts kept for reference (`trainEnpactModel.R`,
  `train_enpact_different_features.R`), replaced by `train_enpact.R`.

## Things worth knowing before you dig further
- Several `.sbatch` files hardcode absolute paths to sibling repos (`TFPred-snakemake`, `Enpact`,
  `TFXcan-snakemake`) rather than calling scripts in this `src/`. Check the actual `python3`/`Rscript` line in
  a wrapper before assuming it runs the local copy of a same-named script.
- Some notebook references are stale (e.g. `reproduceData.qmd` calls `src/trainEnpactModel.R`, which has since
  moved to `src/archive/`). When a script and its notebook usage disagree, trust the script + its sbatch wrapper.
