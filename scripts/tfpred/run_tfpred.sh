
# run enformer
conda activate /lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools
python3 /lus/grand/projects/TFXcan/imlab/shared/pipelines/enformer_pipeline/scripts/enformer_predict.py --parameters /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/metadata/cwas_imputed_config.json  

# aggregate
conda activate /lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools
/lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools/bin/python3 /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/scripts/tfpred/aggregate.py --metadata_file="/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/metadata/aggregation_config_cwas_imputed_AR_Prostate.json" --agg_types="aggByMeanCenter" --hpc="polaris"