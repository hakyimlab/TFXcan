
# run enformer
conda activate /lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools
python3 /lus/grand/projects/TFXcan/imlab/shared/pipelines/enformer_pipeline/scripts/enformer_predict.py --parameters /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/metadata/cwas_imputed_config.json  

# aggregate
conda activate /lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools
/lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools/bin/python3 /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/scripts/tfpred/aggregate.py --metadata_file="/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/metadata/aggregation_config_cwas_imputed_AR_Prostate.json" --agg_types="aggByMeanCenter" --hpc="polaris"


# run enformer: association
conda activate /lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools
/lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools/bin/python3 /lus/grand/projects/TFXcan/imlab/shared/pipelines/enformer_pipeline/scripts/enformer_predict.py --parameters /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/association/config_files/enformer_parameters_geuvadis_AR.json >> sofar.out 2>&1

# aggregate
conda activate /lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools
/lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools/bin/python3 /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/scripts/tfpred/aggregate.py --metadata_file="/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/association/metadata/aggregation_config_geuvadis_association_AR_Prostate.json" --agg_types="aggByMean" --hpc="polaris"

# predict using model
qsub -v 'individuals_data_dir=/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/association/predictions_folder/geuvadis_association_AR_Prostate/predictions_2023-06-02/aggregated_predictions,predict_on=/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/association/metadata/geuvadis_association_AR_Prostate_2023-06-02.successful_predictions.csv,evaluate_rscript=/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/scripts/tfpred/enet_evaluate_individuals.R,output_dir=/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/association/output,agg_method=aggByMean,metainfo=AR_Prostate' enet_evaluate_individuals.pbs

# aggregate freedman foxa1
conda activate /lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools
/lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools/bin/python3 /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/scripts/tfpred/aggregate.py --metadata_file="/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/TFPred_pipeline/config_files/aggregation_config_freedman_FOXA1.json" --agg_types="aggByMeanCenter" --hpc="polaris"


# run enformer
conda activate /lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools
python3 /lus/grand/projects/TFXcan/imlab/shared/pipelines/enformer_pipeline_aggregate/scripts/enformer_predict.py --parameters /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/compare_predictors/config_files/enformer_parameters_1KG_AR.json


conda activate /lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools
python3 /lus/grand/projects/TFXcan/imlab/shared/pipelines/enformer_pipeline_aggregate/scripts/enformer_predict.py --parameters /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/compare_predictors/config_files/cwas_imputed_config.json


# aggregate freedman foxa1
conda activate /lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools
/lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools/bin/python3 /lus/grand/projects/TFXcan/imlab/shared/pipelines/enformer_pipeline_aggregate/scripts/aggregate/aggregate.py --metadata_file="/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/compare_predictors/metadata/aggregation_config_1KG_AR_Prostate.json" --agg_types="aggByCollect" --hpc="polaris"


# predict using model
qsub -v 'individuals_data_dir=/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/compare_predictors/predictions_folder/1KG_AR_Prostate/predictions_2023-06-28/aggregated_predictions,predict_on=/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/compare_predictors/metadata/1000_genome_individuals.txt,evaluate_rscript=/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/scripts/tfpred/enet_evaluate_individuals.R,output_dir=/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/compare_predictors/tfpred_scores,agg_method=aggByCollect,metainfo=AR_Prostate' /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/scripts/tfpred/enet_evaluate_individuals.pbs