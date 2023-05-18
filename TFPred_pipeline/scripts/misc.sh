

# aggregate
/lus/grand/projects/TFXcan/imlab/shared/software/conda_envs/enformer-predict-tools/bin/python3 /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/TFPred_pipeline/scripts/utilities/aggregate.py --metadata_file /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/metadata/aggregation_config_cwas_AR_Prostate.json --agg_types "aggByMeanCenter" --hpc="polaris"


# enet evaluate on individuals
qsub /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/scripts/tfpred/enet_evaluate_individuals_pbs.sh

# # qsub -v 'data_file=/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/data/train-test-val/kawakami/data_2022-12-12/kawakami_aggByCenter_FOXA1_old.csv.gz,metainfo=old' train_enet_model_pbs.sh