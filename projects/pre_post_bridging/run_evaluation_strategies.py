import pickle
import pandas as pd
import json
from src.evaluation_functions_copy import (
    post_normalization_bridging,
    load_sample_metadata,
    get_matching_measureids,
    build_matched_measure_table,
    build_matched_score_table,
    evaluate_strategy,
    evaluate_strategy_distribution,
    evaluate_matched_samples,
    MATRIX_EDTA,
    MATRIX_DS1,
    MATRIX_DS2,
    MATRIX_NDS,
)
from src.data_preparation import prepare_adat_datasets
from src.prophet_scores import get_prophet_score

MADMED_ADATS_PATH = "projects/pre_post_bridging/results/adats/madmed_for_MDACC_comparison"
MADMED_ADATS_PATH = "projects/pre_post_bridging/results/adats/madmed_for_MDACC_comparison"
MED_ADATS_PATH = "projects/pre_post_bridging/results/adats/median_for_MDACC_comparison"
POST_NORMALIZATION_DS_PARAMS_PATH = "data/params/post_normalization/ds_med_08032026.json"
POST_NORMALIZATION_nDS_PARAMS_PATH = "data/params/post_normalization/nds_med_08032026.json"
ALL_MATCHED_SAMPLES_PATH = "data/excels/all_matched_samples_streck.xlsx"
MDACC_PATH = "data/excels/matched MDACC proteomics with scores.xlsx"

"""
Run evaluation of Streck bridging strategies against EDTA reference scores.

This script:
1. Loads sample metadata
2. Loads scored datasets for EDTA and Streck strategies
3. Defines the strategy registry
4. Runs evaluation for each strategy
5. Saves plots and returns summary metrics
"""


# --------------------------------------------------
# Main evaluation entrypoint
# --------------------------------------------------

def main():

    # --------------------------------------------------
    # 1. Load metadata
    # --------------------------------------------------
    print(1)
    info_path = "data/excels/20260309_adat_data_samples_with_sample_info.csv"

    sample_registry = load_sample_metadata(info_path)

    info_df = sample_registry["info"]
    prophet_edta_info_df = sample_registry["PROphet_EDTA"]
    all_edta_info_df = sample_registry['all_EDTA']
    ds_info_df = sample_registry["DS"]
    nds_info_df = sample_registry["nDS"]
    # ds_rfu_df, nds_rfu_df, edta_rfu_df = prepare_base_rfu_datasets("data/adat/base/v_166_OH2026_004.parquet", sample_registry)
    
    # --------------------------------------------------
    # 2. Load Data
    # --------------------------------------------------
    print(2)
    post_anml_prophet_edta = pd.read_csv('data/excels/20260309_adat_data_measurements.csv')
    post_anml_prophet_edta = post_anml_prophet_edta[post_anml_prophet_edta['MeasureId'].isin(prophet_edta_info_df["MeasureId"])]
    madmed_ds_df, madmed_nds_df = prepare_adat_datasets(MADMED_ADATS_PATH, sample_registry['DS'], sample_registry['nDS'], 'anmlSMP.adat') 
    med_ds_df, med_nds_df = prepare_adat_datasets(MED_ADATS_PATH, sample_registry['DS'], sample_registry['nDS'], 'anmlSMP.adat') 
    matched_samples = pd.read_excel(ALL_MATCHED_SAMPLES_PATH)
    matched_samples_ds, matched_samples_nds = get_matching_measureids(matched_samples)
    
    with open(POST_NORMALIZATION_DS_PARAMS_PATH, 'r') as f:
        post_norm_ds_params = json.load(f)
    with open(POST_NORMALIZATION_nDS_PARAMS_PATH, 'r') as f:
        post_norm_nds_params = json.load(f)

    # --------------------------------------------------
    # 3. Score Data
    # --------------------------------------------------
    print(3)
    # Single bridging step strategies
    
    model = pd.read_pickle("data/models/DCB_nosex_gamma1_rescaled.pkl")
    post_anml_prophet_edta["PROphetScore"] = get_prophet_score(post_anml_prophet_edta, model)
    post_anml_prophet_edta.to_csv("data/scores/edta_scores.csv", index=False)     
    madmed_ds_df["PROphetScore"] = get_prophet_score(madmed_ds_df, model)   
    madmed_ds_df.to_csv("data/scores/ds_madmed_pre.csv", index=False)  
    madmed_nds_df["PROphetScore"] = get_prophet_score(madmed_nds_df, model)  
    madmed_nds_df.to_csv("data/scores/nds_madmed_pre.csv", index=False) 
    med_ds_df["PROphetScore"] = get_prophet_score(med_ds_df, model)     
    med_ds_df.to_csv("data/scores/ds_median_pre.csv", index=False)
    med_nds_df["PROphetScore"] = get_prophet_score(med_nds_df, model) 
    med_nds_df.to_csv("data/scores/nds_median_pre.csv", index=False) 
    
    # Double bridging step strategies (post-normalization)
    double_madmed_ds_df = post_normalization_bridging(madmed_ds_df, post_norm_ds_params)
    double_madmed_ds_df['PROphetScore'] = get_prophet_score(double_madmed_ds_df, model)
    double_madmed_ds_df.to_csv("data/scores/ds_madmed_double.csv", index=False)
    double_madmed_nds_df = post_normalization_bridging(madmed_nds_df, post_norm_nds_params)
    double_madmed_nds_df['PROphetScore'] = get_prophet_score(double_madmed_nds_df, model)
    double_madmed_nds_df.to_csv("data/scores/nds_madmed_double.csv", index=False)
    double_med_ds_df = post_normalization_bridging(med_ds_df, post_norm_ds_params)  
    double_med_ds_df['PROphetScore'] = get_prophet_score(double_med_ds_df, model)
    double_med_ds_df.to_csv("data/scores/ds_median_double.csv", index=False)
    double_med_nds_df = post_normalization_bridging(med_nds_df, post_norm_nds_params)
    double_med_nds_df['PROphetScore'] = get_prophet_score(double_med_nds_df, model)
    double_med_nds_df.to_csv("data/scores/nds_median_double.csv", index=False)
    mdacc_matched_samples = pd.read_excel(MDACC_PATH)
    mdacc_edta_ids = all_edta_info_df[all_edta_info_df['MeasureId'].isin(mdacc_matched_samples['PROPHET_EDTA_MeasureId'])] 
    mdacc_edta = post_anml_prophet_edta[post_anml_prophet_edta['MeasureId'].isin(mdacc_edta_ids["MeasureId"])]
    mdacc_edta["PROphetScore"] = get_prophet_score(mdacc_edta, model)
    mdacc_edta.to_csv("data/scores/mdacc_edta_scores.csv", index=False)
    # # --------------------------------------------------
    # 4. Load scored datasets
    # --------------------------------------------------
    print(4)
    # These should contain MeasureId + PROphetScore
    edta_prophet_scores_df = pd.read_csv("data/scores/edta_scores.csv")
    edta_prophet_scores_df = edta_prophet_scores_df[['MeasureId', 'PROphetScore']]
    mdacc_edta_scores_df = pd.read_csv("data/scores/mdacc_edta_scores.csv")
    mdacc_edta_scores_df = mdacc_edta_scores_df[['MeasureId', 'PROphetScore']]
    # matched_edta = pd.read_csv('data/excels/20260309_adat_data_measurements.csv')
    # matched_edta = matched_edta[matched_edta['MeasureId'].str.startswith('OH2026') | matched_edta['MeasureId'].str.startswith('OH2025_05')]
    
    # matched_edta["PROphetScore"] = get_prophet_score(matched_edta, model)
    # matched_edta.to_csv("data/scores/matched_edta_scores.csv", index=False)
    matched_edta = pd.read_csv("data/scores/matched_edta_scores.csv")
    matched_edta = matched_edta[['MeasureId', 'PROphetScore']]
    # edta_prophet_scores_df = pd.concat([edta_prophet_scores_df, matched_edta]).drop_duplicates(subset='MeasureId', keep='first').reset_index(drop=True)
    ds_madmed_pre_scores_df = pd.read_csv("data/scores/ds_madmed_pre.csv")
    ds_median_pre_scores_df = pd.read_csv("data/scores/ds_median_pre.csv")

    nds_madmed_pre_scores_df = pd.read_csv("data/scores/nds_madmed_pre.csv")
    nds_median_pre_scores_df = pd.read_csv("data/scores/nds_median_pre.csv")
    
    double_madmed_ds_df = pd.read_csv("data/scores/ds_madmed_double.csv")
    double_med_ds_df = pd.read_csv("data/scores/ds_median_double.csv")
    double_madmed_nds_df = pd.read_csv("data/scores/nds_madmed_double.csv")
    double_med_nds_df = pd.read_csv("data/scores/nds_median_double.csv")

    # late_plates = 
    
    # Matched sample lookup
    # matched_samples = pd.read_excel("data/excels/all_matched_samples_streck.xlsx")
    # matched_wide_df = build_matched_measure_table(matched_samples)

    # --------------------------------------------------
    # 5. Strategy registry
    # --------------------------------------------------
    print(5)
    strategy_registry = [
        {"name": "DS MADMED pre-normalization single bridging", "matrix": MATRIX_DS2, "scores_df": ds_madmed_pre_scores_df, 'MDACC': True},
        {"name": "DS Median pre-normalization single bridging", "matrix": MATRIX_DS2, "scores_df": ds_median_pre_scores_df, 'MDACC': True},
        {"name": "nDS MADMED pre-normalization single bridging", "matrix": MATRIX_NDS, "scores_df": nds_madmed_pre_scores_df, 'MDACC': False},
        {"name": "nDS Median pre-normalization single bridging", "matrix": MATRIX_NDS, "scores_df": nds_median_pre_scores_df, 'MDACC': False},
        {"name": "DS MADMED-Median double bridging", "matrix": MATRIX_DS2, "scores_df": double_madmed_ds_df, 'MDACC': True},
        {"name": "DS Median-Median double bridging", "matrix": MATRIX_DS2, "scores_df": double_med_ds_df, 'MDACC': True},
        {"name": "nDS MADMED-Median double bridging", "matrix": MATRIX_NDS, "scores_df": double_madmed_nds_df, 'MDACC': False},
        {"name": "nDS Median-Median double bridging", "matrix": MATRIX_NDS, "scores_df": double_med_nds_df, 'MDACC': False},
        # {"name": "MDACC MADMED pre-normalization single bridging", "matrix": MATRIX_DS2, "scores_df": ds_madmed_pre_scores_df, 'MDACC': True},
        # {"name": "MDACC Median pre-normalization single bridging", "matrix": MATRIX_DS2, "scores_df": ds_median_pre_scores_df, 'MDACC': True},
        # {"name": "MDACC MADMED-Median double bridging", "matrix": MATRIX_DS2, "scores_df": double_madmed_ds_df, 'MDACC': True},
        # {"name": "MDACC Median-Median double bridging", "matrix": MATRIX_DS2, "scores_df": double_med_ds_df, 'MDACC': True},

    ]

    # --------------------------------------------------
    # 6. Run evaluation
    # --------------------------------------------------
    print(6)
    # results = []
    evaluation_dict = {}
    for strategy in strategy_registry:

        strategy_name = strategy["name"]
        matrix_type = strategy["matrix"]
        streck_scores_df = strategy["scores_df"]

        print(f"Evaluating {matrix_type} {strategy_name}")
        strategy_distribution_metrics = evaluate_strategy_distribution(strategy_name=strategy_name,
                                                                       matrix_type=matrix_type,
                                                                       edta_scores_df=edta_prophet_scores_df,
                                                                       streck_scores_df=streck_scores_df,
                                                                       save_dir=f"projects/pre_post_bridging/results/analysis/{matrix_type}_{strategy_name}")
        
        evaluation_dict[strategy_name] = {'distribution_metrics': strategy_distribution_metrics}
        if strategy['MDACC']:
            matched_samples_subset = pd.read_excel(MDACC_PATH)
            matched_samples_subset = matched_samples_subset.rename(columns={
            "PROPHET_EDTA_MeasureId": "MeasureId_EDTA",
            "STRECK_MeasureId": "MeasureId_Streck",
            })
            edta_scores_df = mdacc_edta_scores_df
           
          
        else:
            edta_scores_df = edta_prophet_scores_df
            matched_samples_subset = matched_samples_nds
            # else:
            #     matched_samples_subset = matched_samples_ds
        matched_samples_metrics = evaluate_matched_samples(
            strategy_name=strategy_name,
            matrix_type=matrix_type,
            edta_scores_df=edta_scores_df,
            streck_scores_df=streck_scores_df,
            matched_df=matched_samples_subset,
            save_dir=f"projects/pre_post_bridging/results/analysis/{matrix_type}_{strategy_name}"
        )
        evaluation_dict[strategy_name]['matched_samples_metrics'] = matched_samples_metrics
        

        # results.append(metrics)

    # --------------------------------------------------
    # 7. Save summary metrics
    # --------------------------------------------------
    print(7)

    with open("projects/pre_post_bridging/results/analysis/strategy_evaluation_summary.pkl", "wb") as f:
        pickle.dump(evaluation_dict, f)
    # results_df = pd.DataFrame(results)

    # results_df.to_csv("projects/pre_post_bridging/results/strategy_summary.csv", index=False)

    print("\nEvaluation complete. <Insert Austin Powers joke here>")
    # print(results_df.sort_values("agreement", ascending=False))


# --------------------------------------------------
# Script entrypoint
# --------------------------------------------------

if __name__ == "__main__":
    main()