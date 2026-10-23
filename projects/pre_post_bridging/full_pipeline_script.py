
# %%

from plotting_functions import *

# %%
base_rfus_path = 'data/adat/base/v_166_OH2026_004.parquet'
info_path = 'data/excels/20260126_adat_data_samples_with_sample_info.csv'

# %%
base_rfus, info_df, edta_rfu_df, ds_rfu_df, nds_rfu_df = get_base_rfu_data(base_rfus_path, info_path)

# %%
plates = [f"OH2025_05{i}" for i in range(1, 8)] 
plates += [f"OH2026_00{i}" for i in range(1, 4)]

# %%
ds_rfu_df = ds_rfu_df[ds_rfu_df['PlateId'].isin(plates)]

# %%
# nds_rfu_df

# %%
# ds_rfu_df

# %%
aptamers = get_aptamers(ds_rfu_df)

# %%

params_ds_051_003 = fit_mad_log_transform(ds_rfu_df, edta_rfu_df, aptamers)

# %%
params_nds_051_003 = fit_mad_log_transform(nds_rfu_df, edta_rfu_df, aptamers)

# %%

# %%
# import json
# with open('./data/params/pre_hybridization/ds_log_04032026.json', 'w') as fp:
#     json.dump(params_ds_051_003, fp)
# with open('./data/params/pre_hybridization/nds_log_04032026.json', 'w') as fp:
#     json.dump(params_nds_051_003, fp)

# %% [markdown]
# # Transform plates

# %%
from src.altering_plate_multiple_collection_methods import alter_scanner_files_MAD_bridging as alter_multiple_mad
from src.altering_plate_multiple_collection_methods import alter_scanner_files_median_bridging as alter_multiple_median

# %%
plates = [f"OH2026_00{i}" for i in range(4, 10)]

# %%
# for plate in plates:
#     alter_multiple_mad(f'./data/plates/original_with_altered_workbooks/{plate}',
#             f'./data/plates/original_with_altered_workbooks/{plate}/{plate} Workbook.xlsx',
            
#             {'Streck': './data/params/pre_hybridization/ds_log_04032026.json',
#             'NDS': './data/params/pre_hybridization/nds_log_04032026.json'},
#             f'projects/pre_post_bridging/results/plates/madmed_transformed_10032026/{plate}')

# # %%
# for plate in plates:
#     alter_multiple_median(f'./data/plates/original_with_altered_workbooks/{plate}',
#             f'./data/plates/original_with_altered_workbooks/{plate}/{plate} Workbook.xlsx',
            
#             {'Streck': './data/params/pre_hybridization/ds_log_04032026.json',
#             'NDS': './data/params/pre_hybridization/nds_log_04032026.json'},
#             f'projects/pre_post_bridging/results/plates/median_transformed_09032026/{plate}')


# # Flags


# %%
bridged_madmed_adat_path = 'projects/pre_post_bridging/results/adats/madmed_09032026'
bridged_median_adat_path = 'projects/pre_post_bridging/results/adats/median_09032026'

unbridged_adat_path = 'projects/pre_post_bridging/results/adats/no_bridging'

# %%
ds_unbridged_df, ds_madmed_bridged_df, nds_unbridged_df, nds_madmed_bridged_df = flag_report(bridged_madmed_adat_path, unbridged_adat_path, ds_rfu_df, nds_rfu_df)
ds_unbridged_df, ds_median_bridged_df, nds_unbridged_df, nds_median_bridged_df = flag_report(bridged_median_adat_path, unbridged_adat_path, ds_rfu_df, nds_rfu_df)

# %%
model = pd.read_pickle('data/models/DCB_nosex_gamma1_rescaled.pkl')


# %%
ds_unbridged_df['PROphetScore'] = get_prophet_score(ds_unbridged_df, model)
ds_madmed_bridged_df['PROphetScore'] = get_prophet_score(ds_madmed_bridged_df, model)
ds_median_bridged_df['PROphetScore'] = get_prophet_score(ds_median_bridged_df, model)

nds_unbridged_df['PROphetScore'] = get_prophet_score(nds_unbridged_df, model)
nds_madmed_bridged_df['PROphetScore'] = get_prophet_score(nds_madmed_bridged_df, model)
nds_median_bridged_df['PROphetScore'] = get_prophet_score(nds_median_bridged_df, model)

# %%
plot_density_histograms([ds_unbridged_df['PROphetScore'], ds_madmed_bridged_df['PROphetScore'], nds_unbridged_df['PROphetScore'], nds_madmed_bridged_df['PROphetScore']],
                 ['DS Unbridged', 'DS Madmed Bridged', 'nDS Unbridged', 'nDS Madmed Bridged'])
plot_density_histograms([ds_median_bridged_df['PROphetScore'], nds_median_bridged_df['PROphetScore']],
                 ['DS Median Bridged', 'nDS Median Bridged'], title="Median Bridging PROphet Scores", bins=20)

# %%


# %% [markdown]
# # COMPLEMENTARY(!!!) Post-ANML bridging

# %%
final_rfus = pd.read_parquet('data/adat/anmlSMP/v_249_OH2024_036.parquet').reset_index(drop=False)

# %%
final_rfus = add_measure_id(final_rfus)

# %%
final_edta_rfus = final_rfus[final_rfus['MeasureId'].isin(edta_rfu_df['MeasureId'])].drop_duplicates(subset=['MeasureId'])

# %%
ds_bridged_postANML_coefficients = fit_mad_linear_transform(ds_bridged_df, final_edta_rfus, aptamers)
nds_bridged_postANML_coefficients = fit_mad_linear_transform(nds_bridged_df, final_edta_rfus, aptamers)
ds_unbridged_postANML_coefficients = fit_mad_linear_transform(ds_unbridged_df, final_edta_rfus, aptamers)
nds_unbridged_postANML_coefficients = fit_mad_linear_transform(nds_unbridged_df, final_edta_rfus, aptamers)



# %%
ds_bridged_postANML_coefficients['10000-28']

# %%
nds_bridged_postANML_coefficients['10000-28']

# %%
ds_unbridged_postANML_coefficients['10000-28']

# %%
nds_unbridged_postANML_coefficients['10000-28']

# %%
# import json
# with open('./data/params/post_hybridization/ds_med_08032026.json', 'w') as fp:
#     json.dump(ds_bridged_postANML_coefficients, fp)
# with open('./data/params/post_hybridization/nds_med_08032026.json', 'w') as fp:
#     json.dump(nds_bridged_postANML_coefficients, fp)
# with open('./data/params/post_hybridization/ds_unbridged_med_08032026.json', 'w') as fp:
#     json.dump(ds_unbridged_postANML_coefficients, fp)
# with open('./data/params/post_hybridization/nds_unbridged_med_08032026.json', 'w') as fp:
#     json.dump(nds_unbridged_postANML_coefficients, fp)

# %%
from src.bridging import transform_median_post_anml

# %%
ds_bridged_postANML = transform_median_post_anml(ds_bridged_df, ds_bridged_postANML_coefficients)
nds_bridged_postANML = transform_median_post_anml(nds_bridged_df, nds_bridged_postANML_coefficients)
ds_unbridged_postANML = transform_median_post_anml(ds_unbridged_df, ds_unbridged_postANML_coefficients)
nds_unbridged_postANML = transform_median_post_anml(nds_unbridged_df, nds_unbridged_postANML_coefficients)



# %%
nds_bridged_postANML[aptamers]/nds_bridged_df[aptamers]

# %%
ds_bridged_postANML['PROphetScore'] = get_prophet_score(ds_bridged_postANML, model)
nds_bridged_postANML['PROphetScore'] = get_prophet_score(nds_bridged_postANML, model)
ds_unbridged_postANML['PROphetScore'] = get_prophet_score(ds_unbridged_postANML, model)
nds_unbridged_postANML['PROphetScore'] = get_prophet_score(nds_unbridged_postANML, model)

# %%
plot_density_histograms([ds_unbridged_postANML['PROphetScore'], ds_bridged_postANML['PROphetScore'], nds_unbridged_postANML['PROphetScore'], nds_bridged_postANML['PROphetScore']],
                 ['DS Unbridged', 'DS Bridged', 'nDS Unbridged', 'nDS Bridged'], title="Post-ANML median bridging PROphet Scores", bins=20)

# %%
plot_density_histograms([ds_bridged_postANML['PROphetScore'], nds_bridged_postANML['PROphetScore']],
                 ['DS Bridged', 'nDS Bridged'], title="Post-ANML median bridging PROphet Scores", bins=20)

# %%
plot_density_histograms([ds_bridged_postANML['PROphetScore']],
                 ['DS Bridged'], title="Post-ANML median bridging PROphet Scores", bins=20)

# %%
(ds_bridged_postANML['PROphetScore']>=5).sum()/len(ds_bridged_postANML['PROphetScore'])

# %%
(nds_bridged_postANML['PROphetScore']>=5).sum()/len(nds_bridged_postANML['PROphetScore'])

# %% [markdown]
# # Spearman comparisons

# %% [markdown]
# ## 1. Initial

# %%
bridged_madmed_adat_path = 'projects/pre_post_bridging/results/adats/madmed_09032026'
unbridged_adat_path = 'projects/pre_post_bridging/results/adats/no_bridging'

# %%
initial_unbridged = concat_initial_dfs(unbridged_adat_path)
initial_bridged = concat_initial_dfs(bridged_madmed_adat_path)

# %%
initial_unbridged

# %%
initial_bridged

# %%
initial_bridged = add_measure_id(initial_bridged)
initial_unbridged = add_measure_id(initial_unbridged)

# %%
initial_bridged_ds = initial_bridged[initial_bridged['MeasureId'].isin(ds_rfu_df['MeasureId'])]
initial_bridged_nds = initial_bridged[initial_bridged['MeasureId'].isin(nds_rfu_df['MeasureId'])]
initial_unbridged_ds = initial_unbridged[initial_unbridged['MeasureId'].isin(ds_rfu_df['MeasureId'])]
initial_unbridged_nds = initial_unbridged[initial_unbridged['MeasureId'].isin(nds_rfu_df['MeasureId'])]

# %%
initial_bridged_ds[aptamers]

# %%
initial_unbridged_ds[['MeasureId'] + aptamers]

# %%
initial_bridged_ds['10001-7']/initial_unbridged_ds['10001-7']

# %%
spearman = initial_unbridged_ds[aptamers].corrwith(initial_bridged_ds[aptamers], method="spearman")
print(spearman[spearman < 0.999])

# %%


# %%
assert np.allclose(
    initial_unbridged_ds[aptamers].corrwith(initial_bridged_ds[aptamers], method="spearman"),
    1,
    atol=1e-3
)

# %% [markdown]
# ## 2. Final

# %%


# %%


# %%
assert np.allclose(
    ds_bridged_df[aptamers].corrwith(ds_bridged_postANML[aptamers], method="spearman"),
    1,
    atol=1e-12
)

# %%
spearman = ds_bridged_df[aptamers].corrwith(ds_bridged_postANML[aptamers], method="spearman")
print(spearman[spearman < 0.999999])

# %%
ds_bridged_postANML['10000-28'].corr(ds_bridged_df['10000-28'], method="spearman")

# %% [markdown]
# ## 3. Initial vs Final

# %%
spearman = initial_bridged_ds[aptamers].corrwith(ds_bridged_postANML[aptamers], method="spearman")
print(spearman[spearman < 0.999999])

# %%
spearman.hist(bins=100, density=True)
plt.xlabel("Spearman correlation coefficient")
plt.ylabel("Frequency density")
plt.title("Distribution of Spearman correlations between initial and post-ANML, all aptamers")


# %%
raps = [key for key in model.RAPs.keys()]

# %%
spearman[raps].hist(bins=100, density=True)
plt.xlabel("Spearman correlation coefficient, RAPs only")
plt.ylabel("Frequency density")
plt.title("Distribution of Spearman correlations between initial and post-ANML, RAPs only")

# %%


# %% [markdown]
# # Matched Samples

# %%
ds_bridged_df['PROphetScore']

# %%
nds_bridged_df['PROphetScore']

# %%
ds_bridged_postANML['PROphetScore']

# %%
matched_samples = pd.read_excel('data/excels/all_matched_samples_streck.xlsx')

# %%
matched_samples

# %%
matched_edta = matched_samples[matched_samples['SampleRun'] == 'EDTA']
matched_ds = matched_samples[matched_samples['SampleRun'] == 'DS']
matched_nds = matched_samples[matched_samples['SampleRun'] == 'non DS']

# %%
matched_edta

# %% [markdown]
# 

# %%
matched_samples['SubjectId'].value_counts()

# %%
matched_edta_ds = pd.merge(matched_edta, matched_ds, how='inner', on='SubjectId', suffixes=('_edta', '_ds'))

# %%
matched_edta_ds

# %%
matched_edta_nds = pd.merge(matched_edta, matched_nds, how='inner', on='SubjectId', suffixes=('_edta', '_nds'))

# %%
matched_edta_nds.shape

# %%
matched_edta_nds[matched_edta_nds['MeasureId_nds'].isin(initial_unbridged_nds['MeasureId'])]

# %%
matched_edta_ds[matched_edta_ds['MeasureId_ds'].isin(initial_unbridged_ds['MeasureId'])]

# %% [markdown]
# ## Matching EDTA-nDS scores after single and double bridging

# %% [markdown]
# ### Single bridging

# %%
nds_bridged_df[nds_bridged_df['MeasureId'].isin(matched_edta_nds['MeasureId_nds'])]

# %%
nds_bridged_postANML[nds_bridged_postANML['MeasureId'].isin(matched_edta_nds['MeasureId_nds'])]

# %%
nds_bridged_postANML[aptamers]/nds_bridged_df[aptamers]

# %%
nds_bridged_df_short = nds_bridged_df[['MeasureId', 'PROphetScore']]
nds_bridged_postANML_short = nds_bridged_postANML[['MeasureId', 'PROphetScore']]
nds_bridged_df_short.columns = ['MeasureId_nds', 'PROphetScore_nds']
nds_bridged_postANML_short.columns = ['MeasureId_nds', 'PROphetScore_nds']

# %%
nds_bridged_df_short

# %%
nds_bridged_postANML_short

# %%
matched_edta_nds_raw = matched_edta_nds.copy()

# %%
matched_edta_nds_single_bridging = pd.merge(matched_edta_nds, nds_bridged_df_short, how='inner', on='MeasureId_nds')
matched_edta_nds_postANML_double_bridging = pd.merge(matched_edta_nds, nds_bridged_postANML_short, how='inner', on='MeasureId_nds')

# %%
matched_edta_nds_postANML_double_bridging

# %%
nds_bridged_postANML_short

# %%
edta_final_measurements = pd.read_csv('data/excels/20260126_adat_data_measurements.csv')

# %%
edta_final_measurements = edta_final_measurements[edta_final_measurements['MeasureId'].isin(matched_edta_nds['MeasureId_edta'])]

# %%
edta_final_measurements['PROphetScore_edta'] = get_prophet_score(edta_final_measurements, model)

# %%
edta_final_measurements = edta_final_measurements[['MeasureId', 'PROphetScore_edta']]
edta_final_measurements.columns = ['MeasureId_edta', 'PROphetScore_edta']

# %%
matched_edta_nds_single_bridging

# %%
edta_nds = pd.merge(matched_edta_nds_single_bridging, edta_final_measurements, how='inner', on='MeasureId_edta')
edta_nds_postANML = pd.merge(matched_edta_nds_postANML_double_bridging, edta_final_measurements, how='inner', on='MeasureId_edta')

# %%
edta_nds_postANML

# %%
edta_nds[['PROphetScore_edta','PROphetScore_nds']]

# %%
edta_nds_postANML #[['PROphetScore_edta','PROphetScore_nds']]

# %%



# %%
scatter_with_regression(edta_nds['PROphetScore_edta'], edta_nds['PROphetScore_nds'], xlabel="EDTA PROphet Score", ylabel="nDS PROphet Score", 
                        title="EDTA vs nDS PROphet Scores, prehybridization bridging only")

# %%
scatter_with_regression(edta_nds_postANML['PROphetScore_edta'], edta_nds_postANML['PROphetScore_nds'], xlabel="EDTA PROphet Score", ylabel="nDS PROphet Score", 
                        title="EDTA vs nDS PROphet Scores, double bridging")

# %%
confusion_double_bridging = pd.crosstab(edta_nds_postANML['PROphetScore_edta']>=5, edta_nds_postANML['PROphetScore_nds']>=5, rownames=['EDTA'], colnames=['nDS'])

# %%
print(confusion_double_bridging)

# %%
confusion_single_bridging = pd.crosstab(edta_nds['PROphetScore_edta']>=5, edta_nds['PROphetScore_nds']>=5, rownames=['EDTA'], colnames=['nDS'])

# %%
print(confusion_single_bridging)

# %%
import matplotlib.pyplot as plt

# %%
plot_confusion_matrix(confusion_single_bridging, title="Single Bridging Confusion Matrix", xlabel='nDS', ylabel="EDTA", save_path="projects/pre_post_bridging/results/confusion_matrices/single_bridging_confusion.png")

# %%
plot_confusion_matrix(confusion_double_bridging, title="Double Bridging Confusion Matrix",  xlabel='nDS', ylabel="EDTA", save_path="projects/pre_post_bridging/results/confusion_matrices/double_bridging_confusion.png")

# %% [markdown]
# ## Matching EDTA-DS scores after single and double bridging

# %%
matched_ds

# %%
pd.merge(matched_edta, matched_ds, how='inner', on='SubjectId', suffixes=('_edta', '_ds'))['MeasureId_ds'].str.split('_').str[:2].str.join("_").value_counts()

# %%


# %% [markdown]
# # MDACC

# %%
matched_mdacc_df = pd.read_excel('data/excels/matched MDACC proteomics with scores.xlsx')

# %%
matched_mdacc_df

# %%
matched_mdacc_df.columns

# %%
matched_mdacc_df['plate'] = matched_mdacc_df['STRECK_MeasureId'].str.split('_').str[:2].str.join("_")

# %%
matched_mdacc_df['plate'].value_counts()

# %%
streck_matched_plates = matched_mdacc_df[matched_mdacc_df['plate'].str.startswith(
    'OH2026') | matched_mdacc_df['plate'].str.endswith(
        '052') | matched_mdacc_df['plate'].str.endswith(
            '054') | matched_mdacc_df['plate'].str.endswith('055')]['plate'].unique()

# %%
streck_matched_plates

# %%
streck_match_df = matched_mdacc_df[matched_mdacc_df['plate'].str.startswith(
    'OH2026') | matched_mdacc_df['plate'].str.endswith(
        '052') | matched_mdacc_df['plate'].str.endswith(
            '054') | matched_mdacc_df['plate'].str.endswith('055')]

# %%
streck_match_df

# %%
base_rfus.columns[:30]

# %%
base_rfus['SampleType'].value_counts()

# %%
info_df['SampleMatrixType'].value_counts()

# %%
info_df = pd.read_csv('./data/excels/20260309_adat_data_samples_with_sample_info.csv')

# %%
info_df[info_df['PlateId']=='OH2026_004'][['MeasureId', 'SampleMatrixType', 'SubType']]

# %%
info_df['SampleMatrixType'].value_counts()

# %%
info_df[(info_df['PlateId']=='OH2026_009') & (info_df['SampleMatrixType']=='Streck non-DS')][['PlatePosition', 'SampleMatrixType']]

# %%
matched_mdacc_df[matched_mdacc_df['plate'].str.startswith(
    'OH2026') | matched_mdacc_df['plate'].str.endswith(
        '052') | matched_mdacc_df['plate'].str.endswith('054')]

# %%
mdacc_initial_adats_df = concat_initial_dfs('projects/pre_post_bridging/results/adats/madmed_for_MDACC_comparison')
mdacc_final_adats_df = concat_final_dfs('projects/pre_post_bridging/results/adats/madmed_for_MDACC_comparison')

# %%
mdacc_initial_adats_df = add_measure_id(mdacc_initial_adats_df)
mdacc_final_adats_df = add_measure_id(mdacc_final_adats_df)

# %%
mdacc_final_adats_df.head()

# %%
mdacc_streck_final = mdacc_final_adats_df[mdacc_final_adats_df['MeasureId'].isin(streck_match_df['STRECK_MeasureId'])]

# %%
mdacc_streck_final

# %%
ds_bridged_postANML_coefficients

# %%
with open('./data/params/post_hybridization/ds_med_08032026.json', 'r') as fp:
    dummy = json.load(fp)

# %%
for key in dummy.keys():
    if dummy[key] != ds_bridged_postANML_coefficients[key]:
        print(f"Mismatch for {key}: {dummy[key]} vs {ds_bridged_postANML_coefficients[key]}")

# %%
mdacc_streck_final_transformed = transform_median_post_anml(mdacc_streck_final, ds_bridged_postANML_coefficients)

# %%
mdacc_streck_final_transformed['PROphetScore'] = get_prophet_score(mdacc_streck_final_transformed, model)

# %%
mdacc_streck_final['PROphetScore'] = get_prophet_score(mdacc_streck_final, model)

# %%
mdacc_streck_final_transformed_short = mdacc_streck_final_transformed[['MeasureId', 'PROphetScore']]
mdacc_streck_final_short = mdacc_streck_final[['MeasureId', 'PROphetScore']]

# %%
mdacc_streck_final_transformed_short.rename(columns={'MeasureId': 'STRECK_MeasureId', 'PROphetScore': 'PROphetScore_streck'})

# %%
matches = matched_mdacc_df[matched_mdacc_df['STRECK_MeasureId'].isin(mdacc_streck_final_short['MeasureId'])][['STRECK_MeasureId', 'PROPHET_EDTA_MeasureId']]

# %%
matches

# %%
adat_measurements = pd.read_csv('./data/excels/20260309_adat_data_measurements.csv')

# %%
edta_final_measurements = adat_measurements[adat_measurements['MeasureId'].isin(matches['PROPHET_EDTA_MeasureId'])]

# %%
edta_final_measurements['PROphetScore'] = get_prophet_score(edta_final_measurements, model)

# %%
edta_final_measurements = pd.merge(matches, 
                                   edta_final_measurements[['MeasureId', 'PROphetScore']],
                                   how='inner',
                                   left_on='PROPHET_EDTA_MeasureId',
                                   right_on='MeasureId')[['STRECK_MeasureId', 'PROPHET_EDTA_MeasureId', 'PROphetScore']]

# %%
edta_final_measurements = edta_final_measurements.rename(columns={'PROphetScore': 'PROphetScore_edta'})

# %%
# pd.merge(edta_final_measurements,mdacc_streck_final_transformed_short, 

# %%
# edta_final_measurements.head(1)

# %%
# mdacc_streck_final_transformed_short.head(1)

# %%
mdacc_streck_final_transformed_short.rename(columns={'MeasureId': 'STRECK_MeasureId', 'PROphetScore': 'PROphetScore_streck'}, inplace=True)

# %%
mdacc_ds_comparison = pd.merge(edta_final_measurements, mdacc_streck_final_transformed_short, how='inner', on='STRECK_MeasureId')[['PROphetScore_edta', 'PROphetScore_streck']]

# %%
confusion_double_bridging = pd.crosstab(edta_nds_postANML['PROphetScore_edta']>=5, edta_nds_postANML['PROphetScore_nds']>=5, rownames=['EDTA'], colnames=['nDS'])

# %%
# mdacc_ds_comparison.head(1)

# %%
confusion_mdacc_comparison = pd.crosstab(mdacc_ds_comparison['PROphetScore_edta']>=5, mdacc_ds_comparison['PROphetScore_streck']>=5, rownames=['EDTA'], colnames=['STRECK'])

# %%
plot_confusion_matrix(confusion_mdacc_comparison, title="MDACC DS-EDTA Confusion Matrix", xlabel='STRECK', ylabel="EDTA", save_path="projects/pre_post_bridging/results/confusion_matrices/mdacc_comparison_confusion.png")

# %%
scatter_with_regression(mdacc_ds_comparison['PROphetScore_edta'], mdacc_ds_comparison['PROphetScore_streck'], xlabel="EDTA PROphet Score", ylabel="STRECK DS PROphet Score", 
                        title="EDTA vs DS PROphet Scores, double bridging")

# %%


# %%


# %% [markdown]
# ## BOOM

# %%


# # %%


# # %%


# # %%


# # %%
# ds_unbridged_df['PROphetScore'] = get_prophet_score(ds_unbridged_df, model)
# ds_bridged_df['PROphetScore'] = get_prophet_score(ds_bridged_df, model)
# nds_unbridged_df['PROphetScore'] = get_prophet_score(nds_unbridged_df, model)
# nds_bridged_df['PROphetScore'] = get_prophet_score(nds_bridged_df, model)

# # %%
# overlay_histograms([ds_unbridged_df['PROphetScore'], ds_bridged_df['PROphetScore'], nds_unbridged_df['PROphetScore'], nds_bridged_df['PROphetScore']],
#                  ['DS Unbridged', 'DS Bridged', 'nDS Unbridged', 'nDS Bridged']  )

# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%


# # %%
# ds_unbridged_df = get_prophet_score(ds_unbridged_df, model)
# ds_bridged_df = get_prophet_score(ds_bridged_df, model)
# nds_unbridged_df = get_prophet_score(nds_unbridged_df, model)
# nds_bridged_df = get_prophet_score(nds_bridged_df, model)

# # %%


# # %%
# ds_bridged_df = bridged_df[bridged_df['MeasureId'].isin(ds_rfu_df['MeasureId'])]
# ds_unbridged_df = unbridged_df[unbridged_df['MeasureId'].isin(ds_rfu_df['MeasureId'])]

# nds_bridged_df = bridged_df[bridged_df['MeasureId'].isin(nds_rfu_df['MeasureId'])]
# nds_unbridged_df = unbridged_df[unbridged_df['MeasureId'].isin(nds_rfu_df['MeasureId'])]

# # %%
# ds_bridged_df['RowCheck'].value_counts()

# # %%
# ds_unbridged_df['RowCheck'].value_counts()

# # %%
# nds_bridged_df['RowCheck'].value_counts()

# # %%
# nds_unbridged_df['RowCheck'].value_counts()

# # %%
# ds_bridged_df['ProphetScore'] = get_prophet_score(ds_bridged_df, model)
# ds_unbridged_df['ProphetScore'] = get_prophet_score(ds_unbridged_df, model)
# nds_bridged_df['ProphetScore'] = get_prophet_score(nds_bridged_df, model)
# nds_unbridged_df['ProphetScore'] = get_prophet_score(nds_unbridged_df, model)

# # %%
# ds_bridged_df['ProphetScore'].hist(alpha=0.5, label='Bridged DS')

# # %%
# ds_unbridged_df['ProphetScore'].hist(alpha=0.5, label='Unbridged DS')

# # %%
# nds_bridged_df['ProphetScore'].hist(alpha=0.5, label='Bridged nDS')

# # %%
# nds_unbridged_df['ProphetScore'].hist(alpha=0.5, label='Unbridged nDS')

# # %%


# # %%
# plot_overlaid_histograms(ds_bridged_df['ProphetScore'], ds_unbridged_df['ProphetScore'], 'Bridged DS Prophet Scores', 'Unbridged DS Prophet Scores', 'DS Bridged vs Unbridged Prophet Score Distribution')

# # %%


# # %%


# # %%


# # %%


# # %% [markdown]
# # # Assessing initial distributions

# # %%
# ds_rfu_df

# # %%
# scatter_comparison(edta_rfu_df[aptamers], 
#                        ds_rfu_051_003_df[aptamers], 
#                        metric_type='mean',
#                        xlabel='EDTA', 
#                        ylabel='Streck', 
#                        title='Mean comparison of pre-normalization RFUs, no bridging',
#                        trim_percentile=100)

# # %%


# # %%
# scatter_comparison(edta_rfu_df[aptamers], 
#                        ds_rfu_051_003_df[aptamers], 
#                        metric_type='median',
#                        xlabel='EDTA', 
#                        ylabel='Streck', 
#                        title='Median comparison of pre-normalization RFUs, no bridging',
#                        trim_percentile=100)

# # %%
# scatter_comparison(edta_rfu_df[aptamers], 
#                        ds_rfu_051_003_df[aptamers], 
#                        metric_type='median',
#                        xlabel='EDTA', 
#                        ylabel='Streck', 
#                        title='Median comparison of pre-normalization RFUs, no bridging',
#                        trim_percentile=100)

# # %%
# scatter_comparison(edta_rfu_df[aptamers], 
#                        ds_rfu_051_003_df[aptamers], 
#                        metric_type='std',
#                        xlabel='EDTA', 
#                        ylabel='Streck', 
#                        title='Standard deviation comparison of pre-normalization RFUs, no bridging',
#                        trim_percentile=100)

# # %%
# scatter_comparison(edta_rfu_df[aptamers], 
#                        streck_rfu_df[aptamers], 
#                        metric_type='mad',
#                        xlabel='EDTA', 
#                        ylabel='Streck', 
#                        title='Median absolute deviation comparison of pre-normalization RFUs, no bridging',
#                        trim_percentile=100)

# # %%


# # %% [markdown]
# # # README on running SL pipeline

# # %% [markdown]
# # # Reading ADATs and assessing initial bridged values

# # %%
# import os
# adat_path = 'data/adats/'

# initial_adats = []
# final_adats = []
# for folder in os.listdir(adat_path):
#     splits = folder.split('_')
#     initial_file = splits[0] + '_' + splits[1] + '.adat'
#     folder_path = os.path.join(adat_path, folder)
#     files = os.listdir(folder_path)
#     final_file = max(files, key=len)
#     initial_adat, _ = read_adat_file(os.path.join(folder_path, initial_file))
#     final_adat, _ = read_adat_file(os.path.join(folder_path, final_file))
#     initial_adats.append(initial_adat)
#     final_adats.append(final_adat)


# # %%
# initial_bridged_df = add_measure_id(pd.concat(initial_adats, ignore_index=True))
# final_bridged_df = add_measure_id(pd.concat(final_adats, ignore_index=True))

# # %%
# initial_bridged_df

# # %%
# ds_rfu_051_003_df['MeasureId']

# # %%
# initial_bridged_ds_df = initial_bridged_df[initial_bridged_df['MeasureId'].isin(ds_rfu_051_003_df['MeasureId'])]
# initial_bridged_nds_df = initial_bridged_df[initial_bridged_df['MeasureId'].isin(nds_rfu_df['MeasureId'])]
# final_bridged_ds_df = final_bridged_df[final_bridged_df['MeasureId'].isin(ds_rfu_051_003_df['MeasureId'])]
# final_bridged_nds_df = final_bridged_df[final_bridged_df['MeasureId'].isin(nds_rfu_df['MeasureId'])]

# # %%
# final_bridged_ds_df.shape

# # %%
# nds_rfu_df['MeasureId']

# # %%
# scatter_comparison(edta_rfu_df[aptamers], 
#                        initial_bridged_nds_df[aptamers], 
#                        metric_type='median',
#                        xlabel='EDTA', 
#                        ylabel='Non-desalinated Streck', 
#                        title='Standard Median comparison of pre-normalization RFUs, bridged',
#                        trim_percentile=100)

# # %%
# (edta_rfu_df[aptamers].median() / initial_bridged_nds_df[aptamers].median()).mean()

# # %%
# scatter_comparison(np.log(edta_rfu_df[aptamers]), 
#                        np.log(initial_bridged_nds_df[aptamers]), 
#                        metric_type='mad',
#                        xlabel='EDTA', 
#                        ylabel='Non-desalinated Streck', 
#                        title='Standard MAD comparison of pre-normalization RFUs, bridged',
#                        trim_percentile=100)

# # %%
# scatter_comparison(np.log(edta_rfu_df[aptamers]), 
#                        np.log(initial_bridged_ds_df[aptamers]), 
#                        metric_type='median',
#                        xlabel='EDTA', 
#                        ylabel='Desalinated Streck', 
#                        title='Standard Median comparison of pre-normalization RFUs, bridged',
#                        trim_percentile=100)

# # %%
# with open('./data/conversion_coefficients/MAD_median_nds_coefficients_plates_OH2025_054_OH2025_055_OH2026_001_OH2026_002_04022026.json', 'r') as fp:
#     params_nds = json.load(fp)

# # %%
# params_nds

# # %%
# import pandas as pd
# import numpy as np


# # Example usage:
# # df_new = transform_columns(df, your_dict)

# # %%
# nds_rfu_df.shape

# # %%
# transformed_nds = transform_columns(nds_rfu_df, params_nds)

# # %%
# transformed_nds

# # %%
# scatter_comparison(np.log(edta_rfu_df[aptamers]), 
#                        np.log(initial_bridged_nds_df[aptamers]), 
#                        metric_type='median',
#                        xlabel='EDTA', 
#                        ylabel='Desalinated Streck', 
#                        title='Standard Median comparison of pre-normalization RFUs, bridged',
#                        trim_percentile=100)

# # %%
# scatter_comparison(np.log(edta_rfu_df[aptamers]), 
#                        np.log(transformed_nds[aptamers]), 
#                        metric_type='median',
#                        xlabel='EDTA', 
#                        ylabel='Desalinated Streck', 
#                        title='Standard Median comparison of pre-normalization RFUs, bridged',
#                        trim_percentile=100)

# # %% [markdown]
# # # PROphet score comparison

# # %%
# model = pd.read_pickle('./data/models/DCB_nosex_gamma1_rescaled.pkl')

# # %%
# final_bridged_ds_df['PROphetScore'] = get_prophet_score(final_bridged_ds_df, model)
# final_bridged_nds_df['PROphetScore'] = get_prophet_score(final_bridged_nds_df, model)

# # %%
# final_edta_df = pd.read_csv('data/20260208_adat_data_measurements.csv')

# # %%
# final_edta_df = final_edta_df[final_edta_df['MeasureId'].isin(edta_rfu_df['MeasureId'])]

# # %%
# final_edta_df.shape

# # %%
# final_edta_df['PROphetScore'] = get_prophet_score(final_edta_df, model)

# # %%
# plot_overlaid_histograms(final_edta_df['PROphetScore'], final_bridged_ds_df['PROphetScore'], 10, xtitle='EDTA PROphet score', ytitle='Desalted Streck PROphet score', trim_percentile=100)


# # %%
# plot_overlaid_histograms(final_edta_df['PROphetScore'], final_bridged_nds_df['PROphetScore'], 10, xtitle='EDTA PROphet score', ytitle='Non-desalted Streck PROphet score', trim_percentile=100)


# # %%
# (final_bridged_nds_df['PROphetScore'] >= 5).sum()

# # %%
# (final_bridged_ds_df['PROphetScore'] >= 5).sum()

# # %%
# final_bridged_nds_df['PROphetScore'].mean()

# # %%
# final_bridged_ds_df['PROphetScore'].mean()

# # %%
# nds_rfu_df.shape

# # %%
# nds_rfu_df[['MeasureId', 'PlateId', 'PlatePosition', 'PlateRunDate']].to_csv('streck_non_desalted_samples.csv')

# # %%
# nds_rfu_df.shape

# # %%
# ds_rfu_df[['MeasureId', 'PlateId', 'PlatePosition', 'PlateRunDate']].to_csv('streck_desalted_samples.csv')

# # %%
# ds_rfu_df[~ds_rfu_df['PlateId'].isin(['OH2025_049', 'OH2025_050'])][['MeasureId', 'PlateId', 'PlatePosition', 'PlateRunDate']].reset_index(drop=True).to_csv('streck_desalted_samples.csv')

# # %%
# ds_rfu_df[~ds_rfu_df['PlateId'].isin(['OH2025_049', 'OH2025_050'])][['MeasureId', 'PlateId', 'PlatePosition', 'PlateRunDate']]

# # %%
# initial_df = pd.concat(initial_adats, ignore_index=True)
# final_df = pd.concat(final_adats, ignore_index=True)                    

# # %%
# scatter_comparison(edta_rfu_df[aptamers], 
#                        initial_df[aptamers], 
#                        metric_type='median',
#                        xlabel='EDTA', 
#                        ylabel='Streck', 
#                        title='Median comparison of pre-normalization RFUs, after bridging',
#                        trim_percentile=100)

# # %%
# scatter_comparison(edta_rfu_df[aptamers], 
#                        initial_df[aptamers], 
#                        metric_type='mad',
#                        xlabel='EDTA', 
#                        ylabel='Streck', 
#                        title='MAD comparison of pre-normalization RFUs, after bridging',
#                        trim_percentile=100)

# # %%
# scatter_comparison(edta_rfu_df[aptamers], 
#                        final_df[aptamers], 
#                        metric_type='median',
#                        xlabel='EDTA', 
#                        ylabel='Streck', 
#                        title='Median comparison of post-normalization RFUs, after bridging',
#                        trim_percentile=90)

# # %% [markdown]
# # # Comparing distributions

# # %%


# # %%


# # %%


# # %%



