import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import r2_score
from scipy.stats import ks_2samp
from src.dataframe_transformation import add_measure_id
from src.plotting_functions import plot_density_histograms, plot_confusion_matrix

# ----------------------------
# Matrix type definitions
# ----------------------------

MATRIX_EDTA = "EDTA Plasma"
MATRIX_DS1 = "Streck DS1"
MATRIX_DS2 = "Streck DS2"
MATRIX_NDS = "Streck non-DS"


def post_normalization_bridging(df, params):
    aptamers = list(params.keys())
    for aptamer in aptamers:
        df[aptamer] = df[aptamer] * params[aptamer]['med_edta'] / params[aptamer]['med_streck']
    return df

def load_sample_metadata(info_path):
    """
    Load and prepare sample metadata for evaluation.

    Parameters
    ----------
    info_path : str
        Path to the metadata CSV file.

    Returns
    -------
    dict
        Dictionary containing:
        - "info": full cleaned metadata dataframe
        - "EDTA": metadata subset for EDTA samples
        - "DS": metadata subset for Streck DS samples
        - "nDS": metadata subset for Streck non-DS samples
    """
    
    # Load metadata

    info_df = pd.read_csv(info_path)
    available_matrices = set(info_df['SampleMatrixType'].unique())
    missing = {MATRIX_EDTA, MATRIX_DS1, MATRIX_DS2, MATRIX_NDS} - available_matrices
    if missing:
        raise ValueError(f"Missing expected SampleMatrixType values: {missing}")

    # ---- PlateId normalization (historical dataset compatibility) ----

    info_df['PlateId'] = info_df['PlateId'].replace('OH2023_001', 'OH2023_01')
    info_df['PlateId'] = info_df['PlateId'].str.split('_set').str[0]
    info_df['PlateId'] = info_df['PlateId'].str.split('_Lot3').str[0]

    info_df = add_measure_id(info_df)

    # ---- Subset metadata ----

    all_edta_info_df = info_df[
        info_df['SampleMatrixType'] == MATRIX_EDTA].copy()
    
    prophet_edta_info_df = info_df[
        (info_df['SampleMatrixType'] == MATRIX_EDTA) &
        (info_df['SubType'] == 'PROPHET')
    ].copy()
    
    prophet_edta_info_df = prophet_edta_info_df.drop_duplicates(subset=['MeasureId']).reset_index(drop=True)

    ds_info_df = info_df[
        (info_df['SampleMatrixType'].isin([MATRIX_DS1, MATRIX_DS2])) #&
#        (info_df['SubType'] == 'PROPHET')
    ].copy()
    ds_info_df = ds_info_df.drop_duplicates(subset=['MeasureId']).reset_index(drop=True)

    nds_info_df = info_df[
        (info_df['SampleMatrixType'] == MATRIX_NDS)
    ].copy()
    nds_info_df = nds_info_df.drop_duplicates(subset=['MeasureId']).reset_index(drop=True)

    return {
        "info": info_df,
        "all_EDTA": all_edta_info_df,
        "PROphet_EDTA": prophet_edta_info_df,
        "DS": ds_info_df,
        "nDS": nds_info_df
    }
    
    
def build_matched_measure_table(matched_df):
    """
    Convert long-format matched_df into wide-format table expected
    by evaluate_strategy().

    Input columns
    -------------
    MeasureId
    SubjectId
    SampleRun  (values: 'EDTA', 'DS', 'non DS')

    Output columns
    --------------
    SubjectId
    MeasureId_edta
    MeasureId_ds
    MeasureId_nds
    """

    # Pivot the table
    wide = matched_df.pivot(
        index="SubjectId",
        columns="SampleRun",
        values="MeasureId"
    ).reset_index()

    # Rename columns to match evaluation code expectations
    wide = wide.rename(columns={
        "EDTA": "MeasureId_edta",
        "DS": "MeasureId_ds",
        "non DS": "MeasureId_nds"
    })

    return wide    
    
def build_matched_score_table(
    matched_df,
    edta_scores_df,
    streck_scores_df,
    matrix_type="Streck DS2"
):
    """
    Build paired EDTA–Streck PROphet score table for matched samples.

    Parameters
    ----------
    matched_df : pd.DataFrame
        Table containing SubjectId and MeasureIds for EDTA and Streck samples.

    edta_scores_df : pd.DataFrame
        DataFrame containing EDTA MeasureId and PROphetScore.

    streck_scores_df : pd.DataFrame
        DataFrame containing Streck MeasureId and PROphetScore for a strategy.

    matrix_type : str
        Either "DS" or "nDS" to select the correct matched column.

    Returns
    -------
    pd.DataFrame
        Paired score table with columns:
        SubjectId, EDTA_score, Streck_score
    """

    # if matrix_type not in ["DS", "nDS"]:
    #     raise ValueError("matrix_type must be 'DS' or 'nDS'")

    # Select correct matched column
    if "Streck DS" in matrix_type:
        streck_measure_col = "MeasureId_ds"
    elif "Streck nDS" in matrix_type:
        streck_measure_col = "MeasureId_nds"
    # streck_measure_col = f"MeasureId_{matrix_type.lower()}"
    edta_measure_col = "MeasureId_edta"

    # Extract needed columns
    matched_subset = matched_df[
        ["SubjectId", edta_measure_col, streck_measure_col]
    ].copy()

    # Join EDTA scores
    edta_lookup = edta_scores_df[["MeasureId", "PROphetScore"]].rename(
        columns={"MeasureId": edta_measure_col, "PROphetScore": "EDTA_score"}
    )

    matched_subset = matched_subset.merge(
        edta_lookup,
        on=edta_measure_col,
        how="inner"
    )

    # Join Streck scores
    streck_lookup = streck_scores_df[["MeasureId", "PROphetScore"]].rename(
        columns={"MeasureId": streck_measure_col, "PROphetScore": "Streck_score"}
    )

    matched_subset = matched_subset.merge(
        streck_lookup,
        on=streck_measure_col,
        how="inner"
    )

    # Final table
    matched_scores = matched_subset[
        ["SubjectId", "EDTA_score", "Streck_score"]
    ].copy()
    print(f"Matched pairs: {len(matched_scores)}")
    return matched_scores


import pandas as pd


def get_matching_measureids(df: pd.DataFrame):
    """
    Return two dataframes:
    1. Matching EDTA-DS MeasureIds by SubjectId
    2. Matching EDTA-non DS MeasureIds by SubjectId

    Expected columns in df:
    - MeasureId
    - SubjectId
    - SampleRun
    """

    required_cols = {"MeasureId", "SubjectId", "SampleRun"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    tmp = df.copy()
    tmp["SampleRun"] = tmp["SampleRun"].astype(str).str.strip()

    # Adjust these masks if your SampleRun values use slightly different naming
    edta_df = tmp[tmp["SampleRun"].str.fullmatch(r"EDTA", case=False, na=False)]
    ds_df = tmp[tmp["SampleRun"].str.fullmatch(r"DS", case=False, na=False)]
    non_ds_df = tmp[tmp["SampleRun"].str.fullmatch(r"non DS", case=False, na=False)]

    edta_ds_matches = edta_df[["SubjectId", "MeasureId"]].merge(
        ds_df[["SubjectId", "MeasureId"]],
        on="SubjectId",
        suffixes=("_EDTA", "_Streck")
    )

    edta_non_ds_matches = edta_df[["SubjectId", "MeasureId"]].merge(
        non_ds_df[["SubjectId", "MeasureId"]],
        on="SubjectId",
        suffixes=("_EDTA", "_Streck")
    )

    return edta_ds_matches, edta_non_ds_matches




def evaluate_strategy_distribution(
    strategy_name,
    matrix_type,
    edta_scores_df,
    streck_scores_df,
    # matched_df,
    save_dir=None,
    threshold=5
    ):
    """
    Evaluate a Streck bridged distribution against EDTA.

    This function performs a full dataset distribution comparison:

    - Overlay density / histogram plot of EDTA vs Streck PROphet scores
    - KS test statistic and p-value annotation
    """
    edta_scores = edta_scores_df["PROphetScore"]
    streck_scores = streck_scores_df["PROphetScore"]
    ks_stat, ks_pvalue = ks_2samp(edta_scores, streck_scores)
    fig_dist, ax_dist = plot_density_histograms(
        [edta_scores, streck_scores],
        names=["EDTA", f"{matrix_type} {strategy_name}"],
        title=f"PROphet Score Distribution\nEDTA vs {matrix_type} {strategy_name}"
    )
    
    ax_dist.text(
    0.02,
    0.95,
    f"KS = {ks_stat:.3f}\np = {ks_pvalue:.3e}",
    transform=ax_dist.transAxes,
    verticalalignment="top"
    )

    ax_dist.axvline(threshold, linestyle="--", alpha=0.5)
    if save_dir is not None:

        os.makedirs(save_dir, exist_ok=True)

        fig_dist.savefig(
            os.path.join(save_dir, f"{strategy_name}_distribution.png"),
            dpi=300
        )
        plt.close(fig_dist)
    else:
        pass
        

    
    return fig_dist, ax_dist, ks_stat, ks_pvalue
    

def build_matched_scores(streck_df, edta_df, match_df):
    """
    Returns a dataframe with matched scores:
    columns: Score_streck, Score_edta
    """

    # attach Streck scores
    merged = match_df.merge(
        streck_df[["MeasureId", "PROphetScore"]],
        left_on="MeasureId_Streck",
        right_on="MeasureId",
        how="inner"
    ).rename(columns={"PROphetScore": "PROphetScore_streck"}).drop(columns="MeasureId")

    # attach EDTA scores
    merged = merged.merge(
        edta_df[["MeasureId", "PROphetScore"]],
        left_on="MeasureId_EDTA",
        right_on="MeasureId",
        how="inner"
    ).rename(columns={"PROphetScore": "PROphetScore_edta"}).drop(columns="MeasureId")

    return merged[["PROphetScore_streck", "PROphetScore_edta"]]



    
def evaluate_matched_samples(
    strategy_name,
    matrix_type,
    edta_scores_df,
    streck_scores_df,
    matched_df,
    save_dir=None,
    threshold=5
    ):
    
    """
    Matched sample comparison
       - Scatter plot with linear regression (EDTA vs Streck)
       - Confusion matrix using the PROphet threshold
    """
    matches = build_matched_scores(streck_scores_df, edta_scores_df, matched_df)
    x = matches["PROphetScore_edta"]
    y = matches["PROphetScore_streck"]

    # --------------------------------------------------
    # 3. Scatter + regression
    # --------------------------------------------------

    slope, intercept = np.polyfit(x, y, 1)

    y_pred = slope * x + intercept
    r2 = r2_score(y, y_pred)
    fig_scatter, ax = plt.subplots()

    ax.scatter(x, y)

    x_line = np.linspace(min(x), max(x), 100)
    ax.plot(x_line, slope * x_line + intercept, label="Regression")
    ax.plot(
    x_line,
    slope * x_line + intercept,
    label=f"y = {slope:.2f}x + {intercept:.2f}\nR² = {r2:.2f}"
    )   

    ax.plot(
        [min(x), max(x)],
        [min(x), max(x)],
        linestyle=":",
        label="X = Y"
    )

    ax.set_xlabel("EDTA PROphet Score")
    ax.set_ylabel(f"{matrix_type} PROphet Score")
    ax.set_title(f"Matched Samples\n{strategy_name}")
    ax.legend()
    # fig_scatter = plt.figure()
    # plt.scatter(x, y)

    # x_line = np.linspace(min(x), max(x), 100)
    # plt.plot(x_line, slope * x_line + intercept, label="Regression")

    # plt.plot(
    #     [min(x), max(x)],
    #     [min(x), max(x)],
    #     linestyle=":",
    #     label="X = Y"
    # )

    # plt.xlabel("EDTA PROphet Score")
    # plt.ylabel(f"{matrix_type} PROphet Score")
    # plt.title(f"Matched Samples\n{strategy_name}")
    # plt.legend()

    # --------------------------------------------------
    # 4. Confusion matrix
    # --------------------------------------------------

    edta_class = x >= threshold
    streck_class = y >= threshold

    confusion = pd.crosstab(
        edta_class,
        streck_class,
        rownames=["EDTA"],
        colnames=[matrix_type]
    )

    confusion_metrics = plot_confusion_matrix(
        confusion,
        title=f"{strategy_name} Confusion Matrix",
        xlabel=matrix_type,
        ylabel="EDTA"
    )

    agreement = confusion_metrics["accuracy"]
    fig_confusion = confusion_metrics["fig"]

    # --------------------------------------------------
    # 5. Save or show plots
    # --------------------------------------------------

    if save_dir is not None:

        os.makedirs(save_dir, exist_ok=True)

        # fig_dist.savefig(
        #     os.path.join(save_dir, f"{strategy_name}_distribution.png"),
        #     dpi=300
        # )

        fig_scatter.savefig(
            os.path.join(save_dir, f"{strategy_name}_scatter.png"),
            dpi=300
        )
        plt.close(fig_scatter)
        
        fig_confusion.savefig(
            os.path.join(save_dir, f"{strategy_name}_confusion.png"),
            dpi=300
        )
        plt.close(fig_confusion)

    else:
        pass
        # plt.show()

    # --------------------------------------------------
    # 6. Return evaluation metrics
    # --------------------------------------------------

    return {
        "strategy": strategy_name,
        "matrix": matrix_type,
        "n_matched": len(matched_df),
        "slope": slope,
        "r2": r2,
        "scatter_plot": fig_scatter,
        "confusion_plot": fig_confusion,
        "agreement": agreement
    }

    




def evaluate_strategy(
    strategy_name,
    matrix_type,
    edta_scores_df,
    streck_scores_df,
    matched_df,
    save_dir=None,
    threshold=5
):
    """
    Evaluate a Streck bridging strategy against EDTA.

    This function performs two independent evaluations:

    1. Full dataset comparison
       - Overlay density / histogram plot of EDTA vs Streck PROphet scores

    2. Matched sample comparison
       - Scatter plot with linear regression (EDTA vs Streck)
       - Confusion matrix using the PROphet threshold

    Parameters
    ----------
    strategy_name : str
        Name of the strategy (used for plot titles and file names)

    matrix_type : str
        Either "DS" or "nDS". Determines which matched MeasureId column to use.

    edta_scores_df : pd.DataFrame
        DataFrame containing:
            MeasureId
            PROphetScore

    streck_scores_df : pd.DataFrame
        DataFrame containing:
            MeasureId
            PROphetScore
        for the specific strategy.

    matched_df : pd.DataFrame
        Table linking SubjectId to MeasureIds:
            SubjectId
            MeasureId_edta
            MeasureId_ds
            MeasureId_nds

    save_dir : str or None
        Directory where plots will be saved.
        If None, plots will be displayed instead.

    threshold : float
        PROphet classification threshold (default = 5).

    Returns
    -------
    dict
        Summary metrics including:
            n_matched
            slope
            r2
            agreement
    """


    # --------------------------------------------------
    # 1. Distribution comparison (full datasets)
    # --------------------------------------------------

    edta_scores = edta_scores_df["PROphetScore"]
    streck_scores = streck_scores_df["PROphetScore"]
    ks_stat, ks_pvalue = ks_2samp(edta_scores, streck_scores)
    fig_dist, ax_dist = plot_density_histograms(
        [edta_scores, streck_scores],
        names=["EDTA", f"{matrix_type} {strategy_name}"],
        title=f"PROphet Score Distribution\nEDTA vs {matrix_type} {strategy_name}"
    )
    
    ax_dist.text(
    0.02,
    0.95,
    f"KS = {ks_stat:.3f}\np = {ks_pvalue:.3e}",
    transform=ax_dist.transAxes,
    verticalalignment="top"
    )

    ax_dist.axvline(threshold, linestyle="--", alpha=0.5)

    
    # return {
    #     "strategy": [strategy_name]}
    # --------------------------------------------------
    # 2. Build matched score table
    # --------------------------------------------------

    matched_scores = build_matched_score_table(
        matched_df,
        edta_scores_df,
        streck_scores_df,
        matrix_type=matrix_type
    )

    x = matched_scores["EDTA_score"]
    y = matched_scores["Streck_score"]

    # --------------------------------------------------
    # 3. Scatter + regression
    # --------------------------------------------------

    slope, intercept = np.polyfit(x, y, 1)

    y_pred = slope * x + intercept
    r2 = r2_score(y, y_pred)

    fig_scatter = plt.figure()
    plt.scatter(x, y)

    x_line = np.linspace(min(x), max(x), 100)
    plt.plot(x_line, slope * x_line + intercept, label="Regression")

    plt.plot(
        [min(x), max(x)],
        [min(x), max(x)],
        linestyle=":",
        label="X = Y"
    )

    plt.xlabel("EDTA PROphet Score")
    plt.ylabel(f"{matrix_type} PROphet Score")
    plt.title(f"Matched Samples\n{strategy_name}")
    plt.legend()

    # --------------------------------------------------
    # 4. Confusion matrix
    # --------------------------------------------------

    edta_class = x >= threshold
    streck_class = y >= threshold

    confusion = pd.crosstab(
        edta_class,
        streck_class,
        rownames=["EDTA"],
        colnames=[matrix_type]
    )

    metrics = plot_confusion_matrix(
        confusion,
        title=f"{strategy_name} Confusion Matrix",
        xlabel=matrix_type,
        ylabel="EDTA"
    )

    agreement = metrics["accuracy"]

    # --------------------------------------------------
    # 5. Save or show plots
    # --------------------------------------------------

    if save_dir is not None:

        os.makedirs(save_dir, exist_ok=True)

        # fig_dist.savefig(
        #     os.path.join(save_dir, f"{strategy_name}_distribution.png"),
        #     dpi=300
        # )

        fig_scatter.savefig(
            os.path.join(save_dir, f"{strategy_name}_scatter.png"),
            dpi=300
        )
        plt.close(fig_scatter)

    else:
        pass
        # plt.show()

    # --------------------------------------------------
    # 6. Return evaluation metrics
    # --------------------------------------------------

    return {
        "strategy": strategy_name,
        "matrix": matrix_type,
        "n_matched": len(matched_scores),
        "slope": slope,
        "r2": r2,
        "scatter_plot": fig_scatter,
        "agreement": agreement
    }
