
from __future__ import annotations

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd
from scipy import stats
from src.altering_plate_multiple_collection_methods import *
from src.bridging_coefficients import *
from src.distribution_assessment import *
from src.adat_handling import read_adat_file
from src.dataframe_transformation import *
from src.prophet_scores import * 
from modeling import rap_response_predictor
from typing import Sequence, Optional, Union
from src.prophet_scores import get_prophet_score

from sklearn.metrics import r2_score



def get_base_rfu_data(base_rfus_path, info_path):
    
    
    #Reading the data
    base_rfus = pd.read_parquet(base_rfus_path).reset_index(drop=False)
    info_df = pd.read_csv(info_path)

    #Cleaning PlateId values
    base_rfus['PlateId'] = base_rfus['PlateId'].replace('OH2023_001', 'OH2023_01')
    info_df['PlateId'] = info_df['PlateId'].str.split('_set').str[0]
    base_rfus['PlateId'] = base_rfus['PlateId'].str.split('_set').str[0]
    base_rfus['PlateId'] = base_rfus['PlateId'].str.split('_Lot3').str[0]
    base_rfus = add_measure_id(base_rfus)
    
    #Filtering
    edta_info_df = filter_dataframe(info_df, info_df['SampleMatrixType']=='EDTA Plasma',
                                info_df['SubType']=='PROPHET')
    edta_rfu_df = base_rfus[base_rfus['MeasureId'].isin(edta_info_df['MeasureId'])].drop_duplicates(subset=['MeasureId']).reset_index(drop=True)
    
    ds_info_df = filter_dataframe(info_df, info_df['SampleMatrixType']=='Streck',
                                info_df['SubType']=='PROPHET')
    ds_rfu_df = base_rfus[base_rfus['MeasureId'].isin(ds_info_df['MeasureId'])].drop_duplicates(subset=['MeasureId']).reset_index(drop=True)

    nds_info_df = filter_dataframe(info_df, info_df['SampleMatrixType']=='non-DS Streck')
    nds_rfu_df = base_rfus[base_rfus['MeasureId'].isin(nds_info_df['MeasureId'])].drop_duplicates(subset=['MeasureId']).reset_index(drop=True)

    return base_rfus, info_df, edta_rfu_df, ds_rfu_df, nds_rfu_df

def concat_final_dfs(path):
    folder = Path(path)

    dfs = []

    for file in folder.rglob("*.adat"):  # search recursively in subfolders
        if "anmlSMP.adat" in file.name:
            df, _ = read_adat_file(str(file))
            dfs.append(df)
    # print(len(dfs))
    df = pd.concat(dfs, ignore_index=True)
    return df

def concat_initial_dfs(path):
    folder = Path(path)

    dfs = []

    for file in folder.rglob("*.adat"):
        if re.search(r"\d\.adat$", file.name):
            df, _ = read_adat_file(str(file))
            dfs.append(df)

    df = pd.concat(dfs, ignore_index=True)
    return df


def flag_report(bridged_adat_path, unbridged_adat_path, ds_rfu_df, nds_rfu_df):

    bridged_df = concat_final_dfs(bridged_adat_path)
    unbridged_df = concat_final_dfs(unbridged_adat_path)
    
    bridged_df['MeasureId'] = bridged_df['PlateId'] + '_' + bridged_df['PlatePosition']
    unbridged_df['MeasureId'] = unbridged_df['PlateId'] + '_' + unbridged_df['PlatePosition']
    
    ds_bridged_df = bridged_df[bridged_df['MeasureId'].isin(ds_rfu_df['MeasureId'])]
    ds_unbridged_df = unbridged_df[unbridged_df['MeasureId'].isin(ds_rfu_df['MeasureId'])]

    nds_bridged_df = bridged_df[bridged_df['MeasureId'].isin(nds_rfu_df['MeasureId'])]
    nds_unbridged_df = unbridged_df[unbridged_df['MeasureId'].isin(nds_rfu_df['MeasureId'])]
    
    print("DS unbridged flags:")
    print(ds_unbridged_df['RowCheck'].value_counts())
    print("DS bridged flags:")
    print(ds_bridged_df['RowCheck'].value_counts())
    print("nDS unbridged flags:")
    print(nds_unbridged_df['RowCheck'].value_counts())
    print("nDS bridged flags:")
    print(nds_bridged_df['RowCheck'].value_counts())
    print("****************************")
    return ds_unbridged_df, ds_bridged_df, nds_unbridged_df, nds_bridged_df

def overlay_histograms(series_list, names, bins=50, alpha=0.4):
    """
    Plot normalized histograms of multiple pandas Series on the same axes.

    Parameters
    ----------
    series_list : list of pandas.Series
        The data series to plot.
    names : list of str
        Names corresponding to each series.
    bins : int or sequence
        Number of bins or bin edges.
    alpha : float
        Transparency for overlay.
    """

    plt.figure()

    for s, name in zip(series_list, names):
        values = s.dropna().values
        plt.hist(values, bins=bins, density=True, alpha=alpha, label=name)

    plt.xlabel("Value")
    plt.ylabel("Frequency density")
    plt.legend()
    plt.tight_layout()
    plt.show()


ArrayLike = Union[Sequence[float], np.ndarray]

def plot_density_histograms(
    series_list: Sequence[ArrayLike],
    names: Optional[Sequence[str]] = None,
    title: str = "",
    bins: Union[int, Sequence[float], str] = "auto",
    alpha: Optional[float] = None,
    histtype: Optional[str] = None,
    common_range: bool = True,
    dropna: bool = True,
    trim: Optional[float] = None,  # NOW: percent to KEEP (e.g. 0.95 keeps central 95%)
    ax: Optional[plt.Axes] = None,
    linewidth: Optional[float] = None,
    ylabel: str = "Density",
    xlabel: str = "PROphet Score",
    ):
    """
    Plot multiple density histograms on the same axes.

    Parameters
    ----------
    trim : float in (0, 1], optional
        Fraction of data to KEEP from the center (symmetric tails removed).
        Example: trim=0.95 keeps the central 95% => removes 2.5% from each tail.
        trim=1.0 keeps everything (no trimming).
    """

    if not series_list:
        raise ValueError("series_list is empty.")

    n = len(series_list)

    # Auto styling for readability with many overlays
    if histtype is None:
        histtype = "step" if n > 3 else "stepfilled"
    if alpha is None:
        alpha = 1.0 if histtype == "step" else 0.35
    if linewidth is None:
        linewidth = 2.0 if histtype == "step" else 1.0

    # Validate trim (percent to keep)
    if trim is not None:
        if not (0 < trim <= 1):
            raise ValueError("trim must be in (0, 1]. Example: trim=0.95 keeps central 95%.")
        tail = (1.0 - float(trim)) / 2.0
    else:
        tail = None

    # Convert, clean, trim
    cleaned = []
    auto_names = []

    for i, s in enumerate(series_list, start=1):
        arr = np.asarray(s, dtype=float).ravel()

        if dropna:
            arr = arr[np.isfinite(arr)]

        if tail is not None and tail > 0 and arr.size > 0:
            lower = np.quantile(arr, tail)
            upper = np.quantile(arr, 1.0 - tail)
            arr = arr[(arr >= lower) & (arr <= upper)]

        cleaned.append(arr)

        series_name = getattr(s, "name", None)
        auto_names.append(series_name if series_name else f"Series {i}")

    if names is None:
        names = auto_names
    if len(names) != len(cleaned):
        raise ValueError("Length of 'names' must match length of 'series_list'.")

    # Axes setup
    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 5))
    else:
        fig = ax.figure

    # Shared range
    hist_range = None
    if common_range:
        mins = [a.min() for a in cleaned if a.size]
        maxs = [a.max() for a in cleaned if a.size]
        if mins and maxs:
            hist_range = (min(mins), max(maxs))

    # High-contrast colors + linestyles (helps when 4+ series)
    cmap = plt.get_cmap("tab10")
    colors = [cmap(i % 10) for i in range(n)]
    linestyles = ["-", "--", "-.", ":"]  # cycles nicely

    for idx, (arr, label) in enumerate(zip(cleaned, names)):
        if arr.size == 0:
            continue

        color = colors[idx]
        ls = linestyles[idx % len(linestyles)]

        ax.hist(
            arr,
            bins=bins,
            range=hist_range,
            density=True,
            histtype=histtype,
            alpha=alpha,
            linewidth=linewidth,
            linestyle=ls if histtype == "step" else "-",  # linestyle mainly matters for step
            color=color,
            label=label,
        )

    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.legend()
    ax.grid(True, alpha=0.25)

    return fig, ax



def scatter_with_regression(x, y,
                            xlabel="X",
                            ylabel="Y",
                            title="Scatter Plot with Linear Regression"):

    x = np.array(x)
    y = np.array(y)

    # Linear regression
    slope, intercept = np.polyfit(x, y, 1)

    # Regression line
    x_line = np.linspace(min(x), max(x), 100)
    y_line = slope * x_line + intercept

    # R²
    y_pred = slope * x + intercept
    r2 = r2_score(y, y_pred)

    # Plot
    plt.figure()
    plt.scatter(x, y, label="Data")

    # Regression line
    plt.plot(x_line, y_line,
             label=f"Regression: y = {slope:.3f}x + {intercept:.3f}")

    # X = Y reference line
    min_val = min(min(x), min(y))
    max_val = max(max(x), max(y))
    plt.plot([min_val, max_val], [min_val, max_val],
             linestyle=":",
             label="X = Y")

    # Equation text
    equation = f"R² = {r2:.3f}"
    plt.text(min(x), max(y), equation)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    plt.legend()

    plt.show()

    return slope, intercept, r2



def scatter_with_regression(x, y,
                            xlabel="X",
                            ylabel="Y",
                            title="Scatter Plot with Linear Regression"):

    x = np.array(x)
    y = np.array(y)

    # Linear regression
    slope, intercept = np.polyfit(x, y, 1)

    # Regression line
    x_line = np.linspace(min(x), max(x), 100)
    y_line = slope * x_line + intercept

    # R²
    y_pred = slope * x + intercept
    r2 = r2_score(y, y_pred)

    # Plot
    plt.figure()
    plt.scatter(x, y, label="Data")

    # Regression line
    plt.plot(x_line, y_line,
             label=f"Regression: y = {slope:.3f}x + {intercept:.3f}")

    # X = Y reference line
    min_val = min(min(x), min(y))
    max_val = max(max(x), max(y))
    plt.plot([min_val, max_val], [min_val, max_val],
             linestyle=":",
             label="X = Y")

    # Equation text
    equation = f"R² = {r2:.3f}"
    plt.text(min(x), max(y), equation)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    plt.legend()

    plt.show()

    return slope, intercept, r2




def plot_confusion_matrix(conf, title="Confusion Matrix", ylabel="Reference", xlabel="Test", save_path=None):
    """
    Plot a presentation-ready confusion matrix heatmap.

    Parameters
    ----------
    conf : pandas.DataFrame
        2x2 confusion matrix from pd.crosstab
    title : str
        Plot title
    save_path : str or None
        If provided, saves the figure to this path

    Returns
    -------
    metrics : dict
        Dictionary with TP, TN, FP, FN and accuracy
    """

    # Ensure full 2x2 structure
    conf = conf.reindex(index=[False, True], columns=[False, True], fill_value=0)

    TN = conf.loc[False, False]
    FP = conf.loc[False, True]
    FN = conf.loc[True, False]
    TP = conf.loc[True, True]

    total = conf.values.sum()
    accuracy = (TP + TN) / total if total > 0 else 0

    # plt.figure(figsize=(4,4))
    fig, ax = plt.subplots(figsize=(4, 4))

    sns.heatmap(
        conf,
        annot=True,
        fmt="d",
        cmap="Blues",
        cbar=False,
        ax=ax
    )

    ax.set_title(f"{title}\nAgreement = {accuracy:.2%} (n={total})")
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)

    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300)

    plt.show()

    return {
        "TP": TP,
        "TN": TN,
        "FP": FP,
        "FN": FN,
        "accuracy": accuracy,
        "fig": fig
    }


def transform_columns(df, transform_dict):
    """
    Transform specified columns using the formula:
    new_value = (e^b) * (old_value^a)
    where a = med_e/med_s and b = mad_e - (med_e/med_s)*mad_s
    
    Parameters:
    df: DataFrame to transform
    transform_dict: Dictionary with column names as keys and stats as values
    
    Returns:
    DataFrame with transformed columns
    """
    df_transformed = df.copy()
    
    for col_name, stats in transform_dict.items():
        if col_name in df_transformed.columns:
            # Calculate transformation parameters
            a = stats['med_edta'] / stats['med_streck']
            b = stats['mad_edta'] - (stats['med_edta'] / stats['med_streck']) * stats['mad_streck']
            
            # Apply transformation: new_value = (e^b) * (old_value^a)
            df_transformed[col_name] = (np.exp(b)) * (df_transformed[col_name] ** a)
    
    return df_transformed
