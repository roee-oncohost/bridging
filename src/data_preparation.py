import re
import pandas as pd
from pathlib import Path
from src.adat_handling import read_adat_file
from src.dataframe_transformation import add_measure_id

def concat_adats(path, pattern):
    folder = Path(path)
    dfs = []

    for file in folder.rglob("*.adat"):
        if re.search(pattern, file.name):
            df, _ = read_adat_file(str(file))
            dfs.append(df)

    if not dfs:
        return pd.DataFrame()

    return pd.concat(dfs, ignore_index=True)

def split_ds_nds(df, ds_rfu_df, nds_rfu_df):

    ds_df = df[df["MeasureId"].isin(ds_rfu_df["MeasureId"])]
    nds_df = df[df["MeasureId"].isin(nds_rfu_df["MeasureId"])]

    return ds_df, nds_df


def prepare_adat_datasets(adat_path, ds_rfu_df, nds_rfu_df, pattern):

    df = concat_adats(adat_path, pattern)
    df = add_measure_id(df)

    ds_df, nds_df = split_ds_nds(df, ds_rfu_df, nds_rfu_df)

    return ds_df, nds_df


def prepare_base_rfu_datasets(base_rfu_path, info_dfs):
    base_rfu_df = pd.read_parquet(base_rfu_path).reset_index(drop=False)
    ds_info = info_dfs["DS"]
    nds_info = info_dfs["nDS"]
    edta_info = info_dfs["EDTA"]
    base_rfu_df = add_measure_id(base_rfu_df)
    ds_rfu_df = base_rfu_df[base_rfu_df["MeasureId"].isin(ds_info["MeasureId"])]
    nds_rfu_df = base_rfu_df[base_rfu_df["MeasureId"].isin(nds_info["MeasureId"])]
    edta_rfu_df = base_rfu_df[base_rfu_df["MeasureId"].isin(edta_info["MeasureId"])]
    return ds_rfu_df, nds_rfu_df, edta_rfu_df
    