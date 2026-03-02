"""
module docstring

"""

import os
import pandas as pd
import somadata


def read_adat_file(file_path):
    """
    Reads an ADAT file and returns 2 DataFrames.
    The first DataFrame contains the plate data, 
    and the second DataFrame contains the plate targets.
    Args:
        file_path (str): The path to the ADAT file.
    Returns:
        tuple: A tuple containing two DataFrames (plate_data, plate_targets).
    """
    plate_adat = somadata.read_adat(file_path)
    plate_targets_df = plate_adat.columns.to_frame().reset_index(drop=True)
    df = plate_adat.reset_index()
    df.columns = ['_'.join(map(str, col)).strip('_') if isinstance(col, tuple) else col
              for col in df.columns]
    df.columns = [column.split('__')[0] for column in df.columns]
    return pd.DataFrame(df), plate_targets_df

def get_plate_data(adat_file_path):
    """
    Reads an ADAT file and returns a DataFrame with plate data.
    
    Parameters:
    adat_file (str): Path to the ADAT file.
    
    Returns:
    pd.DataFrame: DataFrame containing plate data.
    """
    if not os.path.exists(adat_file_path):
        raise FileNotFoundError(f"The file {adat_file_path} does not exist.")

    # Read the ADAT file using the custom method
    plate_data, general_data = read_adat_file(adat_file_path)

    return plate_data, general_data


def get_adat_files(file_paths):
    """    Reads ADAT files and returns a list of DataFrames 
           with sample data and a DataFrame with protein data.
    """
    # df_list = []
    df_dict = {}
    for file in file_paths:
        sample_df, protein_df = read_adat_file(file)
        df_dict[file] = sample_df
        return df_dict, protein_df
    
# def add_measure_id(df: pd.DataFrame):
#     """adds the MeasureId column to a dataframe

#     Args:
#         df (pd.DataFrame): a proteomics+metadata dataframe
#     Returns:
#         df (pd.DataFrame): the same df but with the extra column
#         success: Boolean
#     """
#     if "MeasureId" in df.columns:
#         return df, True
    
#     if 'PlatePosition' in df.columns and "PlateId" in df.columns:
#         df['PlateId'] = df['PlateId'].str.replace(r'_ABCD$', '', regex=True)
#         df['MeasureId'] = df['PlateId'] + '_' + df['PlatePosition']
#         return df, True
#     return df, False

