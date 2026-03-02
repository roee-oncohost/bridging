import os
import json
import numpy as np
import time
import shutil
import pandas as pd
import warnings

warnings.filterwarnings('ignore')


WORKBOOK_SAMPLE_MATRIX = "Sample Matrix"


def compute_median_mad_stats(source_df, target_df):
    """
    Compute median and MAD (Median Absolute Deviation) for each numerical column
    in both source and target datasets.
    
    Parameters:
    -----------
    source_df : pandas.DataFrame
        Source dataset
    target_df : pandas.DataFrame
        Target dataset
    
    Returns:
    --------
    dict
        Dictionary where each key is a column name and each value is a dict with:
        - 'source_median': median of source column
        - 'source_MAD': MAD of source column
        - 'target_median': median of target column
        - 'target_MAD': MAD of target column
    """
    stats = {}
    
    # Get numerical columns that exist in both datasets
    numerical_cols = source_df.select_dtypes(include=[np.number]).columns
    common_cols = [col for col in numerical_cols if col in target_df.columns]
    
    for col in common_cols:
        # Compute median
        source_median = source_df[col].median()
        target_median = target_df[col].median()
        
        # Compute MAD (Median Absolute Deviation)
        source_mad = (source_df[col] - source_median).abs().median()
        target_mad = (target_df[col] - target_median).abs().median()
        
        stats[col] = {
            'source_median': source_median,
            'source_MAD': source_mad,
            'target_median': target_median,
            'target_MAD': target_mad
        }
        stats["other"] = {
            'source_median': 1,
            'source_MAD': 1,
            'target_median': 1,
            'target_MAD': 1
        }
    
    return stats



def map_to_target_distribution(source_df, stats):
    """
    Map source dataset to target distribution using median and MAD statistics.
    
    For each column, the transformation is:
    1. Standardize using source median and MAD: (x - source_median) / source_MAD
    2. Rescale using target median and MAD: standardized * target_MAD + target_median
    
    This ensures the transformed data has the same median and MAD as the target.
    
    Parameters:
    -----------
    source_df : pandas.DataFrame
        Source dataset to transform
    stats : dict
        Statistics dictionary from compute_median_mad_stats()
    
    Returns:
    --------
    pandas.DataFrame
        Transformed dataset with target distribution characteristics
    """
    # Create a copy to avoid modifying the original
    transformed_df = source_df.copy()
    
    for col, col_stats in stats.items():
        if col in transformed_df.columns:
            source_median = col_stats['source_median']
            source_mad = col_stats['source_MAD']
            target_median = col_stats['target_median']
            target_mad = col_stats['target_MAD']
            
            # Avoid division by zero
            if source_mad == 0:
                # If source MAD is 0, all values are the same, just set to target median
                transformed_df[col] = target_median
            else:
                # Standardize using source stats, then rescale using target stats
                standardized = (transformed_df[col] - source_median) / source_mad
                transformed_df[col] = standardized * target_mad + target_median
    
    return transformed_df

def map_plate_to_distribution(df, normalization_dict):
    df1 = df.copy()
    original_columns = df.columns
    norm_df = pd.DataFrame.from_dict(normalization_dict, orient='index')

# Merge with your main dataframe
    df1 = df1.merge(norm_df, left_on='aptamer', right_index=True, how='left')

    # Vectorized calculation
    df1['gProcessedSignal'] = (
        (df1['gProcessedSignal'] - df1['source_median']) * 
        df1['target_MAD'] / df1['source_MAD'] + 
        df1['target_median']
    )

    # Drop the temporary columns if you don't need them
    df1 = df1[original_columns]
    return df1
    
def map_plate_to_distribution_log(df, normalization_dict):
    df1 = df.copy()
    original_columns = df.columns
    norm_df = pd.DataFrame.from_dict(normalization_dict, orient='index')

# Merge with your main dataframe
    df1 = df1.merge(norm_df, left_on='aptamer', right_index=True, how='left')

    # Vectorized calculation
    df1['gProcessedSignalLog'] = (
        (np.log(df1['gProcessedSignal']) - df1['source_median']) * 
        df1['target_MAD'] / df1['source_MAD'] + 
        df1['target_median']
    )
    
    df1['gProcessedSignal'] = np.exp(df1['gProcessedSignalLog'])

    # Drop the temporary columns if you don't need them
    df1 = df1[original_columns]
    return df1








def replace_with_conditional_mean(df):
    """
    For each aptamer (except "other"), replace all gProcessedSignal values with the mean
    of gProcessedSignal values where gIsFeatPopnOL < 0.5 for that aptamer.
    
    This vectorized version is much faster than looping, especially for large datasets.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with columns: 'aptamer', 'gProcessedSignal', 'gIsFeatPopnOL'
    
    Returns:
    --------
    pandas.DataFrame
        DataFrame with updated gProcessedSignal values
    """
    # Create a copy to avoid modifying the original
    df_copy = df.copy()
    
    # Create a mask for rows we want to process (not "other" and gIsFeatPopnOL < 0.5)
    process_mask = (df_copy['aptamer'] != 'other') & (df_copy['gIsFeatPopnOL'] < 0.5)
    
    # Calculate mean for each aptamer using only rows where gIsFeatPopnOL < 0.5
    # This creates a Series with aptamer as index and mean as value
    aptamer_means = df_copy[process_mask].groupby('aptamer')['gProcessedSignal'].mean()
    
    # Map the means back to all rows (not just those with gIsFeatPopnOL < 0.5)
    # For rows with aptamer="other" or aptamers not in aptamer_means, this will be NaN
    df_copy.loc[df_copy['aptamer'] != 'other', 'gProcessedSignal'] = \
        df_copy.loc[df_copy['aptamer'] != 'other', 'aptamer'].map(aptamer_means)
    
    # If an aptamer has no rows with gIsFeatPopnOL < 0.5, restore original values
    # (the map will have set them to NaN)
    nan_mask = df_copy['gProcessedSignal'].isna()
    if nan_mask.any():
        df_copy.loc[nan_mask, 'gProcessedSignal'] = df.loc[nan_mask, 'gProcessedSignal']
    
    return df_copy


def read_text_file(path_to_text):
    """
    Read a scanner text file 

    Args:
        path_to_text (str): path to the scanner text file

    Returns:
        str: the text
    """
    with open(path_to_text, 'r') as f:
        text_data = f.read()
    return text_data


def get_sections(text):
    """
    Each scanner file includes 3 dataframes (sections).
    We split the text into 3 texts, each depicting a single dataframe 

    Args:
        text (str): scanner text file

    Returns:
        list: list of text strings
    """
    sections = text.strip().split('*\n')
    return sections



def analyze_lines(lines):
    """
    Every scanner file is parsed into 3 sections (lists of text lines).
    These may include a line for types, header line, header name line,
    and lines of actual data.
    This method parses the given section line by line

    Args:
        lines (list): a list of lines which form a section (dataframe) of the scanner file

    Returns:
        tuple: a tuple of the different kinds of lines
    """
    type_row = None
    header_row = None
    header_name = ''
    data_rows = []
    for line in lines:
        if line.startswith('TYPE\t'):
            type_row = line.split('\t')[1:] 
        elif line.startswith('FEPARAMS\t') or line.startswith('STATS\t') or line.startswith('FEATURES\t'):
            parts = line.split('\t')
            header_name = parts[0]
            header_row = parts[1:]
        elif line.startswith('DATA\t'):
            data_rows.append(line.split('\t')[1:])
    return type_row, header_row, header_name, data_rows


def create_df(type_row, header_row, data_rows, header_name, type_mappings):
    """
    Convert parsed rows to a dataframe
    
    Args: the aforementioned parsed lines, as well as type_mappings (which are later used to rebuild the text file)
        
    Returns:
        pandas DataFrame based on the parsed lines
    """
    
    type_mappings[header_name] = type_row
    df = pd.DataFrame(data_rows, columns=header_row)
    for col, dtype in zip(header_row, type_row):
        if col in df.columns:
            try:
                if dtype == 'integer':
                    df[col] = pd.to_numeric(df[col], errors='coerce').astype('Int64')
                elif dtype == 'float':
                    df[col] = pd.to_numeric(df[col], errors='coerce')
                elif dtype == 'boolean':
                    df[col] = df[col].map({'0': False, '1': True, 'False': False, 'True': True})
                # 'text' stays as string
            except Exception as e:
                print(f"Warning: Could not convert column {col} to {dtype}: {e}")

    
    return df 


def dataframes_to_text(dataframes, type_mappings):
    """
    Convert dataframes back to the original text format.
    
    Args:
        dataframes: Dictionary of dataframes (e.g., {'FEPARAMS': df1, 'STATS': df2})
        type_mappings: Dictionary mapping dataframe names to their type lists
                      (e.g., {'FEPARAMS': ['text', 'integer', ...], ...})
    
    Returns:
        String in the original text format
    """
    sections = []
    
    for name, df in dataframes.items():
        lines = []
        
        # Get the type mapping for this dataframe
        types = type_mappings.get(name, ['text'] * len(df.columns))
        
        # TYPE row
        type_line = 'TYPE\t' + '\t'.join(types)
        lines.append(type_line)
        
        # Header row
        header_line = name + '\t' + '\t'.join(df.columns)
        lines.append(header_line)
        
        # DATA rows
        for _, row in df.iterrows():
            # Convert values back to strings, handling special cases
            values = []
            for val in row:
                if pd.isna(val):
                    values.append('')
                elif isinstance(val, bool):
                    values.append('1' if val else '0')
                else:
                    values.append(str(val))
            
            data_line = 'DATA\t' + '\t'.join(values)
            lines.append(data_line)
        
        sections.append('\n'.join(lines))
    
    # Join sections with asterisk separator
    return '\n*\n'.join(sections) + '\n'


def make_text_file(file_name, text):
    """
    write a text file from the reconstituted text

    Args:
        file_name (str): name of file to be written
        text (str): text to be written to file
    """
    with open(file_name, 'w') as f:
        f.write(text)


def read_workbook(workbook_path) -> pd.DataFrame:
    ### Create DF from Sheet-3 
    sheet_name_tab3 = 'Tab 3 - Assayed Sample List'
    col_smb_tab3 = 'A, F, J'  # v7    
    col_names_tab3 = ['Sample Number', 'SampleId', WORKBOOK_SAMPLE_MATRIX]
    
    # Read the excel with convert the value of the column SampleID to string
    df_sheet3 = pd.read_excel(workbook_path, sheet_name=sheet_name_tab3, 
                              usecols=col_smb_tab3, dtype={'SampleId': str})
    
    # Skip IF One is EMPTY
    df_data_tab3 = df_sheet3[col_names_tab3].dropna()
    
    ### Create DF from Sheet-4 
    sheet_name_tab4 = 'Tab 4 - Plate Map'
    col_smb_tab4 = 'A, E, H, K'
    # Read the excel with convert the value of column NAME to string
    df_sheet4 = pd.read_excel(workbook_path, sheet_name=sheet_name_tab4, 
                              usecols=col_smb_tab4, skiprows=6, dtype={'NAME': str})
    df_sheet4['Slide #'] = df_sheet4['Slide #'].ffill()
    df_data_tab4 = df_sheet4[df_sheet4['NAME'].notna()]
    df_data_tab4 = df_data_tab4[df_data_tab4['NAME'].astype(str) != "0"]
    # Create a DICT of the Aliquots-Barcodes
    excel_df = pd.merge(df_data_tab3,
                    df_data_tab4,
                    how='inner',
                    left_on='SampleId', right_on='NAME').drop('NAME', axis=1)
    excel_df['slide'] = excel_df['Slide #'].ffill()
    excel_df.rename(columns={'PDF Subarray': 'pdf_subarray'}, inplace=True)
    return excel_df


def match_file(row, files):
    """
    Matches each well in the workbook to the corresponding scanner file 

    Args:
        row (pd.DataFrame row): a row in tab 4 of the workbook
        files (list): list of the scanner files in the folder

    Returns:
        str: name of the file corresponding to the row
    """
    for f in files:
        if str(row['slide']) in f and f.endswith(row['pdf_subarray'] + '.txt'):
            return f
    return None

def match_probe_aptamer(df):
    """
    matches the probe names in the "FEATURES" dataframe to the corresponding aptamers

    Args:
        df (pd.DataFrame): the "FEATURES" dataframe

    Returns:
        dict: a mapping from probe names to aptamers
    """
    probe_names = list(set(df['ProbeName'].to_list()))
    probe_names = [probe_name for probe_name in probe_names if probe_name.startswith('anti-')]
    probe_to_aptamer_dict = {probe_name: probe_name.split('anti-')[1].split('_')[0] for probe_name in probe_names}
    return probe_to_aptamer_dict


def transform_feature_df_MAD_bridging(df1, params):
    df = df1.copy()
    probe_to_aptamer_dict = match_probe_aptamer(df)
    original_columns = df.columns
    df['aptamer'] = df['ProbeName'].map(probe_to_aptamer_dict).fillna('other')
    df = replace_with_conditional_mean(df)
    # aptamers = list(probe_to_aptamer_dict.values())
    # df = map_to_target_distribution(df, params)
    # df = map_plate_to_distribution(df, params)
    df = map_plate_to_distribution_log(df, params)
    # df['conversion_coefficient'] = df['aptamer'].map(params)
    # if file in streck_files:
    # df['gProcessedSignal'] = df['gProcessedSignal'] * df['conversion_coefficient']
    df = df[original_columns]
    return df


def alter_scanner_file_MAD_bridging(input_path, file, output_path, conversion_coefficients, already_handled_files=[]):
    file_path = os.path.join(input_path, file)
    # if file_path.endswith('.txt'):
    text = read_text_file(file_path)
    # if is_type_file:
    print(f"Altered: {file}")
    sections = get_sections(text)
    dataframes_dict = {}
    type_mappings = {}
    for section in sections:
        lines = section.strip().split('\n')
        type_row, header_row, header_name, data_rows = analyze_lines(lines)
        df = create_df(type_row, header_row, data_rows, header_name, type_mappings)
        if header_name=='FEATURES':
            df = transform_feature_df_MAD_bridging(df, conversion_coefficients)
            # In alter_scanner_files, after creating the dataframe:
# if header_name == 'FEATURES':
            # print(f"File: {file}")
            # print(f"gProcessedSignal dtype: {df['gProcessedSignal'].dtype}")
            # print(f"Sample values: {df['gProcessedSignal'].head()}")
            # print(f"Any NaN values: {df['gProcessedSignal'].isna().any()}")
            # original_columns = df.columns
        dataframes_dict[header_name] = df
    reconstituted_text = dataframes_to_text(dataframes_dict, type_mappings)
# else:
    # reconstituted_text = text
    output_file_path = os.path.join(output_path, file)
    # if file not in already_handled_files:
    make_text_file(output_file_path, reconstituted_text)
    already_handled_files.append(file)
    return already_handled_files


def alter_scanner_files_MAD_bridging(text_path, workbook_path, params_paths, output_path): #, streck_wells_list=[]):
    """
    This method uses the previous ones to alter the scanner files

    Args:
        text_path (str): path to the original scanner files
        workbook_path (str): path to the workbook
        params_paths (dict): dictionary mapping sample types to their corresponding params file paths
        output_path (str): path to the output folder
        
    """
    

    workbook_df = read_workbook(workbook_path)
    
    os.makedirs(output_path, exist_ok=True)
    shutil.copy(workbook_path, output_path)


    files = [
        f for f in os.listdir(text_path)
        if os.path.isfile(os.path.join(text_path, f))
    ]
    text_files = [file for file in files if file.endswith('.txt')]
    already_handled_files = []
    for sample_type, params_path in params_paths.items():
        with open(params_path, 'r') as fp:
                params = json.load(fp)
        samples_to_alter  = workbook_df[workbook_df[WORKBOOK_SAMPLE_MATRIX].isin([sample_type])]

        samples_to_alter['filename'] = samples_to_alter.apply(match_file, axis=1, args=(text_files,))
        type_files = samples_to_alter['filename'].to_list()
        reconstituted_texts = {}
        for file in type_files:
            # output_file_path = os.path.join(output_path, file)
            already_handled_files.append(alter_scanner_file_MAD_bridging(text_path, file, output_path, params, already_handled_files))
            # reconstituted_texts[file] = reconstituted_text

        # reconstituted_texts["streck_files"] = type_files 
    files_to_copy = [text_file for text_file in text_files if text_file not in already_handled_files]
    for file in files_to_copy:
        shutil.copy(os.path.join(text_path, file), output_path)
    print('Done!')
    return True


if __name__ == '__main__':
    plates = [f'OH2025_05{i}' for i in range(1, 8)]
    plates += [f'OH2026_00{i}' for i in range(1, 4)]

    for plate in plates:
        start = time.time()
            # plate = f"OH2026_00{i}"
            
    # i = 2
    # if i == 3:
            # if i > -5:
        alter_scanner_files_MAD_bridging(f'projects/report/data/plates/original/original_with_altered_workbooks/{plate}',
                    f'projects/report/data/plates/original/original_with_altered_workbooks/{plate}/{plate} Workbook.xlsx',
                    
                    {'Streck': 'data/conversion_coefficients/ds_log_16022026.json',
                    'NDS': 'data/conversion_coefficients/nds_log_16022026.json'},
                    f'projects/report/data/plates/madmed_transformed_23022026/{plate}',
                    
        )
                # alter_scanner_files_MAD_bridging(f'data/plates/multiple_types_experiment/OH2025_05{i}_altered_workbook',
                #                     f'data/plates/multiple_types_experiment/OH2025_05{i}_altered_workbook/OH2025_05{i} Workbook.xlsx',
                #                     {'Streck': 'data/conversion_coefficients/MAD_median_streck_coefficients.json',
                #                     'nds': 'data/conversion_coefficients/nds_params.json'},
                #                     f'data/plates/multiple_types_experiment/results/testing_transformation/OH2025_05{i}_for testing',
                                
                # )
        print(f"Time for plate {plate}: {time.time() - start} seconds")

        


