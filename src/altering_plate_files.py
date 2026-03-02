"""_summary_

    Returns:
    _type_: _description_
"""

import os
import json
import shutil
import pandas as pd
from bridging import *

WORKBOOK_SAMPLE_MATRIX = "Sample Matrix"

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



# def read_workbook(file_path,
#                 tab_name=None,
#                 columns=['sample', 'well', 'clogged', 'low volume', 'name', 'lot number', 
#                 'sample type', 'slide', 'initial_chamber', 'subarray', "pdf_subarray", 'c2 aspiration',
#                  'leak', 'sample notes', 'assay notes']):
#     """
#     Read a workbook excel file (tab 4)

#     Args:
#         file_path (str): path to excel file
#         tab_name (str, optional): The tab to be read. Defaults to None.
#         columns (list, optional): the columns in the tab. Defaults to ['sample', 'well', 'clogged', 'low volume', 'name', 'lot number', 'sample type', 'slide', 'initial_chamber', 'subarray', "pdf_subarray", 'c2 aspiration', 'leak', 'sample notes', 'assay notes'].

#     Returns:
#         pd.DataFrame: a dataframe of the workbook tab
#     """
#     if not tab_name:
#         xls = pd.ExcelFile(file_path)
#         tab_name = xls.sheet_names[3]
#     workbook_df = pd.read_excel(file_path, sheet_name=tab_name, skiprows=7)
#     workbook_df = workbook_df.iloc[:, :len(columns)]
#     workbook_df.columns = columns
#     workbook_df['slide'] = workbook_df['slide'].ffill()
#     return workbook_df

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
    df_data_tab4 = df_sheet4[df_sheet4['NAME'].notna()]
    
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

def transform_feature_df_MAD_bridging(df, params):
    probe_to_aptamer_dict = match_probe_aptamer(df)
    original_columns = df.columns
    df['aptamer'] = df['ProbeName'].map(probe_to_aptamer_dict).fillna('other')
    aptamers = df['aptamer']
    transform_mad_log(df, params, aptamers)
    df['conversion_coefficient'] = df['aptamer'].map(params)
    # if file in streck_files:
    df['gProcessedSignal'] = df['gProcessedSignal'] * df['conversion_coefficient']
    df = df[original_columns]
    return df


def aggregation_function(df):
    
    return round(s.mean(), 1)

def switch_to_aggregate(df1):
    df = df1.copy()
    for aptamer in df['aptamer'].unique():
        if aptamer[0].isdigit():
            idx = df[df['aptamer']==aptamer].index
            df.loc[idx, 'gProcessedSignal'] = aggregation_function(
                df.loc[idx, 'gProcessedSignal'])
    return df


def transform_feature_df_normalization_bridging(df, a_b_dict):
    
    a_b_dict['other'] = {'a': 1,
                         'b': 0}
    original_columns = df.columns
    probe_to_aptamer_dict = match_probe_aptamer(df)

    # transform_streck_df(streck_df, a_b_dict, aptamers)
    df['aptamer'] = df['ProbeName'].map(probe_to_aptamer_dict).fillna('other')
    df = switch_to_aggregate(df)
    a_dict = {key: value['a'] for key, value in a_b_dict.items()}
    b_dict = {key: value['b'] for key, value in a_b_dict.items()}

    
    df['A'] = df['aptamer'].map(a_dict)
    df['B'] = df['aptamer'].map(b_dict)

    # if file in streck_files:
    for i, row in df.iterrows():
        row['gProcessedSignal'] = normalize_value(row['gProcessedSignal'], row['A'], row['B'])
    df = df[original_columns]
    return df

def transform_feature_df_mad_bridging(df, params):
    

    original_columns = df.columns
    probe_to_aptamer_dict = match_probe_aptamer(df)

    # transform_streck_df(streck_df, a_b_dict, aptamers)
    df['aptamer'] = df['ProbeName'].map(probe_to_aptamer_dict).fillna('other')
    df = switch_to_aggregate(df)
    a_dict = {key: value['a'] for key, value in params.items()}
    b_dict = {key: value['b'] for key, value in params.items()}

    
    df['A'] = df['aptamer'].map(a_dict)
    df['B'] = df['aptamer'].map(b_dict)

    # if file in streck_files:
    for i, row in df.iterrows():
        row['gProcessedSignal'] = normalize_value(row['gProcessedSignal'], row['A'], row['B'])
    df = df[original_columns]
    return df



def alter_scanner_file_MAD_bridging(input_path, file, output_path, is_streck_file, conversion_coefficients):
    file_path = os.path.join(input_path, file)
    # if file_path.endswith('.txt'):
    text = read_text_file(file_path)
    if is_streck_file:
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
                print(f"File: {file}")
                print(f"gProcessedSignal dtype: {df['gProcessedSignal'].dtype}")
                print(f"Sample values: {df['gProcessedSignal'].head()}")
                print(f"Any NaN values: {df['gProcessedSignal'].isna().any()}")
                # original_columns = df.columns
            dataframes_dict[header_name] = df
        reconstituted_text = dataframes_to_text(dataframes_dict, type_mappings)
        
    else:
        reconstituted_text = text
    output_file_path = os.path.join(output_path, file)
    make_text_file(output_file_path, reconstituted_text)
    return reconstituted_text


def alter_scanner_file_normalization_bridging(input_path, file, output_path, is_streck_file, conversion_coefficients):
    file_path = os.path.join(input_path, file)
    # if file_path.endswith('.txt'):
    text = read_text_file(file_path)
    if is_streck_file:
        print(f"Altered: {file}")
        sections = get_sections(text)
        dataframes_dict = {}
        type_mappings = {}
        for section in sections:
            lines = section.strip().split('\n')
            type_row, header_row, header_name, data_rows = analyze_lines(lines)
            df = create_df(type_row, header_row, data_rows, header_name, type_mappings)
            if header_name=='FEATURES':
                df = transform_feature_df_normalization_bridging(df, conversion_coefficients)
                # In alter_scanner_files, after creating the dataframe:
# if header_name == 'FEATURES':
                print(f"File: {file}")
                print(f"gProcessedSignal dtype: {df['gProcessedSignal'].dtype}")
                print(f"Sample values: {df['gProcessedSignal'].head()}")
                print(f"Any NaN values: {df['gProcessedSignal'].isna().any()}")
                # original_columns = df.columns
            dataframes_dict[header_name] = df
        reconstituted_text = dataframes_to_text(dataframes_dict, type_mappings)
        
    else:
        reconstituted_text = text
    output_file_path = os.path.join(output_path, file)
    make_text_file(output_file_path, reconstituted_text)
    return reconstituted_text


def alter_scanner_file_MAD_bridging(input_path, file, output_path, is_streck_file, params):
    file_path = os.path.join(input_path, file)
    # if file_path.endswith('.txt'):
    text = read_text_file(file_path)
    if is_streck_file:
        print(f"Altered: {file}")
        sections = get_sections(text)
        dataframes_dict = {}
        type_mappings = {}
        for section in sections:
            lines = section.strip().split('\n')
            type_row, header_row, header_name, data_rows = analyze_lines(lines)
            df = create_df(type_row, header_row, data_rows, header_name, type_mappings)
            if header_name=='FEATURES':
                df = transform_feature_df_mad_bridging(df, params)
                # In alter_scanner_files, after creating the dataframe:
# if header_name == 'FEATURES':
                print(f"File: {file}")
                print(f"gProcessedSignal dtype: {df['gProcessedSignal'].dtype}")
                print(f"Sample values: {df['gProcessedSignal'].head()}")
                print(f"Any NaN values: {df['gProcessedSignal'].isna().any()}")
                # original_columns = df.columns
            dataframes_dict[header_name] = df
        reconstituted_text = dataframes_to_text(dataframes_dict, type_mappings)
        
    else:
        reconstituted_text = text
    output_file_path = os.path.join(output_path, file)
    make_text_file(output_file_path, reconstituted_text)
    return reconstituted_text


def alter_scanner_files_MAD_bridging(text_path, workbook_path, params_path, output_path): #, streck_wells_list=[]):
    """
    This method uses the previous ones to alter the scanner files

    Args:
        text_path (str): path to the original scanner files
        workbook_path (str): path to the workbook
        streck_conversion_coefficients_path (str): path to the bridging factors file
        streck_wells_list (list, optional): List of wells to be bridged (currently from Streck to EDTA). Defaults to [].
    """

    workbook_df = read_workbook(workbook_path)
      
    params = pd.read_csv(params_path)
    # with open(params_path, 'r') as fp:
    #     streck_conversion_coefficients = json.load(fp)
    STRECK_VALUE = 'Streck'

    os.makedirs(output_path, exist_ok=True)
    shutil.copy(workbook_path, output_path)

    
    streck_workbook_df  = workbook_df[workbook_df[WORKBOOK_SAMPLE_MATRIX].isin([STRECK_VALUE])]

    # streck_workbook_df  = workbook_df[workbook_df['well'].isin(streck_wells_list)]
    files = [
        f for f in os.listdir(text_path)
        if os.path.isfile(os.path.join(text_path, f))
    ]
    text_files = [file for file in files if file.endswith('.txt')]
    streck_workbook_df['filename'] = streck_workbook_df.apply(match_file, axis=1, args=(text_files,))
    streck_files = streck_workbook_df['filename'].to_list()
    reconstituted_texts = {}
    for file in text_files:
        # output_file_path = os.path.join(output_path, file)
        reconstituted_text = alter_scanner_file_MAD_bridging(text_path, file, output_path, file in streck_files, params)
        reconstituted_texts[file] = reconstituted_text

    reconstituted_texts["streck_files"] = streck_files 
    print('Done!')
    return reconstituted_texts

def alter_scanner_files_normalization_bridging(text_path, workbook_path, streck_conversion_coefficients_path, output_path): #, streck_wells_list=[]):
    """
    This method uses the previous ones to alter the scanner files

    Args:
        text_path (str): path to the original scanner files
        workbook_path (str): path to the workbook
        streck_conversion_coefficients_path (str): path to the bridging factors file
        streck_wells_list (list, optional): List of wells to be bridged (currently from Streck to EDTA). Defaults to [].
    """

    workbook_df = read_workbook(workbook_path)
    

    with open(streck_conversion_coefficients_path, 'r') as fp:
        streck_conversion_coefficients = json.load(fp)
    STRECK_VALUE = 'Streck'

    os.makedirs(output_path, exist_ok=True)
    shutil.copy(workbook_path, output_path)
    
    streck_workbook_df  = workbook_df[workbook_df[WORKBOOK_SAMPLE_MATRIX].isin([STRECK_VALUE])]

    # streck_workbook_df  = workbook_df[workbook_df['well'].isin(streck_wells_list)]
    files = [
        f for f in os.listdir(text_path)
        if os.path.isfile(os.path.join(text_path, f))
    ]
    text_files = [file for file in files if file.endswith('.txt')]
    streck_workbook_df['filename'] = streck_workbook_df.apply(match_file, axis=1, args=(text_files,))
    streck_files = streck_workbook_df['filename'].to_list()
    reconstituted_texts = {}
    for file in text_files:
        # output_file_path = os.path.join(output_path, file)
        reconstituted_text = alter_scanner_file_normalization_bridging(text_path, file, output_path, file in streck_files, streck_conversion_coefficients)
        reconstituted_texts[file] = reconstituted_text

    reconstituted_texts["streck_files"] = streck_files 
    print('Done!')
    return reconstituted_texts


# def alter_scanner_files_MAD_bridging(text_path, workbook_path, streck_conversion_coefficients_path, output_path): #, streck_wells_list=[]):
#     """
#     This method uses the previous ones to alter the scanner files

#     Args:
#         text_path (str): path to the original scanner files
#         workbook_path (str): path to the workbook
#         streck_conversion_coefficients_path (str): path to the bridging factors file
#         streck_wells_list (list, optional): List of wells to be bridged (currently from Streck to EDTA). Defaults to [].
#     """

#     workbook_df = read_workbook(workbook_path)
    

#     with open(streck_conversion_coefficients_path, 'r') as fp:
#         streck_conversion_coefficients = json.load(fp)
#     STRECK_VALUE = 'Streck'

#     os.makedirs(output_path, exist_ok=True)
#     shutil.copy(workbook_path, output_path)
    
#     streck_workbook_df  = workbook_df[workbook_df[WORKBOOK_SAMPLE_MATRIX].isin([STRECK_VALUE])]

#     # streck_workbook_df  = workbook_df[workbook_df['well'].isin(streck_wells_list)]
#     files = [
#         f for f in os.listdir(text_path)
#         if os.path.isfile(os.path.join(text_path, f))
#     ]
#     text_files = [file for file in files if file.endswith('.txt')]
#     streck_workbook_df['filename'] = streck_workbook_df.apply(match_file, axis=1, args=(text_files,))
#     streck_files = streck_workbook_df['filename'].to_list()
#     reconstituted_texts = {}
#     for file in text_files:
#         # output_file_path = os.path.join(output_path, file)
#         reconstituted_text = alter_scanner_file_MAD_bridging(text_path, file, output_path, file in streck_files, streck_conversion_coefficients)
#         reconstituted_texts[file] = reconstituted_text

#     reconstituted_texts["streck_files"] = streck_files 
#     print('Done!')
#     return reconstituted_texts




if __name__ == '__main__':
    print(os.getcwd())
    # for i in range(1, 6):
    i = 1
    alter_scanner_files_MAD_bridging(f'data/plates/OH2025_05{i}',
                        f'data/plates/OH2025_05{i}/OH2025_05{i} Workbook.xlsx',
                        'data/conversion_coefficients/MAD_median_streck_coefficients.csv',
                        f'data/plates/normalization_bridging/OH2025_05{i}_normalized',
                        # ['A6', 'B4', 'C1', 'C5', 'D3', 'E5', 'F1', 'H2', 'H4']
    )
    # i = 1
        # alter_scanner_files_normalization_bridging(f'data/plates/OH2025_05{i}',
        #                 f'data/plates/OH2025_05{i}/OH2025_05{i} Workbook.xlsx',
        #                 'data/conversion_coefficients/prehyb_normalization_edta_streck.json',
        #                 f'data/plates/normalization_bridging/OH2025_05{i}_normalized',
        #                 # ['A6', 'B4', 'C1', 'C5', 'D3', 'E5', 'F1', 'H2', 'H4']
        # )
    
    # alter_scanner_files_normalization_bridging(f'data/plates/OH2025_05{i}',
    #                 f'data/plates/OH2025_05{i}/OH2025_05{i} Workbook.xlsx',
    #                 'data/conversion_coefficients/prehyb_normalization_edta_streck.json',
    #                 f'data/plates/normalization_bridging/OH2025_05{i}_testing_value_processing',
    #                 # ['A6', 'B4', 'C1', 'C5', 'D3', 'E5', 'F1', 'H2', 'H4']
    # )
    # read_workbook('./data/test_files/test Workbook.xlsx')
    # ['A1'], 'A3', 'B2', 'B5', 'C2', 'C3', 'C4', 'D4', 'E1',
#                         'E4', 'E5', 'F2', 'F3', 'G1', 'G2', 'H3', 'H4', 'H5'])
