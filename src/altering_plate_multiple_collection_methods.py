"""_summary_

    Returns:
    _type_: _description_
"""

import os
import json
import time
import shutil
import pandas as pd
from src.bridging import *
from src.bridging_coefficients import *

import warnings
warnings.filterwarnings('ignore')

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


def transform_feature_df_MAD_bridging(df1, params):
    df = df1.copy()
    probe_to_aptamer_dict = match_probe_aptamer(df)
    original_columns = df.columns
    df['aptamer'] = df['ProbeName'].map(probe_to_aptamer_dict).fillna('other')
    df = aggregate(df)
    aptamers = list(probe_to_aptamer_dict.values())
    df = transform_mad_log(df, params, aptamers)
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
    # return reconstituted_texts


def test_aggregation_only(text_path, workbook_path, params_path, output_path, sample_type='Streck'): #, streck_wells_list=[]):
    """
    This method uses the previous ones to alter the scanner files

    Args:
        text_path (str): path to the original scanner files
        workbook_path (str): path to the workbook
        streck_conversion_coefficients_path (str): path to the bridging factors file
        streck_wells_list (list, optional): List of wells to be bridged (currently from Streck to EDTA). Defaults to [].
    """

    workbook_df = read_workbook(workbook_path)
    with open(params_path, 'r') as fp:
        params = json.load(fp)
          
    # params = pd.read_csv(params_path)
    # with open(params_path, 'r') as fp:
    #     streck_conversion_coefficients = json.load(fp)
    STRECK_VALUE = sample_type

    os.makedirs(output_path, exist_ok=True)
    shutil.copy(workbook_path, output_path)

    
    samples_to_alter  = workbook_df

    # streck_workbook_df  = workbook_df[workbook_df['well'].isin(streck_wells_list)]
    files = [
        f for f in os.listdir(text_path)
        if os.path.isfile(os.path.join(text_path, f))
    ]
    text_files = [file for file in files if file.endswith('.txt')]
    samples_to_alter['filename'] = samples_to_alter.apply(match_file, axis=1, args=(text_files,))
    streck_files = samples_to_alter['filename'].to_list()
    reconstituted_texts = {}
    for file in text_files:
        # output_file_path = os.path.join(output_path, file)
        reconstituted_text = alter_scanner_file_aggregation_only(text_path, file, output_path, file in streck_files, params)
        reconstituted_texts[file] = reconstituted_text

    reconstituted_texts["streck_files"] = streck_files 
    print('Done!')
    return reconstituted_texts

def alter_scanner_file_aggregation_only(input_path, file, output_path, is_streck_file, conversion_coefficients):
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
                df = transform_feature_df_aggregation_only(df, conversion_coefficients)
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

def transform_feature_df_aggregation_only(df1, params):
    df = df1.copy()
    probe_to_aptamer_dict = match_probe_aptamer(df)
    original_columns = df.columns
    df['aptamer'] = df['ProbeName'].map(probe_to_aptamer_dict).fillna('other')
    df = aggregate(df)
    # aptamers = list(probe_to_aptamer_dict.values())
    # df = transform_mad_log(df, params, aptamers)
    # df['conversion_coefficient'] = df['aptamer'].map(params)
    # if file in streck_files:
    # df['gProcessedSignal'] = df['gProcessedSignal'] * df['conversion_coefficient']
    df = df[original_columns]
    return df


if __name__ == '__main__':
    print(os.getcwd())
    
    # i = 4    # i = 1
    
    
    for i in range(1, 7):
        start = time.time()
            
        # i = 2
        if i == 3:
            alter_scanner_files_MAD_bridging(f'data/plates/OH2025_05{i}',
                    f'data/plates/OH2025_05{i}/OH2025_05{i} Workbook.xlsx',
                    
                    {'Streck': 'data/conversion_coefficients/MAD_median_ds_coefficients_plates_OH2025_051_to_OH_2026_003_04022026.json',
                    'nds': 'data/conversion_coefficients/MAD_median_nds_coefficients_plates_OH2025_054_OH2025_055_OH_2026_001_OH_2026_002_04022026.json'},
                    f'data/plates/debugging/OH2025_05{i}_only_streck_bridging',
                    
        )
            # alter_scanner_files_MAD_bridging(f'data/plates/multiple_types_experiment/OH2025_05{i}_altered_workbook',
            #                     f'data/plates/multiple_types_experiment/OH2025_05{i}_altered_workbook/OH2025_05{i} Workbook.xlsx',
            #                     {'Streck': 'data/conversion_coefficients/MAD_median_streck_coefficients.json',
            #                     'nds': 'data/conversion_coefficients/nds_params.json'},
            #                     f'data/plates/multiple_types_experiment/results/testing_transformation/OH2025_05{i}_for testing',
                                
            # )
            print(f"Time for plate OH2025_05{i}: {time.time() - start} seconds")
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
