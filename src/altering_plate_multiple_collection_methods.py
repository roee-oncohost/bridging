"""_summary_

    Returns:
    _type_: _description_
"""

import os
import json
import time
import shutil
from collections.abc import Mapping
from numbers import Real
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd
from src.bridging import aggregate, transform_mad_log, transform_median_only

WORKBOOK_SAMPLE_MATRIX = "Sample Matrix"

DEFAULT_HCE_APTAMERS = frozenset({
    '2171-12', '2178-55', '2194-91', '2229-54', '2249-25', '2273-34',
    '2288-7', '2305-52', '2312-13', '2359-65', '2430-52', '2513-7',
})

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
                elif isinstance(val, (bool, np.bool_)):
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
    missing_tab3_columns = sorted(set(col_names_tab3) - set(df_sheet3.columns))
    if missing_tab3_columns:
        raise ValueError(
            f"{sheet_name_tab3!r} is missing required columns: "
            f"{missing_tab3_columns}"
        )
    
    # Skip IF One is EMPTY
    df_data_tab3 = df_sheet3[col_names_tab3].dropna().copy()
    duplicate_sample_ids = df_data_tab3.loc[
        df_data_tab3['SampleId'].duplicated(keep=False), 'SampleId'
    ].astype(str).unique().tolist()
    if duplicate_sample_ids:
        raise ValueError(
            "Tab 3 contains duplicate SampleId values and cannot be joined "
            f"unambiguously: {duplicate_sample_ids[:10]}"
        )
    
    ### Create DF from Sheet-4 
    sheet_name_tab4 = 'Tab 4 - Plate Map'
    col_smb_tab4 = 'A, E, H, K'
    # Read the excel with convert the value of column NAME to string
    df_sheet4 = pd.read_excel(workbook_path, sheet_name=sheet_name_tab4, 
                              usecols=col_smb_tab4, skiprows=6, dtype={'NAME': str})
    required_tab4_columns = {'NAME', 'Slide #', 'PDF Subarray'}
    missing_tab4_columns = sorted(required_tab4_columns - set(df_sheet4.columns))
    if missing_tab4_columns:
        raise ValueError(
            f"{sheet_name_tab4!r} is missing required columns: "
            f"{missing_tab4_columns}"
        )
    df_data_tab4 = df_sheet4[df_sheet4['NAME'].notna()].copy()
    df_data_tab4['Slide #'] = df_data_tab4['Slide #'].ffill()
    relevant_tab4 = df_data_tab4[
        df_data_tab4['NAME'].isin(df_data_tab3['SampleId'])
    ].copy()
    duplicate_names = relevant_tab4.loc[
        relevant_tab4['NAME'].duplicated(keep=False), 'NAME'
    ].astype(str).unique().tolist()
    if duplicate_names:
        raise ValueError(
            "Tab 4 contains duplicate NAME values and cannot be joined "
            f"unambiguously: {duplicate_names[:10]}"
        )
    # Create a DICT of the Aliquots-Barcodes
    excel_df = pd.merge(df_data_tab3,
                    relevant_tab4,
                    how='left',
                    left_on='SampleId', right_on='NAME',
                    validate='one_to_one', indicator=True)
    missing_plate_map = excel_df[excel_df['_merge'] != 'both']
    if not missing_plate_map.empty:
        examples = missing_plate_map['SampleId'].astype(str).head(10).tolist()
        raise ValueError(
            f"Tab 4 is missing plate-map rows for {len(missing_plate_map)} "
            f"Tab 3 samples: {examples}"
        )
    excel_df = excel_df.drop(columns=['NAME', '_merge'])
    excel_df['slide'] = excel_df['Slide #']
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
    def scanner_token(value):
        if pd.isna(value):
            return ''
        if isinstance(value, Real) and not isinstance(value, (bool, np.bool_)):
            numeric_value = float(value)
            if np.isfinite(numeric_value) and numeric_value.is_integer():
                return str(int(numeric_value))
        return str(value).strip()

    slide = scanner_token(row['slide'])
    pdf_subarray = scanner_token(row['pdf_subarray'])
    if not slide or not pdf_subarray:
        return None
    expected_suffix = f"_{slide}_{pdf_subarray}.txt"
    matches = [
        filename for filename in files
        if filename.endswith(expected_suffix)
    ]
    if len(matches) > 1:
        raise ValueError(
            f"Workbook row matched multiple scanner files for slide={slide!r}, "
            f"pdf_subarray={pdf_subarray!r}: {matches}"
        )
    return matches[0] if matches else None

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


def _aptamers_to_transform(probe_to_aptamer_dict, params, excluded_aptamers, strict_parameters):
    """Return transformable aptamers after exclusions and parameter validation."""
    excluded = set(excluded_aptamers or ())
    candidates = sorted(set(probe_to_aptamer_dict.values()) - excluded)
    missing = [aptamer for aptamer in candidates if aptamer not in params]
    if missing and strict_parameters:
        preview = ', '.join(missing[:10])
        suffix = ' ...' if len(missing) > 10 else ''
        raise ValueError(
            f"Missing bridging parameters for {len(missing)} aptamers: "
            f"{preview}{suffix}"
        )
    return [aptamer for aptamer in candidates if aptamer in params]


def _validate_aggregated_signals(df, aptamers, *, require_positive):
    """Fail before a scanner transform can create NaN or infinite RFUs."""
    selected = df['aptamer'].isin(aptamers)
    if not selected.any():
        return
    signals = pd.to_numeric(df.loc[selected, 'gProcessedSignal'], errors='coerce')
    invalid = signals.isna() | ~np.isfinite(signals)
    if require_positive:
        invalid |= signals <= 0
    if invalid.any():
        bad_aptamers = sorted(
            set(df.loc[signals.index[invalid], 'aptamer'].astype(str))
        )
        requirement = 'finite and positive' if require_positive else 'finite'
        raise ValueError(
            f"Aggregated scanner signals must be {requirement}; affected "
            f"aptamers: {bad_aptamers[:10]}"
        )


def transform_feature_df_MAD_bridging(
    df1,
    params,
    excluded_aptamers=DEFAULT_HCE_APTAMERS,
    strict_parameters=True,
):
    df = df1.copy()
    probe_to_aptamer_dict = match_probe_aptamer(df)
    original_columns = df.columns
    df['aptamer'] = df['ProbeName'].map(probe_to_aptamer_dict).fillna('other')
    df = aggregate(df)
    aptamers_to_alter = _aptamers_to_transform(
        probe_to_aptamer_dict,
        params,
        excluded_aptamers,
        strict_parameters,
    )
    _validate_aggregated_signals(
        df,
        aptamers_to_alter,
        require_positive=True,
    )
    df = transform_mad_log(df, params, aptamers_to_alter)
    _validate_aggregated_signals(
        df,
        aptamers_to_alter,
        require_positive=True,
    )
    # df['conversion_coefficient'] = df['aptamer'].map(params)
    # if file in streck_files:
    # df['gProcessedSignal'] = df['gProcessedSignal'] * df['conversion_coefficient']
    df = df[original_columns]
    return df

def transform_feature_df_median_bridging(
    df1,
    params,
    excluded_aptamers=(),
    strict_parameters=True,
):
    df = df1.copy()
    probe_to_aptamer_dict = match_probe_aptamer(df)
    original_columns = df.columns
    df['aptamer'] = df['ProbeName'].map(probe_to_aptamer_dict).fillna('other')
    df = aggregate(df)
    aptamers = _aptamers_to_transform(
        probe_to_aptamer_dict,
        params,
        excluded_aptamers,
        strict_parameters,
    )
    _validate_aggregated_signals(df, aptamers, require_positive=False)
    df = transform_median_only(df, params, aptamers)
    _validate_aggregated_signals(df, aptamers, require_positive=False)
    # df['conversion_coefficient'] = df['aptamer'].map(params)
    # if file in streck_files:
    # df['gProcessedSignal'] = df['gProcessedSignal'] * df['conversion_coefficient']
    df = df[original_columns]
    return df



def alter_scanner_file_MAD_bridging(
    input_path,
    file,
    output_path,
    conversion_coefficients,
    already_handled_files=None,
    excluded_aptamers=DEFAULT_HCE_APTAMERS,
    strict_parameters=True,
):
    if not file:
        raise ValueError("A scanner filename is required")
    if already_handled_files is None:
        already_handled_files = []
    file_path = os.path.join(input_path, file)
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Scanner file does not exist: {file_path}")
    # if file_path.endswith('.txt'):
    text = read_text_file(file_path)
    # if is_type_file:
    print(f"Altered: {file}")
    sections = get_sections(text)
    dataframes_dict = {}
    type_mappings = {}
    feature_section_count = 0
    for section in sections:
        lines = section.strip().split('\n')
        type_row, header_row, header_name, data_rows = analyze_lines(lines)
        df = create_df(type_row, header_row, data_rows, header_name, type_mappings)
        if header_name=='FEATURES':
            feature_section_count += 1
            df = transform_feature_df_MAD_bridging(
                df,
                conversion_coefficients,
                excluded_aptamers=excluded_aptamers,
                strict_parameters=strict_parameters,
            )
            # In alter_scanner_files, after creating the dataframe:
# if header_name == 'FEATURES':
            # print(f"File: {file}")
            # print(f"gProcessedSignal dtype: {df['gProcessedSignal'].dtype}")
            # print(f"Sample values: {df['gProcessedSignal'].head()}")
            # print(f"Any NaN values: {df['gProcessedSignal'].isna().any()}")
            # original_columns = df.columns
        dataframes_dict[header_name] = df
    if feature_section_count != 1:
        raise ValueError(
            f"Scanner file {file!r} must contain exactly one FEATURES section; "
            f"found {feature_section_count}"
        )
    reconstituted_text = dataframes_to_text(dataframes_dict, type_mappings)
# else:
    # reconstituted_text = text
    output_file_path = os.path.join(output_path, file)
    # if file not in already_handled_files:
    os.makedirs(output_path, exist_ok=True)
    make_text_file(output_file_path, reconstituted_text)
    if file not in already_handled_files:
        already_handled_files.append(file)
    return already_handled_files


def alter_scanner_file_median_bridging(
    input_path,
    file,
    output_path,
    conversion_coefficients,
    already_handled_files=None,
    excluded_aptamers=(),
    strict_parameters=True,
):
    if not file:
        raise ValueError("A scanner filename is required")
    if already_handled_files is None:
        already_handled_files = []
    file_path = os.path.join(input_path, file)
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Scanner file does not exist: {file_path}")
    # if file_path.endswith('.txt'):
    text = read_text_file(file_path)
    # if is_type_file:
    print(f"Altered: {file}")
    sections = get_sections(text)
    dataframes_dict = {}
    type_mappings = {}
    feature_section_count = 0
    for section in sections:
        lines = section.strip().split('\n')
        type_row, header_row, header_name, data_rows = analyze_lines(lines)
        df = create_df(type_row, header_row, data_rows, header_name, type_mappings)
        if header_name=='FEATURES':
            feature_section_count += 1
            df = transform_feature_df_median_bridging(
                df,
                conversion_coefficients,
                excluded_aptamers=excluded_aptamers,
                strict_parameters=strict_parameters,
            )
            # In alter_scanner_files, after creating the dataframe:
# if header_name == 'FEATURES':
            # print(f"File: {file}")
            # print(f"gProcessedSignal dtype: {df['gProcessedSignal'].dtype}")
            # print(f"Sample values: {df['gProcessedSignal'].head()}")
            # print(f"Any NaN values: {df['gProcessedSignal'].isna().any()}")
            # original_columns = df.columns
        dataframes_dict[header_name] = df
    if feature_section_count != 1:
        raise ValueError(
            f"Scanner file {file!r} must contain exactly one FEATURES section; "
            f"found {feature_section_count}"
        )
    reconstituted_text = dataframes_to_text(dataframes_dict, type_mappings)
# else:
    # reconstituted_text = text
    output_file_path = os.path.join(output_path, file)
    # if file not in already_handled_files:
    os.makedirs(output_path, exist_ok=True)
    make_text_file(output_file_path, reconstituted_text)
    if file not in already_handled_files:
        already_handled_files.append(file)
    return already_handled_files



def _manifest_factors_to_legacy(data, expected_methods):
    from src.bridging_factors import validate_bridging_factors

    validate_bridging_factors(data)
    method = data.get('method')
    if method not in expected_methods:
        raise ValueError(
            f"Factor manifest method {method!r} is incompatible with this scanner "
            f"transform; expected one of {sorted(expected_methods)}"
        )
    factors = data.get('factors')
    if not isinstance(factors, Mapping) or not factors:
        raise ValueError("Factor manifest must contain a non-empty 'factors' mapping")

    if method == 'mad_log':
        required = {'source_median', 'source_mad', 'target_median', 'target_mad'}
        converted = {}
        for aptamer, factor in factors.items():
            missing = required - set(factor)
            if missing:
                raise ValueError(
                    f"Factor manifest aptamer {aptamer!r} is missing {sorted(missing)}"
                )
            converted[aptamer] = {
                'med_streck': factor['source_median'],
                'mad_streck': factor['source_mad'],
                'med_edta': factor['target_median'],
                'mad_edta': factor['target_mad'],
            }
        return converted

    converted = {}
    for aptamer, factor in factors.items():
        if not isinstance(factor, Mapping) or 'factor' not in factor:
            raise ValueError(
                f"Factor manifest aptamer {aptamer!r} is missing its multiplier"
            )
        # transform_median_only multiplies med_edta / med_streck. Encoding the
        # fitted multiplier directly avoids changing the public factor schema.
        converted[aptamer] = {
            'med_streck': 1.0,
            'med_edta': factor['factor'],
        }
    return converted


def _load_bridging_parameters(parameter_source, expected_methods):
    """Load a legacy factor JSON, an in-memory factor mapping, or a factor manifest."""
    if isinstance(parameter_source, (str, bytes, os.PathLike)):
        with open(parameter_source, 'r', encoding='utf-8') as fp:
            data = json.load(fp)
    elif isinstance(parameter_source, Mapping):
        data = dict(parameter_source)
    else:
        raise TypeError(
            "Each params_paths value must be a JSON path or an in-memory mapping"
        )

    if 'factors' in data:
        return _manifest_factors_to_legacy(data, expected_methods)
    for key in ('parameters', 'per_aptamer'):
        if key in data:
            data = data[key]
            break
    if not isinstance(data, Mapping):
        raise ValueError("Bridging parameters must resolve to an aptamer mapping")
    return dict(data)


def _matched_scanner_files(samples, text_files, sample_type):
    matched = samples.copy()
    matched['filename'] = matched.apply(match_file, axis=1, args=(text_files,))
    missing = matched[matched['filename'].isna()]
    if not missing.empty:
        identifying_columns = [
            column for column in ('SampleId', 'slide', 'pdf_subarray')
            if column in missing.columns
        ]
        examples = missing[identifying_columns].head(5).to_dict('records')
        raise FileNotFoundError(
            f"No scanner file matched {len(missing)} {sample_type!r} workbook rows. "
            f"Examples: {examples}"
        )
    duplicated_files = matched[matched['filename'].duplicated(keep=False)]
    if not duplicated_files.empty:
        identifying_columns = [
            column for column in ('SampleId', 'slide', 'pdf_subarray', 'filename')
            if column in duplicated_files.columns
        ]
        examples = duplicated_files[identifying_columns].head(10).to_dict('records')
        raise ValueError(
            f"Multiple {sample_type!r} workbook rows resolve to the same scanner "
            f"file. Examples: {examples}"
        )
    return matched['filename'].to_list()


def _alter_scanner_files(
    text_path,
    workbook_path,
    params_paths,
    output_path,
    *,
    method,
    factor_methods,
    alter_file,
    excluded_aptamers,
    strict_parameters,
):
    text_directory = Path(text_path)
    workbook_file = Path(workbook_path)
    destination = Path(output_path)
    if not text_directory.is_dir():
        raise NotADirectoryError(f"Scanner input directory does not exist: {text_path}")
    if not workbook_file.is_file():
        raise FileNotFoundError(f"Workbook does not exist: {workbook_path}")
    if not isinstance(params_paths, Mapping) or not params_paths:
        raise ValueError("params_paths must be a non-empty sample-type mapping")
    if not isinstance(strict_parameters, bool):
        raise TypeError("strict_parameters must be a boolean")
    if isinstance(excluded_aptamers, (str, bytes)):
        raise TypeError("excluded_aptamers must be a collection, not a string")
    excluded_aptamers = tuple(excluded_aptamers or ())
    if text_directory.resolve() == destination.resolve():
        raise ValueError("output_path must differ from the scanner input directory")
    if (destination / workbook_file.name).resolve() == workbook_file.resolve():
        raise ValueError(
            "output_path would overwrite the source workbook; choose a different "
            "directory"
        )
    if destination.exists() and not destination.is_dir():
        raise NotADirectoryError(f"output_path is not a directory: {destination}")

    workbook_df = read_workbook(workbook_file)
    if WORKBOOK_SAMPLE_MATRIX not in workbook_df.columns:
        raise ValueError(
            f"Workbook data is missing the {WORKBOOK_SAMPLE_MATRIX!r} column"
        )

    files = [
        filename for filename in os.listdir(text_directory)
        if os.path.isfile(os.path.join(text_directory, filename))
    ]
    text_files = sorted(filename for filename in files if filename.endswith('.txt'))
    if not text_files:
        raise FileNotFoundError(f"No scanner .txt files found in {text_path}")

    file_owners = {}
    sample_type_details = {}
    transform_jobs = []
    for sample_type, parameter_source in params_paths.items():
        params = _load_bridging_parameters(parameter_source, factor_methods)
        samples_to_alter = workbook_df[
            workbook_df[WORKBOOK_SAMPLE_MATRIX].isin([sample_type])
        ]
        if samples_to_alter.empty:
            available = sorted(
                str(value)
                for value in workbook_df[WORKBOOK_SAMPLE_MATRIX].dropna().unique()
            )
            raise ValueError(
                f"No workbook rows have sample matrix {sample_type!r}. "
                f"Available values: {available}"
            )
        type_files = _matched_scanner_files(samples_to_alter, text_files, sample_type)
        for filename in type_files:
            previous_owner = file_owners.get(filename)
            if previous_owner is not None and previous_owner != sample_type:
                raise ValueError(
                    f"Scanner file {filename!r} matched both {previous_owner!r} "
                    f"and {sample_type!r}"
                )
            file_owners[filename] = sample_type
            transform_jobs.append((filename, params))
        sample_type_details[str(sample_type)] = {
            'sample_count': int(len(samples_to_alter)),
            'files': type_files,
        }

    files_to_copy = [
        filename for filename in text_files
        if filename not in file_owners
    ]

    destination.parent.mkdir(parents=True, exist_ok=True)
    already_handled_files = []
    with tempfile.TemporaryDirectory(
        dir=destination.parent,
        prefix=f".{destination.name}.staging.",
    ) as staging_name:
        staging = Path(staging_name)
        shutil.copy2(workbook_file, staging / workbook_file.name)
        for filename, params in transform_jobs:
            already_handled_files = alter_file(
                text_directory,
                filename,
                staging,
                params,
                already_handled_files=already_handled_files,
                excluded_aptamers=excluded_aptamers,
                strict_parameters=strict_parameters,
            )
        for filename in files_to_copy:
            shutil.copy2(text_directory / filename, staging / filename)

        if destination.exists():
            for staged_item in staging.iterdir():
                staged_item.replace(destination / staged_item.name)
        else:
            staging.replace(destination)

    manifest = {
        'method': method,
        'input_path': os.fspath(text_directory),
        'workbook_path': os.fspath(workbook_file),
        'output_path': os.fspath(destination),
        'altered_files': list(already_handled_files),
        'copied_files': files_to_copy,
        'sample_types': sample_type_details,
        'excluded_aptamers': sorted(set(excluded_aptamers or ())),
    }
    print('Done!')
    return manifest


def alter_scanner_files_MAD_bridging(
    text_path,
    workbook_path,
    params_paths,
    output_path,
    *,
    excluded_aptamers=DEFAULT_HCE_APTAMERS,
    strict_parameters=True,
):
    """Apply log-MAD bridging to selected scanner files and return a run manifest."""
    return _alter_scanner_files(
        text_path,
        workbook_path,
        params_paths,
        output_path,
        method='mad_log',
        factor_methods={'mad_log'},
        alter_file=alter_scanner_file_MAD_bridging,
        excluded_aptamers=excluded_aptamers,
        strict_parameters=strict_parameters,
    )


def alter_scanner_files_median_bridging(
    text_path,
    workbook_path,
    params_paths,
    output_path,
    *,
    excluded_aptamers=(),
    strict_parameters=True,
):
    """Apply median-only bridging to selected scanner files and return a run manifest."""
    return _alter_scanner_files(
        text_path,
        workbook_path,
        params_paths,
        output_path,
        method='median_multiplier',
        factor_methods={'median_ratio', 'log_median_shift'},
        alter_file=alter_scanner_file_median_bridging,
        excluded_aptamers=excluded_aptamers,
        strict_parameters=strict_parameters,
    )



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
    
    # plates = os.listdir('data/plates/original_')
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
