import os
import json
import shutil
import pandas as pd

from src.bridging import *
from src.bridging_coefficients import *

WORKBOOK_SAMPLE_MATRIX = "Sample Matrix"


def read_text_file(path_to_text):
    with open(path_to_text, 'r') as f:
        return f.read()


def get_sections(text):
    return text.strip().split('*\n')


def analyze_lines(lines):
    type_row = None
    header_row = None
    header_name = ''
    data_rows = []

    for line in lines:
        if line.startswith('TYPE\t'):
            type_row = line.split('\t')[1:]
        elif line.startswith(('FEPARAMS\t', 'STATS\t', 'FEATURES\t')):
            parts = line.split('\t')
            header_name = parts[0]
            header_row = parts[1:]
        elif line.startswith('DATA\t'):
            data_rows.append(line.split('\t')[1:])

    return type_row, header_row, header_name, data_rows


def create_df(type_row, header_row, data_rows, header_name, type_mappings):
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
            except Exception as e:
                print(f"Warning: Could not convert column {col} to {dtype}: {e}")

    return df


def dataframes_to_text(dataframes, type_mappings):
    sections = []

    for name, df in dataframes.items():
        lines = []

        types = type_mappings.get(name, ['text'] * len(df.columns))
        lines.append('TYPE\t' + '\t'.join(types))
        lines.append(name + '\t' + '\t'.join(df.columns))

        for _, row in df.iterrows():
            values = []
            for val in row:
                if pd.isna(val):
                    values.append('')
                elif isinstance(val, bool):
                    values.append('1' if val else '0')
                else:
                    values.append(str(val))

            lines.append('DATA\t' + '\t'.join(values))

        sections.append('\n'.join(lines))

    return '\n*\n'.join(sections) + '\n'


def make_text_file(file_name, text):
    with open(file_name, 'w') as f:
        f.write(text)


def read_workbook(workbook_path):
    sheet_name_tab3 = 'Tab 3 - Assayed Sample List'
    col_smb_tab3 = 'A, F, J'
    col_names_tab3 = ['Sample Number', 'SampleId', WORKBOOK_SAMPLE_MATRIX]

    df_sheet3 = pd.read_excel(workbook_path, sheet_name=sheet_name_tab3,
                              usecols=col_smb_tab3, dtype={'SampleId': str})

    df_data_tab3 = df_sheet3[col_names_tab3].dropna()

    sheet_name_tab4 = 'Tab 4 - Plate Map'
    col_smb_tab4 = 'A, E, H, K'

    df_sheet4 = pd.read_excel(workbook_path, sheet_name=sheet_name_tab4,
                              usecols=col_smb_tab4, skiprows=6, dtype={'NAME': str})

    df_data_tab4 = df_sheet4[df_sheet4['NAME'].notna()]
    df_data_tab4['Slide #'] = df_data_tab4['Slide #'].ffill()

    excel_df = pd.merge(df_data_tab3,
                        df_data_tab4,
                        how='inner',
                        left_on='SampleId', right_on='NAME').drop('NAME', axis=1)

    excel_df['slide'] = excel_df['Slide #']
    excel_df.rename(columns={'PDF Subarray': 'pdf_subarray'}, inplace=True)

    return excel_df


def match_file(row, files):
    for f in files:
        if str(row['slide']) in f and f.endswith(row['pdf_subarray'] + '.txt'):
            return f
    return None


def match_probe_aptamer(df):
    probe_names = list(set(df['ProbeName'].to_list()))
    probe_names = [p for p in probe_names if p.startswith('anti-')]

    return {
        p: p.split('anti-')[1].split('_')[0]
        for p in probe_names
    }


def transform_feature_df_MAD_bridging(df1, params):
    hce_list = [
        '2171-12','2178-55','2194-91','2229-54','2249-25','2273-34',
        '2288-7','2305-52','2312-13','2359-65','2430-52','2513-7'
    ]

    df = df1.copy()
    probe_map = match_probe_aptamer(df)

    original_columns = df.columns
    df['aptamer'] = df['ProbeName'].map(probe_map).fillna('other')

    df = aggregate(df)

    aptamers = list(probe_map.values())
    aptamers_to_alter = [a for a in aptamers if a not in hce_list]

    df = transform_mad_log(df, params, aptamers_to_alter)

    return df[original_columns]


def alter_scanner_file_MAD_bridging(input_path, file, output_path, params, already_handled_files=[]):
    file_path = os.path.join(input_path, file)

    text = read_text_file(file_path)
    print(f"Altered: {file}")

    sections = get_sections(text)
    dataframes_dict = {}
    type_mappings = {}

    for section in sections:
        lines = section.strip().split('\n')
        type_row, header_row, header_name, data_rows = analyze_lines(lines)

        df = create_df(type_row, header_row, data_rows, header_name, type_mappings)

        if header_name == 'FEATURES':
            df = transform_feature_df_MAD_bridging(df, params)

        dataframes_dict[header_name] = df

    text_out = dataframes_to_text(dataframes_dict, type_mappings)

    output_file_path = os.path.join(output_path, file)
    make_text_file(output_file_path, text_out)

    already_handled_files.append(file)
    return already_handled_files


def alter_scanner_files_MAD_bridging(text_path, workbook_path, params_paths, output_path):
    workbook_df = read_workbook(workbook_path)

    os.makedirs(output_path, exist_ok=True)
    shutil.copy(workbook_path, output_path)

    files = [
        f for f in os.listdir(text_path)
        if os.path.isfile(os.path.join(text_path, f))
    ]

    text_files = [f for f in files if f.endswith('.txt')]
    already_handled_files = []

    for sample_type, params_path in params_paths.items():
        with open(params_path, 'r') as fp:
            params = json.load(fp)

        samples = workbook_df[
            workbook_df[WORKBOOK_SAMPLE_MATRIX].isin([sample_type])
        ]

        samples['filename'] = samples.apply(match_file, axis=1, args=(text_files,))
        type_files = samples['filename'].to_list()

        for file in type_files:
            already_handled_files.append(
                alter_scanner_file_MAD_bridging(
                    text_path, file, output_path, params, already_handled_files
                )
            )

    files_to_copy = [f for f in text_files if f not in already_handled_files]

    for file in files_to_copy:
        shutil.copy(os.path.join(text_path, file), output_path)

    print('Done!')
    
################################################################################################################
################################################# EXAMPLE USE ##################################################
################################################################################################################
# if __name__ == '__main__':
#     pass
    
    # for plate in plates:
            

        # alter_scanner_files_MAD_bridging(<plate folder>,
        #         <plate folder>/plate Workbook.xlsx',
                
        #         {'DS': <Streck (DS) bridging coefficients file path>,
        #         'nDS': <Streck (nDS) bridging coefficients file path>},
        #         f'output folder path (will be created if does not exist)',
                
        # )