import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from src.dataframe_transformation import get_aptamers


################################################################################################################
############################################ GENERAL FUNCTIONS #################################################
################################################################################################################

def aggregate_(df1):
    import time
    # start = time.time()
    df = df1.copy()
    for aptamer in df['aptamer'].unique():
        if aptamer[0].isdigit():
            df2 = df[df['aptamer']==aptamer]
            idx = df2.index
            relevant_df = df2[~df2['gIsFeatPopnOL']]
            mean = relevant_df['gProcessedSignal'].mean()
            df.loc[idx, 'gProcessedSignal'] = mean
    # print(f"Aggregation took {time.time() - start} seconds")
    return df

def aggregate(df1):
    import time
    # start = time.time()
    df = df1.copy()
    
    # Filter to aptamers starting with a digit
    mask = df['aptamer'].str[0].str.isdigit()
    
    # Calculate means for each aptamer where gIsFeatPopnOL == 0
    means = (df[mask & (df['gIsFeatPopnOL'] == 0)]
             .groupby('aptamer')['gProcessedSignal']
             .mean())
    
    # Map the means back to all rows with those aptamers
    df.loc[mask, 'gProcessedSignal'] = df.loc[mask, 'aptamer'].map(means)
    # print(f"Aggregation took {time.time() - start} seconds")
    return df

################################################################################################################
########################################## MEDIAN/MEAN MAPPING #################################################
################################################################################################################

def bridge_dataframes(old_df: pd.DataFrame, new_df: pd.DataFrame)-> dict:
    
    aptamers = get_aptamers(old_df)
    numerators = old_df[aptamers].median()
    denominators = new_df[aptamers].median()
    ratios = numerators/denominators
    ratios = ratios.reset_index()
    ratios.columns = ['Aptamer', 'Ratio']
    ratios_dict = {row['Aptamer']: row['Ratio'] for index, row in ratios.iterrows()}
    return ratios_dict

def bridge_dataframes_mean(old_df: pd.DataFrame, new_df: pd.DataFrame)-> dict:
    
    aptamers = get_aptamers(old_df)
    numerators = old_df[aptamers].mean()
    denominators = new_df[aptamers].mean()
    ratios = numerators/denominators
    ratios = ratios.reset_index()
    ratios.columns = ['Aptamer', 'Ratio']
    ratios_dict = {row['Aptamer']: row['Ratio'] for index, row in ratios.iterrows()}
    return ratios_dict






################################################################################################################
########################################### QUANTILE MAPPING ###################################################
################################################################################################################

def calculate_quantiles(data, n_quantiles=100):
    """
    Calculate quantiles for a dataset
    
    Parameters:
    - data: numpy array or list of values
    - n_quantiles: number of quantiles to calculate (default 100)
    
    Returns:
    - numpy array of quantile values
    """
    return np.percentile(data, np.linspace(0, 100, n_quantiles + 1))

def quantile_map(values, source_quantiles, target_quantiles):
    """
    Map values from source distribution to target distribution using interpolation
    
    Parameters:
    - values: numpy array of values to map
    - source_quantiles: quantiles of the source distribution (streck)
    - target_quantiles: quantiles of the target distribution (edta)
    
    Returns:
    - numpy array of mapped values
    """
    # Use numpy's interp for efficient linear interpolation
    mapped = np.interp(values, source_quantiles, target_quantiles)
    return mapped

def quantile_bridge_streck_to_edta(edta_data, streck_data, n_quantiles=10):
    """
    Bridge streck data to edta distribution using quantile mapping
    
    Parameters:
    - edta_data: numpy array of edta values (reference distribution)
    - streck_data: numpy array of streck values to be bridged
    - n_quantiles: number of quantiles to use (default 100)
    
    Returns:
    - bridged_data: streck data mapped to edta distribution
    - streck_quantiles: quantiles of streck distribution (for inspection)
    - edta_quantiles: quantiles of edta distribution (for inspection)
    """
    # Calculate quantiles for both distributions
    edta_quantiles = calculate_quantiles(edta_data, n_quantiles)
    streck_quantiles = calculate_quantiles(streck_data, n_quantiles)
    
    # Map streck values to edta space
    bridged_data = quantile_map(streck_data, streck_quantiles, edta_quantiles)
    
    return bridged_data, streck_quantiles, edta_quantiles

def learn_quantile_mapping(edta_data, streck_data, n_quantiles=100):
    """
    Learn the quantile mapping from edta and streck training data
    This creates a "mapper" that can be applied to new streck data
    
    Parameters:
    - edta_data: pandas DataFrame or numpy array (n_samples_edta, n_variables) for edta training
    - streck_data: pandas DataFrame or numpy array (n_samples_streck, n_variables) for streck training
    - n_quantiles: number of quantiles to use
    
    Returns:
    - mapper: dictionary containing the learned mapping for each variable
    """
    # Handle DataFrame inputs
    is_dataframe = isinstance(edta_data, pd.DataFrame)
    column_names = None
    
    if is_dataframe:
        column_names = edta_data.columns.tolist()
        edta_array = edta_data.values
        streck_array = streck_data.values
    else:
        edta_array = edta_data
        streck_array = streck_data
        column_names = [f'Var_{i}' for i in range(edta_array.shape[1])]
    
    n_vars = edta_array.shape[1]
    mapper = {
        'column_names': column_names,
        'n_quantiles': n_quantiles,
        'mappings': {}
    }
    
    for i in range(n_vars):
        var_name = column_names[i]
        edta_quantiles = calculate_quantiles(edta_array[:, i], n_quantiles)
        streck_quantiles = calculate_quantiles(streck_array[:, i], n_quantiles)
        
        mapper['mappings'][var_name] = {
            'edta_quantiles': edta_quantiles,
            'streck_quantiles': streck_quantiles
        }
    
    return mapper

def apply_quantile_mapping(streck_new_data, mapper):
    """
    Apply a learned quantile mapping to new streck data (test set)
    
    Parameters:
    - streck_new_data: pandas DataFrame or numpy array of new streck data to bridge
    - mapper: dictionary returned by learn_quantile_mapping()
    
    Returns:
    - bridged_data: pandas DataFrame or numpy array with bridged values
    """
    # Handle DataFrame inputs
    is_dataframe = isinstance(streck_new_data, pd.DataFrame)
    
    if is_dataframe:
        streck_array = streck_new_data.values
        column_names = streck_new_data.columns.tolist()
    else:
        streck_array = streck_new_data
        column_names = mapper['column_names']
    
    # Verify columns match
    if column_names != mapper['column_names']:
        raise ValueError(f"Column mismatch! Expected {mapper['column_names']}, got {column_names}")
    
    n_vars = streck_array.shape[1]
    bridged_array = np.zeros_like(streck_array)
    
    for i in range(n_vars):
        var_name = column_names[i]
        streck_quantiles = mapper['mappings'][var_name]['streck_quantiles']
        edta_quantiles = mapper['mappings'][var_name]['edta_quantiles']
        
        # Apply the mapping
        bridged_array[:, i] = quantile_map(streck_array[:, i], streck_quantiles, edta_quantiles)
    
    # Return DataFrame if input was DataFrame
    if is_dataframe:
        bridged_data = pd.DataFrame(bridged_array, columns=column_names, index=streck_new_data.index)
    else:
        bridged_data = bridged_array
    
    return bridged_data

def save_mapper(mapper, filename):
    """Save the learned mapping to a file"""
    import pickle
    with open(filename, 'wb') as f:
        pickle.dump(mapper, f)
    print(f"Mapper saved to {filename}")

def load_mapper(filename):
    """Load a saved mapping from a file"""
    import pickle
    with open(filename, 'rb') as f:
        mapper = pickle.load(f)
    print(f"Mapper loaded from {filename}")
    return mapper

def bridge_multiple_variables(edta_data, streck_data, n_quantiles=100):
    """
    Bridge multiple variables at once (convenience function for training data)
    
    NOTE: This is for exploring/testing. For production use:
    1. Use learn_quantile_mapping() on training data
    2. Use apply_quantile_mapping() on test data
    
    Parameters:
    - edta_data: pandas DataFrame or numpy array of shape (n_samples_edta, n_variables) for edta
    - streck_data: pandas DataFrame or numpy array of shape (n_samples_streck, n_variables) for streck
              Note: can have different number of rows than edta_data
    - n_quantiles: number of quantiles to use
    
    Returns:
    - bridged_data: pandas DataFrame (if input was DataFrame) or numpy array with bridged values
                    Same shape as streck_data
    - quantile_info: dictionary containing quantile information for each variable
    """
    mapper = learn_quantile_mapping(edta_data, streck_data, n_quantiles)
    bridged_data = apply_quantile_mapping(streck_data, mapper)
    return bridged_data, mapper['mappings']

def plot_mapping_analysis(edta_data, streck_data, bridged_data, variable_name="Variable"):
    """
    Create diagnostic plots for the quantile mapping
    
    Parameters:
    - edta_data: original edta data
    - streck_data: original streck data
    - bridged_data: bridged streck data
    - variable_name: name for plot titles
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: Distributions comparison
    axes[0, 0].hist(edta_data, bins=30, alpha=0.5, label='edta', density=True)
    axes[0, 0].hist(streck_data, bins=30, alpha=0.5, label='streck Original', density=True)
    axes[0, 0].hist(bridged_data, bins=30, alpha=0.5, label='streck Bridged', density=True)
    axes[0, 0].set_xlabel('Value')
    axes[0, 0].set_ylabel('Density')
    axes[0, 0].set_title(f'{variable_name}: Distribution Comparison')
    axes[0, 0].legend()
    
    # Plot 2: Q-Q plot (quantile-quantile)
    edta_sorted = np.sort(edta_data)
    streck_sorted = np.sort(streck_data)
    min_len = min(len(edta_sorted), len(streck_sorted))
    edta_sample = edta_sorted[np.linspace(0, len(edta_sorted)-1, min_len, dtype=int)]
    streck_sample = streck_sorted[np.linspace(0, len(streck_sorted)-1, min_len, dtype=int)]
    
    axes[0, 1].scatter(streck_sample, edta_sample, alpha=0.5, s=10)
    axes[0, 1].plot([streck_sample.min(), streck_sample.max()], 
                     [streck_sample.min(), streck_sample.max()], 
                     'r--', label='y=x (perfect match)')
    axes[0, 1].set_xlabel('streck Quantiles')
    axes[0, 1].set_ylabel('edta Quantiles')
    axes[0, 1].set_title(f'{variable_name}: Q-Q Plot (streck vs edta)')
    axes[0, 1].legend()
    
    # Plot 3: Mapping function
    streck_q = calculate_quantiles(streck_data, 100)
    edta_q = calculate_quantiles(edta_data, 100)
    axes[1, 0].plot(streck_q, edta_q, 'b-', linewidth=2)
    axes[1, 0].set_xlabel('streck Value')
    axes[1, 0].set_ylabel('edta Value (Mapped)')
    axes[1, 0].set_title(f'{variable_name}: Mapping Function')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Plot 4: Before/After scatter
    axes[1, 1].scatter(streck_data, bridged_data, alpha=0.3, s=10)
    axes[1, 1].set_xlabel('streck Original')
    axes[1, 1].set_ylabel('streck Bridged')
    axes[1, 1].set_title(f'{variable_name}: Original vs Bridged streck')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig

def compare_statistics(edta_data, streck_data, bridged_data, variable_names=None):
    """
    Compare statistics between edta, streck original, and streck bridged
    
    Parameters:
    - edta_data: pandas DataFrame or numpy array (n_samples_edta, n_variables) for edta
    - streck_data: pandas DataFrame or numpy array (n_samples_streck, n_variables) for streck
    - bridged_data: pandas DataFrame or numpy array (n_samples_streck, n_variables) for bridged streck
    - variable_names: list of variable names (optional, auto-detected from DataFrame)
    
    Returns:
    - pandas DataFrame with comparison statistics
    """
    # Handle DataFrame inputs
    if isinstance(edta_data, pd.DataFrame):
        variable_names = edta_data.columns.tolist()
        edta_array = edta_data.values
        streck_array = streck_data.values
        bridged_array = bridged_data.values
    else:
        edta_array = edta_data
        streck_array = streck_data
        bridged_array = bridged_data
    
    n_vars = edta_array.shape[1] if len(edta_array.shape) > 1 else 1
    
    if variable_names is None:
        variable_names = [f'Var_{i}' for i in range(n_vars)]
    
    stats_list = []
    
    if n_vars == 1:
        edta_array = edta_array.reshape(-1, 1)
        streck_array = streck_array.reshape(-1, 1)
        bridged_array = bridged_array.reshape(-1, 1)
    
    for i in range(n_vars):
        stats_list.append({
            'Variable': variable_names[i],
            'edta_n_samples': len(edta_array[:, i]),
            'streck_n_samples': len(streck_array[:, i]),
            'edta_mean': np.mean(edta_array[:, i]),
            'edta_median': np.median(edta_array[:, i]),
            'edta_std': np.std(edta_array[:, i]),
            'streck_mean': np.mean(streck_array[:, i]),
            'streck_median': np.median(streck_array[:, i]),
            'streck_std': np.std(streck_array[:, i]),
            'Bridged_mean': np.mean(bridged_array[:, i]),
            'Bridged_median': np.median(bridged_array[:, i]),
            'Bridged_std': np.std(bridged_array[:, i]),
        })
    
    return pd.DataFrame(stats_list)

# Example usage
"""
import pandas as pd
import numpy as np

# ===== TRAINING PHASE =====
# Load your TRAINING data
edta_train = pd.read_csv('edta_train.csv')  # 1000 rows x 30 columns
streck_train = pd.read_csv('streck_train.csv')  # 30 rows x 30 columns

print(f"edta training shape: {edta_train.shape}")
print(f"streck training shape: {streck_train.shape}")

# Learn the mapping from training data
mapper = learn_quantile_mapping(edta_train, streck_train, n_quantiles=100)

# Save the mapper for future use
save_mapper(mapper, 'quantile_mapper.pkl')

# Optional: bridge the training data to check quality
streck_train_bridged = apply_quantile_mapping(streck_train, mapper)
stats = compare_statistics(edta_train, streck_train, streck_train_bridged)
print(stats)

# ===== TEST/PRODUCTION PHASE =====
# Load the saved mapper
mapper = load_mapper('quantile_mapper.pkl')

# Load your NEW streck test data (any number of rows)
streck_test = pd.read_csv('streck_test.csv')  # e.g., 100 rows x 30 columns

# Apply the learned mapping to test data
streck_test_bridged = apply_quantile_mapping(streck_test, mapper)

# Now use streck_test_bridged with your model
# predictions = model.predict(streck_test_bridged)

# Save bridged test data
streck_test_bridged.to_csv('streck_test_bridged.csv', index=False)
"""


################################################################################################################
########################################### LOG-NORMALIZED SD-BASED MAPPING ####################################
################################################################################################################


def get_a_b(target_distribution, source_distribution):
    log_e = target_distribution.apply(np.log)
    log_s = source_distribution.apply(np.log)
    sigma_log_e = log_e.std()
    sigma_log_s = log_s.std()
    mu_log_e = log_e.median()
    mu_log_s = log_s.median()
    a = sigma_log_e/sigma_log_s
    b = mu_log_e - a * mu_log_s
    return a, b


def get_a_b_dict(aptamers, edta_df, streck_df):
    a_b_dict = {}
    for aptamer in aptamers:
        a, b = get_a_b(edta_df[aptamer], streck_df[aptamer])
        a_b_dict[aptamer] = {'a': a,
                             'b': b}
    return a_b_dict

def normalize_value(x, a, b):
    return np.exp(b) * (x ** a)
    

def transform_streck_df(streck_df, a_b_dict, aptamers):
    df = streck_df.copy()
    
    for aptamer in aptamers:
        a, b = a_b_dict[aptamer]['a'], a_b_dict[aptamer]['b'] 
        df[aptamer] = df[aptamer].apply(lambda x: np.exp(b) * (x ** a))
    return df




################################################################################################################
########################################## LOG-NORMALIZED MAD-BASED MAPPING ####################################
################################################################################################################

# def mad(x):
#     """Unscaled median absolute deviation."""
#     med = np.nanmedian(x)
#     return np.nanmedian(np.abs(x - med))

# def fit_mad_log_transform(streck_train, edta_train, aptamers=[]):
#     """
#     Fit MAD-based log-scale alignment parameters.
    
#     Parameters:
#     - streck_train: pandas dataframe of Streck signals (or generally any other method's signals), to be transformed to fit a target distribution
#     - edta_train: pandas DataFrame of EDTA signals (any other method's signals), to serve as target distribution 
#     - aptamers: list of column names. Defaults to empty list, in which case we use all columns whose names begin with a digit

    
#     Returns:
#     - dict of parameters to reuse on test data
#     """
#     if not aptamers:
#         aptamers = [column for column in streck_train.columns if column[0].isdigit()]
        
#     log_streck = np.log(streck_train[aptamers])
#     log_edta = np.log(edta_train[aptamers])

#     params = pd.DataFrame({
#         "med_streck": log_streck.median(),
#         "mad_streck": log_streck.apply(mad),
#         "med_edta": log_edta.median(),
#         "mad_edta": log_edta.apply(mad),
#     })

#     # avoid division by zero later
#     params["mad_streck"].replace(0, np.nan, inplace=True)

#     return params


def transform_mad_log_(streck_df, params, aptamers): #, return_log=False):
    """
    Apply a fitted MAD-based log transform to new S data.
    """
    # if not aptamers:
    #     aptamers = [column for column in streck_df.columns if column[0].isdigit()]

    for aptamer in aptamers:
        # if aptamer != 'other':
        param = params[aptamer]
        temp_df = streck_df[streck_df['aptamer']==aptamer].copy()
        idx = temp_df.index 
        temp_df['gProcessedSignal'] = temp_df['gProcessedSignal'].apply(np.log)
        a = param['mad_edta']/param['mad_streck']
        b = param['med_edta'] - a * param['med_streck']
        log_streck_star = temp_df['gProcessedSignal'] * a + b
        streck_df.loc[idx, 'gProcessedSignal'] = np.exp(log_streck_star) 
        
        # streck_df.loc[idx, 'gProcessedSignal'] = temp_df['gpro']
        # temp_df.loc[idx, 'value'] = temp_df[temp_df['key'] == 1]['value'].apply(np.log)
    # log_streck = np.log(streck_df[aptamers])


    return streck_df

def transform_mad_log(streck_df, params, aptamers):
    """
    Apply a fitted MAD-based log transform to new S data.
    """
    df = streck_df.copy()
    mask = df['aptamer'].isin(aptamers)
    
    # Create parameter lookup DataFrames
    param_df = pd.DataFrame.from_dict(params, orient='index')
    a_map = (param_df['mad_edta'] / param_df['mad_streck']).to_dict()
    b_map = (param_df['med_edta'] - 
             (param_df['mad_edta'] / param_df['mad_streck']) * param_df['med_streck']).to_dict()
    
    # Vectorized transformation
    a = df.loc[mask, 'aptamer'].map(a_map)
    b = df.loc[mask, 'aptamer'].map(b_map)
    
    log_signal = np.log(df.loc[mask, 'gProcessedSignal'])
    log_streck_star = log_signal * a + b
    df.loc[mask, 'gProcessedSignal'] = np.exp(log_streck_star)
    
    return df

def transform_mad_linear(streck_df, params, aptamers):
    """
    Apply a fitted MAD-based log transform to new S data.
    """
    df = streck_df.copy()
    mask = df['aptamer'].isin(aptamers)
    
    # Create parameter lookup DataFrames
    param_df = pd.DataFrame.from_dict(params, orient='index')
    a_map = (param_df['mad_edta'] / param_df['mad_streck']).to_dict()
    b_map = (param_df['med_edta'] - 
             ((param_df['mad_edta'] / param_df['mad_streck']) * param_df['med_streck'])).to_dict()
    
    # Vectorized transformation
    a = df.loc[mask, 'aptamer'].map(a_map)
    b = df.loc[mask, 'aptamer'].map(b_map)
    
    log_signal = np.log(df.loc[mask, 'gProcessedSignal'])
    streck_star = df.loc[mask, 'gProcessedSignal'] * a + b
    df.loc[mask, 'gProcessedSignal'] = streck_star # np.exp(streck_star)
    
    return df