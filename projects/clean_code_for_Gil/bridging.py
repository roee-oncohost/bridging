import numpy as np
import pandas as pd


def aggregate(df1):
    df = df1.copy()
    
    mask = df['aptamer'].str[0].str.isdigit()
    
    means = (df[mask & (df['gIsFeatPopnOL'] == 0)]
             .groupby('aptamer')['gProcessedSignal']
             .mean())
    
    df.loc[mask, 'gProcessedSignal'] = df.loc[mask, 'aptamer'].map(means)
    return df


def transform_mad_log(streck_df, params, aptamers):
    df = streck_df.copy()
    mask = df['aptamer'].isin(aptamers)
    
    param_df = pd.DataFrame.from_dict(params, orient='index')
    
    a_map = (param_df['mad_edta'] / param_df['mad_streck']).to_dict()
    b_map = (
        param_df['med_edta'] -
        (param_df['mad_edta'] / param_df['mad_streck']) * param_df['med_streck']
    ).to_dict()
    
    a = df.loc[mask, 'aptamer'].map(a_map)
    b = df.loc[mask, 'aptamer'].map(b_map)
    
    log_signal = np.log(df.loc[mask, 'gProcessedSignal'])
    log_streck_star = log_signal * a + b
    
    df.loc[mask, 'gProcessedSignal'] = np.exp(log_streck_star)
    
    return df