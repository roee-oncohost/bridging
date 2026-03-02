import numpy as np
import pandas as pd

def mad(x):
    """Unscaled median absolute deviation."""
    med = np.nanmedian(x)
    return np.nanmedian(np.abs(x - med))

def fit_mad_log_transform(streck_train, edta_train, aptamers):
    """
    Fit MAD-based log-scale alignment parameters.
    
    Returns a dict of parameters to reuse on test data.
    """

    log_streck = np.log(streck_train[aptamers])
    log_edta = np.log(edta_train[aptamers])

    params = pd.DataFrame({
        "med_streck": log_streck.median(),
        "mad_streck": log_streck.apply(mad),
        "med_edta": log_edta.median(),
        "mad_edta": log_edta.apply(mad),
    })

    # avoid division by zero later
    params["mad_streck"].replace(0, np.nan, inplace=True)
    params_dict = params.to_dict('index')

    return params_dict

def fit_mad_linear_transform(streck_train, edta_train, aptamers):
    """
    Fit MAD-based log-scale alignment parameters.
    
    Returns a dict of parameters to reuse on test data.
    """

    linear_streck = streck_train[aptamers]
    linear_edta = edta_train[aptamers]

    params = pd.DataFrame({
        "med_streck": linear_streck.median(),
        "mad_streck": linear_streck.apply(mad),
        "med_edta": linear_edta.median(),
        "mad_edta": linear_edta.apply(mad),
    })

    # avoid division by zero later
    params["mad_streck"].replace(0, np.nan, inplace=True)
    params_dict = params.to_dict('index')

    return params_dict


def fit_multiple_log_transforms(streck_train, edta_train, aptamers, sample_types):
    params = {}
    for sample_type in sample_types:
        streck_subset = streck_train[streck_train['SampleType'] == sample_type]
        # edta_subset = edta_train[edta_train['SampleType'] == sample_type]
        params[sample_type] = fit_mad_log_transform(streck_subset, edta_train, aptamers)
        # yield sample_type, params
    return params

# def mad(x):
#     """Unscaled median absolute deviation."""
#     med = np.nanmedian(x)
#     return np.nanmedian(np.abs(x - med))

def fit_sd_med_log_transform(streck_train, edta_train, aptamers):
    """
    Fit sd/mean-based log-scale alignment parameters.
    
    Returns a dict of parameters to reuse on test data.
    """

    log_streck = np.log(streck_train[aptamers])
    log_edta = np.log(edta_train[aptamers])

    params = pd.DataFrame({
        "mean_streck": log_streck.mean(),
        "std_streck": log_streck.std(),
        "mean_edta": log_edta.mean(),
        "std_edta": log_edta.std(),
    })

    # avoid division by zero later
    params["std_streck"].replace(0, np.nan, inplace=True)
    params_dict = params.to_dict('index')

    return params_dict
