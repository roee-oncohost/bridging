import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


def plot_overlaid_histograms(source, target, bins=50, alpha=0.6, 
                             figsize=(10, 6), xtitle='Transformed', ytitle='Target', title="Distribution Comparison",
                             trim_percentile=None):
    """
    Plot overlaid histograms of transformed and target distributions.
    
    Parameters:
    -----------
    source : array-like
        Transformed distribution samples
    target : array-like
        Target distribution samples
    bins : int or str, default=50
        Number of bins for histogram (or 'auto', 'sturges', etc.)
    alpha : float, default=0.6
        Transparency level (0-1)
    figsize : tuple, default=(10, 6)
        Figure size (width, height)
    title : str
        Plot title
    trim_percentile : float, optional
        If specified, trim to middle percentile (e.g., 95 for middle 95%)
        
    Returns:
    --------
    fig, ax : matplotlib figure and axis objects
    """
    source = np.asarray(source)
    target = np.asarray(target)
    
    # Trim to middle percentile if requested
    if trim_percentile is not None:
        lower = (100 - trim_percentile) / 2
        upper = 100 - lower
        
        transformed_lower, transformed_upper = np.percentile(source, [lower, upper])
        target_lower, target_upper = np.percentile(target, [lower, upper])
        
        source = source[(source >= transformed_lower) & 
                                  (source <= transformed_upper)]
        target = target[(target >= target_lower) & (target <= target_upper)]
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot histograms
    ax.hist(source, bins=bins, alpha=alpha, label=xtitle, 
            color='blue', density=True, edgecolor='black', linewidth=0.5)
    ax.hist(target, bins=bins, alpha=alpha, label=ytitle, 
            color='red', density=True, edgecolor='black', linewidth=0.5)
    
    # Add labels and legend
    ax.set_xlabel('Value', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    return fig, ax


def plot_qq(source, target, figsize=(8, 8), title="Q-Q Plot",
            trim_percentile=None):
    """
    Create a Q-Q plot comparing transformed and target distributions.
    
    Parameters:
    -----------
    source : array-like
        Transformed distribution samples
    target : array-like
        Target distribution samples
    figsize : tuple, default=(8, 8)
        Figure size (width, height)
    title : str
        Plot title
    trim_percentile : float, optional
        If specified, trim to middle percentile (e.g., 95 for middle 95%)
        
    Returns:
    --------
    fig, ax : matplotlib figure and axis objects
    """
    # Calculate quantiles
    source = np.asarray(source)
    target = np.asarray(target)
    
    # Trim to middle percentile if requested
    if trim_percentile is not None:
        lower = (100 - trim_percentile) / 2
        upper = 100 - lower
        
        transformed_lower, transformed_upper = np.percentile(source, [lower, upper])
        target_lower, target_upper = np.percentile(target, [lower, upper])
        
        source = source[(source >= transformed_lower) & 
                                  (source <= transformed_upper)]
        target = target[(target >= target_lower) & (target <= target_upper)]
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Sort both arrays
    transformed_sorted = np.sort(source)
    target_sorted = np.sort(target)
    
    # Use the minimum length for comparison
    n = min(len(transformed_sorted), len(target_sorted))
    
    # Interpolate to get matching quantiles
    quantiles = np.linspace(0, 1, n)
    transformed_quantiles = np.quantile(source, quantiles)
    target_quantiles = np.quantile(target, quantiles)
    
    # Plot Q-Q
    ax.scatter(target_quantiles, transformed_quantiles, alpha=0.6, s=20)
    
    # Add diagonal reference line (perfect match)
    min_val = min(target_quantiles.min(), transformed_quantiles.min())
    max_val = max(target_quantiles.max(), transformed_quantiles.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2, 
            label='Perfect Match', zorder=3)
    
    # Add labels
    ax.set_xlabel('Target Quantiles', fontsize=12)
    ax.set_ylabel('Transformed Quantiles', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_aspect('equal', adjustable='box')
    
    plt.tight_layout()
    return fig, ax


def compare_distributions_stats(source, target):
    """
    Compute statistical tests comparing calibrated and target distributions.
    
    Parameters:
    -----------
    calibrated : array-like
        Calibrated values
    target : array-like
        Target distribution values
        
    Returns:
    --------
    results : dict
        Dictionary containing test results:
        - 'ks_statistic': Kolmogorov-Smirnov test . Good similarity: <= 0.05
        - 'ks_pvalue': K-S test p-value. Good similarity: >= 0.05
        - 'wasserstein': Wasserstein distance The lower the better
        - 'mean_diff': Difference in means
        - 'std_diff': Difference in standard deviations
        - 'log_mean_diff': Difference in log-space means
        - 'log_std_diff': Difference in log-space standard deviations
    """
    source = np.asarray(source)
    target = np.asarray(target)
    
    # Kolmogorov-Smirnov test
    ks_stat, ks_pval = stats.ks_2samp(source, target)
    
    # Wasserstein distance
    wasserstein = stats.wasserstein_distance(source, target)
    
    # Basic statistics
    mean_diff = np.mean(source) - np.mean(target)
    std_diff = np.std(source, ddof=1) - np.std(target, ddof=1)
    
    # Log-space statistics (if positive)
    log_mean_diff = np.nan
    log_std_diff = np.nan
    if np.all(source > 0) and np.all(target > 0):
        log_calibrated = np.log(source)
        log_target = np.log(target)
        log_mean_diff = np.mean(log_calibrated) - np.mean(log_target)
        log_std_diff = np.std(log_calibrated, ddof=1) - np.std(log_target, ddof=1)
    
    results = {
        'ks_statistic': ks_stat,
        'ks_pvalue': ks_pval,
        'wasserstein': wasserstein,
        'mean_diff': mean_diff,
        'std_diff': std_diff,
        'log_mean_diff': log_mean_diff,
        'log_std_diff': log_std_diff
    }
    
    return results


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def scatter_comparison(df1, 
                       df2, 
                       metric_type='mean',
                       xlabel='EDTA', 
                       ylabel='Streck', 
                       title='Comparison Across Proteins',
                       trim_percentile=100,
                       add_abels=False):
    """
    Create a scatter plot comparing specified metrics ('mean', 'std', 'median', 'mad') of columns between two dataframes.
    
    Parameters:
    -----------
    df1 : pd.DataFrame
        First dataframe with numerical columns
    df2 : pd.DataFrame
        Second dataframe with numerical columns (should have same columns as df1)
    metric_type : str
        Type of metric to compute for comparison ('mean', 'std', 'median', 'mad')
    xlabel : str
        Label for x-axis
    ylabel : str
        Label for y-axis
    title : str
        Title for the plot
    trim_percent : float
        Percentage (0-100) of extreme points to trim from each end based on 
        distance from origin. Default is 0 (no trimming).
    
    Returns:
    --------
    fig, ax : matplotlib figure and axis objects
    """
    # Calculate standard deviations for each column
    if metric_type == 'mean':
        metric_df1 = df1.mean()
        metric_df2 = df2.mean()
    elif metric_type == 'std':
        metric_df1 = df1.std()
        metric_df2 = df2.std()
    elif metric_type == 'median':
        metric_df1 = df1.median()
        metric_df2 = df2.median()
    elif metric_type == 'mad':
        metric_df1 = (df1 - df1.median()).abs().median()
        metric_df2 = (df2 - df2.median()).abs().median()
    else:
        raise ValueError("Invalid metric_type. Choose from 'mean', 'std', 'median', 'mad'.")
    
    # Create dataframe for easy filtering
    metric_data = pd.DataFrame({
        'column': metric_df1.index,
        'metric1': metric_df1.values,
        'metric2': metric_df2.values
    })
    
    # Apply trimming if specified
    if trim_percentile < 100:
        # Calculate distance from origin for each point
        metric_data['distance'] = np.sqrt(metric_data['metric1']**2 + metric_data['metric2']**2)
        
        # Calculate percentiles
        lower_percentile = (100 - trim_percentile) / 2
        upper_percentile = 100 - lower_percentile
        
        lower_bound = np.percentile(metric_data['distance'], lower_percentile)
        upper_bound = np.percentile(metric_data['distance'], upper_percentile)
        
        # Filter data
        metric_data = metric_data[(metric_data['distance'] >= lower_bound) & 
                           (metric_data['distance'] <= upper_bound)]
    
    # Create scatter plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(metric_data['metric1'], metric_data['metric2'], s=100, alpha=0.6)
    
    # Add labels for each point (column names)
    if add_abels:
        for _, row in metric_data.iterrows():
            ax.annotate(row['column'], (row['metric1'], row['metric2']), 
                    xytext=(5, 5), textcoords='offset points')
    
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.grid(True, alpha=0.3)
    
    # Add a diagonal reference line (y=x)
    if len(metric_data) > 0:
        max_val = max(metric_data['metric1'].max(), metric_data['metric2'].max())
        min_val = min(metric_data['metric1'].min(), metric_data['metric2'].min())
        ax.plot([min_val, max_val], [min_val, max_val], 
                'r--', alpha=0.5, label='y=x reference')
        ax.legend()
    
    plt.tight_layout()
    return fig, ax


# # Example usage:
# if __name__ == "__main__":
#     # Example dataframes
#     df1 = pd.DataFrame({
#         'col_A': [1, 2, 3, 4, 5],
#         'col_B': [10, 20, 30, 40, 50],
#         'col_C': [100, 110, 105, 115, 120],
#         'col_D': [1000, 1100, 1050, 1150, 1200]  # Potential outlier
#     })
    
#     df2 = pd.DataFrame({
#         'col_A': [2, 4, 6, 8, 10],
#         'col_B': [5, 15, 25, 35, 45],
#         'col_C': [95, 100, 105, 110, 115],
#         'col_D': [900, 1000, 1050, 1100, 1150]  # Potential outlier
#     })
    
#     # Plot without trimming
#     fig1, ax1 = scatter_comparison(df1, df2)
#     plt.show()
    
#     # Plot with 10% trimming
#     fig2, ax2 = scatter_comparison(df1, df2, trim_percent=90,
#                                     title='Comparison with 10% Trimming')
#     plt.show()