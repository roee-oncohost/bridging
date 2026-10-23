import pandas as pd
import numpy as np


def update_reference_from_matrix(df_ref, df_values):
    """
    Update reference statistics in df_ref using distributions from df_values.

    Parameters
    ----------
    df_ref : pd.DataFrame
        Must contain columns:
        ['SeqId', 'RefMedian', 'RefSD', 'RefLog10Median', 'RefLog10MAD']
        plus any additional metadata columns.

    df_values : pd.DataFrame
        Wide matrix where:
        - columns = SeqId
        - rows = samples (e.g., patients)
        - values = measurements

    Returns
    -------
    pd.DataFrame
        Updated df_ref with recomputed reference statistics.
    """

    # =========================
    # 1. Prepare values (handle log safety)
    # =========================

    # Replace non-positive values with NaN (to avoid log10 issues)
    df_positive = df_values.where(df_values > 0)

    # =========================
    # 2. Compute statistics per SeqId (column-wise)
    # =========================

    # Linear stats
    ref_median = df_positive.median(axis=0)
    ref_sd = df_positive.std(axis=0, ddof=1)

    # Log stats
    log_df = np.log10(df_positive)
    ref_log_median = log_df.median(axis=0)
    ref_log_mad = (log_df.sub(ref_log_median, axis=1).abs()).median(axis=0)

    # =========================
    # 3. Build new reference table
    # =========================

    df_new_ref = pd.DataFrame({
        "SeqId": df_values.columns.astype(str),
        "RefMedian": ref_median.values,
        "RefSD": ref_sd.values,
        "RefLog10Median": ref_log_median.values,
        "RefLog10MAD": ref_log_mad.values
    })

    # =========================
    # 4. Update original reference dataframe
    # =========================

    df_ref = df_ref.copy()
    df_ref_columns = df_ref.columns
    # Ensure matching types
    df_ref["SeqId"] = df_ref["SeqId"].astype(str)
    df_new_ref["SeqId"] = df_new_ref["SeqId"].astype(str)

    # Align on SeqId and update only relevant columns
    df_ref = df_ref.set_index("SeqId")
    df_new_ref = df_new_ref.set_index("SeqId")

    df_ref.update(df_new_ref)

    return df_ref.reset_index()[df_ref_columns]
