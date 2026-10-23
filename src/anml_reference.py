"""Create SomaLogic ANML reference tables from a sample RFU matrix.

The public entry point, :func:`create_anml_reference`, replaces the four
distribution columns used by ANML while preserving every other field and the
original row/column order of the supplied reference table.
"""

from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import tempfile
from collections.abc import Mapping
from typing import Any, Literal

import numpy as np
import pandas as pd


ANML_STAT_COLUMNS = (
    "RefMedian",
    "RefSD",
    "RefLog10Median",
    "RefLog10MAD",
)
ANML_REFERENCE_MANIFEST_SCHEMA = "somalogic-anml-reference-manifest/v1"

MissingPolicy = Literal["error", "keep"]
ExtraPolicy = Literal["error", "ignore"]
InsufficientPolicy = Literal["error", "keep"]


def _read_reference(
    reference: pd.DataFrame | str | Path,
    seq_id_column: str,
) -> tuple[pd.DataFrame, str]:
    if isinstance(reference, pd.DataFrame):
        return reference.copy(), "<dataframe>"
    if isinstance(reference, (str, Path)):
        reference_path = Path(reference)
        if not reference_path.is_file():
            raise FileNotFoundError(f"ANML reference file does not exist: {reference_path}")
        return (
            pd.read_csv(
                reference_path,
                sep="\t",
                dtype=str,
                keep_default_na=False,
            ),
            str(reference_path),
        )
    raise TypeError("reference must be a pandas DataFrame or a path to a tab-separated file")


def _normalise_identifier(value: Any) -> str:
    return str(value).strip()


def _validate_policy(name: str, value: str, allowed: tuple[str, ...]) -> None:
    if value not in allowed:
        choices = ", ".join(repr(choice) for choice in allowed)
        raise ValueError(f"{name} must be one of {choices}; got {value!r}")


def _validate_reference(reference_df: pd.DataFrame, seq_id_column: str) -> list[str]:
    duplicate_columns = reference_df.columns[reference_df.columns.duplicated()].tolist()
    if duplicate_columns:
        raise ValueError(f"ANML reference contains duplicate columns: {duplicate_columns[:10]}")
    required = {seq_id_column, *ANML_STAT_COLUMNS}
    missing_columns = sorted(required - set(reference_df.columns))
    if missing_columns:
        raise ValueError(f"ANML reference is missing required columns: {missing_columns}")
    if reference_df.empty:
        raise ValueError("ANML reference contains no rows")
    if reference_df[seq_id_column].isna().any():
        raise ValueError(f"ANML reference contains null {seq_id_column} values")

    seq_ids = [_normalise_identifier(value) for value in reference_df[seq_id_column]]
    empty_positions = [index for index, value in enumerate(seq_ids) if not value]
    if empty_positions:
        raise ValueError(
            f"ANML reference contains blank {seq_id_column} values at row positions "
            f"{empty_positions[:10]}"
        )

    duplicates = pd.Series(seq_ids)[pd.Series(seq_ids).duplicated(keep=False)].unique().tolist()
    if duplicates:
        raise ValueError(
            f"ANML reference contains duplicate {seq_id_column} values: {duplicates[:10]}"
        )
    return seq_ids


def _validate_retained_statistics(
    reference_df: pd.DataFrame,
    reference_seq_ids: list[str],
    retained_seq_ids: list[str],
) -> None:
    """Validate existing statistics that will be retained by a keep policy."""
    if not retained_seq_ids:
        return
    row_by_seq_id = {
        seq_id: row_position for row_position, seq_id in enumerate(reference_seq_ids)
    }
    positions = [row_by_seq_id[seq_id] for seq_id in retained_seq_ids]
    for column in ANML_STAT_COLUMNS:
        original = reference_df.iloc[positions][column]
        numeric = pd.to_numeric(original, errors="coerce")
        invalid = numeric.isna() | ~np.isfinite(numeric)
        if column == "RefMedian":
            invalid |= numeric <= 0
        elif column in {"RefSD", "RefLog10MAD"}:
            invalid |= numeric < 0
        if invalid.any():
            invalid_seq_ids = [
                retained_seq_ids[position]
                for position, is_invalid in enumerate(invalid.to_numpy())
                if is_invalid
            ]
            raise ValueError(
                f"Retained ANML statistics contain invalid {column} values for "
                f"SeqIds: {invalid_seq_ids[:10]}"
            )


def _normalised_value_columns(values: pd.DataFrame) -> dict[str, Any]:
    if not isinstance(values, pd.DataFrame):
        raise TypeError(
            "values must be a pandas DataFrame with samples in rows and SeqIds in columns"
        )
    if len(values.index) == 0:
        raise ValueError("values contains no samples")
    if len(values.columns) == 0:
        raise ValueError("values contains no columns")

    column_lookup: dict[str, Any] = {}
    duplicate_ids: list[str] = []
    for original_column in values.columns:
        normalised = _normalise_identifier(original_column)
        if not normalised:
            raise ValueError("values contains a blank column name")
        if normalised in column_lookup:
            duplicate_ids.append(normalised)
        else:
            column_lookup[normalised] = original_column
    if duplicate_ids:
        raise ValueError(
            "values contains duplicate SeqId columns after string normalisation: "
            f"{sorted(set(duplicate_ids))[:10]}"
        )
    return column_lookup


def _numeric_values(
    values: pd.DataFrame,
    column_lookup: Mapping[str, Any],
    seq_ids: list[str],
) -> pd.DataFrame:
    original = values.loc[:, [column_lookup[seq_id] for seq_id in seq_ids]].copy()
    original.columns = seq_ids
    numeric = original.apply(pd.to_numeric, errors="coerce")
    invalid = original.notna() & numeric.isna()
    infinite = numeric.notna() & ~np.isfinite(numeric)
    invalid_by_seq_id = (invalid | infinite).sum(axis=0)
    invalid_counts = {
        seq_id: int(count) for seq_id, count in invalid_by_seq_id.items() if count
    }
    numeric = numeric.replace([np.inf, -np.inf], np.nan)

    if invalid_counts:
        preview = dict(list(invalid_counts.items())[:10])
        raise ValueError(
            "values contains non-numeric or infinite measurements for matched SeqIds; "
            f"invalid cell counts: {preview}"
        )
    return numeric


def create_anml_reference(
    reference: pd.DataFrame | str | Path,
    values: pd.DataFrame,
    output_path: str | Path | None = None,
    *,
    seq_id_column: str = "SeqId",
    min_positive_samples: int = 3,
    missing_policy: MissingPolicy = "error",
    extra_policy: ExtraPolicy = "ignore",
    insufficient_policy: InsufficientPolicy = "error",
    metadata: Mapping[str, Any] | None = None,
    overwrite: bool = False,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Create an ANML reference from a wide sample-by-SeqId RFU matrix.

    Non-positive measurements are excluded because the ANML reference contains
    log10 statistics. Missing reference SeqIds and SeqIds with too few positive
    observations can either raise an error (the safe default) or retain their
    original reference statistics. Columns in ``values`` that are not reference
    SeqIds, such as sample metadata, can be ignored or rejected.

    Parameters
    ----------
    reference:
        Existing ANML reference DataFrame or path to a tab-separated reference.
    values:
        Wide DataFrame with samples in rows and SeqIds in columns.
    output_path:
        Optional path for the new tab-separated ANML reference file.
    seq_id_column:
        Name of the reference identifier column. Values are matched to the
        string-normalized columns of ``values``.
    min_positive_samples:
        Minimum finite values greater than zero required for each updated SeqId.
        Must be at least two so ``RefSD`` with ``ddof=1`` is defined.
    missing_policy:
        ``"error"`` rejects reference SeqIds absent from ``values``; ``"keep"``
        retains their original statistics.
    extra_policy:
        ``"error"`` rejects non-reference columns; ``"ignore"`` treats them as
        metadata and records them in the manifest.
    insufficient_policy:
        ``"error"`` rejects matched SeqIds below ``min_positive_samples``;
        ``"keep"`` retains their original statistics.
    metadata:
        Optional JSON-native provenance mapping copied to the manifest. Object
        keys must be strings and values must round-trip through JSON unchanged.
    overwrite:
        Replace an existing ``output_path`` when true. The source reference can
        never be overwritten. Writes are atomic within the destination folder.

    Returns
    -------
    tuple[pandas.DataFrame, dict]
        The updated reference and a JSON-serializable audit manifest.
    """

    _validate_policy("missing_policy", missing_policy, ("error", "keep"))
    _validate_policy("extra_policy", extra_policy, ("error", "ignore"))
    _validate_policy("insufficient_policy", insufficient_policy, ("error", "keep"))
    if not isinstance(seq_id_column, str) or not seq_id_column.strip():
        raise ValueError("seq_id_column must be a non-empty string")
    if isinstance(min_positive_samples, bool) or not isinstance(min_positive_samples, int):
        raise TypeError("min_positive_samples must be an integer")
    if min_positive_samples < 2:
        raise ValueError("min_positive_samples must be at least 2 so RefSD is defined")
    if not isinstance(overwrite, bool):
        raise TypeError("overwrite must be a boolean")

    if metadata is not None and not isinstance(metadata, Mapping):
        raise TypeError("metadata must be a mapping")
    metadata_dict = dict(metadata or {})
    try:
        encoded_metadata = json.dumps(metadata_dict, allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "metadata must contain JSON-serializable finite values"
        ) from exc
    if json.loads(encoded_metadata) != metadata_dict:
        raise ValueError(
            "metadata must use JSON-native values and string object keys"
        )

    reference_df, reference_source = _read_reference(reference, seq_id_column)
    reference_columns = reference_df.columns.tolist()
    reference_seq_ids = _validate_reference(reference_df, seq_id_column)
    reference_set = set(reference_seq_ids)

    value_columns = _normalised_value_columns(values)
    value_seq_ids = list(value_columns)
    matched_seq_ids = [seq_id for seq_id in reference_seq_ids if seq_id in value_columns]
    missing_seq_ids = [seq_id for seq_id in reference_seq_ids if seq_id not in value_columns]
    extra_columns = [seq_id for seq_id in value_seq_ids if seq_id not in reference_set]

    if missing_seq_ids and missing_policy == "error":
        raise ValueError(
            f"values is missing {len(missing_seq_ids)} reference SeqIds: {missing_seq_ids[:10]}"
        )
    if extra_columns and extra_policy == "error":
        raise ValueError(
            f"values contains {len(extra_columns)} columns absent from the reference: "
            f"{extra_columns[:10]}"
        )
    if not matched_seq_ids:
        raise ValueError("values and the ANML reference have no SeqIds in common")

    numeric = _numeric_values(values, value_columns, matched_seq_ids)
    positive = numeric.where(numeric > 0)
    positive_counts = positive.count().astype(int)
    insufficient_seq_ids = [
        seq_id for seq_id in matched_seq_ids if positive_counts[seq_id] < min_positive_samples
    ]
    if insufficient_seq_ids and insufficient_policy == "error":
        preview = {seq_id: int(positive_counts[seq_id]) for seq_id in insufficient_seq_ids[:10]}
        raise ValueError(
            f"{len(insufficient_seq_ids)} SeqIds have fewer than {min_positive_samples} "
            f"positive samples: {preview}"
        )

    updated_seq_ids = [seq_id for seq_id in matched_seq_ids if seq_id not in insufficient_seq_ids]
    retained_seq_ids = missing_seq_ids + insufficient_seq_ids
    _validate_retained_statistics(
        reference_df,
        reference_seq_ids,
        retained_seq_ids,
    )
    positive_for_update = positive.loc[:, updated_seq_ids]
    log_values = np.log10(positive_for_update)
    log_medians = log_values.median(axis=0)

    statistics = pd.DataFrame(
        {
            "RefMedian": positive_for_update.median(axis=0),
            "RefSD": positive_for_update.std(axis=0, ddof=1),
            "RefLog10Median": log_medians,
            "RefLog10MAD": log_values.sub(log_medians, axis=1).abs().median(axis=0),
        }
    )
    statistics_array = statistics.to_numpy(dtype=float)
    if not np.isfinite(statistics_array).all():
        invalid_columns = statistics.columns[
            ~np.isfinite(statistics_array).all(axis=0)
        ].tolist()
        raise ValueError(
            "Computed ANML statistics contain non-finite values; affected "
            f"statistics: {invalid_columns}"
        )
    if (statistics["RefMedian"] <= 0).any():
        raise ValueError("Computed RefMedian values must be positive")
    if (statistics[["RefSD", "RefLog10MAD"]] < 0).any().any():
        raise ValueError("Computed RefSD and RefLog10MAD values must be non-negative")

    result = reference_df.copy()
    row_by_seq_id = {seq_id: row_index for row_index, seq_id in enumerate(reference_seq_ids)}
    row_positions = [row_by_seq_id[seq_id] for seq_id in updated_seq_ids]
    column_positions = [result.columns.get_loc(column) for column in ANML_STAT_COLUMNS]
    result.iloc[row_positions, column_positions] = statistics.loc[
        updated_seq_ids, list(ANML_STAT_COLUMNS)
    ].to_numpy()
    result = result.loc[:, reference_columns]

    written_path: str | None = None
    if output_path is not None:
        destination = Path(output_path)
        overwrites_source = (
            reference_source != "<dataframe>"
            and destination.resolve() == Path(reference_source).resolve()
        )
        if overwrites_source:
            raise ValueError("output_path must not overwrite the source ANML reference")
        if destination.is_dir():
            raise IsADirectoryError(f"output_path is a directory: {destination}")
        if destination.exists() and not overwrite:
            raise FileExistsError(
                f"Refusing to overwrite existing ANML reference: {destination}"
            )
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary_path: Path | None = None
        reserved_destination = False
        try:
            if not overwrite:
                try:
                    descriptor = os.open(
                        destination,
                        os.O_CREAT | os.O_EXCL | os.O_WRONLY,
                        0o600,
                    )
                except FileExistsError as exc:
                    raise FileExistsError(
                        "Refusing to overwrite ANML reference created by another "
                        f"process: {destination}"
                    ) from exc
                else:
                    os.close(descriptor)
                    reserved_destination = True
            with tempfile.NamedTemporaryFile(
                "w",
                encoding="utf-8",
                newline="",
                dir=destination.parent,
                prefix=f".{destination.name}.",
                suffix=".tmp",
                delete=False,
            ) as handle:
                temporary_path = Path(handle.name)
                result.to_csv(handle, sep="\t", index=False)
                handle.flush()
                os.fsync(handle.fileno())
            temporary_path.replace(destination)
        except Exception:
            if temporary_path is not None:
                temporary_path.unlink(missing_ok=True)
            if reserved_destination:
                destination.unlink(missing_ok=True)
            raise
        written_path = str(destination)

    positive_count_manifest = {
        seq_id: int(positive_counts[seq_id]) for seq_id in matched_seq_ids
    }
    nonpositive_count_manifest = {
        seq_id: int((numeric[seq_id].notna() & (numeric[seq_id] <= 0)).sum())
        for seq_id in matched_seq_ids
    }
    manifest: dict[str, Any] = {
        "schema": ANML_REFERENCE_MANIFEST_SCHEMA,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "method": "positive_rfu_anml_reference_statistics",
        "reference_source": reference_source,
        "output_path": written_path,
        "seq_id_column": seq_id_column,
        "statistics": {
            "RefMedian": "median of positive RFU",
            "RefSD": "sample standard deviation of positive RFU (ddof=1)",
            "RefLog10Median": "median of log10 positive RFU",
            "RefLog10MAD": "unscaled median absolute deviation of log10 positive RFU",
        },
        "policies": {
            "missing": missing_policy,
            "extra": extra_policy,
            "insufficient": insufficient_policy,
        },
        "min_positive_samples": min_positive_samples,
        "sample_count": int(len(values)),
        "reference_seq_id_count": int(len(reference_seq_ids)),
        "matched_seq_id_count": int(len(matched_seq_ids)),
        "updated_seq_id_count": int(len(updated_seq_ids)),
        "coverage_fraction": float(len(matched_seq_ids) / len(reference_seq_ids)),
        "updated_seq_ids": updated_seq_ids,
        "missing_seq_ids": missing_seq_ids,
        "extra_columns": extra_columns,
        "insufficient_seq_ids": insufficient_seq_ids,
        "positive_sample_counts": positive_count_manifest,
        "nonpositive_sample_counts": nonpositive_count_manifest,
        "metadata": metadata_dict,
    }
    json.dumps(manifest, allow_nan=False)
    return result, manifest


__all__ = [
    "ANML_REFERENCE_MANIFEST_SCHEMA",
    "ANML_STAT_COLUMNS",
    "create_anml_reference",
]
