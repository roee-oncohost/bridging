"""Fit, validate, save, and load SomaLogic bridging-factor manifests.

Every method maps a ``source`` distribution to a ``target`` distribution.
For the current project that usually means Streck -> EDTA, but neutral names
make the transformation direction explicit in every saved artifact.
"""

from __future__ import annotations

import json
import math
import os
import tempfile
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd


FACTOR_SCHEMA = "somalogic-bridging-factors"
FACTOR_SCHEMA_VERSION = 1
SUPPORTED_METHODS = frozenset(
    {"mad_log", "median_ratio", "log_median_shift", "mad_linear"}
)


def _mad(values: np.ndarray) -> float:
    median = float(np.median(values))
    return float(np.median(np.abs(values - median)))


def _numeric_values(
    series: pd.Series,
    *,
    dataframe_name: str,
    aptamer: str,
    require_positive: bool,
) -> np.ndarray:
    """Return finite, non-null numeric RFUs and reject invalid input."""
    if pd.api.types.is_bool_dtype(series.dtype):
        raise TypeError(
            f"{dataframe_name}[{aptamer!r}] is boolean; RFUs must be numeric."
        )
    non_null = series.dropna()
    if non_null.map(lambda value: isinstance(value, (bool, np.bool_))).any():
        raise TypeError(
            f"{dataframe_name}[{aptamer!r}] contains boolean RFUs; "
            "RFUs must be numeric."
        )
    numeric = pd.to_numeric(series, errors="coerce")
    invalid_numeric = series.notna() & numeric.isna()
    if invalid_numeric.any():
        examples = series.loc[invalid_numeric].astype(str).head(3).tolist()
        raise TypeError(
            f"{dataframe_name}[{aptamer!r}] contains non-numeric RFUs; "
            f"examples: {examples}."
        )
    values = numeric.dropna().to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(
            f"{dataframe_name}[{aptamer!r}] contains infinite RFUs."
        )
    if require_positive and np.any(values <= 0):
        count = int(np.count_nonzero(values <= 0))
        raise ValueError(
            f"{dataframe_name}[{aptamer!r}] contains {count} non-positive "
            "RFU value(s); this method requires strictly positive RFUs."
        )
    return values


def _resolve_aptamers(
    source_df: pd.DataFrame,
    target_df: pd.DataFrame,
    aptamers: Sequence[str] | None,
) -> list[str]:
    for name, dataframe in (("source_df", source_df), ("target_df", target_df)):
        if not isinstance(dataframe, pd.DataFrame):
            raise TypeError(f"{name} must be a pandas DataFrame.")
        if dataframe.columns.has_duplicates:
            duplicated = dataframe.columns[dataframe.columns.duplicated()].tolist()
            raise ValueError(f"{name} contains duplicate columns: {duplicated[:5]}.")

    if aptamers is None:
        target_columns = set(target_df.columns)
        resolved = [
            column
            for column in source_df.columns
            if isinstance(column, str)
            and column
            and column[0].isdigit()
            and column in target_columns
        ]
        if not resolved:
            raise ValueError(
                "No shared digit-prefixed aptamer columns were found. Pass "
                "aptamers explicitly if another column convention is used."
            )
        return resolved

    if isinstance(aptamers, (str, bytes)):
        raise TypeError("aptamers must be a sequence of names, not a string.")
    resolved = list(aptamers)
    if not resolved:
        raise ValueError("aptamers must not be empty.")
    if not all(isinstance(aptamer, str) and aptamer for aptamer in resolved):
        raise TypeError("Every aptamer name must be a non-empty string.")
    if len(resolved) != len(set(resolved)):
        raise ValueError("aptamers contains duplicate names.")
    missing_source = [aptamer for aptamer in resolved if aptamer not in source_df]
    missing_target = [aptamer for aptamer in resolved if aptamer not in target_df]
    if missing_source or missing_target:
        raise ValueError(
            "Requested aptamers are missing. "
            f"source missing={missing_source[:10]}; "
            f"target missing={missing_target[:10]}."
        )
    return resolved


def _label(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string.")
    return value.strip()


def _metadata(value: Mapping[str, Any] | None) -> dict[str, Any]:
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise TypeError("metadata must be a mapping.")
    result = dict(value)
    try:
        encoded = json.dumps(result, allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise ValueError("metadata must contain JSON-serializable values.") from exc
    if json.loads(encoded) != result:
        raise ValueError(
            "metadata must use JSON-native values and string object keys."
        )
    return result


def _fit_one(
    source_values: np.ndarray,
    target_values: np.ndarray,
    method: str,
) -> dict[str, Any]:
    counts = {
        "source_n": int(source_values.size),
        "target_n": int(target_values.size),
    }
    if method in {"mad_log", "log_median_shift"}:
        source_fit = np.log(source_values)
        target_fit = np.log(target_values)
    else:
        source_fit = source_values
        target_fit = target_values

    source_median = float(np.median(source_fit))
    target_median = float(np.median(target_fit))

    if method in {"mad_log", "mad_linear"}:
        source_mad = _mad(source_fit)
        target_mad = _mad(target_fit)
        if source_mad == 0:
            raise ValueError("source MAD is zero, so a scale transform cannot be fit")
        a = target_mad / source_mad
        b = target_median - a * source_median
        return {
            **counts,
            "source_median": source_median,
            "source_mad": source_mad,
            "target_median": target_median,
            "target_mad": target_mad,
            "a": float(a),
            "b": float(b),
        }
    if method == "median_ratio":
        if source_median == 0:
            raise ValueError("source median is zero, so a ratio cannot be fit")
        return {
            **counts,
            "source_median": source_median,
            "target_median": target_median,
            "factor": float(target_median / source_median),
        }
    if method == "log_median_shift":
        log_shift = target_median - source_median
        try:
            multiplier = math.exp(log_shift)
        except OverflowError as exc:
            raise ValueError(
                "log-median shift is too large to produce a finite multiplier"
            ) from exc
        return {
            **counts,
            "source_log_median": source_median,
            "target_log_median": target_median,
            "log_shift": float(log_shift),
            "factor": float(multiplier),
        }
    raise AssertionError(f"Unhandled method: {method}")


def create_bridging_factors(
    source_df: pd.DataFrame,
    target_df: pd.DataFrame,
    *,
    method: str = "mad_log",
    aptamers: Sequence[str] | None = None,
    source_label: str = "source",
    target_label: str = "target",
    min_valid_samples: int = 3,
    metadata: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Fit per-aptamer factors that map ``source_df`` to ``target_df``.

    The inputs are wide matrices with samples in rows and aptamer SeqIds in
    columns. Their sample counts do not need to match. Log-domain methods
    reject non-positive RFUs rather than silently excluding observations.
    """
    if not isinstance(method, str) or method not in SUPPORTED_METHODS:
        raise ValueError(
            f"Unsupported method {method!r}; expected one of "
            f"{sorted(SUPPORTED_METHODS)}."
        )
    if isinstance(min_valid_samples, bool) or not isinstance(min_valid_samples, int):
        raise TypeError("min_valid_samples must be an integer.")
    if min_valid_samples < 1:
        raise ValueError("min_valid_samples must be at least 1.")

    aptamer_list = _resolve_aptamers(source_df, target_df, aptamers)
    require_positive = method in {"mad_log", "log_median_shift"}
    factors: dict[str, dict[str, Any]] = {}
    for aptamer in aptamer_list:
        source_values = _numeric_values(
            source_df[aptamer],
            dataframe_name="source_df",
            aptamer=aptamer,
            require_positive=require_positive,
        )
        target_values = _numeric_values(
            target_df[aptamer],
            dataframe_name="target_df",
            aptamer=aptamer,
            require_positive=require_positive,
        )
        if source_values.size < min_valid_samples:
            raise ValueError(
                f"source_df[{aptamer!r}] has {source_values.size} valid RFUs; "
                f"at least {min_valid_samples} are required."
            )
        if target_values.size < min_valid_samples:
            raise ValueError(
                f"target_df[{aptamer!r}] has {target_values.size} valid RFUs; "
                f"at least {min_valid_samples} are required."
            )
        try:
            factors[aptamer] = _fit_one(source_values, target_values, method)
        except ValueError as exc:
            raise ValueError(f"Cannot fit aptamer {aptamer!r}: {exc}.") from exc

    manifest: dict[str, Any] = {
        "schema": FACTOR_SCHEMA,
        "schema_version": FACTOR_SCHEMA_VERSION,
        "method": method,
        "direction": {
            "source": _label(source_label, "source_label"),
            "target": _label(target_label, "target_label"),
        },
        "fit": {
            "source_sample_count": int(len(source_df)),
            "target_sample_count": int(len(target_df)),
            "min_valid_samples": min_valid_samples,
        },
        "aptamers": aptamer_list,
        "aptamer_count": len(aptamer_list),
        "factors": factors,
        "metadata": _metadata(metadata),
    }
    validate_bridging_factors(manifest)
    return manifest


def _finite(value: Any, location: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{location} must be a number.")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{location} must be finite.")
    return result


def _positive_int(value: Any, location: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{location} must be a positive integer.")
    return value


def _validate_mad_factor(
    aptamer: str,
    factor: Mapping[str, Any],
) -> None:
    source_median = _finite(
        factor.get("source_median"), f"factors[{aptamer!r}].source_median"
    )
    source_mad = _finite(
        factor.get("source_mad"), f"factors[{aptamer!r}].source_mad"
    )
    target_median = _finite(
        factor.get("target_median"), f"factors[{aptamer!r}].target_median"
    )
    target_mad = _finite(
        factor.get("target_mad"), f"factors[{aptamer!r}].target_mad"
    )
    a = _finite(factor.get("a"), f"factors[{aptamer!r}].a")
    b = _finite(factor.get("b"), f"factors[{aptamer!r}].b")
    if source_mad <= 0 or target_mad < 0:
        raise ValueError(f"factors[{aptamer!r}] has invalid MAD values.")
    expected_a = target_mad / source_mad
    expected_b = target_median - expected_a * source_median
    if not math.isclose(a, expected_a, rel_tol=1e-12, abs_tol=1e-15):
        raise ValueError(f"factors[{aptamer!r}].a is inconsistent with its MADs.")
    if not math.isclose(b, expected_b, rel_tol=1e-12, abs_tol=1e-15):
        raise ValueError(
            f"factors[{aptamer!r}].b is inconsistent with its medians and MADs."
        )


def _validate_ratio_factor(
    aptamer: str,
    factor: Mapping[str, Any],
) -> None:
    source_median = _finite(
        factor.get("source_median"), f"factors[{aptamer!r}].source_median"
    )
    target_median = _finite(
        factor.get("target_median"), f"factors[{aptamer!r}].target_median"
    )
    multiplier = _finite(factor.get("factor"), f"factors[{aptamer!r}].factor")
    if source_median <= 0 or target_median <= 0:
        raise ValueError(
            f"factors[{aptamer!r}] medians must be positive for median_ratio."
        )
    if not math.isclose(
        multiplier,
        target_median / source_median,
        rel_tol=1e-12,
        abs_tol=1e-15,
    ):
        raise ValueError(f"factors[{aptamer!r}].factor is inconsistent with medians.")


def _validate_log_shift_factor(
    aptamer: str,
    factor: Mapping[str, Any],
) -> None:
    source_median = _finite(
        factor.get("source_log_median"),
        f"factors[{aptamer!r}].source_log_median",
    )
    target_median = _finite(
        factor.get("target_log_median"),
        f"factors[{aptamer!r}].target_log_median",
    )
    shift = _finite(
        factor.get("log_shift"), f"factors[{aptamer!r}].log_shift"
    )
    multiplier = _finite(factor.get("factor"), f"factors[{aptamer!r}].factor")
    expected_shift = target_median - source_median
    if not math.isclose(shift, expected_shift, rel_tol=1e-12, abs_tol=1e-15):
        raise ValueError(f"factors[{aptamer!r}].log_shift is inconsistent.")
    try:
        expected_multiplier = math.exp(expected_shift)
    except OverflowError as exc:
        raise ValueError(
            f"factors[{aptamer!r}].log_shift produces an infinite multiplier."
        ) from exc
    if not math.isclose(
        multiplier, expected_multiplier, rel_tol=1e-12, abs_tol=1e-15
    ):
        raise ValueError(f"factors[{aptamer!r}].factor is inconsistent.")


def validate_bridging_factors(manifest: Mapping[str, Any]) -> None:
    """Validate the public factor-manifest schema, raising ``ValueError``."""
    if not isinstance(manifest, Mapping):
        raise ValueError("Factor manifest must be a mapping.")
    required = {
        "schema",
        "schema_version",
        "method",
        "direction",
        "fit",
        "aptamers",
        "aptamer_count",
        "factors",
        "metadata",
    }
    missing = sorted(required - set(manifest))
    if missing:
        raise ValueError(f"Factor manifest is missing keys: {missing}.")
    if manifest["schema"] != FACTOR_SCHEMA:
        raise ValueError(f"Unsupported factor schema: {manifest['schema']!r}.")
    if (
        isinstance(manifest["schema_version"], bool)
        or not isinstance(manifest["schema_version"], int)
        or manifest["schema_version"] != FACTOR_SCHEMA_VERSION
    ):
        raise ValueError(
            f"Unsupported factor schema version: {manifest['schema_version']!r}."
        )

    method = manifest["method"]
    if not isinstance(method, str) or method not in SUPPORTED_METHODS:
        raise ValueError(f"Unsupported factor method: {method!r}.")
    direction = manifest["direction"]
    if not isinstance(direction, Mapping):
        raise ValueError("direction must be a mapping.")
    source_label = _label(direction.get("source"), "direction.source")
    target_label = _label(direction.get("target"), "direction.target")
    if source_label == target_label:
        raise ValueError("direction.source and direction.target must differ.")

    fit = manifest["fit"]
    if not isinstance(fit, Mapping):
        raise ValueError("fit must be a mapping.")
    for key in ("source_sample_count", "target_sample_count", "min_valid_samples"):
        _positive_int(fit.get(key), f"fit.{key}")

    aptamers = manifest["aptamers"]
    factors = manifest["factors"]
    if not isinstance(aptamers, list) or not aptamers:
        raise ValueError("aptamers must be a non-empty list.")
    if not all(isinstance(aptamer, str) and aptamer for aptamer in aptamers):
        raise ValueError("aptamers must contain non-empty strings.")
    if len(aptamers) != len(set(aptamers)):
        raise ValueError("aptamers contains duplicate names.")
    if (
        isinstance(manifest["aptamer_count"], bool)
        or not isinstance(manifest["aptamer_count"], int)
        or manifest["aptamer_count"] != len(aptamers)
    ):
        raise ValueError("aptamer_count does not match aptamers.")
    if not isinstance(factors, Mapping) or set(factors) != set(aptamers):
        raise ValueError("factors keys must exactly match aptamers.")

    for aptamer in aptamers:
        factor = factors[aptamer]
        if not isinstance(factor, Mapping):
            raise ValueError(f"factors[{aptamer!r}] must be a mapping.")
        source_n = _positive_int(
            factor.get("source_n"), f"factors[{aptamer!r}].source_n"
        )
        target_n = _positive_int(
            factor.get("target_n"), f"factors[{aptamer!r}].target_n"
        )
        if not fit["min_valid_samples"] <= source_n <= fit["source_sample_count"]:
            raise ValueError(
                f"factors[{aptamer!r}].source_n is inconsistent with fit counts."
            )
        if not fit["min_valid_samples"] <= target_n <= fit["target_sample_count"]:
            raise ValueError(
                f"factors[{aptamer!r}].target_n is inconsistent with fit counts."
            )
        if method in {"mad_log", "mad_linear"}:
            _validate_mad_factor(aptamer, factor)
        elif method == "median_ratio":
            _validate_ratio_factor(aptamer, factor)
        else:
            _validate_log_shift_factor(aptamer, factor)

    try:
        _metadata(manifest["metadata"])
        json.dumps(dict(manifest), allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise ValueError("Factor manifest is not JSON serializable.") from exc


def save_bridging_factors(
    manifest: Mapping[str, Any],
    path: str | os.PathLike[str],
    *,
    overwrite: bool = False,
    indent: int = 2,
) -> Path:
    """Validate and atomically write a factor manifest as UTF-8 JSON."""
    validate_bridging_factors(manifest)
    if not isinstance(overwrite, bool):
        raise TypeError("overwrite must be a boolean.")
    if isinstance(indent, bool) or not isinstance(indent, int) or indent < 0:
        raise ValueError("indent must be a non-negative integer.")
    destination = Path(path)
    if destination.is_dir():
        raise IsADirectoryError(f"Factor path is a directory: {destination}")
    if destination.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing file: {destination}")
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
                    f"Refusing to overwrite file created by another process: "
                    f"{destination}"
                ) from exc
            else:
                os.close(descriptor)
                reserved_destination = True
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            newline="\n",
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            json.dump(dict(manifest), handle, indent=indent, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        temporary_path.replace(destination)
    except Exception:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
        if reserved_destination:
            destination.unlink(missing_ok=True)
        raise
    return destination


def load_bridging_factors(path: str | os.PathLike[str]) -> dict[str, Any]:
    """Load and validate a UTF-8 JSON factor manifest."""
    source = Path(path)
    with source.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)
    validate_bridging_factors(manifest)
    return manifest


__all__ = [
    "FACTOR_SCHEMA",
    "FACTOR_SCHEMA_VERSION",
    "SUPPORTED_METHODS",
    "create_bridging_factors",
    "load_bridging_factors",
    "save_bridging_factors",
    "validate_bridging_factors",
]

