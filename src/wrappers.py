"""Stable, high-level entry points for the SomaLogic bridging toolkit."""

from __future__ import annotations

from collections.abc import Collection, Mapping
from os import PathLike
from pathlib import Path
from typing import Any

from .anml_reference import create_anml_reference
from .bridging_factors import (
    create_bridging_factors,
    load_bridging_factors,
    save_bridging_factors,
    validate_bridging_factors,
)


def _load_factor_sources(
    factors_by_sample_matrix: Mapping[str, Mapping[str, Any] | str | PathLike[str]],
) -> dict[str, dict[str, Any]]:
    if not isinstance(factors_by_sample_matrix, Mapping) or not factors_by_sample_matrix:
        raise ValueError(
            "factors_by_sample_matrix must be a non-empty mapping from workbook "
            "sample-matrix values to factor manifests or JSON paths."
        )

    loaded: dict[str, dict[str, Any]] = {}
    for sample_matrix, factor_source in factors_by_sample_matrix.items():
        if not isinstance(sample_matrix, str) or not sample_matrix.strip():
            raise ValueError("Every sample-matrix key must be a non-empty string.")
        if isinstance(factor_source, (str, PathLike)):
            manifest = load_bridging_factors(factor_source)
        elif isinstance(factor_source, Mapping):
            manifest = dict(factor_source)
            validate_bridging_factors(manifest)
        else:
            raise TypeError(
                f"Factor source for {sample_matrix!r} must be a manifest or JSON path."
            )
        loaded[sample_matrix] = manifest
    return loaded


def apply_bridging_to_plate(
    scanner_input_path: str | PathLike[str],
    workbook_path: str | PathLike[str],
    factors_by_sample_matrix: Mapping[
        str, Mapping[str, Any] | str | PathLike[str]
    ],
    output_path: str | PathLike[str],
    *,
    excluded_aptamers: Collection[str] | None = None,
    strict_parameters: bool = True,
    overwrite: bool = False,
    verify_source_labels: bool = True,
) -> dict[str, Any]:
    """Apply fitted factors to selected samples in a SomaLogic scanner plate.

    ``factors_by_sample_matrix`` maps exact ``Sample Matrix`` values from the
    workbook (for example ``"Streck"`` and ``"NDS"``) to manifests returned by
    :func:`create_bridging_factors`, or to paths written by
    ``save_bridging_factors``.

    All manifests in one call must use the same scanner-compatible method:
    ``mad_log``, ``median_ratio``, or ``log_median_shift``. ``mad_linear`` is a
    post-normalization method and is deliberately rejected here.

    When ``excluded_aptamers`` is omitted, the evaluated project behavior is
    retained: MADMED excludes the 12 HCE aptamers and median methods exclude
    none. Pass an explicit collection to make another policy reproducible.
    Existing non-empty output directories are rejected unless ``overwrite`` is
    true; even then, unrelated files are never deleted. By default each
    manifest's source label must exactly equal its ``Sample Matrix`` mapping
    key, which prevents applying a factor in the reverse direction. Set
    ``verify_source_labels=False`` only when legacy labels cannot be aligned.
    """
    for name, value in (
        ("strict_parameters", strict_parameters),
        ("overwrite", overwrite),
        ("verify_source_labels", verify_source_labels),
    ):
        if not isinstance(value, bool):
            raise TypeError(f"{name} must be a boolean.")

    manifests = _load_factor_sources(factors_by_sample_matrix)
    if verify_source_labels:
        mismatches = {
            sample_matrix: manifest["direction"]["source"]
            for sample_matrix, manifest in manifests.items()
            if manifest["direction"]["source"] != sample_matrix
        }
        if mismatches:
            raise ValueError(
                "Factor source labels must exactly match their workbook Sample "
                f"Matrix keys; mismatches: {mismatches}. Set "
                "verify_source_labels=False only for intentionally different "
                "legacy labels."
            )
    methods = {manifest["method"] for manifest in manifests.values()}
    if methods == {"mad_log"}:
        mode = "mad_log"
    elif methods and methods <= {"median_ratio", "log_median_shift"}:
        mode = "median"
    elif "mad_linear" in methods:
        raise ValueError(
            "mad_linear factors are intended for already-normalized values, not "
            "scanner plate files."
        )
    else:
        raise ValueError(
            "A plate call cannot mix MADMED and median-only factor methods; "
            f"received {sorted(methods)}."
        )

    scanner_directory = Path(scanner_input_path)
    workbook = Path(workbook_path)
    destination = Path(output_path)
    if scanner_directory.resolve() == destination.resolve():
        raise ValueError("output_path must differ from scanner_input_path.")
    if workbook.resolve() == destination.resolve():
        raise ValueError("output_path must be a directory, not the workbook path.")
    if (destination / workbook.name).resolve() == workbook.resolve():
        raise ValueError(
            "output_path would overwrite the source workbook; choose a different "
            "directory."
        )
    if destination.exists() and not destination.is_dir():
        raise NotADirectoryError(f"output_path is not a directory: {destination}")
    if destination.exists() and any(destination.iterdir()) and not overwrite:
        raise FileExistsError(
            f"Refusing to write into non-empty output directory: {destination}"
        )

    if isinstance(excluded_aptamers, (str, bytes)):
        raise TypeError(
            "excluded_aptamers must be a collection of names, not a string."
        )
    if excluded_aptamers is not None:
        exclusions = tuple(excluded_aptamers)
        if not all(isinstance(value, str) and value for value in exclusions):
            raise ValueError("excluded_aptamers must contain non-empty strings.")
    else:
        exclusions = None

    # Import lazily so factor/reference creation does not require plotting or
    # scanner-file dependencies at package import time.
    from .altering_plate_multiple_collection_methods import (
        DEFAULT_HCE_APTAMERS,
        alter_scanner_files_MAD_bridging,
        alter_scanner_files_median_bridging,
    )

    if exclusions is None:
        exclusions = tuple(DEFAULT_HCE_APTAMERS) if mode == "mad_log" else ()

    if mode == "mad_log":
        run_manifest = alter_scanner_files_MAD_bridging(
            scanner_directory,
            workbook,
            manifests,
            destination,
            excluded_aptamers=exclusions,
            strict_parameters=strict_parameters,
        )
    else:
        run_manifest = alter_scanner_files_median_bridging(
            scanner_directory,
            workbook,
            manifests,
            destination,
            excluded_aptamers=exclusions,
            strict_parameters=strict_parameters,
        )

    run_manifest["factor_manifests"] = {
        sample_matrix: {
            "method": manifest["method"],
            "direction": manifest["direction"],
            "aptamer_count": manifest["aptamer_count"],
            "metadata": manifest["metadata"],
        }
        for sample_matrix, manifest in manifests.items()
    }
    run_manifest["factor_methods"] = sorted(methods)
    run_manifest["source_labels_verified"] = verify_source_labels
    return run_manifest


__all__ = [
    "apply_bridging_to_plate",
    "create_anml_reference",
    "create_bridging_factors",
    "load_bridging_factors",
    "save_bridging_factors",
    "validate_bridging_factors",
]
