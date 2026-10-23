import math

import pandas as pd
import pytest

from src.bridging_factors import (
    create_bridging_factors,
    load_bridging_factors,
    save_bridging_factors,
    validate_bridging_factors,
)


APTAMER = "10000-28"


@pytest.mark.parametrize(
    ("method", "source", "target", "expected"),
    [
        (
            "mad_log",
            [1.0, math.e, math.e**2],
            [math.e**2, math.e**4, math.e**6],
            {"a": 2.0, "b": 2.0},
        ),
        (
            "median_ratio",
            [1.0, 2.0, 3.0],
            [2.0, 4.0, 6.0],
            {"factor": 2.0},
        ),
        (
            "log_median_shift",
            [1.0, math.e, math.e**2],
            [math.e, math.e**2, math.e**3],
            {"factor": math.e, "log_shift": 1.0},
        ),
        (
            "mad_linear",
            [1.0, 2.0, 3.0],
            [5.0, 7.0, 9.0],
            {"a": 2.0, "b": 3.0},
        ),
    ],
)
def test_create_bridging_factors(method, source, target, expected):
    manifest = create_bridging_factors(
        pd.DataFrame({APTAMER: source}),
        pd.DataFrame({APTAMER: target}),
        method=method,
        source_label="Streck DS2",
        target_label="EDTA Plasma",
        metadata={"cohort": "unit-test"},
    )

    assert manifest["method"] == method
    assert manifest["direction"] == {
        "source": "Streck DS2",
        "target": "EDTA Plasma",
    }
    assert manifest["aptamers"] == [APTAMER]
    assert manifest["aptamer_count"] == 1
    for key, value in expected.items():
        assert manifest["factors"][APTAMER][key] == pytest.approx(value)
    validate_bridging_factors(manifest)


def test_log_method_rejects_non_positive_rfus():
    source = pd.DataFrame({APTAMER: [0.0, 1.0, 2.0]})
    target = pd.DataFrame({APTAMER: [1.0, 2.0, 3.0]})

    with pytest.raises(ValueError, match="non-positive"):
        create_bridging_factors(source, target, method="mad_log")


def test_manifest_round_trip_and_safe_overwrite(tmp_path):
    manifest = create_bridging_factors(
        pd.DataFrame({APTAMER: [1.0, 2.0, 3.0]}),
        pd.DataFrame({APTAMER: [2.0, 4.0, 6.0]}),
        method="median_ratio",
    )
    path = tmp_path / "factors.json"

    assert save_bridging_factors(manifest, path) == path
    assert load_bridging_factors(path) == manifest
    with pytest.raises(FileExistsError):
        save_bridging_factors(manifest, path)


def test_validator_rejects_tampered_factor():
    manifest = create_bridging_factors(
        pd.DataFrame({APTAMER: [1.0, 2.0, 3.0]}),
        pd.DataFrame({APTAMER: [2.0, 4.0, 6.0]}),
        method="median_ratio",
    )
    manifest["factors"][APTAMER]["factor"] = 99.0

    with pytest.raises(ValueError, match="inconsistent"):
        validate_bridging_factors(manifest)


@pytest.mark.parametrize("field", ["schema_version", "aptamer_count"])
def test_validator_rejects_boolean_integer_fields(field):
    manifest = create_bridging_factors(
        pd.DataFrame({APTAMER: [1.0, 2.0, 3.0]}),
        pd.DataFrame({APTAMER: [2.0, 4.0, 6.0]}),
        method="median_ratio",
    )
    manifest[field] = True

    with pytest.raises(ValueError):
        validate_bridging_factors(manifest)


def test_validator_rejects_inconsistent_per_aptamer_counts():
    manifest = create_bridging_factors(
        pd.DataFrame({APTAMER: [1.0, 2.0, 3.0]}),
        pd.DataFrame({APTAMER: [2.0, 4.0, 6.0]}),
        method="median_ratio",
    )
    manifest["factors"][APTAMER]["source_n"] = 4

    with pytest.raises(ValueError, match="source_n is inconsistent"):
        validate_bridging_factors(manifest)


def test_boolean_rfus_are_rejected_even_in_object_columns():
    source = pd.DataFrame({APTAMER: pd.Series([1.0, True, 3.0], dtype=object)})
    target = pd.DataFrame({APTAMER: [2.0, 4.0, 6.0]})

    with pytest.raises(TypeError, match="boolean RFUs"):
        create_bridging_factors(source, target, method="median_ratio")


def test_overwrite_true_replaces_existing_manifest(tmp_path):
    manifest = create_bridging_factors(
        pd.DataFrame({APTAMER: [1.0, 2.0, 3.0]}),
        pd.DataFrame({APTAMER: [2.0, 4.0, 6.0]}),
        method="median_ratio",
    )
    path = tmp_path / "factors.json"
    path.write_text("user-owned", encoding="utf-8")

    save_bridging_factors(manifest, path, overwrite=True)

    assert load_bridging_factors(path) == manifest


def test_validator_rejects_non_string_method_with_value_error():
    manifest = create_bridging_factors(
        pd.DataFrame({APTAMER: [1.0, 2.0, 3.0]}),
        pd.DataFrame({APTAMER: [2.0, 4.0, 6.0]}),
        method="median_ratio",
    )
    manifest["method"] = ["median_ratio"]

    with pytest.raises(ValueError, match="Unsupported factor method"):
        validate_bridging_factors(manifest)
