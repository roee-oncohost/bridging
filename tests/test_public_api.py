from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import somalogic_bridging as sb
from src.altering_plate_multiple_collection_methods import (
    DEFAULT_HCE_APTAMERS,
    _manifest_factors_to_legacy,
    transform_feature_df_MAD_bridging,
    transform_feature_df_median_bridging,
)


APTAMER = "10000-28"


def _manifest(method: str = "median_ratio"):
    source = pd.DataFrame({APTAMER: [1.0, 2.0, 3.0]})
    if method == "mad_log":
        target = pd.DataFrame({APTAMER: [2.0, 4.0, 8.0]})
    else:
        target = pd.DataFrame({APTAMER: [2.0, 4.0, 6.0]})
    return sb.create_bridging_factors(
        source,
        target,
        method=method,
        source_label="Streck DS2",
        target_label="EDTA Plasma",
    )


def _feature_frame():
    return pd.DataFrame(
        {
            "ProbeName": ["anti-10000-28_x", "anti-10000-28_y", "control"],
            "gProcessedSignal": [1.0, 3.0, 99.0],
            "gIsFeatPopnOL": [0, 0, 0],
        }
    )


def test_public_functions_are_exposed():
    assert callable(sb.create_bridging_factors)
    assert callable(sb.apply_bridging_to_plate)
    assert callable(sb.create_anml_reference)


def test_median_manifest_can_transform_scanner_features():
    manifest = _manifest("median_ratio")
    params = _manifest_factors_to_legacy(
        manifest, {"median_ratio", "log_median_shift"}
    )

    transformed = transform_feature_df_median_bridging(_feature_frame(), params)

    # Replicate signals aggregate to 2.0, then the target/source factor is 2.
    assert transformed.loc[:1, "gProcessedSignal"].tolist() == pytest.approx([4.0, 4.0])
    assert transformed.loc[2, "gProcessedSignal"] == 99.0


def test_mad_manifest_can_transform_scanner_features():
    manifest = _manifest("mad_log")
    params = _manifest_factors_to_legacy(manifest, {"mad_log"})

    transformed = transform_feature_df_MAD_bridging(_feature_frame(), params)

    expected = manifest["factors"][APTAMER]
    value = 2.0
    expected_value = float(np.exp(expected["a"] * np.log(value) + expected["b"]))
    assert transformed.loc[:1, "gProcessedSignal"].tolist() == pytest.approx(
        [expected_value, expected_value]
    )


def test_high_level_plate_api_selects_median_runner(monkeypatch, tmp_path):
    import src.altering_plate_multiple_collection_methods as plate_module

    scanner_path = tmp_path / "scanner"
    scanner_path.mkdir()
    workbook_path = tmp_path / "plate.xlsx"
    workbook_path.write_bytes(b"placeholder")
    output_path = tmp_path / "output"
    captured = {}

    def fake_runner(input_path, workbook, manifests, output, **kwargs):
        captured.update(
            input_path=Path(input_path),
            workbook=Path(workbook),
            manifests=manifests,
            output=Path(output),
            kwargs=kwargs,
        )
        return {"method": "median_ratio", "altered_files": [], "copied_files": []}

    monkeypatch.setattr(plate_module, "alter_scanner_files_median_bridging", fake_runner)

    result = sb.apply_bridging_to_plate(
        scanner_path,
        workbook_path,
        {"Streck DS2": _manifest("median_ratio")},
        output_path,
    )

    assert captured["kwargs"]["excluded_aptamers"] == ()
    assert result["factor_manifests"]["Streck DS2"]["method"] == "median_ratio"


def test_high_level_plate_api_uses_hce_default_for_mad(monkeypatch, tmp_path):
    import src.altering_plate_multiple_collection_methods as plate_module

    scanner_path = tmp_path / "scanner"
    scanner_path.mkdir()
    workbook_path = tmp_path / "plate.xlsx"
    workbook_path.write_bytes(b"placeholder")
    captured = {}

    def fake_runner(*args, **kwargs):
        captured.update(kwargs)
        return {"method": "mad_log", "altered_files": [], "copied_files": []}

    monkeypatch.setattr(plate_module, "alter_scanner_files_MAD_bridging", fake_runner)
    sb.apply_bridging_to_plate(
        scanner_path,
        workbook_path,
        {"Streck DS2": _manifest("mad_log")},
        tmp_path / "output",
    )

    assert set(captured["excluded_aptamers"]) == set(DEFAULT_HCE_APTAMERS)


def test_high_level_plate_api_rejects_nonempty_output(tmp_path):
    scanner_path = tmp_path / "scanner"
    scanner_path.mkdir()
    workbook_path = tmp_path / "plate.xlsx"
    workbook_path.write_bytes(b"placeholder")
    output_path = tmp_path / "output"
    output_path.mkdir()
    (output_path / "owned.txt").write_text("do not overwrite", encoding="utf-8")

    with pytest.raises(FileExistsError, match="non-empty"):
        sb.apply_bridging_to_plate(
            scanner_path,
            workbook_path,
            {"Streck DS2": _manifest("median_ratio")},
            output_path,
        )


def test_high_level_plate_api_never_overwrites_source_workbook(tmp_path):
    scanner_path = tmp_path / "scanner"
    scanner_path.mkdir()
    workbook_path = tmp_path / "plate.xlsx"
    original = b"user-owned-workbook"
    workbook_path.write_bytes(original)

    with pytest.raises(ValueError, match="source workbook"):
        sb.apply_bridging_to_plate(
            scanner_path,
            workbook_path,
            {"Streck DS2": _manifest("median_ratio")},
            tmp_path,
            overwrite=True,
        )
    assert workbook_path.read_bytes() == original


def test_high_level_plate_api_rejects_linear_factors(tmp_path):
    scanner_path = tmp_path / "scanner"
    scanner_path.mkdir()
    workbook_path = tmp_path / "plate.xlsx"
    workbook_path.write_bytes(b"placeholder")

    with pytest.raises(ValueError, match="already-normalized"):
        sb.apply_bridging_to_plate(
            scanner_path,
            workbook_path,
            {"Streck DS2": _manifest("mad_linear")},
            tmp_path / "output",
        )


def test_high_level_plate_api_rejects_reversed_or_mismatched_direction(tmp_path):
    scanner_path = tmp_path / "scanner"
    scanner_path.mkdir()
    workbook_path = tmp_path / "plate.xlsx"
    workbook_path.write_bytes(b"placeholder")
    manifest = _manifest("median_ratio")
    manifest["direction"] = {"source": "EDTA Plasma", "target": "Streck DS2"}

    with pytest.raises(ValueError, match="source labels"):
        sb.apply_bridging_to_plate(
            scanner_path,
            workbook_path,
            {"Streck DS2": manifest},
            tmp_path / "output",
        )


def test_high_level_plate_api_rejects_string_exclusions(tmp_path):
    scanner_path = tmp_path / "scanner"
    scanner_path.mkdir()
    workbook_path = tmp_path / "plate.xlsx"
    workbook_path.write_bytes(b"placeholder")

    with pytest.raises(TypeError, match="not a string"):
        sb.apply_bridging_to_plate(
            scanner_path,
            workbook_path,
            {"Streck DS2": _manifest("median_ratio")},
            tmp_path / "output",
            excluded_aptamers="2171-12",
        )


def test_mad_transform_rejects_nonpositive_or_all_outlier_signals():
    manifest = _manifest("mad_log")
    params = _manifest_factors_to_legacy(manifest, {"mad_log"})
    feature_frame = _feature_frame()
    feature_frame.loc[:1, "gIsFeatPopnOL"] = 1

    with pytest.raises(ValueError, match="finite and positive"):
        transform_feature_df_MAD_bridging(feature_frame, params)
