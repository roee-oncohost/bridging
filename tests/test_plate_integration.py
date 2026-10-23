from pathlib import Path

from openpyxl import Workbook
import pandas as pd
import pytest

import somalogic_bridging as sb


APTAMER = "10000-28"


def _write_workbook(
    path: Path,
    *,
    include_plate_row: bool = True,
    duplicate_coordinate: bool = False,
) -> None:
    workbook = Workbook()
    samples = workbook.active
    samples.title = "Tab 3 - Assayed Sample List"
    samples["A1"] = "Sample Number"
    samples["F1"] = "SampleId"
    samples["J1"] = "Sample Matrix"
    samples["A2"] = 1
    samples["F2"] = "SAMPLE-1"
    samples["J2"] = "Streck DS2"

    plate_map = workbook.create_sheet("Tab 4 - Plate Map")
    plate_map["A7"] = "Plate Position"
    plate_map["E7"] = "NAME"
    plate_map["H7"] = "Slide #"
    plate_map["K7"] = "PDF Subarray"
    if include_plate_row:
        plate_map["A8"] = "A1"
        plate_map["E8"] = "SAMPLE-1"
        plate_map["H8"] = 1
        plate_map["K8"] = 1
    if duplicate_coordinate:
        samples["A3"] = 2
        samples["F3"] = "SAMPLE-2"
        samples["J3"] = "Streck DS2"
        plate_map["A9"] = "A2"
        plate_map["E9"] = "SAMPLE-2"
        plate_map["H9"] = 1
        plate_map["K9"] = 1
    workbook.save(path)


def _scanner_text(first_signal: float = 1.0, second_signal: float = 3.0) -> str:
    return (
        "TYPE\ttext\n"
        "FEPARAMS\tKey\n"
        "DATA\tvalue\n"
        "*\n"
        "TYPE\tinteger\n"
        "STATS\tCount\n"
        "DATA\t1\n"
        "*\n"
        "TYPE\ttext\tfloat\tboolean\n"
        "FEATURES\tProbeName\tgProcessedSignal\tgIsFeatPopnOL\n"
        f"DATA\tanti-{APTAMER}_x\t{first_signal}\t0\n"
        f"DATA\tanti-{APTAMER}_y\t{second_signal}\t0\n"
        "DATA\tcontrol\t99\t0\n"
    )


def test_apply_bridging_to_plate_round_trip(tmp_path):
    scanner_directory = tmp_path / "scanner"
    scanner_directory.mkdir()
    selected_name = "run_1_1.txt"
    untouched_name = "run_2_1.txt"
    (scanner_directory / selected_name).write_text(_scanner_text(), encoding="utf-8")
    untouched_bytes = _scanner_text(7.0, 9.0).encode("utf-8")
    (scanner_directory / untouched_name).write_bytes(untouched_bytes)

    workbook_path = tmp_path / "plate.xlsx"
    _write_workbook(workbook_path)
    original_workbook = workbook_path.read_bytes()
    factors = sb.create_bridging_factors(
        pd.DataFrame({APTAMER: [1.0, 2.0, 3.0]}),
        pd.DataFrame({APTAMER: [2.0, 4.0, 6.0]}),
        method="median_ratio",
        source_label="Streck DS2",
        target_label="EDTA Plasma",
    )
    output_directory = tmp_path / "bridged"

    manifest = sb.apply_bridging_to_plate(
        scanner_directory,
        workbook_path,
        {"Streck DS2": factors},
        output_directory,
    )

    bridged_text = (output_directory / selected_name).read_text(encoding="utf-8")
    assert bridged_text.count(f"DATA\tanti-{APTAMER}") == 2
    assert bridged_text.count("\t4.0\t0") == 2
    assert (output_directory / untouched_name).read_bytes() == untouched_bytes
    assert (output_directory / workbook_path.name).read_bytes() == original_workbook
    assert manifest["method"] == "median_multiplier"
    assert manifest["altered_files"] == [selected_name]
    assert manifest["copied_files"] == [untouched_name]
    assert manifest["sample_types"]["Streck DS2"]["sample_count"] == 1
    assert manifest["factor_methods"] == ["median_ratio"]
    assert manifest["source_labels_verified"] is True


def test_missing_plate_map_sample_fails_before_output_is_published(tmp_path):
    scanner_directory = tmp_path / "scanner"
    scanner_directory.mkdir()
    (scanner_directory / "run_1_1.txt").write_text(
        _scanner_text(), encoding="utf-8"
    )
    workbook_path = tmp_path / "plate.xlsx"
    _write_workbook(workbook_path, include_plate_row=False)
    factors = sb.create_bridging_factors(
        pd.DataFrame({APTAMER: [1.0, 2.0, 3.0]}),
        pd.DataFrame({APTAMER: [2.0, 4.0, 6.0]}),
        method="median_ratio",
        source_label="Streck DS2",
        target_label="EDTA Plasma",
    )
    output_directory = tmp_path / "bridged"

    with pytest.raises(ValueError, match="missing plate-map rows"):
        sb.apply_bridging_to_plate(
            scanner_directory,
            workbook_path,
            {"Streck DS2": factors},
            output_directory,
        )

    assert not output_directory.exists()


def test_duplicate_scanner_coordinate_is_rejected(tmp_path):
    scanner_directory = tmp_path / "scanner"
    scanner_directory.mkdir()
    (scanner_directory / "run_1_1.txt").write_text(
        _scanner_text(), encoding="utf-8"
    )
    workbook_path = tmp_path / "plate.xlsx"
    _write_workbook(workbook_path, duplicate_coordinate=True)
    factors = sb.create_bridging_factors(
        pd.DataFrame({APTAMER: [1.0, 2.0, 3.0]}),
        pd.DataFrame({APTAMER: [2.0, 4.0, 6.0]}),
        method="median_ratio",
        source_label="Streck DS2",
        target_label="EDTA Plasma",
    )

    with pytest.raises(ValueError, match="same scanner file"):
        sb.apply_bridging_to_plate(
            scanner_directory,
            workbook_path,
            {"Streck DS2": factors},
            tmp_path / "bridged",
        )
