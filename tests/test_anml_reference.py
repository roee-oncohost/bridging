import json

import numpy as np
import pandas as pd
import pytest

from src.anml_reference import ANML_STAT_COLUMNS, create_anml_reference


def _reference() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Target": ["protein-b", "protein-a"],
            "SeqId": ["B", "A"],
            "Type": ["Protein", "Protein"],
            "RefLog10Median": [1.0, 1.3],
            "RefLog10MAD": [0.1, 0.2],
            "RefMedian": [10.0, 20.0],
            "RefSD": [1.0, 2.0],
        },
        index=[10, 10],
    )


def test_create_anml_reference_updates_by_seqid_and_writes_tsv(tmp_path):
    values = pd.DataFrame(
        {
            "A": [1.0, 10.0, 100.0, 0.0],
            "B": [2.0, 2.0, 4.0, 8.0],
            "PlateId": ["P1", "P1", "P2", "P2"],
        }
    )
    output_path = tmp_path / "new_reference.txt"

    result, manifest = create_anml_reference(
        _reference(),
        values,
        output_path,
        metadata={"cohort": "unit-test"},
    )

    assert result.columns.tolist() == _reference().columns.tolist()
    assert result["SeqId"].tolist() == ["B", "A"]
    a = result.set_index("SeqId").loc["A"]
    assert a["RefMedian"] == pytest.approx(10.0)
    assert a["RefSD"] == pytest.approx(np.std([1.0, 10.0, 100.0], ddof=1))
    assert a["RefLog10Median"] == pytest.approx(1.0)
    assert a["RefLog10MAD"] == pytest.approx(1.0)

    written = pd.read_csv(output_path, sep="\t")
    pd.testing.assert_frame_equal(written, result.reset_index(drop=True), check_dtype=False)
    assert manifest["updated_seq_ids"] == ["B", "A"]
    assert manifest["extra_columns"] == ["PlateId"]
    assert manifest["positive_sample_counts"]["A"] == 3
    assert manifest["nonpositive_sample_counts"]["A"] == 1
    assert manifest["metadata"] == {"cohort": "unit-test"}
    json.dumps(manifest, allow_nan=False)


def test_keep_policies_leave_missing_and_insufficient_seqids_unchanged():
    reference = _reference()
    values = pd.DataFrame({"A": [0.0, 5.0], "metadata": [1, 2]})

    result, manifest = create_anml_reference(
        reference,
        values,
        min_positive_samples=2,
        missing_policy="keep",
        insufficient_policy="keep",
    )

    pd.testing.assert_frame_equal(result, reference)
    assert manifest["updated_seq_id_count"] == 0
    assert manifest["missing_seq_ids"] == ["B"]
    assert manifest["insufficient_seq_ids"] == ["A"]


@pytest.mark.parametrize(
    ("values", "message"),
    [
        (pd.DataFrame({"A": [1.0, 2.0, 3.0]}), "missing 1 reference SeqIds"),
        (
            pd.DataFrame({"A": [1.0, 2.0, 3.0], "B": [1.0, 2.0, 3.0], "extra": [1, 2, 3]}),
            "columns absent from the reference",
        ),
    ],
)
def test_strict_coverage_policies_raise(values, message):
    kwargs = {"extra_policy": "error"} if "extra" in values.columns else {}
    with pytest.raises(ValueError, match=message):
        create_anml_reference(_reference(), values, **kwargs)


def test_duplicate_reference_seqids_are_rejected():
    reference = _reference()
    reference.loc[:, "SeqId"] = ["A", "A"]
    with pytest.raises(ValueError, match="duplicate SeqId"):
        create_anml_reference(reference, pd.DataFrame({"A": [1.0, 2.0, 3.0]}))


def test_non_numeric_matched_measurements_are_rejected():
    values = pd.DataFrame({"A": [1.0, "bad", 3.0], "B": [1.0, 2.0, 3.0]})
    with pytest.raises(ValueError, match="non-numeric or infinite"):
        create_anml_reference(_reference(), values)


def test_source_reference_cannot_be_overwritten(tmp_path):
    source = tmp_path / "reference.txt"
    _reference().to_csv(source, sep="\t", index=False)
    values = pd.DataFrame({"A": [1.0, 2.0, 3.0], "B": [1.0, 2.0, 3.0]})

    with pytest.raises(ValueError, match="must not overwrite"):
        create_anml_reference(source, values, source)


def test_existing_output_requires_explicit_overwrite(tmp_path):
    destination = tmp_path / "existing.txt"
    destination.write_text("user-owned", encoding="utf-8")
    values = pd.DataFrame({"A": [1.0, 2.0, 3.0], "B": [1.0, 2.0, 3.0]})

    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        create_anml_reference(_reference(), values, destination)
    assert destination.read_text(encoding="utf-8") == "user-owned"

    create_anml_reference(_reference(), values, destination, overwrite=True)
    assert pd.read_csv(destination, sep="\t")["SeqId"].tolist() == ["B", "A"]


def test_path_reference_preserves_literal_metadata_and_textual_seqid(tmp_path):
    source = tmp_path / "reference.txt"
    reference = _reference().iloc[[0]].copy()
    reference.loc[:, "SeqId"] = "00123"
    reference.loc[:, "Target"] = "NA"
    reference["Comment"] = ""
    reference["BatchCode"] = "001"
    reference.to_csv(source, sep="\t", index=False)
    values = pd.DataFrame({"00123": [1.0, 2.0, 3.0]})

    result, _ = create_anml_reference(source, values)

    assert result.loc[result.index[0], "SeqId"] == "00123"
    assert result.loc[result.index[0], "Target"] == "NA"
    assert result.loc[result.index[0], "Comment"] == ""
    assert result.loc[result.index[0], "BatchCode"] == "001"


@pytest.mark.parametrize("metadata", [{"bad": float("nan")}, [("pair", 1)]])
def test_metadata_must_be_a_finite_json_mapping(metadata):
    values = pd.DataFrame({"A": [1.0, 2.0, 3.0], "B": [1.0, 2.0, 3.0]})
    with pytest.raises((TypeError, ValueError), match="metadata"):
        create_anml_reference(_reference(), values, metadata=metadata)


def test_invalid_retained_reference_statistics_are_rejected():
    reference = _reference()
    reference.loc[reference["SeqId"] == "B", "RefSD"] = float("inf")
    values = pd.DataFrame({"A": [1.0, 2.0, 3.0]})

    with pytest.raises(ValueError, match="Retained ANML statistics"):
        create_anml_reference(reference, values, missing_policy="keep")


def test_nonfinite_computed_statistics_are_rejected():
    reference = _reference().iloc[[1]].copy()
    values = pd.DataFrame({"A": [1.0e308, 9.0e307, 8.0e307]})

    with np.errstate(over="ignore", invalid="ignore"):
        with pytest.raises(ValueError, match="non-finite"):
            create_anml_reference(reference, values)


def test_failed_write_leaves_no_destination_or_temporary_file(monkeypatch, tmp_path):
    destination = tmp_path / "failed.txt"
    values = pd.DataFrame({"A": [1.0, 2.0, 3.0], "B": [1.0, 2.0, 3.0]})

    def fail_to_csv(*args, **kwargs):
        raise RuntimeError("simulated write failure")

    monkeypatch.setattr(pd.DataFrame, "to_csv", fail_to_csv)
    with pytest.raises(RuntimeError, match="simulated"):
        create_anml_reference(_reference(), values, destination)

    assert not destination.exists()
    assert list(tmp_path.glob(".failed.txt.*.tmp")) == []


def test_required_reference_columns_are_validated():
    reference = _reference().drop(columns=[ANML_STAT_COLUMNS[0]])
    with pytest.raises(ValueError, match="missing required columns"):
        create_anml_reference(reference, pd.DataFrame({"A": [1.0, 2.0, 3.0]}))
