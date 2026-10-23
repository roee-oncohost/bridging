# SomaLogic bridging toolkit

This repository provides a reusable, auditable API for fitting SomaLogic assay
bridging factors, applying scanner-level factors to a plate, and creating a new
ANML reference file.

The supported public workflows are:

```python
from somalogic_bridging import (
    apply_bridging_to_plate,
    create_anml_reference,
    create_bridging_factors,
)
```

The scripts and notebooks under `projects/` are retained as analysis history.
New code should use the public package above rather than importing a project
script.

## Installation

Python 3.10 or newer is required. From the repository root:

```bash
python -m venv .venv
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

Install the optional exploratory-analysis stack with:

```bash
python -m pip install -e ".[analysis,dev]"
```

`requirements.txt` is the historical full environment snapshot. The editable
install above is the supported installation route for this library.

## 1. Create bridging factors

Inputs are wide pandas DataFrames: samples in rows and SomaLogic SeqIds in
columns. By default, shared columns whose names begin with a digit are treated
as aptamers. Pass `aptamers=[...]` to select them explicitly.

```python
import pandas as pd

from somalogic_bridging import (
    create_bridging_factors,
    save_bridging_factors,
)

source_rfu = pd.read_csv("streck_training.csv")
target_rfu = pd.read_csv("edta_training.csv")

factors = create_bridging_factors(
    source_rfu,
    target_rfu,
    method="mad_log",
    source_label="Streck DS2",  # exact workbook Sample Matrix value
    target_label="EDTA Plasma",
    min_valid_samples=3,
    metadata={"study": "bridging-training-2026"},
)

save_bridging_factors(factors, "factors/streck_ds2_to_edta.json")
```

Every manifest records its schema version, source-to-target direction, sample
counts, aptamer coverage, fit parameters, and user-supplied provenance. Saving
is atomic and refuses to replace a file unless `overwrite=True` is explicit.

### Factor methods

| Method | Transformation from source `x` to target scale | Intended use |
|---|---|---|
| `mad_log` | `exp(a * log(x) + b)`; `a` is the target/source log-MAD ratio | Raw scanner RFUs; aligns log median and unscaled log MAD |
| `median_ratio` | `x * target_median / source_median` | Multiplicative raw-RFU bridging |
| `log_median_shift` | `x * exp(target_log_median - source_log_median)` | Multiplicative raw-RFU bridging fitted in log space |
| `mad_linear` | `a * x + b`; `a` is the target/source linear-MAD ratio | Already-normalized values; not accepted by the scanner-plate API |

Log-domain methods reject non-positive training RFUs. All methods reject
non-numeric or infinite values, insufficient per-aptamer coverage, and invalid
or internally inconsistent manifests.

Use `load_bridging_factors(...)` and `validate_bridging_factors(...)` when
consuming a stored factor artifact.

## 2. Apply factors to scanner files

Map each exact `Sample Matrix` value in the standard SomaLogic workbook to its
factor manifest or saved JSON path:

```python
from somalogic_bridging import apply_bridging_to_plate

run = apply_bridging_to_plate(
    scanner_input_path="plates/OH2026_001/raw",
    workbook_path="plates/OH2026_001/OH2026_001 Workbook.xlsx",
    factors_by_sample_matrix={
        "Streck DS2": "factors/streck_ds2_to_edta.json",
        "Streck nDS": "factors/streck_nds_to_edta.json",
    },
    output_path="plates/OH2026_001/bridged",
)
```

The source label in each manifest must match its mapping key. This guards
against accidentally applying an EDTA-to-Streck artifact in the reverse
direction. For older artifacts whose labels cannot be changed, the check can be
explicitly disabled with `verify_source_labels=False`.

The plate workflow:

- reads the standard `Tab 3 - Assayed Sample List` and `Tab 4 - Plate Map`
  workbook sheets;
- matches scanner filenames by the exact `_<Slide #>_<PDF Subarray>.txt`
  suffix;
- aggregates non-outlier replicate probes before bridging;
- transforms only samples whose `Sample Matrix` has a supplied manifest;
- copies unselected scanner files and the workbook unchanged; and
- returns a JSON-serializable run manifest listing altered and copied files.

For `mad_log`, the 12 HCE aptamers used by the historical project workflow are
excluded by default. Median methods exclude none. Make a different policy
explicit with `excluded_aptamers=[...]`; use an empty tuple to disable the MAD
default. `strict_parameters=True` (the default) rejects a scanner aptamer with
no fitted factor.

All manifests in one plate call must use `mad_log`, or they must all use the
multiplicative `median_ratio`/`log_median_shift` family. Existing non-empty
output directories are rejected unless `overwrite=True`; unrelated files are
not deleted. Processing occurs in a staging directory so validation failures do
not publish a partial new plate. For an atomic handoff, use a new output
directory: publishing over an existing directory updates its files individually
and can leave unrelated files from an earlier run in place.

## 3. Create an ANML reference

Start from an existing reference table so all SomaLogic metadata columns and
their order are preserved. Supply positive pre-normalization RFUs in a wide
sample-by-SeqId DataFrame:

```python
import pandas as pd

from somalogic_bridging import create_anml_reference

rfu = pd.read_parquet("reference_cohort_rfu.parquet")

new_reference, audit = create_anml_reference(
    reference="Reference_V4.1_Plasma_ANML.txt",
    values=rfu,
    output_path="Reference_V4.1_Plasma_ANML_new-cohort.txt",
    min_positive_samples=3,
    metadata={"cohort": "new-cohort", "assay_version": "V4.1"},
)
```

For every updated SeqId, the function calculates:

- `RefMedian`: median of positive RFUs;
- `RefSD`: sample standard deviation of positive RFUs (`ddof=1`);
- `RefLog10Median`: median of `log10` positive RFUs; and
- `RefLog10MAD`: unscaled median absolute deviation in `log10` space.

Non-positive values are excluded from these ANML statistics and counted in the
returned audit manifest. Missing reference SeqIds, extra value columns, and
SeqIds with insufficient positive observations have explicit `error`, `keep`,
or `ignore` policies. Safe defaults reject missing and insufficient SeqIds and
ignore extra metadata columns.

The source reference can never be overwritten. Any other existing destination
also requires `overwrite=True`, and the tab-separated UTF-8 output is written
atomically.

## Reproducibility guidance

- Fit factors on a designated training/reference cohort, then reuse the saved
  manifest; do not refit on the plate being evaluated.
- Treat source and target direction as part of the artifact identity.
- Fit separate manifests for biologically or procedurally distinct sample
  matrices such as desalted and non-desalted Streck.
- Keep cohort, assay version, plate range, code version, and exclusions in the
  manifest metadata or alongside the returned run manifest.
- Validate a newly generated ANML reference in the intended downstream
  normalization workflow before production use.

## Tests

```bash
python -m pytest
```

The test suite covers all factor formulas and schema checks, safe serialization,
ANML statistics and alignment, overwrite protection, scanner-feature
transformation, and a workbook-plus-scanner plate round trip.
