# CSDA V3 replication package: REG-HD-TFKM / REG-HD-TRKM

Public replication package accompanying the CSDA submission **Regularized high-dimensional Trimmed Factorial and Reduced K-means: stability, recovery invariants and structural identification** by Cinzia Di Nuzzo and Matteo Farnè.

**Release:** CSDA V3 final, 23 October 2026.

## Quick start

CHECK and REPORT reconstruct the archived numerical results; they do **not** refit the estimators. On Windows x64, `SETUP_PYTHON_WINDOWS.cmd` creates a local Python runtime/environment without administrator rights. From MATLAB in the repository root:

```matlab
RUN_REPLICATION("SETUP", "Confirm", true); % first setup only
RUN_REPLICATION("CHECK");
RUN_REPLICATION("REPORT");
RUN_REPLICATION("SELFTEST");
```

Without MATLAB, `RUN_REPORT_WINDOWS.cmd` regenerates the saved-result report after Python setup. macOS/Linux users can use `setup_python.sh` with Python 3.13.

## What is reproduced

- Main Table 3, Main Table 4 mechanism checks, Supplement Tables S1-S3 and Figures 2-3 from the archived 41-cell / 1,230-data-set study.
- Supplement Table S4 from the FINAL32 known-parameter benchmark: 32 data sets, six methods, 192 valid fits; REG-HD-TFKM uses HIER4 and REG-HD-TRKM uses RV--CV.
- Banking tables from non-identifying saved aggregates only.
- Figure 1 from the deterministic conceptual-figure script.

`REPORT` validates printed manuscript values against full-precision archived inputs. `SELFTEST` runs small MATLAB contract tests. `FULL` is an explicit, expensive replay route and is never called automatically. A passing REPORT is a saved-result reproduction, not a certification that a new full estimator replay has been completed.

## Confidential banking data

The euro-area supervisory banking microdata are confidential and are **not included**. No bank names, bank identifiers, bank-level observations or bank-indexed partitions are distributed. The public package contains only non-identifying aggregate outputs needed to reproduce the manuscript tables/diagnostics. See `documentation/PRIVATE_BANKING.md`.

## Repository contents

- `results/reference/` archived non-confidential result inputs used by REPORT
- `reproduce/` Python saved-result reconstruction and validation
- `code/matlab/` MATLAB fitting, selection, diagnostic and application code
- `documentation/` protocol, provenance and data dictionaries
- `tests/` regression and package-contract tests
- `validation/` example outputs and build/author receipts

## Citation

Please cite the accompanying manuscript. A journal DOI is not yet available at submission time; `CITATION.cff` will be updated after publication.

## Release status

See `LICENSE_STATUS.md` before redistribution of third-party dependencies. This package does not grant access to the confidential supervisory source data.
