# Banking replication with authorized private inputs

## Public-route calculation

REPORT uses only dataset/group-level aggregates: selections, group counts, trimming cross-counts, contingency tables and saved subspace overlaps. ARIs are calculated from the contingency counts, and angles from the overlaps. It does not claim to estimate the overlap or score correlations from the confidential raw matrix.

No bank names, bank-indexed membership files, raw matrix rows or saved bank-specific scores are included. Source code for the original application/analyzers is retained under `code/matlab/banking/`; it is not automatically run by REPORT or the synthetic FULL pipeline. The authors must confirm authorization before distributing banking-derived aggregates publicly.

## Stronger local read-only calculation

With authorized original inputs outside the package:

```text
python reproduce/private_banking.py --inputs "D:\CSDA_BANKING_EXPORTS" --sources "D:\CSDA_PRIVATE_SOURCES" --output "D:\CSDA_BANKING_CHECK"
```

The output directory must be new. `--inputs` must contain:

```text
ANOR_POSTSELECTION_ANALYSIS_V1_output.zip
analysis_descriptive.zip
ANOR_TRKM_final_V1.zip
ECB_TFKM_final_V1(1).zip
ANOR_TFKM_DESCRIPTIVE_REFIT.mat
ECB_TRKM_DESCRIPTIVE_ANALYSIS.mat
```

`--sources` must contain `data_ECB_ANOR.mat` and either `ECB_cov_input.mat` or the original `R2_ECB_APPLICATION_V12_CEIL_package*.zip` containing it. Each requested source is required to be unique. These matrices are not copied into the replication package.

The reader reconstructs distances/assignments from saved projectors and centroids, rechecks prescribed duplicate-row removal, contingencies, angles, residuals, row-space support, size associations and variable profiles. It does **not** estimate a new clustering fit. It uses temporary extraction directories that are cleaned on normal completion; generated output tables are aggregate-only. A forcibly terminated Python process can leave its temporary directory behind; check the OS temporary area before redistributing any work folder.

62 checks were executed and passed in the package-building environment. The earlier audit had a 63rd check against a redundant historical anonymous label export; that comparison is now optional through `--previous-labels`, rather than requiring redistribution of bank-indexed labels. All other numerical operations are retained. `PRIVATE_READER_ADAPTATION.diff` records these I/O and optional-comparison changes.

Run results do not authorize data release. Keep the private input directory and all working extracts outside any public repository. The package's default ZIP and example REPORT contain no such inputs.
