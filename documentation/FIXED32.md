# Table S4: FINAL32 fixed-parameter experiment

Manuscript target: **CSDA submission dated 23 October 2026, Supplement Table S4**.

The source is `results/reference/fixed32/FIXED_all_results.csv`. It contains
192 valid fits: four scenario IDs (33001, 33004, 33005, 33007), eight
replications each, and six methods. The eight worker CSVs are matched row by
row against this merged file and are never counted as additional observations.
The three supplied summaries are recomputed and checked column by column.

The original method IDs are `TFKM_V2_HIER4`, `TRKM_RVCV`, `TKM`, `PCA_TKM`,
`RPCA_TKM`, and `TCLUST`. The manuscript labels REG-HD-TFKM and REG-HD-TRKM
denote the first two respectively. REG-HD-TFKM uses HIER4 start arbitration
and REG-HD-TRKM uses RV--CV. TKM and TCLUST receive `(g0, alpha0)` because they
do not have a latent-rank parameter. `FIXED4_PROTOCOL.txt` is included
unmodified.

Mean and median ARI in Table S4 use `ari_inliers_full`, not
`ari_inliers_retained`. Runtime uses `runtime_seconds`. Every saved fit records
100 starts.

## Original producer source versus a new estimator run

`code/matlab/experiments/fixed32/TFKM_V2_FINAL32_COMPETITORS_V1.zip` is the
recovered original source archive, included without modification and isolated
from the FRESH MATLAB paths. Its merger writes `summary/FIXED_all_results.csv`.
The archive also contains development stages not used by REPORT and not used as
evidence in Table S4. Its saved-context runners require the original local
study/merged-Units/template MAT files and the listed MATLAB/FSDA dependencies.
Their execution is not certified here. REPORT requires none of those files and
never refits an estimator. Runtime columns are historical recorded values, not
timing predictions for another machine.
