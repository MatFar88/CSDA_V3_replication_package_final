# PATCH1 — TRKM `NATIVE010_ONE_STEP`

This patch changes **only downstream structural selection for TRKM**. Candidate fitting is unchanged.

## Why existing fresh candidate fits remain valid

Each fresh session computes the complete `(r,g,alpha)` cube before `select_final_pipeline_v1` is called by the merge step. The patch does not change the DGP, validation seeds, multistart seeds, number of starts, HD fit, candidate metrics, or checkpoint keys. Therefore candidate rows already completed under the original package are exactly the rows required by PATCH1 and must **not** be recomputed.

The safeguard decision was based on the consumed Confirm4 offline audit. At the time PATCH1 was adopted, only about 15 fresh candidate rows had been produced, and all were REG_TFKM_HD. PATCH1 leaves TFKM unchanged. No TRKM fresh candidate had been used, and the count was far below the 150 rows required for one complete method/scenario/replicate cube, so no fresh final structural selection was available to drive the change.

## Frozen PATCH1 rule

TFKM is unchanged.

TRKM pipeline:

`CURRENT -> G_STRICT_W -> R_STRICT_NR2 -> NATIVE010_ONE_STEP -> frozen alpha`

`NATIVE010_ONE_STEP` attempts **at most one** `g -> g-1` move. Let

`r_candidate = min(r_current, g_candidate-1)`.

Accept only when both hold:

- `R2_candidate >= R2_reference - 0.01`;
- `RV_CV_candidate >= RV_CV_reference`.

The step is never iterated. The frozen alpha selector is then re-applied on the resulting `(r,g)` path.

## Existing output root

Do **not** rerun setup or DGP certification. Once this patched package is on the MATLAB path, register the metadata once:

```matlab
apply_trkm_native010_patch_v1(output_root)
```

This edits only `study.mat` metadata and writes a patch marker. Existing result/checkpoint files are untouched.
