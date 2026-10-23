TFKM HIER_STABLE4 FULL-PIPELINE CHALLENGE v1
============================================

PURPOSE
-------
Development-only full-pipeline test of the prospectively frozen challenger
HIER_STABLE4 discovered after the v11 fixed-triplet multistart basin audit.

ONLY CHANGE
-----------
TFKM multistart arbitration at every fixed (r,g,alpha):

  frozen v11:  SIGN_ALL_STABLE (SIGN family only)
  challenger:  HIER_STABLE4

Everything else remains unchanged:
- 100 starts and 50/35/10/5 family allocation
- deterministic start seeds
- robust initialization definitions
- corrected alternating TFKM fit
- ridge scales tauW/tauT
- validity contract
- candidate diagnostics E, R2, W_ratio, W_raw
- RANK_EN + GMAX_R_E
- G_STRICT_W
- R_STRICT_NR2
- GAIN_OBJ_L0800
- RG vote -> alpha composition
- TRKM is not refitted or modified by this challenge

HIER_STABLE4
------------
1. Within WINSOR, SIGN, TKM_ORIGINAL and RP_TKM separately, choose the
   terminal fit with maximum TFKM stability:
      ARI(mutually retained labels) + 0.2*Jaccard(trim sets).
   Exact ties -> lower trimmed objective.
2. Compare at most four family winners with the SAME similarity.
3. Choose the most stable family representative; exact ties -> lower objective.
4. The second-stage winner stability is stored as selector_score = N and is
   therefore the native evidence supplied to the unchanged RANK_EN pipeline.

WHY THIS IS FAMILY-NEUTRAL
--------------------------
The 50/35/10/5 allocation controls exploration within a family. The final
competition has at most one representative per family, so allocation size does
not determine final voting weight.

FIRST RUN
---------
Open RUN_HIER_STABLE4_FIRST.m and set output_root to the ORIGINAL v11 output.
Start with CLEAN_BRIDGE. Four workers can be run in parallel in four MATLAB
sessions using worker_id 1..4, n_workers=4.

OUTPUTS
-------
Written only under:
  <original output_root>\hier_stable4_full_pipeline_v1\

The original frozen v11 files are never overwritten.

PRIMARY INTERPRETATION
----------------------
After merge, inspect:
- mean/min ARI and delta vs current
- structural_over_rate
- over_r_rate and over_g_rate separately
- any-under rate (safety)
- exact (r,g)
- stage1 -> stage2 changes
- final initialization-family usage

IMPORTANT SCIENTIFIC STATUS
---------------------------
The v11 datasets have already been inspected and were used to choose
HIER_STABLE4 over other challengers. Therefore this full-pipeline run is
DEVELOPMENT evidence, not independent validation. If HIER_STABLE4 is adopted,
a new untouched validation panel is required for confirmatory performance.
