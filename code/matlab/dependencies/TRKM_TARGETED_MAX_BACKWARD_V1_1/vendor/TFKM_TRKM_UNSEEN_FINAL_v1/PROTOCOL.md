# Final fresh validation protocol v1

This is a new experiment. No Confirm4 data, replicate IDs, or seeds are reused.
The old calibrated registry is used only to define the four DGP parameter cells; the clean cell is pre-specified at delta=10, rho=2.

## Independence

- DGP certification seed: `271828 + scenario_id*100000 + cert_rep`.
- Validation seed: `314159 + scenario_id*100000 + validation_rep`.
- New scenario IDs: 23001, 23004, 23005, 23007.
- Certification reps and validation reps are both numbered 1:8 but live in disjoint seed namespaces.

The validation runners assert that the independent DGP certification passed before any validation fit can be launched.

## Frozen operational pipeline

No truth is used.

1. Run the frozen 100-start HD fits for the complete 15 x 10 `(r,g,alpha)` cube.
2. Stage 1:
   - TFKM: frozen CURRENT joint selector plus the pre-frozen `GMAX_R_E` boundary rank guard.
   - TRKM: frozen CURRENT joint selector.
3. Stage 2:
   - `G_STRICT_W`: lower g only under strict non-worsening of W_raw and W_ratio relative to the Stage-1 reference.
   - `R_STRICT_NR2`: lower r only under strict non-worsening of R2 and method-native selector evidence.
4. PATCH1, TRKM only: apply `NATIVE010_ONE_STEP` after Stage 2. Attempt at most one `g -> g-1` move with `r_candidate=min(r_current,g_candidate-1)`; accept only if `R2_candidate >= R2_reference - .01` and `RV_CV_candidate >= RV_CV_reference`. Never iterate.
5. Re-apply the already-frozen alpha selector along the resulting `(r,g)` path.

## Fresh validation endpoints

- ARI on genuine observations.
- Exact `(r,g)` and no-under `(r,g)`.
- r/g overshoot.
- Exact alpha and exact joint recovery.
- Positive ARI regret vs the true-fixed cube row.
- Trim precision/recall/F1 where contamination exists.
- Conditional quality of structural over-selection:
  - mean ARI among overselected fits;
  - minimum ARI among overselected fits;
  - fraction of overselected fits with ARI >= .95.

The .95 value is an interpretive validation threshold only. It does not alter model selection.

## Post-selection forensic layer

Truth is opened only after final selection. By default every structurally overselected unit or ARI<.95 is deterministically refitted at the already-selected `(r,g,alpha)` using the same production seeds. Reproduction against the stored cube ARI is asserted at 1e-8.

The forensic output is designed to distinguish:

- extra clusters that merely absorb contamination;
- pure or impure fragmentation of a true group;
- mixed cross-truth overclustering;
- extra directions that are random-noise-like;
- directions associated with outlier status;
- sample-adaptive directions associated with labels despite being orthogonal to the population truth span;
- genuine estimated-subspace distortion;
- low ARI recoverable by merging versus residual cross-truth misallocation.

No forensic diagnostic is permitted to change the frozen procedure in this experiment.


## PATCH1 provenance and resumability

PATCH1 was accepted using only the consumed Confirm4 offline safeguard audit. Candidate fitting is unchanged, so candidate rows/checkpoints already produced by the fresh sessions remain valid and are reused. At adoption, only about 15 fresh candidate rows existed and all were REG_TFKM_HD. TFKM is unchanged by PATCH1; no TRKM fresh candidate had been used, and no complete 150-row unit/final structural selection existed.
