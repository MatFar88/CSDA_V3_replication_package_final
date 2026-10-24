# Data definitions and reporting contract

## Keys and independent units

`scenario_id,replicate` identify the 1,230 synthetic datasets. `method` expands these to 2,460 method-dataset units; the methods use the same generated dataset. The scenario catalogue has 27 MAIN, 5 WEAK_DELTA7 and 9 STRUCTURE cells. Each has 30 replications.

`data_seed =314159265+100*scenario_id+replicate`; TFKM pair seed adds110000, TRKM adds220000. Original candidates reuse the prescribed pair seed; the wrapper does not add a new candidate offset. The original MAT study/config is supplied for MATLAB, while CSVs permit independent inspection without a MATLAB table decoder.

## File roles

- `FRESH_TRUTH_UNIT_RESULTS.csv`: final frozen operational outputs evaluated after truth was opened; one row per method/dataset.
- `POSTHOC_UNIT_DIAGNOSTIC.csv`: additional evaluations of already fitted candidates and decomposition of operational outcomes. It does not redefine the selector.
- `THEORY_AUDITS_ALL.csv`: one-step genuine-label/mask events, empirical-reference spectral diagnostics and in-sample projected classification.
- `CLOSURE_LOCAL_RETURN_ALL.csv`: 24,600 baseline/perturbation records. Figure2 uses BASELINE plus LABEL only, never TRIM or JOINT.
- `CLOSURE_GRANK_ALL.csv`: 1,230 component-respecting distortion paths and empirical-reference midpoint rank diagnostics.
- `FRESH_INVARIANTS_ALL.csv`: saved invariant diagnostics; the radial branch is cross-checked against its separate saved decomposition.
- `A2_TRKM_RADIAL_ARCHIVED_UNIT.csv`: separates radial-only, archived-fallback and final decisions. Table4's510/510 and2/720 concern **radial-only**, not finalTRKM after fallback.
- `WEAK_TRUESTATE_POOL_ALL.csv`: 300 method/dataset records for 150 WEAK datasets; each describes 100 valid starts at generating parameters. Available means some saved start ARI>=.90; Selected means the arbitrator's chosen start ARI>=.90, unconditionally over datasets.

## Evaluation targets

`structure_exact_POSTHOC`: unique terminal structure equals `(r0,g0)`.

`action_correct_POSTHOC`: positive action is the target for leverage only. Clean and bridge target no-positive-trimming; bridge nevertheless has nominal generative contamination .10. Action is a dataset-level decision, not an observation-level contamination indicator.

`joint_correct_POSTHOC`: exact structure AND correct action AND, for leverage, exact trimming level. It is **not** only the conjunction of structure and action.

`truth_state_inlier_ARI_POSTHOC` and `best_alpha0_inlier_ARI_POSTHOC`: state-conditional/reference partition diagnostics, not aggregate selected-output ARIs. Rows without a materialized terminal fit have no invented ARI. No-action descriptive completion is excluded from the default package.

## Empirical reference versus population theory

Legacy column names such as `population_target_correct`, `beta_min`, `projected_oracle_exact` are preserved in source CSVs for provenance. In these archived evaluations, `centroid_source=REALIZED_TRUE_GROUP_CENTROIDS`: the centroids and between reference are constructed from empirical true-group means in the generating subspace, with sample weights. They are not saved population centroids. The projected classification is in-sample and truth-assisted.

The rank statistic uses two stratified halves, but its comparison target and midpoint threshold use the empirical reference. In V15 terminology the threshold is `beta_min_ref/2`, and the diagnostic Weyl ratio is `beta_min_ref/(2*cross_error_ref)`. These diagnostics do not validate the independent cross-fitted construction of the theorem end-to-end.

`one_step_truth_state_recovered` records genuine-label plus mask recovery; it does not certify every contaminant's pre-trimming label. High ARI, exact genuine-label recovery, exact mask recovery and latent-rank identification are different events.

## Group-number path

The reporter recalculates `gain[g+1]=Rhat_g-Rhat_g+1`, checks the minimum underfit gain and first overfit gain, and reconstructs their ratio. The original diagnostic enumerates true-group merges and uses a numerical multistart first split. The reporter does not turn that computed split into a certified unrestricted global optimum.

## Banking

ANOR files denote ratios to assets (363x382); ECB files denote log-level FINREP (363x361). No-action partitions are the explicitly separate descriptive alpha0 fits. The comparisons restrict the descriptive partition to the selected method's retained sample. On log levels the retained TRKM margins1/272 sharply constrain ARI; the package does not reinterpret near-zeroARI as disagreement between two balanced splits.

## Rounding and failures

Full precision is kept in output CSVs. Integer counts and structural keys are checked exactly; float consistency comparisons use stated absolute tolerances. Published displays are compared at their stated rounding precision.

The four Supplement Table S2 median differences at the last printed digit, previously described as `DISPLAY_ROUNDING_DIFFERENCE`, belong to the V1 generation. They were corrected in the current manuscript, as documented in the [author's V1 receipt](../validation/AUTHOR_V1_RECEIPT.md). For the V3 release, all 10 `TableS2` checks in the [saved-result report](../validation/example_report/checks.csv) are `PASS`; the [V3 build validation](../validation/BUILD_VALIDATION_V3.json) records 262 passed checks overall and confirms verification of manuscript display values.

Missing input, duplicate key, mismatched seeds, invalid contingency or a meaningful numerical discrepancy raises an error and writes failure information, not a substituted result.
