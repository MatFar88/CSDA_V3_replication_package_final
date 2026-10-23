TFKM/TRKM MULTISTART BASIN AUDIT v1
===================================

PURPOSE
-------
This is a post-freeze forensic audit of the v11 multistart layer.
It does NOT modify, tune, or relabel the frozen v11 estimator.

Primary scientific question for TFKM:
Did a non-SIGN initialization already reach a materially better terminal
basin at the SAME fixed (r,g,alpha), but the frozen SIGN-only eligibility
rule prevented it from competing?

Secondary question:
If such a basin exists, can a truth-blind, family-neutral challenger select
it more reliably?

TRKM is NOT changed. Its RV_CV selector already scores all valid terminal
fits. TRKM is included only for a truth-aware audit of how close RV_CV is to
the best available basin.

PRIMARY TEST BED
----------------
Scenario 23001: clean
Scenario 23004: bridge
Validation reps 1:8

FROZEN TFKM CURRENT
-------------------
SIGN_ALL_STABLE among valid SIGN_PCA_TKM terminal fits only.
The package first reproduces the stored v11 current winner. If reproduction
fails, the audit stops before interpretation.

TFKM TRUTH-BLIND CHALLENGERS
----------------------------
1. RV_CV_ALL
   Existing robust-view CV score applied to every valid terminal fit.

2. FB_STABLE_ALL
   Family-balanced stability. For each candidate fit, pairwise
   ARI(mutually retained) + 0.2*Jaccard(trim sets) is averaged within each
   initialization family and then equally across families. This removes the
   50/35/10/5 allocation imbalance.

3. HIER_STABLE4
   User-proposed hierarchical challenger:
   - choose one winner inside each family using the frozen within-family
     stability idea;
   - compare the resulting family winners by stability with equal family
     representation;
   - objective is only an exact tie-break.

4. HIER_RVCV4
   Same first stage (one within-family stability winner per family), then
   compare the family winners by RV_CV (geometry then objective tie-breaks).

DIAGNOSTIC-ONLY BASELINES
-------------------------
- ALL_STABLE_NAIVE_DIAGNOSTIC: all valid starts pooled with no family
  balancing. It is deliberately NOT recommended because the family
  allocations are 50/35/10/5 for TFKM.
- MINOBJ_ALL_DIAGNOSTIC: minimum residual objective among all valid starts.

TRUTH-OPEN REFERENCES
---------------------
- ORACLE_BEST_ALL
- ORACLE_BEST_SIGN
These are diagnostics only and are never used by a truth-blind challenger.

KEY QUANTITIES
--------------
best_all_ari       = best inlier ARI among all valid terminal fits
best_sign_ari      = best inlier ARI among valid SIGN fits
basin_gap          = best_all_ari - best_sign_ari
selector_regret    = best_all_ari - ARI selected by the rule

A 0.05 gap is flagged in the text report only as a forensic magnitude. It
is NOT a tuning threshold and is NOT part of any estimator.

RUN
---
Unzip the package anywhere and in MATLAB run:

  cd('<package_folder>')
  RUN_BASIN_AUDIT_FIRST

Then point to your EXISTING v11 output root, e.g.

  output_root = "C:\\Users\\...\\TFKM_TRKM_FINAL_FRESH_VALIDATION_v1_output";

Fast decisive TFKM clean/bridge pass (16 pools):

  report = run_multistart_basin_audit_v1(output_root,"TFKM_SELECTED");

Primary both-method pass (32 pools):

  report = run_multistart_basin_audit_v1(output_root,"PRIMARY");

Add TFKM true-fixed controls:

  report = run_multistart_basin_audit_v1(output_root,"TFKM_SELECTED_TRUE");

The runner is resumable: a completed pool is not refitted.

OUTPUTS
-------
<output_root>/multistart_basin_audit_v1/summary/
  audit_targets.csv
  all_pool_rules.csv
  all_pool_starts.csv
  all_family_summary.csv
  current_pool_diagnostics.csv
  rule_summary.csv
  tfkm_challenger_vs_current.csv
  BASIN_AUDIT_REPORT.txt

Each audited pool also has its own starts/rules/family CSV and a compact MAT
file under multistart_basin_audit_v1/pools/.

INTERPRETATION
--------------
1. best_all >> best_SIGN:
   a better non-SIGN basin was available but excluded by frozen eligibility.

2. best_all approximately equals best_SIGN but current << best_all:
   the issue is the stability arbitration inside SIGN, not family exclusion.

3. best_all itself is poor at the current selected triplet:
   the structural (r,g,alpha) choice is the dominant failure; changing the
   basin selector at that fixed triplet cannot repair it.

4. A challenger improves fixed-triplet ARI:
   this is forensic evidence only. It does NOT retroactively validate a new
   full automatic selector, because changing the terminal fit also changes
   E/N and potentially the later RANK_EN path. Any adopted new estimator
   would need independent validation.
