TFKM/TRKM FINAL UNSEEN CONFIRMATION v1
======================================

ONE ROOT
--------
Use one brand-new output root for both methods, e.g.
  C:\Users\matteo.farne2\Documents\TFKM_TRKM_UNSEEN_FINAL_v1_output

WHAT IS FROZEN
--------------
DGP definitions: same four development cells.
New scenario IDs: 33001, 33004, 33005, 33007.
New validation seed namespace: 1618033 + sid*100000 + replicate.
Replicates: 8.
Grid: 15 (r,g) pairs x 10 alpha values.
100 starts/candidate.

TFKM candidate fit:
  HIER_STABLE4, exactly as frozen after development.
  Extra bookkeeping saves cross-family ARI and J components; it does NOT alter selection.
TFKM structure:
  SUPPORTED_WEAK.
TFKM regime/alpha:
  tau* < .01               -> LEVERAGE -> alpha = unique argmax J(alpha>0)
  .01 <= tau* <= .05       -> CLEAN    -> alpha = 0
  tau* > .05               -> BRIDGE_AMBIGUOUS -> alpha not identified

TRKM candidate fit:
  unchanged RV_CV multistart core.
TRKM structure:
  RV_R2_LEX_CONSERVATIVE.
TRKM regime/alpha:
  no RV break              -> CLEAN -> alpha = 0
  0 < first RV break < .01 -> BRIDGE_AMBIGUOUS -> alpha not identified
  first RV break >= .01    -> LEVERAGE -> alpha = FIRST_RV_LOSS high alpha

TRUTH FIREWALL
--------------
Operational candidate CSVs contain no r0/g0/alpha0/contamination/ARI/trim-F1.
Truth metrics are written to:
  _SEALED_TRUTH_DO_NOT_OPEN_BEFORE_SELECTION
The operational selector never reads this folder.
Do not inspect it manually before operational_selection/SELECTION_FROZEN.mat exists.

RUN
---
1) Open RUN_UNSEEN_FIRST.m and run the setup section.
2) Run sessions A:H, ideally in parallel MATLAB windows:
     run_unseen_session_v1("A",output_root)
     ...
     run_unseen_session_v1("H",output_root)
3) Check:
     inspect_unseen_progress_v1(output_root)
4) When complete:
     M = merge_unseen_operational_v1(output_root);
     R_UNSEEN = select_unseen_operational_v1(output_root);
5) Only then:
     E_UNSEEN = evaluate_unseen_truth_v1(output_root);
     E_UNSEEN.summary

DO NOT RETUNE AFTER UNSEEN RESULTS.
