TFKM/TRKM FINAL FRESH VALIDATION v1
===================================

PURPOSE
- Completely new datasets: new scenario IDs and new seed namespaces.
- The clean TFKM-favourable moderate cell is pre-recalibrated to delta=10,rho=2.
- DGP certification and validation use DIFFERENT seed namespaces.
- Frozen method pipeline:
    TFKM Stage1 = CURRENT + GMAX_R_E boundary rank guard
    TRKM Stage1 = CURRENT
    Stage2 g    = G_STRICT_W
    Stage2 r    = R_STRICT_NR2
    TRKM PATCH1 = NATIVE010_ONE_STEP (one g-down step max; R2 loss <= .01 and RV_CV non-decrease)
    alpha       = frozen alpha selector re-applied after structural selection
- No truth enters fitting or selection.
- Truth-open forensics are run only AFTER final selection and never feed back.

COMPUTE SIZE
- 4 scenarios x 8 untouched validation replicates x 2 methods x 15 (r,g) pairs x 10 alpha values
  = 9600 candidate fits, each with the frozen 100-start internal multistart.
- Sessions A:H each contain one replicate = 1200 candidate fits and are resumable.

RUN ORDER
1) In the package folder:
   RUN_ME_FIRST

2) Create a NEW output directory, e.g.
   output_root = "C:\Users\matteo.farne2\Documents\TFKM_TRKM_FINAL_FRESH_VALIDATION_v1_output";
   study = setup_final_fresh_validation_v1(output_root);

3) Cheap independent DGP gate. This does NOT touch validation seeds:
   cert = certify_fresh_dgp_catalog_v1(output_root);
   All four scenarios must show cluster_cert_rate=1 and rank_cert_rate=1.
   Validation sessions are programmatically locked until this passes.

4) Run A:H. They may run concurrently in separate MATLAB sessions:
   run_final_fresh_session_v1("A",output_root);
   ...
   run_final_fresh_session_v1("H",output_root);

5) If this output_root was created with the pre-PATCH1 package, register PATCH1 ONCE (no refits):
   apply_trkm_native010_patch_v1(output_root)

6) Merge and apply the frozen PATCH1 final selector:
   report = merge_final_fresh_validation_v1(output_root);

7) Forensic explanation layer. Default is TARGETED: every structural over-selection or ARI<.95.
   These sessions may again run concurrently:
   run_final_fresh_forensics_session_v1("A",output_root);
   ...
   run_final_fresh_forensics_session_v1("H",output_root);
   forensic = merge_final_fresh_forensics_v1(output_root);

WHAT THE FORENSICS EXPLAIN
For extra clusters:
- maximum-overlap one-to-one core matching between true and fitted clusters;
- unmatched clusters are explicitly marked as extra;
- counts of genuine inliers, contamination points, retained points;
- dominant true group and inlier purity;
- paired core cluster for the same true group;
- centroid separation in true U0 and estimated Uhat;
- merge WCSS cost and within-true-group split silhouette;
- signatures: OUTLIER_ONLY, CONTAMINATION_DOMINANT, PURE_FRAGMENT,
  IMPURE_FRAGMENT, MIXED_CROSS_TRUTH.

For extra directions:
- exact directions inside span(Uhat) orthogonal to U0;
- R2 with true labels and with fitted labels;
- R2 with outlier status;
- random U0-perp true-label R2 q95 reference;
- alignment with leading orthogonal residual/noise PCs;
- dominant external association: OUTLIER_STATUS, TRUE_LABEL_SAMPLE_ASSOC, RANDOMLIKE.

For low ARI:
- truth-open merge ARI tells whether loss is largely recoverable by merging fragments;
- subspace coverage/min principal cosine quantify latent-space distortion;
- validation-draw DGP cluster/rank warnings distinguish an unlucky generated draw from estimator failure;
- failure_signature combines these mechanisms in one interpretable field.

PRIMARY OVER-SELECTION INTERPRETATION
The report explicitly gives, by method:
- mean and minimum ARI conditional on structural over-selection;
- fraction of overselected fits with ARI >= .95;
- low-ARI rate and positive regret.
This is a validation endpoint, not a selector rule.
