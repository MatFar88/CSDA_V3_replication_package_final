# Native MATLAB reproduction: explicit, staged, resumable

## What is and is not tested

The Python REPORT path was executed in the build environment. MATLAB was not installed there. Original `.m` bytes are preserved; the new orchestration functions were inspected statically, not executed. Begin with SELFTEST and FIT_SMOKE in MATLAB, not with the full computation. The original fitting/generator source calls `kmeans` and `pdist`, so its Statistics and Machine Learning Toolbox dependency is checked. No exact historical MATLAB release/hardware is inferred; each new execution logs the actual release and toolboxes.

Use a short local path (for example `D:\CSDA_REPLAY`) and a normal JVM-enabled MATLAB session. Do not add every source directory recursively. Do not point Workspace to the original historical result directory. PREPARE refuses an unrecognized nonempty directory, stages original configuration/dependencies, and checks unchanged staged files on resume.

## Fresh 41-cell study

```matlab
cd("C:\CSDA_REPLICATION_PACKAGE_V1");
RUN_REPLICATION("SELFTEST");
RUN_REPLICATION("FULL","Workspace","D:\CSDA_REPLAY", ...
                "Stage","PREPARE","Confirm",true);
RUN_REPLICATION("FULL","Workspace","D:\CSDA_REPLAY", ...
                "Stage","DGP_SMOKE","Confirm",true);
RUN_REPLICATION("FULL","Workspace","D:\CSDA_REPLAY", ...
                "Stage","FIT_SMOKE","Confirm",true);
```

DGP_SMOKE generates the 41 first-replicate synthetic datasets but does not fit them. FIT_SMOKE performs exactly one fixed candidate `(r,g,alpha)=(1,2,0)` for each method using the original 100-start allocations, writes to `replication_smoke/`, and does not populate production candidate cubes. It is a runtime test, not another scientific experiment.

### Parallel fit and selection

After preparation/smoke succeeds, use six independent MATLAB sessions, with `Worker` from 1 to 6. For example in session 1:

```matlab
RUN_REPLICATION("FULL","Workspace","D:\CSDA_REPLAY", ...
  "Stage","FIT","Worker",1,"Workers",6,"Confirm",true);
```

Each assigned dataset runs TFKM and TRKM. Change only Worker in the other sessions. Do not run the same worker twice. Keep Workers unchanged when resuming a stage. Checkpoints are the original runners' candidate CSV/cube MAT and COMPLETE markers. Interrupted sessions resume without recomputing stored candidates. A leftover `.lock` file after a process crash must be removed manually only after confirming that its session is no longer active.

After **all** FIT sessions finish, run SELECT with the same six worker IDs:

```matlab
RUN_REPLICATION("FULL","Workspace","D:\CSDA_REPLAY", ...
  "Stage","SELECT","Worker",1,"Workers",6,"Confirm",true);
```

SELECT uses the original method-specific finalizers. Backward propagation may regenerate the relevant start pools, so SELECT is not a no-refit reporting operation.

After all SELECT workers finish, run once:

```matlab
RUN_REPLICATION("FULL","Workspace","D:\CSDA_REPLAY", ...
  "Stage","SEAL","Confirm",true);
RUN_REPLICATION("FULL","Workspace","D:\CSDA_REPLAY", ...
  "Stage","EVALUATE","Confirm",true);
```

SEAL is collection only. The adapter collects TRKM per-unit results without requiring obsolete SIGN-only TFKM results, then calls the original corrected-TFKM merger. It checks 1,230 prediction keys per method and rejects truth columns. It snapshots the regenerated predictions before evaluation. The legacy directory name `FROZEN_PRETRUTH_20261010` is kept only because original evaluators expect it; a new `REPLICATION_PREDICTIONS_SEALED.json` states the current timestamp. This is replication of an already opened design, not a new historical seal or an independent holdout experiment.

EVALUATE uses the original truth-open evaluator, posthoc decomposition and no-refit diagnostics. The original positive-action evaluation reads the actual `FinalLabels`/`FinalTrim` in `blind_final.mat`. It does not substitute the initial cube representative at the same `(r,g,alpha)`.

### Diagnostic reproduction

Run six DIAGNOSTICS workers after EVALUATE has completed:

```matlab
RUN_REPLICATION("FULL","Workspace","D:\CSDA_REPLAY", ...
  "Stage","DIAGNOSTICS","Worker",1,"Workers",6,"Confirm",true);
```

This includes the original theory, closure, invariant and WEAK100-start-pool calculations. The WEAK component does refit/replay pools; it is explicitly a diagnostic stage. After all six finish:

```matlab
RUN_REPLICATION("FULL","Workspace","D:\CSDA_REPLAY", ...
  "Stage","MERGE_DIAGNOSTICS","Confirm",true);
```

Then compare reconstructed outputs to the paper's numerical targets:

```text
python run_replication.py REPORT --fresh-root "D:\CSDA_REPLAY\TFKM_TRKM_FRESH_GENERALIZATION_V1_output" --output "D:\CSDA_REPLAY_REPORT"
```

REPORT is intentionally strict about expected counts and relevant values. Inspect discrepancies; do not retune thresholds to remove them. Plot/PDF binary hashes need not match because plotting metadata and renderers differ.

For a fully serial run, `Stage="ALL"` executes the above production stages in order with Workers=1. FIT_SMOKE is deliberately separate. ALL requires explicit confirmation and may be very expensive: the fitting grid has 369,000 candidate-method-dataset combinations with 100 starts each, before backward and diagnostic replays. No runtime estimate is inferred from incomplete historical timing information.

A limited production trial can use `Stage="FIT", "MaxUnits",1`; unlike FIT_SMOKE, that trial computes both full 150-candidate grids for one assigned dataset. It is not a fast smoke test. MaxUnits is accepted only for FIT/SELECT, and partial outputs cannot be sealed as a complete study.

## Configuration and source isolation

The distributed original `study.mat`, `STUDY_CATALOG.csv` and `UNIT_MANIFEST.csv` are copied without changes. The fresh catalogue/CFG/manifest determine n=120,41cells,30replicates. Unused legacy defaults inside the options template are not used to redefine this design. Original pair seeds and all candidate rules remain unchanged.

The source storage includes both original method dependency trees; runtime staging recreates their expected paths. Common files were found byte-identical across the recovered TFKM and TRKM trees. The original path-lock functions still enforce which tree resolves each method. The wrapper restores the caller's MATLAB path on return. It does not invoke an external mutable installation or download code.


## Separate FINAL32 fixed-parameter comparison

Saved-result Table S4 is fully covered by the bundled 192-row
`results/reference/fixed32/FIXED_all_results.csv` and its cross-check files.
It is independent of the 41-cell FRESH run. The original producer source is
included as an isolated ZIP under code/matlab/experiments/fixed32; consult
`documentation/FIXED32.md` for its saved-context prerequisites. The FRESH FULL
wrapper does not recreate this different experiment, execute its earlier
PS/BQS stages or turn an unrelated catalogue into fixed-parameter evidence.

