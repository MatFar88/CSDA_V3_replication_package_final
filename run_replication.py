#!/usr/bin/env python3
"""Single Python entry point: CHECK and REPORT never fit an estimator."""
from __future__ import annotations
import argparse,datetime,json,platform,shutil,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',type=str.upper,choices=['CHECK','REPORT','FULL'])
    p.add_argument('--output',type=Path,help='New empty output directory; references are never overwritten.')
    p.add_argument('--fresh-root',type=Path,help='Optional complete replayed fresh output tree (REPORT).')
    p.add_argument('--no-figures',action='store_true',help='Recompute all table/figure values without plotting.')
    p.add_argument('--strict-complete',action='store_true',help='Return nonzero while any published experiment is unavailable.')
    p.add_argument('--fixed-results',type=Path,help='FINAL32 FIXED_all_results.csv (192 rows); bundled by default.')
    p.add_argument('--matlab',default='matlab',help='MATLAB executable for FULL.')
    p.add_argument('--confirm-expensive-run',action='store_true')
    args=p.parse_args(argv)
    try:
        from reproduce.common import ReplicationError,verify_integrity
        verified=verify_integrity(ROOT)
        from reproduce.optional import resolve_inputs
        args.fixed_results=resolve_inputs(ROOT,args.fixed_results)
        if args.mode=='CHECK':
            import numpy,pandas,scipy,matplotlib
            info={'integrity_files_verified':verified,'python':platform.python_version(),
                  'numpy':numpy.__version__,'pandas':pandas.__version__,'scipy':scipy.__version__,
                  'matplotlib':matplotlib.__version__,'matlab_on_path':bool(shutil.which(args.matlab)),
                  'REPORT_available':True,'FULL_MATLAB_runtime_tested_in_build':False,
                  'V1_SELFTEST_and_REPORT_passed_on_author_PC':True,
                  'manuscript_target':'CSDA submission 2026-10-23',
                  'fixed_benchmark':'FINAL32: 32 datasets / 192 fits',
                  'fixed32_original_integrated':args.fixed_results is not None,
                  'missing_original_result_blocks':(['FINAL32 fixed results'] if args.fixed_results is None else []),
                  'confidential_banking_microdata_included':False}
            print(json.dumps(info,indent=2))
            return 2 if args.strict_complete and info['missing_original_result_blocks'] else 0
        if args.mode=='FULL':
            if not args.confirm_expensive_run:
                raise ReplicationError('FULL requires --confirm-expensive-run. It may run 369000 candidate fits with 100 starts each. Use MATLAB stage commands for resumable batches.')
            if args.output is None:raise ReplicationError('FULL requires an explicit --output workspace.')
            executable=shutil.which(args.matlab)
            if not executable:raise ReplicationError('MATLAB executable not found. CHECK/REPORT do not need MATLAB.')
            work=args.output.resolve()
            q=lambda s:str(s).replace("'","''")
            command=f"addpath('{q(ROOT)}'); RUN_REPLICATION('FULL','Workspace','{q(work)}','Confirm',true);"
            return subprocess.run([executable,'-batch',command],check=False).returncode
        from reproduce.report import run_report
        dest=args.output or ROOT/'output'/('report_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
        dest=dest.resolve()
        for protected in [ROOT/'results',ROOT/'code',ROOT/'reproduce',ROOT/'documentation',ROOT/'config',ROOT/'tests',ROOT/'validation']:
            if dest==protected or protected in dest.parents:raise ReplicationError('Output points inside a protected input/code directory.')
        if dest==ROOT or dest in ROOT.parents:raise ReplicationError('Output cannot be the package directory or an ancestor.')
        if args.fresh_root:
            fresh=args.fresh_root.resolve()
            if dest==fresh or dest in fresh.parents or fresh in dest.parents:raise ReplicationError('REPORT output must be separate from the replayed fresh inputs.')
        info=run_report(ROOT,dest,fresh=args.fresh_root,figures=not args.no_figures,fixed_results=args.fixed_results)
        print(f"REPORT: {info['checks_passed']} checks passed. Manuscript display values verified.")
        print('No estimator refitted. Banking inputs: aggregates only.')
        if info['pending']:print('Pending original results: '+', '.join(info['pending']))
        else:print('All numerical table blocks reconstructed from supplied results for the CSDA submission manuscript. No FULL rerun certified.')
        print('Open: '+str(dest/'REPORT.html'))
        return 2 if args.strict_complete and not info['all_paper_experiments_complete'] else 0
    except (OSError,ValueError,RuntimeError,ImportError,subprocess.SubprocessError) as exc:
        print(f'ERROR: {exc}',file=sys.stderr)
        return 1

if __name__=='__main__':sys.exit(main())
