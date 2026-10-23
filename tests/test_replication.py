"""Fast tests. No MATLAB calls, estimator fitting, or confidential data needed."""
from __future__ import annotations
import hashlib,json,subprocess,sys,tempfile,unittest
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from reproduce.common import ReplicationError,Checks,ari_from_contingency,auc_rank,load_csv,sha256,verify_integrity
from reproduce.report import run_report

class Statistics(unittest.TestCase):
    def test_ari_same_partition(self):self.assertEqual(ari_from_contingency([[3,0],[0,4]]),1.)
    def test_ari_label_permutation(self):self.assertEqual(ari_from_contingency([[0,3],[4,0]]),1.)
    def test_ari_identical_single_group(self):self.assertEqual(ari_from_contingency([[7]]),1.)
    def test_ari_single_group_vs_split(self):self.assertEqual(ari_from_contingency([[3,4]]),0.)
    def test_ari_ratios(self):self.assertAlmostEqual(ari_from_contingency([[0,284],[5,0],[0,2]]),.8274906143503161)
    def test_ari_log(self):self.assertAlmostEqual(ari_from_contingency([[1,168],[0,104]]),-.0028153763)
    def test_invalid_contingency(self):
        for x in [[[1.5,2]],[[1,-1]],[[float('nan'),2]],[[1]]]:
            with self.assertRaises(ReplicationError):ari_from_contingency(x)
    def test_auc_perfect(self):self.assertEqual(auc_rank([0,0,1,1],[0,1,2,3]),1.)
    def test_auc_reversed(self):self.assertEqual(auc_rank([0,0,1,1],[3,2,1,0]),0.)
    def test_auc_ties(self):self.assertEqual(auc_rank([0,1,0,1],[1,1,1,1]),.5)
    def test_auc_bad_class(self):
        with self.assertRaises(ReplicationError):auc_rank([1,1],[2,3])
    def test_source_integrity(self):self.assertGreater(verify_integrity(ROOT),300)
    def test_tamper_detection(self):
        with tempfile.TemporaryDirectory() as d:
            r=Path(d);(r/'documentation').mkdir();(r/'x').write_text('reference')
            (r/'documentation/SHA256SUMS.json').write_text(json.dumps({'x':sha256(r/'x')}))
            self.assertEqual(verify_integrity(r),1);(r/'x').write_text('changed')
            with self.assertRaises(ReplicationError):verify_integrity(r)
    def test_schema_rejects_duplicate_keys(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x.csv';pd.DataFrame({'k':[1,1],'x':[2,3]}).to_csv(p,index=False)
            with self.assertRaises(ReplicationError):load_csv(p,['k'],2,['k'])
    def test_schema_rejects_missing_columns(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x.csv';pd.DataFrame({'k':[1]}).to_csv(p,index=False)
            with self.assertRaises(ReplicationError):load_csv(p,['absent'])
    def test_schema_rejects_missing_file(self):
        with self.assertRaises(ReplicationError):load_csv(ROOT/'not_a_file.csv')
    def test_rounding_is_not_silently_passed(self):
        c=Checks()
        with self.assertRaises(ReplicationError):c.published('test',40.7848954412,40.79,2)
        c.published('corrected',40.7848954412,40.78,2)
        self.assertEqual(c.rows[-1]['status'],'PASS')

class PackageContract(unittest.TestCase):
    def test_original_source_hashes(self):
        t=pd.read_csv(ROOT/'documentation/SOURCE_FILES.csv')
        for r in t.itertuples(index=False):self.assertEqual(sha256(ROOT/r.package_path),r.sha256,r.package_path)
    def test_banking_allowlist_has_no_person_ids(self):
        forbidden={'bank_index','bank_name','bank_id','original_365_row','original_364_row','original_row_index'}
        for f in (ROOT/'results/reference/banking').rglob('*'):
            if f.is_file():
                self.assertEqual(f.suffix,'.csv')
                self.assertFalse(forbidden&set(pd.read_csv(f,nrows=0).columns),f.name)
    def test_no_action_completion_not_distributed_as_stage(self):
        self.assertEqual(list((ROOT/'code').rglob('RECOVER_SELECTED_ARI*')),[])
    def test_correct_tfkm_stage_in_wrapper(self):
        s=(ROOT/'RUN_REPLICATION.m').read_text()
        self.assertIn('RUN_TFKM_FRESH_HIER4_RERUN_UNIT_V1(sid,rep,workspace)',s)
        self.assertNotIn('RUN_TFKM_TRKM_FRESH_UNIT_FIT_V1(sid,rep,"REG_TFKM_HD"',s)
    def test_positive_terminal_source_is_original_finalizer(self):
        s=(ROOT/'code/matlab/stages/evaluate/RUN_TFKM_TRKM_FRESH_TRUTH_OPEN_V1.m').read_text()
        self.assertIn('FinalLabels',s);self.assertIn('blind_final.mat',s)
    def test_full_requires_confirmation(self):
        p=subprocess.run([sys.executable,str(ROOT/'run_replication.py'),'FULL'],capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0);self.assertIn('confirm-expensive-run',p.stderr)
    def test_existing_output_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p=subprocess.run([sys.executable,str(ROOT/'run_replication.py'),'REPORT','--output',d,'--no-figures'],capture_output=True,text=True)
            self.assertNotEqual(p.returncode,0);self.assertIn('already exists',p.stderr)
    def test_protected_output_is_refused(self):
        p=subprocess.run([sys.executable,str(ROOT/'run_replication.py'),'REPORT','--output',str(ROOT/'results/reference/do_not_write'),'--no-figures'],capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0);self.assertFalse((ROOT/'results/reference/do_not_write').exists())

class NumericalRegression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory();cls.a=Path(cls.temp.name)/'a';cls.b=Path(cls.temp.name)/'b'
        cls.info=run_report(ROOT,cls.a,figures=False);run_report(ROOT,cls.b,figures=False)
    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()
    def test_zero_refits(self):self.assertEqual(self.info['estimator_refits'],0)
    def test_all_required_available_checks(self):self.assertGreaterEqual(self.info['checks_passed'],150)
    def test_manuscript_values_verified(self):self.assertTrue(self.info['manuscript_values_verified'])
    def test_all_saved_result_blocks_available(self):
        self.assertTrue(self.info['all_paper_experiments_complete'])
        self.assertFalse(self.info['full_replay_executed'])
    def test_repeatability_csv(self):
        for f in (self.a/'tables').glob('*.csv'):self.assertEqual(f.read_bytes(),(self.b/'tables'/f.name).read_bytes(),f.name)
    def test_joint_counts(self):
        t=pd.read_csv(self.a/'tables/Main_Table3.csv');self.assertEqual(list(t.Joint_correct_count),[662,866])
    def test_no_aggregate_final_ari(self):
        t=pd.read_csv(self.a/'tables/Main_Table3.csv');self.assertNotIn('final_inlier_ARI_POSTHOC',t.columns)
    def test_weak_zero_availability(self):
        t=pd.read_csv(self.a/'tables/Figure3_values.csv');q=t[(t.p_over_n==2)&(t.rho==2)]
        self.assertTrue(q.Available.eq(0).all());self.assertTrue(q.Projected_truth_assisted_exact.eq(1).all())
    def test_all_current_blocks_complete(self):
        self.assertEqual(self.info['pending'],[])
    def test_bank_source_level_is_explicit(self):
        t=pd.read_csv(self.a/'tables/Main_Table6.csv');self.assertTrue(t.Input_level.str.contains('saved aggregate overlap').all())


if __name__=='__main__':unittest.main()
