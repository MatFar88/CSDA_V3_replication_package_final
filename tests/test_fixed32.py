"""Tests for FINAL32 integration. Modified fixtures are never reported as evidence."""
import json,sys,tempfile,unittest,shutil
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from reproduce.fixed32 import fixed_tables
from reproduce.common import Checks,ReplicationError

class Fixed32Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.file=ROOT/'results/reference/fixed32/FIXED_all_results.csv'
        cls.expected=json.loads((ROOT/'config/expected_paper.json').read_text())
        cls.checks=Checks();cls.tables=fixed_tables(ROOT,cls.file,cls.checks,cls.expected)
    def test_all_worker_and_summary_checks_pass(self):
        self.assertGreater(len(self.checks.rows),90)
        self.assertTrue(self.checks.frame().status.eq('PASS').all())
    def test_six_methods_and_correct_denominators(self):
        t=self.tables['Supplement_TableS4'];self.assertEqual(len(t),6);self.assertTrue(t.N.eq(32).all())
    def test_not_relabelled_as_old_experiment(self):
        self.assertFalse((ROOT/'results/additional').exists())
        self.assertFalse((ROOT/'reproduce/received_catalogues.py').exists())
        self.assertEqual(json.loads((ROOT/'config/result_inputs.json').read_text())['fixed32'], 'results/reference/fixed32/FIXED_all_results.csv')
    def test_all_inlier_not_retained_only_ari(self):
        t=self.tables['Supplement_TableS4'].set_index('method')
        self.assertAlmostEqual(t.loc['REG-HD-TFKM','Mean_ARI'],.983640987189653,12)
        self.assertLess(t.loc['REG-HD-TFKM','Mean_ARI'],1)
    def test_medians_are_calculated_from_individual_fits(self):
        t=self.tables['Supplement_TableS4'];self.assertTrue(t.Median_ARI.eq(1).all())
    def bad_input(self,mutator,regex):
        t=pd.read_csv(self.file);mutator(t)
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/'changed.csv';t.to_csv(f,index=False)
            with self.assertRaisesRegex(ReplicationError,regex):fixed_tables(ROOT,f,Checks(),self.expected)
    def test_wrong_experiment_rejected(self):
        self.bad_input(lambda t:t.__setitem__('scenario_id',t.scenario_id+1),'four declared scenarios')
    def test_changed_original_metric_rejected(self):
        self.bad_input(lambda t:t.loc.__setitem__((0,'ari_inliers_full'),.5),'matches merged records')
    def test_invalid_fit_rejected(self):
        self.bad_input(lambda t:t.loc.__setitem__((0,'valid'),0),'all fits valid')
    def test_changed_method_identity_rejected(self):
        self.bad_input(lambda t:t.loc.__setitem__((0,'method'),'WRONG_SELECTOR'),'six original method')
    def test_display_targets_checked(self):
        exp=json.loads(json.dumps(self.expected));exp['table_s4']['TKM'][0]=.1845
        with self.assertRaises(ReplicationError):fixed_tables(ROOT,self.file,Checks(),exp)

if __name__=='__main__':unittest.main()
