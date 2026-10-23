"""Preflight tests only; no Windows installation or network is simulated as real."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('csda_preflight', ROOT/'tools/check_python_environment.py')
pf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pf)

class Preflight(unittest.TestCase):
    def test_actual_reference_environment_passes(self):
        self.assertTrue(pf.assess(ROOT/'requirements.txt')['ready'])

    def test_wrong_python_minor_is_rejected(self):
        with mock.patch.object(pf.sys, 'version_info', (3, 12, 0)):
            r = pf.assess(ROOT/'requirements.txt', base_only=True)
        self.assertFalse(r['ready'])
        self.assertTrue(any('3.13' in x for x in r['errors']))

    def test_nonvenv_installation_is_refused(self):
        with mock.patch.object(pf.sys, 'prefix', pf.sys.base_prefix):
            r = pf.assess(ROOT/'requirements.txt', base_only=True, require_venv=True)
        self.assertFalse(r['ready'])

    def test_free_threaded_variant_is_not_claimed_tested(self):
        with mock.patch.object(pf.sysconfig, 'get_config_var', return_value=1):
            r = pf.assess(ROOT/'requirements.txt', base_only=True)
        self.assertFalse(r['ready'])

    def test_missing_dependency_has_actionable_message(self):
        with mock.patch.object(pf.importlib.metadata, 'version', side_effect=pf.importlib.metadata.PackageNotFoundError('missing')):
            r = pf.assess(ROOT/'requirements.txt')
        self.assertFalse(r['ready'])
        self.assertTrue(any('matplotlib' in x for x in r['errors']))

    def test_wrong_version_is_not_silently_accepted(self):
        with mock.patch.object(pf.importlib.metadata, 'version', return_value='0.0.0'):
            r = pf.assess(ROOT/'requirements.txt')
        self.assertFalse(r['ready'])
        self.assertTrue(any('expected' in x for x in r['errors']))

    def test_malformed_requirements_are_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'requirements.txt'
            p.write_text('numpy>=1\n')
            r = pf.assess(p)
        self.assertFalse(r['ready'])

    def test_setup_contract_is_explicit_in_source_not_runtime_test(self):
        s = (ROOT/'tools/setup_windows.ps1').read_text()
        self.assertIn("'--no-bin', '--no-registry'", s)
        self.assertIn('Get-FileHash -LiteralPath $archive -Algorithm SHA256', s)
        self.assertIn("'UV_PYTHON_INSTALL_REGISTRY' = '0'", s)
        self.assertIn('if (-not $Yes)', s)
        self.assertNotIn('Set-ExecutionPolicy', s)
        self.assertNotIn('winget install', s)
        self.assertNotIn('RUN_REPLICATION("FULL"', s)
        self.assertNotIn('Invoke-Expression', s)
        self.assertNotIn('Remove-Item -LiteralPath $Venv', s)

    def test_matlab_setup_is_separate_and_requires_consent(self):
        s = (ROOT/'RUN_REPLICATION.m').read_text()
        self.assertIn('if mode=="SETUP"', s)
        self.assertIn("assert(o.Confirm,'SETUP downloads", s)
        self.assertIn('ReportingSetupMissing', s)

if __name__ == '__main__':
    unittest.main()
