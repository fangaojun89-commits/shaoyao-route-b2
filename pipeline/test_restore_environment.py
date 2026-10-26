"""Generated command stubs test installation control flow, not a native R restore."""
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

spec = importlib.util.spec_from_file_location('restore_environment', Path(__file__).with_name('restore_environment.py'))
restore = importlib.util.module_from_spec(spec)
spec.loader.exec_module(restore)


class RestoreControlFlow(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.output = Path(self.temporary.name) / 'new-library-root'
        self.plan = {'output': str(self.output), 'r': '/provided/bin/R',
                     'packages': ['first', 'second'], 'commands': [['install-first'], ['install-second']],
                     'acceptance_command': ['check-load'], 'scientific_model_runs': 0}

    def test_existing_output_preserved(self):
        self.output.mkdir()
        marker = self.output / 'retained.txt'
        marker.write_text('keep')
        with self.assertRaises(FileExistsError):
            restore.install(self.plan, run=lambda *a, **k: self.fail('Must not start'))
        self.assertEqual(marker.read_text(), 'keep')

    def test_build_failure_stops_before_later_packages(self):
        calls = []
        def failing(command, **kwargs):
            calls.append(command)
            return SimpleNamespace(returncode=3)
        with self.assertRaisesRegex(RuntimeError, 'Stage failed'):
            restore.install(self.plan, run=failing)
        self.assertEqual(calls, [['install-first']])
        receipt = json.loads((self.output / 'receipt.json').read_text())
        self.assertEqual(receipt['status'], 'INCOMPLETE_FILES_RETAINED')
        self.assertEqual(receipt['completed_packages'], [])
        self.assertTrue((self.output / 'logs/first.stderr').exists())

    def test_success_requires_acceptance_artifact(self):
        with self.assertRaises(FileNotFoundError):
            restore.install(self.plan, run=lambda *a, **k: SimpleNamespace(returncode=0))
        self.assertEqual(json.loads((self.output / 'receipt.json').read_text())['status'], 'INCOMPLETE_FILES_RETAINED')

    def test_failed_load_check_is_not_success(self):
        def fake(command, **kwargs):
            return SimpleNamespace(returncode=2 if command == ['check-load'] else 0)
        with self.assertRaisesRegex(RuntimeError, 'verify_library'):
            restore.install(self.plan, run=fake)
        receipt = json.loads((self.output / 'receipt.json').read_text())
        self.assertEqual(receipt['completed_packages'], ['first', 'second'])
        self.assertEqual(receipt['status'], 'INCOMPLETE_FILES_RETAINED')

    def test_serial_build_isolated_environment_and_final_load(self):
        calls = []
        def fake(command, **kwargs):
            calls.append(command)
            self.assertEqual(kwargs['env']['R_LIBS_USER'], str(self.output / 'library'))
            self.assertEqual(kwargs['env']['R_LIBS_SITE'], 'NULL')
            self.assertNotIn('R_HOME', kwargs['env'])
            self.assertEqual(kwargs['env']['MAKEFLAGS'], '-j1')
            self.assertEqual(kwargs['env']['TMPDIR'], str(self.output / 'tmp'))
            if command == ['check-load']:
                (self.output / 'acceptance').mkdir()
                (self.output / 'acceptance/RESULT.json').write_text(json.dumps({
                    'status': '70_INSTALLED_VERSIONS_AND_EIGHT_ROOT_LOADS_PASS'}))
            return SimpleNamespace(returncode=0)
        receipt = restore.install(self.plan, run=fake)
        self.assertEqual(calls, [['install-first'], ['install-second'], ['check-load']])
        self.assertEqual(receipt['status'], 'RESTORED_CORE_LOAD_CHECK_PASS')
        self.assertEqual(receipt['scientific_model_runs'], 0)


if __name__ == '__main__':
    unittest.main()
