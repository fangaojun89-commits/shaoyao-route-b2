"""Small synthetic IO-contract tests; never launch R or read scientific inputs."""
import contextlib
import importlib.util
import io
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('matched_entry', ROOT / 'pipeline/g205_matched_model.py')
entry = importlib.util.module_from_spec(spec)
spec.loader.exec_module(entry)


class EntryTests(unittest.TestCase):
    def setUp(self):
        # Explicit data-disk test scratch; fixtures retained, no recursive deletion.
        self.base = Path(tempfile.mkdtemp(prefix='g205_entry_', dir=os.environ['G205_TEST_SCRATCH']))
        self.repo = self.base / 'repo'
        self.repo.mkdir()
        self.prepared = self.base / 'synthetic_STAR_s0.rds'
        self.prepared.write_bytes(b'synthetic test only; not a biological object')
        self.resource_root = self.base / 'resources'
        self.resource_root.mkdir()
        resource = self.resource_root / 'synthetic-resource'
        resource.write_bytes(b'synthetic resource')
        self.library = self.base / 'library'
        self.library.mkdir()
        self.rscript = self.base / 'NEVER_EXECUTE'
        self.rscript.write_bytes(b'not an executable')
        source = 'source.txt'
        (self.repo / source).write_bytes(b'original source fixture')
        bridge_name = 'bridge.json'
        entry.write_json(self.repo / bridge_name, {'prepared_outputs': {
            'prepared/G205_STAR_s0_inputs.rds': entry.identity(self.prepared)}})
        contract = dict(sources={source: entry.identity(self.repo / source)},
                        prepared_identity=entry.identity(self.prepared), prepared_identity_source=bridge_name,
                        prepared_identity_source_pin=entry.identity(self.repo / bridge_name),
                        resources={'test': {'relative_path': resource.name, 'identity': entry.identity(resource)}})
        for name in [entry.CONTRACT, entry.ADAPTER, entry.ENTRY, 'pipeline/g205_matched_model.py']:
            path = self.repo / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b'synthetic entry fixture')
        entry.write_json(self.repo / entry.CONTRACT, contract)
        self.output = self.base / 'new_output'
        self.args = ['--prepared', str(self.prepared), '--resource-root', str(self.resource_root),
                     '--rscript', str(self.rscript), '--r-library', str(self.library), '--output', str(self.output)]

    def call(self, extra=()):
        with mock.patch.object(entry, 'REPO', self.repo), mock.patch.object(entry.subprocess, 'run',
                side_effect=AssertionError('No child process is allowed in these tests')), \
                contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            entry.main(self.args + list(extra))

    def test_default_is_plan_without_R(self):
        self.call()
        result = json.loads((self.output / 'receipts/EXECUTION_PARENT.json').read_text())
        self.assertEqual(result['status'], 'PLAN_ONLY_NO_R_PROCESS')
        self.assertFalse(result['historical_native_custody_recreated'])
        self.assertFalse(result['numerical_equivalence_to_historical_results_established'])

    def test_wrong_prepared_rejected_before_output(self):
        self.prepared.write_bytes(b'other quantifier')
        with self.assertRaises(ValueError):
            self.call()
        self.assertFalse(self.output.exists())

    def test_wrong_resource_rejected_before_output(self):
        (self.resource_root / 'synthetic-resource').write_bytes(b'changed resource')
        with self.assertRaises(ValueError):
            self.call()
        self.assertFalse(self.output.exists())

    def test_wrong_source_rejected_before_output(self):
        (self.repo / 'source.txt').write_bytes(b'changed source')
        with self.assertRaises(ValueError):
            self.call()
        self.assertFalse(self.output.exists())

    def test_existing_output_refused(self):
        self.output.mkdir()
        sentinel = self.output / 'untouched'
        sentinel.write_bytes(b'preserve')
        with self.assertRaises(SystemExit):
            self.call()
        self.assertEqual(sentinel.read_bytes(), b'preserve')

    def test_modes_exclusive(self):
        with self.assertRaises(SystemExit):
            self.call(['--execute', '--prepare-only'])
        self.assertFalse(self.output.exists())


if __name__ == '__main__':
    unittest.main()
