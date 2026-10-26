"""Small control-flow tests with synthetic files and a mocked R subprocess."""
import argparse
import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch


def run_tests(repo, output):
    spec = importlib.util.spec_from_file_location('g205_model', repo / 'pipeline/g205_model.py')
    launcher = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(launcher)
    tests = []
    with tempfile.TemporaryDirectory(prefix='g205_synthetic_', dir=output) as temp:
        root = Path(temp)
        synthetic = root / 'repo'
        sources = ['analysis/g205_preparation/prepare.R',
                   'analysis/study_sources/G205_SOURCE_ROLES.json',
                   'pipeline/g205_model.py', 'analysis/g205_execution/execute.R',
                   'analysis/g205_execution/write_result.R', 'worker.R', 'full.json', 'baseline.json']
        for name in sources:
            path = synthetic / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('SYNTHETIC_NO_SCIENTIFIC_CONTENT\n', encoding='utf-8')
        prepared = root / 'input.rds'
        prepared.write_bytes(b'SYNTHETIC_INPUT')
        library = root / 'library'
        library.mkdir()
        executable = root / 'Rscript.mock'
        executable.write_bytes(b'NOT_EXECUTABLE_ONLY_MOCK')
        contract_name = 'analysis/g205_preparation/CONTRACT.json'
        contract = dict(sources={x: launcher.identity(synthetic / x)
                                for x in ['worker.R', 'full.json', 'baseline.json']},
                        prepared_identity=launcher.identity(prepared), resources={},
                        worker='worker.R', full_config='full.json', baseline_config='baseline.json')
        launcher.write_json(synthetic / contract_name, contract)
        launcher.write_json(synthetic / 'MANIFEST.json', {'files': {
            x: launcher.identity(synthetic / x) for x in
            [contract_name, 'analysis/g205_preparation/prepare.R',
             'analysis/study_sources/G205_SOURCE_ROLES.json']}})
        launcher.REPO = synthetic

        def invoke(name, execute=False):
            args = ['g205_model.py', '--prepared', str(prepared), '--resource-root', str(root),
                    '--rscript', str(executable), '--r-library', str(library),
                    '--output', str(root / name)] + (['--execute'] if execute else [])
            with patch.object(sys, 'argv', args):
                launcher.main()

        with patch.object(launcher.subprocess, 'run') as child:
            invoke('plan')
            child.assert_not_called()
        assert json.loads((root / 'plan/receipts/EXECUTION_PARENT.json').read_text())['status'] == 'PLAN_ONLY_NO_R_PROCESS'
        tests.append('default_plan_starts_no_R_process')

        def fake_success(command, **kwargs):
            assert kwargs['env']['OMP_NUM_THREADS'] == '1'
            assert kwargs['env']['R_LIBS_USER'] == library.as_posix()
            assert command[1] == '--vanilla'
            launcher.write_json(kwargs['cwd'] / 'receipts/EXECUTION_QA.json',
                                {'status': 'PASS_ORIGINAL_CONSUMER_AND_OUTPUT_READBACK'})
            return subprocess.CompletedProcess(command, 0)

        with patch.object(launcher.subprocess, 'run', side_effect=fake_success) as child:
            invoke('execution', execute=True)
            assert child.call_count == 1
        result = json.loads((root / 'execution/receipts/EXECUTION_PARENT.json').read_text())
        assert result['status'] == 'PASS' and result['inputs_unchanged'] and result['sources_unchanged']
        tests.append('explicit_execution_mock_command_environment_and_receipt')

        with patch.object(launcher.subprocess, 'run', return_value=subprocess.CompletedProcess([], 2)):
            try:
                invoke('child_failure', execute=True)
                raise AssertionError('Failed child was accepted')
            except SystemExit as error:
                assert error.code == 1
        assert json.loads((root / 'child_failure/receipts/EXECUTION_PARENT.json').read_text())['status'] == 'FAIL'
        tests.append('failed_child_rejected')

        with patch.object(launcher.subprocess, 'run', side_effect=OSError('synthetic launch error')):
            try:
                invoke('launch_failure', execute=True)
                raise AssertionError('Launch failure was accepted')
            except OSError:
                pass
        assert json.loads((root / 'launch_failure/receipts/EXECUTION_PARENT.json').read_text())['status'] == 'FAIL_TO_LAUNCH_R'
        tests.append('launch_failure_recorded')

        try:
            invoke('plan')
            raise AssertionError('Existing output was accepted')
        except SystemExit as error:
            assert error.code == 2
        tests.append('existing_output_rejected')
        prepared.write_bytes(b'ALTERED_SYNTHETIC_INPUT')
        with patch.object(launcher.subprocess, 'run') as child:
            try:
                invoke('bad_input', execute=True)
                raise AssertionError('Wrong input was accepted')
            except ValueError:
                pass
            child.assert_not_called()
        assert not (root / 'bad_input').exists()
        tests.append('incorrect_input_rejected_before_output_or_process')
    launcher.write_json(output / 'RESULT.json', dict(status='PASS_SYNTHETIC_LAUNCHER',
                        tests=tests, real_R_subprocesses=0, scientific_fits=0))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('repo', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    run_tests(args.repo.resolve(), args.output.resolve())
