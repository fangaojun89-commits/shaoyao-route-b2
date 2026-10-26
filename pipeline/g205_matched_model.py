"""Plan, prepare, or explicitly fit the original G205 matched STAR s0 kernel."""
import argparse
import datetime
import hashlib
import json
import os
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CONTRACT = 'analysis/g205_matched_execution/CONTRACT.json'
ADAPTER = 'analysis/g205_matched_execution/adapter.R'
ENTRY = 'analysis/g205_matched_execution/execute.R'


def identity(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            digest.update(block)
    return {'bytes': path.stat().st_size, 'sha256': digest.hexdigest()}


def require_identity(path, expected, role):
    if identity(path) != expected:
        raise ValueError('Identity mismatch: ' + role)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--resource-root', type=Path, required=True)
    parser.add_argument('--rscript', type=Path, required=True)
    parser.add_argument('--r-library', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--prepare-only', action='store_true', help='Check runtime and original preparation; no fits or random draws')
    mode.add_argument('--execute', action='store_true', help='Explicitly run the unchanged matched STAR kernel once')
    args = parser.parse_args(argv)
    if args.output.exists():
        parser.error('output must be a new directory')
    contract = json.loads((REPO / CONTRACT).read_text(encoding='utf-8'))
    source_pins = dict(contract['sources'])
    source_pins[contract['prepared_identity_source']] = contract['prepared_identity_source_pin']
    for relative, pin in source_pins.items():
        require_identity(REPO / relative, pin, relative)
    bridge = json.loads((REPO / contract['prepared_identity_source']).read_text(encoding='utf-8'))
    if bridge['prepared_outputs']['prepared/G205_STAR_s0_inputs.rds'] != contract['prepared_identity']:
        raise ValueError('Prepared identity differs from accepted public input bridge')
    for relative in [CONTRACT, ADAPTER, ENTRY, 'pipeline/g205_matched_model.py']:
        source_pins[relative] = identity(REPO / relative)
    prepared = args.prepared.resolve(strict=True)
    require_identity(prepared, contract['prepared_identity'], 'STAR s0 prepared input')
    inputs = {prepared.as_posix(): identity(prepared)}
    resources = {}
    for role, item in contract['resources'].items():
        path = (args.resource_root / item['relative_path']).resolve(strict=True)
        require_identity(path, item['identity'], role)
        resources[role] = path.as_posix()
        inputs[path.as_posix()] = identity(path)
    rscript = args.rscript.resolve(strict=True)
    library = args.r_library.resolve(strict=True)
    if not rscript.is_file() or not library.is_dir():
        raise ValueError('Existing Rscript file and R library directory required')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output / 'receipts').mkdir()
    (output / 'results').mkdir()
    selected = 'execute' if args.execute else 'prepare' if args.prepare_only else 'plan'
    spec = dict(schema='G205_MATCHED_PORTABLE_RUN/v1', repo=REPO.as_posix(), output=output.as_posix(),
                mode=selected, prepared=prepared.as_posix(), prepared_identity=contract['prepared_identity'],
                resource_paths=resources, input_pins=inputs, source_pins=source_pins,
                contract=contract, historical_native_custody_recreated=False)
    run = output / 'receipts/RUN.json'
    write_json(run, spec)
    command = [str(rscript), '--vanilla', str(REPO / ENTRY), str(run)]
    record = dict(status='PLAN_ONLY_NO_R_PROCESS', mode=selected, command=command,
                  input_pins=inputs, source_pins=source_pins,
                  historical_native_custody_recreated=False,
                  numerical_equivalence_to_historical_results_established=False)
    parent = output / 'receipts/EXECUTION_PARENT.json'
    if selected == 'plan':
        write_json(parent, record)
        print(json.dumps({'status': record['status'], 'output': str(output)}))
        return
    env = os.environ.copy()
    env.pop('R_HOME', None)
    env.update(R_LIBS_USER=library.as_posix(), RENV_CONFIG_AUTOLOADER_ENABLED='FALSE',
               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', R_DATATABLE_NUM_THREADS='1')
    locale = 'English_United States.utf8' if os.name == 'nt' else 'C.UTF-8'
    env.update(LANG=locale, LC_ALL=locale, LC_CTYPE=locale)
    record.update(status='RUNNING_' + selected.upper(), started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    write_json(parent, record)
    try:
        # Original conditions stay in memory; retain only output fingerprints.
        process = subprocess.run(command, env=env, cwd=output, capture_output=True, check=False)
        unchanged = all(identity(Path(p)) == pin for p, pin in inputs.items())
        sources_unchanged = all(identity(REPO / p) == pin for p, pin in source_pins.items())
        qa_path = output / 'receipts/EXECUTION_QA.json'
        qa = json.loads(qa_path.read_text(encoding='utf-8')) if qa_path.exists() else {}
        expected = 'PREPARED_ORIGINAL_FUNCTIONS_NO_FIT' if selected == 'prepare' else 'NEW_NATIVE_RESULT_SAVED_NOT_NUMERICALLY_ACCEPTED'
        success = process.returncode == 0 and unchanged and sources_unchanged and qa.get('status') == expected
        record.update(status='PASS' if success else 'FAIL', returncode=process.returncode,
                      inputs_unchanged=unchanged, sources_unchanged=sources_unchanged,
                      stdout_bytes=len(process.stdout), stdout_sha256=hashlib.sha256(process.stdout).hexdigest(),
                      stderr_bytes=len(process.stderr), stderr_sha256=hashlib.sha256(process.stderr).hexdigest(),
                      raw_output_text_saved=False,
                      output_pins={p.relative_to(output).as_posix(): identity(p) for p in output.rglob('*')
                                   if p.is_file() and p != parent})
    except Exception as error:
        success = False
        record.update(status='FAIL', error_type=type(error).__name__)
    record['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    write_json(parent, record)
    print(json.dumps({'status': record['status'], 'output': str(output)}))
    if not success:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
