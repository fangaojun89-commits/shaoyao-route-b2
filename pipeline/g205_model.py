"""Plan or explicitly execute the unchanged G205 primary checked consumer."""
import argparse
import datetime
import hashlib
import json
import os
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def identity(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            digest.update(block)
    return {'bytes': path.stat().st_size, 'sha256': digest.hexdigest()}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--resource-root', type=Path, required=True)
    parser.add_argument('--rscript', required=True)
    parser.add_argument('--r-library', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--execute', action='store_true', help='Run the original real-data consumer and its independent numerical references; omitted by default')
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output must be a new directory')
    manifest = json.loads((REPO / 'MANIFEST.json').read_text(encoding='utf-8'))
    contract_path = 'analysis/g205_preparation/CONTRACT.json'
    sources = {}
    for name in [contract_path, 'analysis/g205_preparation/prepare.R', 'analysis/study_sources/G205_SOURCE_ROLES.json']:
        actual = identity(REPO / name)
        if actual != manifest['files'][name]:
            raise ValueError('Released preparation identity mismatch: ' + name)
        sources[name] = actual
    contract = json.loads((REPO / contract_path).read_text(encoding='utf-8'))
    for name, expected in contract['sources'].items():
        if identity(REPO / name) != expected:
            raise ValueError('Original source/configuration mismatch: ' + name)
        sources[name] = expected
    for name in ['pipeline/g205_model.py', 'analysis/g205_execution/execute.R', 'analysis/g205_execution/write_result.R']:
        sources[name] = identity(REPO / name)
    prepared = args.prepared.resolve(strict=True)
    if identity(prepared) != contract['prepared_identity']:
        raise ValueError('Prepared HISAT2 s0 object mismatch')
    inputs = {str(prepared): identity(prepared)}
    resources = {}
    for role, item in contract['resources'].items():
        path = ((REPO / 'analysis/g205_preparation/gene_set_manifest.tsv') if role == 'gene_set_manifest'
                else (args.resource_root / item['relative_path'])).resolve(strict=True)
        if identity(path) != item['identity']:
            raise ValueError('Frozen resource mismatch: ' + role)
        resources[role] = path.as_posix()
        inputs[str(path)] = identity(path)
    library = args.r_library.resolve(strict=True)
    if not library.is_dir():
        raise ValueError('R library must be an existing directory')
    rscript = Path(args.rscript).resolve(strict=True)
    if not rscript.is_file():
        raise ValueError('Rscript must be an existing executable file')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output / 'receipts').mkdir()
    (output / 'derived').mkdir()
    spec = dict(repo=REPO.as_posix(), output=output.as_posix(), prepared=prepared.as_posix(),
                prepared_identity=contract['prepared_identity'], resource_paths=resources,
                source_pins=sources, worker=contract['worker'], full_config=contract['full_config'],
                baseline_config=contract['baseline_config'])
    spec_path = output / 'receipts/RUN.json'
    write_json(spec_path, spec)
    command = [str(rscript), '--vanilla', str(REPO / 'analysis/g205_execution/execute.R'), str(spec_path)]
    record = dict(status='PLAN_ONLY_NO_R_PROCESS', mode='execute' if args.execute else 'plan',
                  command=command, input_pins=inputs, source_pins=sources,
                  historical_native_custody_recreated=False,
                  numerical_equivalence_to_historical_results_established=False)
    if not args.execute:
        write_json(output / 'receipts/EXECUTION_PARENT.json', record)
        print(json.dumps({'status': record['status'], 'output': str(output)}))
        return
    env = os.environ.copy()
    env.pop('R_HOME', None)
    env.update(R_LIBS_USER=library.as_posix(), RENV_CONFIG_AUTOLOADER_ENABLED='FALSE',
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    locale = 'English_United States.utf8' if os.name == 'nt' else 'C.UTF-8'
    env.update(LANG=locale, LC_ALL=locale, LC_CTYPE=locale)
    record['started_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    record['status'] = 'EXECUTING_ORIGINAL_CONSUMER'
    write_json(output / 'receipts/EXECUTION_PARENT.json', record)
    try:
        with (output / 'receipts/execution.stdout.txt').open('wb') as stdout, (output / 'receipts/execution.stderr.txt').open('wb') as stderr:
            process = subprocess.run(command, env=env, cwd=output, stdout=stdout, stderr=stderr, check=False)
    except OSError as error:
        record.update(status='FAIL_TO_LAUNCH_R', error_type=type(error).__name__,
                      finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        write_json(output / 'receipts/EXECUTION_PARENT.json', record)
        raise
    unchanged = all(identity(Path(path)) == pin for path, pin in inputs.items())
    sources_unchanged = all(identity(REPO / path) == pin for path, pin in sources.items())
    qa_path = output / 'receipts/EXECUTION_QA.json'
    qa = json.loads(qa_path.read_text(encoding='utf-8')) if qa_path.exists() else {}
    success = process.returncode == 0 and unchanged and sources_unchanged and qa.get('status') == 'PASS_ORIGINAL_CONSUMER_AND_OUTPUT_READBACK'
    record.update(status='PASS' if success else 'FAIL', returncode=process.returncode,
                  finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  inputs_unchanged=unchanged, sources_unchanged=sources_unchanged,
                  output_pins={p.relative_to(output).as_posix(): identity(p) for p in output.rglob('*') if p.is_file() and p.name != 'EXECUTION_PARENT.json'})
    write_json(output / 'receipts/EXECUTION_PARENT.json', record)
    print(json.dumps({'status': record['status'], 'output': str(output)}))
    if not success:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
