"""Prepare the original G205 design/TMM/membership/consumer without fitting models."""
import argparse
import datetime
import hashlib
import json
import os
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def identity(path):
    sha = hashlib.sha256()
    size = 0
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            size += len(block)
            sha.update(block)
    return {'bytes': size, 'sha256': sha.hexdigest()}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--resource-root', type=Path, required=True,
                        help='Root containing the exact rebuilt data/processed/gene_sets files')
    parser.add_argument('--rscript', required=True)
    parser.add_argument('--r-library', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output must be a new directory')
    manifest = json.loads((REPO / 'MANIFEST.json').read_text(encoding='utf-8'))
    contract_path = 'analysis/g205_preparation/CONTRACT.json'
    dependencies = ['pipeline/prepare_g205_model.py', contract_path,
                    'analysis/g205_preparation/prepare.R', 'analysis/study_sources/G205_SOURCE_ROLES.json']
    source_pins = {}
    for name in dependencies:
        actual = identity(REPO / name)
        if actual != manifest['files'][name]:
            raise ValueError('Release identity mismatch: ' + name)
        source_pins[name] = actual
    contract = json.loads((REPO / contract_path).read_text(encoding='utf-8'))
    for name, expected in contract['sources'].items():
        if identity(REPO / name) != expected:
            raise ValueError('Original source/config mismatch: ' + name)
        source_pins[name] = expected
    prepared = args.prepared.resolve(strict=True)
    if identity(prepared) != contract['prepared_identity']:
        raise ValueError('Prepared HISAT2 s0 object mismatch; regenerate with the documented input preparation entry')
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
        raise ValueError('R library is not a directory')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output / 'receipts').mkdir()
    (output / 'derived').mkdir()
    spec = dict(repo=REPO.as_posix(), output=output.as_posix(),
                prepared=prepared.as_posix(), prepared_identity=contract['prepared_identity'],
                resource_paths=resources, source_pins=source_pins,
                worker=contract['worker'], full_config=contract['full_config'],
                baseline_config=contract['baseline_config'])
    spec_path = output / 'receipts/RUN.json'
    write_json(spec_path, spec)
    command = [args.rscript, '--vanilla', (REPO / 'analysis/g205_preparation/prepare.R').as_posix(), spec_path.as_posix()]
    env = os.environ.copy()
    env.update(R_LIBS_USER=library.as_posix(), RENV_CONFIG_AUTOLOADER_ENABLED='FALSE',
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    locale = 'English_United States.utf8' if os.name == 'nt' else 'C.UTF-8'
    env.update(LANG=locale, LC_ALL=locale, LC_CTYPE=locale)
    start = datetime.datetime.now(datetime.timezone.utc).isoformat()
    with (output / 'receipts/stdout.txt').open('wb') as stdout, (output / 'receipts/stderr.txt').open('wb') as stderr:
        process = subprocess.run(command, env=env, cwd=output, stdout=stdout, stderr=stderr, check=False)
    unchanged = all(identity(Path(p)) == expected for p, expected in inputs.items())
    sources_unchanged = all(identity(REPO / p) == expected for p, expected in source_pins.items())
    qa_path = output / 'receipts/PREPARATION_QA.json'
    qa = json.loads(qa_path.read_text(encoding='utf-8')) if qa_path.exists() else {}
    success = process.returncode == 0 and unchanged and sources_unchanged and qa.get('status') == 'PASS_ORIGINAL_G205_PREPARATION_WITHOUT_MODEL_FITS'
    outputs = {p.relative_to(output).as_posix(): identity(p) for p in output.rglob('*') if p.is_file()}
    write_json(output / 'receipts/PARENT.json', dict(status='PASS' if success else 'FAIL',
               started_utc=start, finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
               command=command, returncode=process.returncode, inputs=inputs, source_pins=source_pins,
               inputs_unchanged=unchanged, sources_unchanged=sources_unchanged,
               outputs=outputs, new_real_data_models_run=0, historical_native_custody_recreated=False))
    print(json.dumps({'status': 'PASS' if success else 'FAIL', 'output': str(output)}))
    if not success:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
