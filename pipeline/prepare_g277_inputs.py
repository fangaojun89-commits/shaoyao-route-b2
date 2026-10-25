"""Prepare unchanged GSE277064 count inputs and original prefit design definitions."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

REPO = Path(__file__).resolve().parents[1]


def identity(path):
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'bytes': path.stat().st_size, 'sha256': digest}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rscript', required=True)
    parser.add_argument('--r-library', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Use a new output directory')
    library = args.r_library.resolve(strict=True)
    if not library.is_dir():
        parser.error('R library must be a directory')
    contract_name = 'data/g277_inputs/PROVENANCE.json'
    source_name = 'analysis/study_sources/code/R/analyze_GSE277064.R'
    names = ['pipeline/prepare_g277_inputs.py', 'analysis/g277_preparation/prepare.R',
             source_name, contract_name]
    manifest = json.loads((REPO / 'MANIFEST.json').read_text(encoding='utf-8'))
    pins = {name: identity(REPO / name) for name in names}
    for name, actual in pins.items():
        if actual != manifest['files'][name]:
            raise ValueError('Release source identity mismatch: ' + name)
    contract = json.loads((REPO / contract_name).read_text(encoding='utf-8'))
    for item in contract['files']:
        if identity(REPO / item['public_path']) != item['identity']:
            raise ValueError('GSE277064 input identity mismatch: ' + item['public_path'])
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    for name in ('prepared', 'receipts', 'project'):
        (output / name).mkdir()
    for item in contract['files']:
        path = output / 'project' / item['runtime_path']
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / item['public_path'], path)
        if identity(path) != item['identity']:
            raise ValueError('Staged input mismatch: ' + item['runtime_path'])
    spec = dict(repo=str(REPO), project=str(output / 'project'), output=str(output),
                source=source_name, source_pins=pins)
    spec_path = output / 'receipts/INVOCATION.json'
    spec_path.write_text(json.dumps(spec, indent=2) + '\n', encoding='utf-8')
    env = os.environ.copy()
    env['R_LIBS_USER'] = str(library)
    env['R_LIBS'] = str(library)
    env['LANG'] = 'en_US.UTF-8'
    env['LC_ALL'] = 'en_US.UTF-8'
    cmd = [args.rscript, '--vanilla', str(REPO / 'analysis/g277_preparation/prepare.R'), str(spec_path)]
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    result = subprocess.run(cmd, cwd=output, env=env, capture_output=True)
    (output / 'receipts/stdout.txt').write_bytes(result.stdout)
    (output / 'receipts/stderr.txt').write_bytes(result.stderr)
    receipt = dict(command=cmd, started_utc=started, native_exit=result.returncode,
                   source_pins=pins, scientific_model_fits=0,
                   normalization_or_scoring=False)
    if result.returncode == 0:
        receipt['outputs'] = {str(p.relative_to(output)): identity(p)
                              for p in (output / 'prepared').glob('*') if p.is_file()}
    (output / 'receipts/EXECUTION.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    if result.returncode:
        raise RuntimeError('Original input/design preparation failed; retained receipts contain diagnostics')
    print(json.dumps(dict(status='PASS_G277_PREFIT_INPUT_PREPARATION',
                          outputs=len(receipt['outputs']), scientific_model_fits=0)))


if __name__ == '__main__':
    main()
