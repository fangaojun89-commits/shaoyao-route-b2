"""Rebuild frozen resources with unchanged original code and authenticated local inputs."""
import argparse
import datetime
import hashlib
import json
import os
import platform
import shutil
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def identity(path):
    digest = hashlib.sha256()
    size = 0
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            size += len(block)
            digest.update(block)
    return {'bytes': size, 'sha256': digest.hexdigest()}


def relative(value):
    path = Path(value)
    if path.is_absolute() or '..' in path.parts or ':' in value or '\\' in value:
        raise ValueError('Expected a relative resource path')
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True,
                        help='Root containing the eight original data/raw/gene_sets files')
    parser.add_argument('--rscript', required=True)
    parser.add_argument('--r-library', default='')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output must be a new directory')
    manifest = json.loads((REPO / 'MANIFEST.json').read_text(encoding='utf-8'))
    contract_path = 'analysis/resource_build/CONTRACT.json'
    for name in (contract_path, 'analysis/resource_build/preflight.R',
                 'analysis/resource_build/check_consumers.R', 'analysis/resource_build/original_run_r.ps1',
                 'analysis/GSE16710/code/R/freeze_gene_sets.R', 'analysis/GSE16710/renv.lock',
                 'analysis/GSE16710/metadata/gene_set_manifest.tsv', 'pipeline/build_resources.py'):
        if identity(REPO / name) != manifest['files'][name]:
            raise ValueError('Release dependency identity mismatch: ' + name)
    contract = json.loads((REPO / contract_path).read_text(encoding='utf-8'))
    source_root = args.source_root.resolve(strict=True)
    for record in contract['inputs']:
        if identity(source_root / relative(record['path'])) != record['identity']:
            raise ValueError('Frozen input mismatch: ' + record['path'])
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    receipts = output / 'receipts'
    receipts.mkdir()
    for directory in ('code/R', 'metadata', 'results/diagnostics'):
        (output / directory).mkdir(parents=True, exist_ok=True)
    for record in contract['inputs']:
        target = output / relative(record['path'])
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source_root / record['path'], target)
        if identity(target) != record['identity']:
            raise ValueError('Staged input mismatch')
    builder = output / 'code/R/freeze_gene_sets.R'
    shutil.copyfile(REPO / 'analysis/GSE16710/code/R/freeze_gene_sets.R', builder)
    env = os.environ.copy()
    # Match the original launcher, which explicitly selected UTF-8 English
    # collation. Inheriting C changes GMT order and RDS string representation.
    locale = contract['locales']['Windows' if os.name == 'nt' else 'other']
    for key in ('LANG', 'LC_ALL', 'LC_CTYPE'):
        env[key] = locale
    env['RENV_CONFIG_AUTOLOADER_ENABLED'] = 'FALSE'
    if args.r_library:
        env['R_LIBS_USER'] = str(Path(args.r_library).resolve(strict=True))
    commands = [
        [args.rscript, '--vanilla', str(REPO / 'analysis/resource_build/preflight.R'), str(REPO), str(receipts)],
        [args.rscript, '--vanilla', str(builder)],
    ]
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    (receipts / 'COMMANDS.json').write_text(json.dumps(commands, indent=2) + '\n', encoding='utf-8')
    for index, command in enumerate(commands):
        with (receipts / f'{index}.stdout.txt').open('wb') as stdout, (receipts / f'{index}.stderr.txt').open('wb') as stderr:
            completed = subprocess.run(command, cwd=output, env=env, stdout=stdout, stderr=stderr, check=False)
        if completed.returncode:
            raise RuntimeError(f'Stage {index} failed; retained logs in {receipts}')
    observed = {}
    for record in contract['outputs']:
        actual = identity(output / relative(record['path']))
        if actual != record['identity']:
            raise ValueError('Reconstructed resource differs from accepted bytes: ' + record['path'])
        observed[record['path']] = actual
    consumer_command = [args.rscript, '--vanilla', str(REPO / 'analysis/resource_build/check_consumers.R'),
                        str(REPO), str(output)]
    commands.append(consumer_command)
    (receipts / 'COMMANDS.json').write_text(json.dumps(commands, indent=2) + '\n', encoding='utf-8')
    with (receipts / '2.stdout.txt').open('wb') as stdout, (receipts / '2.stderr.txt').open('wb') as stderr:
        consumer = subprocess.run(consumer_command, cwd=output, env=env, stdout=stdout, stderr=stderr, check=False)
    if consumer.returncode:
        raise RuntimeError('Original resource consumer checks failed; inspect retained stage 2 logs')
    consumer_qa = json.loads((receipts / 'RESOURCE_CONSUMER_QA.json').read_text(encoding='utf-8'))
    if consumer_qa['status'] != 'PASS_ORIGINAL_RESOURCE_CONSUMERS_WITHOUT_MODELS':
        raise ValueError('Unexpected consumer check result')
    for record in contract['inputs']:
        if identity(source_root / record['path']) != record['identity'] or identity(output / record['path']) != record['identity']:
            raise ValueError('Source or staged input changed during construction')
    receipt = {'status': 'PASS_ORIGINAL_RESOURCE_BUILD_EXACT_OUTPUTS',
               'started_utc': started, 'finished_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'python': platform.python_version(), 'system': platform.platform(),
               'requested_locale': locale,
               'contract': identity(REPO / contract_path), 'builder': identity(builder),
               'source_inputs_unchanged': True, 'input_files': len(contract['inputs']),
               'outputs': observed, 'scientific_model_fits': 0,
               'resource_consumer_validation': identity(receipts / 'RESOURCE_CONSUMER_QA.json'),
               'new_environment_installation': False, 'whole_study_restoration_verified': False}
    (receipts / 'RESOURCE_BUILD_COMPLETE.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
