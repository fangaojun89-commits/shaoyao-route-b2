"""Check inputs, test synthetic data, or explicitly run the original G128 sequence."""
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
    sha = hashlib.sha256()
    size = 0
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            sha.update(block)
            size += len(block)
    return {'bytes': size, 'sha256': sha.hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('check', 'synthetic', 'run'))
    parser.add_argument('--counts-root', type=Path)
    parser.add_argument('--resources', type=Path)
    parser.add_argument('--cohort', choices=('raw', 'adapter'))
    parser.add_argument('--rscript', required=True)
    parser.add_argument('--r-library', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output must be a new directory; existing results are never overwritten')
    if args.mode == 'run' and platform.system() != 'Linux':
        parser.error('formal model execution requires Linux; check and synthetic accept Windows')
    if args.mode == 'synthetic':
        if any(x is not None for x in (args.counts_root, args.resources, args.cohort)):
            parser.error('synthetic mode accepts no real scientific inputs')
    elif args.counts_root is None or args.resources is None or args.cohort is None:
        parser.error('check/run require --counts-root, --resources and explicit --cohort')
    manifest = json.loads((REPO / 'MANIFEST.json').read_text(encoding='utf-8'))
    contract_name = 'analysis/g128_execution/CONTRACT.json'
    dependencies = ['pipeline/g128_model.py', contract_name,
                    'analysis/g128_execution/entry.R', 'analysis/g128_execution/synthetic.R',
                    'analysis/g128_execution/compressed_io.R',
                    'analysis/GSE16710/expected/S1_GATE.tsv', 'data/count_inputs/PROVENANCE.json']
    for name in dependencies:
        if identity(REPO / name) != manifest['files'][name]:
            raise ValueError('Release identity mismatch: ' + name)
    contract = json.loads((REPO / contract_name).read_text(encoding='utf-8'))
    for record in contract['sources']:
        if identity(REPO / record['path']) != record['identity']:
            raise ValueError('Original source mismatch: ' + record['role'])
    source_pins = {name: identity(REPO / name) for name in dependencies}
    source_pins.update({r['path']: r['identity'] for r in contract['sources']})
    inputs = {}
    paths = {}
    if args.mode != 'synthetic':
        provenance = json.loads((REPO / 'data/count_inputs/PROVENANCE.json').read_text(encoding='utf-8'))
        roles = {'count_matrix': args.cohort + '/primary_all30/counts',
                 'feature_annotation': 'annotation', 'sample_ledger': 'ledger'}
        for role, source_role in roles.items():
            record, = [r for r in provenance['files'] if r['dataset'] == 'GSE128623' and r['role'] == source_role]
            path = (args.counts_root / record['runtime_relative_path']).resolve(strict=True)
            actual = identity(path)
            if actual != {'bytes': record['original_bytes'], 'sha256': record['original_sha256']}:
                raise ValueError('Accepted input mismatch: ' + role)
            paths[role] = path.as_posix()
            inputs[str(path)] = actual
        resource = args.resources.resolve(strict=True)
        if identity(resource) != contract['frozen_resource']:
            raise ValueError('Frozen resource identity mismatch')
        paths['gene_set_resource'] = resource.as_posix()
        inputs[str(resource)] = identity(resource)
    output = args.output.resolve()
    library = args.r_library.resolve(strict=True)
    if not library.is_dir():
        raise ValueError('R library must be a directory')
    output.mkdir(parents=True, exist_ok=False)
    receipts = output / 'receipts'
    receipts.mkdir()
    if args.mode != 'synthetic':
        frozen = output / 'frozen_design.tsv'
        shutil.copyfile(paths['sample_ledger'], frozen)
        paths['frozen_design'] = frozen.as_posix()
        paths['s1_zero'] = (REPO / 'analysis/GSE16710/expected/S1_GATE.tsv').as_posix()
        inputs[str(frozen)] = identity(frozen)
        inputs[paths['s1_zero']] = identity(Path(paths['s1_zero']))
    paths['staging'] = (output / 'model_output').as_posix()
    paths['log_directory'] = receipts.as_posix()
    spec = {'mode': args.mode, 'cohort': args.cohort, 'paths': paths,
            'repo': REPO.as_posix(), 'output': output.as_posix()}
    spec_path = receipts / 'RUN.json'
    spec_path.write_text(json.dumps(spec, indent=2) + '\n', encoding='utf-8')
    env = os.environ.copy()
    for key in ('LANG', 'LC_ALL', 'LC_CTYPE'):
        env[key] = 'English_United States.utf8' if os.name == 'nt' else 'en_US.UTF-8'
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'R_DATATABLE_NUM_THREADS'):
        env[key] = '1'
    env['R_LIBS_USER'] = str(library)
    env['RENV_CONFIG_AUTOLOADER_ENABLED'] = 'FALSE'
    command = [args.rscript, '--vanilla', str(REPO / 'analysis/g128_execution/entry.R'), str(spec_path)]
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    (receipts / 'COMMAND.json').write_text(json.dumps(command, indent=2) + '\n', encoding='utf-8')
    with (receipts / 'stdout.txt').open('wb') as stdout, (receipts / 'stderr.txt').open('wb') as stderr:
        completed = subprocess.run(command, cwd=output, env=env, stdout=stdout, stderr=stderr, check=False)
    unchanged = all(identity(Path(path)) == pin for path, pin in inputs.items())
    sources_unchanged = all(identity(REPO / path) == pin for path, pin in source_pins.items())
    result = {'exit_code': completed.returncode, 'mode': args.mode, 'cohort': args.cohort,
              'started_utc': started, 'finished_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'python': platform.python_version(), 'system': platform.platform(),
              'inputs': inputs, 'inputs_unchanged': unchanged,
              'sources_unchanged': sources_unchanged,
              'release_manifest': identity(REPO / 'MANIFEST.json'), 'contract': identity(REPO / contract_name),
              'whole_study_restoration_verified': False}
    if completed.returncode == 0:
        result['validation'] = json.loads((receipts / 'R_VALIDATION.json').read_text(encoding='utf-8'))
        result['outputs'] = {p.relative_to(output).as_posix(): identity(p)
                             for p in sorted(output.rglob('*')) if p.is_file()}
    (receipts / 'COMPLETE.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    if completed.returncode != 0 or not unchanged or not sources_unchanged:
        raise RuntimeError('Execution or input-preservation check failed; retained logs: ' + str(receipts))
    print(json.dumps({'status': result['validation']['status'], 'mode': args.mode,
                      'inputs_unchanged': unchanged, 'output': str(output)}))


if __name__ == '__main__':
    main()
