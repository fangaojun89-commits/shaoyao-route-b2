"""JGW03 original primary sequence, with a separate no-fit source-input check."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

REPO = Path(__file__).resolve().parents[1]
BASE = 'analysis/jgw_execution/'


def identity(path):
    path = Path(path)
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return {'bytes': path.stat().st_size, 'sha256': h.hexdigest()}


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def write(path, obj):
    Path(path).write_bytes((json.dumps(obj, indent=2, ensure_ascii=False)+'\n').encode('utf-8'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['check', 'run'])
    parser.add_argument('--input-root', type=Path, required=True,
                        help='Root containing the three pinned data/raw/<accession>/ GEO files')
    parser.add_argument('--resource-root', type=Path, required=True,
                        help='Completed build_resources.py output root')
    parser.add_argument('--rscript', required=True)
    parser.add_argument('--r-library', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error('Use a new output directory; completed results are never overwritten')
    library = args.r_library.resolve(strict=True)
    if not library.is_dir():
        parser.error('R library must be a directory')
    contract = read(REPO/BASE/'CONTRACT.json')
    manifest = read(REPO/'MANIFEST.json')
    dependencies = ['pipeline/jgw_model.py', BASE+'check_inputs.R', BASE+'CONTRACT.json']
    dependencies += [r['public_path'] for r in contract['public_inputs']]
    dependencies += [r['public_path'] for r in contract['expected_inputs']]
    source_pins = {}
    for name in dependencies:
        actual = identity(REPO/name)
        if actual != manifest['files'][name]:
            raise ValueError('Release file changed: '+name)
        source_pins[name] = actual
    staging = []
    for r in contract['public_inputs']:
        staging.append((REPO/r['public_path'], r['runtime_path'], r['identity']))
    for r in contract['source_inputs']:
        staging.append((args.input_root/r['runtime_path'], r['runtime_path'], r['identity']))
    for r in contract['resources']:
        staging.append((args.resource_root/r['runtime_path'], r['runtime_path'], r['identity']))
    # Validate all paths before writing or invoking R. Keep parent inputs immutable.
    input_pins = {}
    for source, relative, pin in staging:
        source = source.resolve(strict=True)
        if Path(relative).is_absolute() or '..' in Path(relative).parts:
            raise ValueError('Invalid relative staging path')
        if identity(source) != pin:
            raise ValueError('Pinned source input changed: '+str(source))
        input_pins[str(source)] = pin
    output.mkdir(parents=True, exist_ok=False)
    receipts = output/'receipts'; receipts.mkdir()
    project = output/'project'; project.mkdir()
    for source, relative, pin in staging:
        dest = project/relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, dest)
        if identity(dest) != pin:
            raise ValueError('Staged file changed: '+relative)
    spec = {'repo': REPO.as_posix(), 'project': project.as_posix(),
            'output': output.as_posix(), 'contract': contract,
            'mode': args.mode, 'source_pins': source_pins}
    write(receipts/'RUN.json', spec)
    env = os.environ.copy()
    locale = 'English_United States.utf8' if os.name == 'nt' else 'en_US.UTF-8'
    env.update(R_LIBS='', R_LIBS_USER=library.as_posix(), RENV_CONFIG_AUTOLOADER_ENABLED='FALSE',
               ROUTE_B2_JGW03_RESUME='false', LANG=locale, LC_ALL=locale,
               LC_CTYPE=locale, TZ='Asia/Shanghai', OMP_NUM_THREADS='1',
               OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    commands = [[args.rscript, '--vanilla', str(REPO/BASE/'check_inputs.R'), str(receipts/'RUN.json')]]
    # run is an explicit future reconstruction, never invoked by check.
    if args.mode == 'run':
        commands.extend([[args.rscript, '--vanilla', str(project/source), str(project)]
                         for source in contract['sequence']])
    write(receipts/'COMMANDS.json', commands)
    start = datetime.datetime.now(datetime.timezone.utc).isoformat()
    exits = []
    for i, command in enumerate(commands):
        with (receipts/f'{i:02d}_stdout.txt').open('wb') as stdout, (receipts/f'{i:02d}_stderr.txt').open('wb') as stderr:
            result = subprocess.run(command, cwd=output, env=env, stdout=stdout, stderr=stderr, check=False)
        exits.append(result.returncode)
        if result.returncode:
            break
    unchanged = all(identity(Path(path)) == pin for path, pin in input_pins.items())
    sources_same = all(identity(REPO/path) == pin for path, pin in source_pins.items())
    staged_same = all(identity(project/relative) == pin for _, relative, pin in staging)
    success = len(exits) == len(commands) and not any(exits) and unchanged and sources_same and staged_same
    qa = read(receipts/'INPUT_CHECK.json') if (receipts/'INPUT_CHECK.json').exists() else None
    if not qa or qa.get('status') != 'PASS_ORIGINAL_SOURCE_READERS_NO_FITS':
        success = False
    products = {p.relative_to(output).as_posix(): identity(p)
                for p in sorted(output.rglob('*')) if p.is_file()}
    report = {'status': 'PASS' if success else 'FAIL', 'mode': args.mode,
              'started_utc': start, 'finished_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'exit_codes': exits, 'input_pins': input_pins, 'source_pins': source_pins,
              'external_inputs_unchanged': unchanged, 'release_sources_unchanged': sources_same,
              'staged_sources_inputs_resources_unchanged': staged_same,
              'input_check': qa, 'outputs': products,
              'primary_analysis_invoked': args.mode == 'run' and len(exits) > 1,
              'whole_study_restoration_verified': False,
              'scope': 'Original primary and dose-summary sequence; separate rotations, conditional references and later reporting revisions are not executed'}
    write(receipts/'COMPLETE.json', report)
    print(json.dumps({'status': report['status'], 'mode': args.mode, 'output': str(output)}))
    if not success:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
