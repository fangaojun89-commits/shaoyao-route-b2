"""Inspect pinned acquisition records or check user-supplied files; never download or install."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import tarfile


REPO = Path(__file__).resolve().parents[1]


def digest(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def catalog():
    inputs = json.loads((REPO / 'data/acquisition/INPUTS.json').read_text(encoding='utf-8'))
    packages = json.loads((REPO / 'environment/G205_core/SOURCES.json').read_text(encoding='utf-8'))
    accepted = json.loads((REPO / 'environment/G205_core/ACCEPTED_ENVIRONMENT.json').read_text(encoding='utf-8'))
    if len(inputs['records']) != 71 or len({r['path'] for r in inputs['records']}) != 71:
        raise ValueError('Input selection changed')
    rows = packages['packages']
    if len(rows) != 70 or len({r['package'] for r in rows}) != 70:
        raise ValueError('Package selection changed')
    if [r['build_position'] for r in rows] != list(range(1, 71)):
        raise ValueError('Package build order changed')
    versions = {r['package']: r['version'] for r in accepted['installed_versions']}
    if versions != {r['package']: r['version'] for r in rows}:
        raise ValueError('Package catalog differs from accepted installation')
    for row in inputs['records']:
        p = Path(row['path'])
        if p.is_absolute() or '..' in p.parts:
            raise ValueError('Invalid relative input path')
    for row in rows:
        if row['filename'] != row['package'] + '_' + row['version'] + '.tar.gz':
            raise ValueError('Invalid source filename')
        if Path(row['filename']).name != row['filename']:
            raise ValueError('Non-flat source filename')
    return inputs, packages


def check_files(root, rows, package=False):
    root = root.resolve(strict=True)
    checked = []
    for row in rows:
        relative = row['filename'] if package else row['path']
        path = root / relative
        result = {'path': relative, 'status': 'MISSING'}
        if path.is_file():
            # Symlinks are permitted: scientific storage may be mounted elsewhere.
            if path.stat().st_size != row['bytes']:
                result['status'] = 'SIZE_MISMATCH'
            else:
                algorithms = [a for a in ('sha256', 'md5') if a in row]
                if not algorithms:
                    raise ValueError('No recorded checksum')
                result['status'] = 'MATCH' if all(digest(path, a) == row[a] for a in algorithms) else 'CHECKSUM_MISMATCH'
                if package and result['status'] == 'MATCH':
                    with tarfile.open(path, 'r:gz') as archive:
                        members = [m for m in archive.getmembers() if len(m.name.strip('/').split('/')) == 2 and m.name.endswith('/DESCRIPTION')]
                        if len(members) != 1 or not members[0].isfile():
                            raise ValueError('Ambiguous package DESCRIPTION: ' + relative)
                        raw = archive.extractfile(members[0]).read()
                    if hashlib.sha256(raw).hexdigest() != row['description_sha256']:
                        result['status'] = 'DESCRIPTION_MISMATCH'
        checked.append(result)
    return {'status': 'ALL_FILES_MATCH' if all(r['status'] == 'MATCH' for r in checked) else 'INCOMPLETE_OR_MISMATCH',
            'checked': len(checked), 'matched': sum(r['status'] == 'MATCH' for r in checked),
            'files': checked, 'downloads': 0, 'installations': 0, 'scientific_model_runs': 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['inspect', 'check-inputs', 'check-packages', 'install-plan'])
    parser.add_argument('--root', type=Path, help='Input catalog root, or flat directory of 70 package archives')
    parser.add_argument('--r', type=Path, help='Explicit separately provisioned Linux R executable, for plan only')
    parser.add_argument('--library', type=Path, help='Fresh isolated Linux library path, for plan only')
    args = parser.parse_args()
    inputs, packages = catalog()
    if args.command == 'inspect':
        result = {'status': 'CATALOG_CONSISTENCY_PASS', 'input_records': 71,
                  'FASTQ_files': sum(r['use'] == 'raw_counting_input' for r in inputs['records']),
                  'reference_files_without_recovered_URL': sum(r['url'] is None for r in inputs['records']),
                  'G205_core_packages': 70, 'local_input_bytes_checked': False,
                  'whole_study_restore_proven': False}
    else:
        if args.root is None:
            parser.error('--root required')
        is_package = args.command in ('check-packages', 'install-plan')
        result = check_files(args.root, packages['packages'] if is_package else inputs['records'], is_package)
        if args.command == 'install-plan' and result['status'] == 'ALL_FILES_MATCH':
            if args.r is None or args.library is None:
                parser.error('--r and --library required')
            if not args.r.is_absolute() or not args.library.is_absolute():
                parser.error('Plan paths must be absolute')
            if args.library.exists() or args.library.is_symlink():
                parser.error('Library must be a fresh path')
            result['status'] = 'SOURCE_CHECKED_INSTALL_PLAN_ONLY'
            result['required_platform'] = 'Linux with independently provisioned R 4.6.1 and compatible compiler/native libraries'
            result['runtime_checked'] = False
            result['commands'] = [[str(args.r), 'CMD', 'INSTALL', '--library=' + str(args.library),
                                   str(args.root.resolve() / r['filename'])] for r in packages['packages']]
            result['environment_overrides'] = {
                'R_LIBS': '', 'R_LIBS_SITE': 'NULL', 'R_LIBS_USER': str(args.library),
                'R_PROFILE': '/dev/null', 'R_PROFILE_USER': '/dev/null',
                'R_ENVIRON': '/dev/null', 'R_ENVIRON_USER': '/dev/null',
                'R_MAKEVARS_USER': '/dev/null', 'R_MAKEVARS_SITE': '/dev/null',
                'MAKEFLAGS': '-j1', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1',
                'MKL_NUM_THREADS': '1', 'LANG': 'C.UTF-8', 'LC_ALL': 'C.UTF-8', 'TZ': 'UTC',
                'RENV_CONFIG_AUTOLOADER_ENABLED': 'FALSE'}
    print(json.dumps(result, indent=2))
    return int(result['status'] == 'INCOMPLETE_OR_MISMATCH')


if __name__ == '__main__':
    sys.exit(main())
