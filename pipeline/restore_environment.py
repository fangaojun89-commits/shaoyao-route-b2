"""Restore the catalogued G205 R core into a fresh Linux library; no scientific fits.

Defaults to source/runtime validation and a printed plan. Installation requires
--execute and never modifies an existing output directory or installs the R runtime.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('acquisition', REPO / 'pipeline/acquisition.py')
acquisition = importlib.util.module_from_spec(spec)
spec.loader.exec_module(acquisition)


def isolated_environment(r, output):
    env = os.environ.copy()
    env.pop('R_HOME', None)
    env.update({
        'R_LIBS': '', 'R_LIBS_SITE': 'NULL', 'R_LIBS_USER': str(output / 'library'),
        'R_PROFILE': '/dev/null', 'R_PROFILE_USER': '/dev/null',
        'R_ENVIRON': '/dev/null', 'R_ENVIRON_USER': '/dev/null',
        'R_MAKEVARS_USER': '/dev/null', 'R_MAKEVARS_SITE': '/dev/null',
        'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
        'MAKEFLAGS': '-j1', 'LANG': 'C.UTF-8', 'LC_ALL': 'C.UTF-8', 'TZ': 'UTC',
        'RENV_CONFIG_AUTOLOADER_ENABLED': 'FALSE',
        'PATH': str(r.parent) + ':/usr/bin:/bin'})
    return env


def make_plan(sources, r, output):
    if sys.platform != 'linux':
        raise RuntimeError('Use a separately provisioned Linux R/toolchain')
    sources = sources.resolve(strict=True)
    r = r.absolute()
    output = output.absolute()
    if output.exists() or output.is_symlink():
        raise ValueError('Output must be a fresh directory; existing files are preserved')
    if not output.parent.is_dir() or output.resolve().is_relative_to(sources):
        raise ValueError('Output parent must exist outside the source archive directory')
    if not r.is_file() or not (r.parent / 'Rscript').is_file():
        raise ValueError('Explicit existing R and adjacent Rscript are required')
    _, catalog = acquisition.catalog()
    checked = acquisition.check_files(sources, catalog['packages'], package=True)
    if checked['status'] != 'ALL_FILES_MATCH':
        raise ValueError('Exact source archive or DESCRIPTION check failed')
    env = isolated_environment(r, output)
    version = subprocess.run([str(r), '--vanilla', '--slave', '-e',
                              'cat(as.character(getRversion()))'], env=env,
                             capture_output=True, text=True, timeout=30, check=True)
    if version.stdout.strip() != '4.6.1':
        raise ValueError('R 4.6.1 is required by this catalog')
    return {
        'status': 'SOURCE_AND_RUNTIME_CHECKED_PLAN', 'r': str(r),
        'r_sha256': acquisition.digest(r), 'r_version': '4.6.1',
        'output': str(output), 'source_root': str(sources),
        'catalog_sha256': acquisition.digest(REPO / 'environment/G205_core/SOURCES.json'),
        'packages': [row['package'] for row in catalog['packages']],
        'commands': [[str(r), 'CMD', 'INSTALL', '--library=' + str(output / 'library'),
                      str(sources / row['filename'])] for row in catalog['packages']],
        'acceptance_command': [str(r.parent / 'Rscript'), '--vanilla',
                               str(REPO / 'environment/G205_core/verify_library.R'),
                               str(output / 'library'),
                               str(REPO / 'environment/G205_core/SOURCES.json'),
                               str(output / 'acceptance')],
        'free_bytes': shutil.disk_usage(output.parent).free,
        'scientific_model_runs': 0, 'whole_study_cold_restore': False}


def install(plan, run=subprocess.run):
    output = Path(plan['output'])
    # mkdir is exclusive even if another caller appeared after preflight.
    output.mkdir(exist_ok=False)
    for name in ('library', 'logs', 'tmp'):
        (output / name).mkdir()
    env = isolated_environment(Path(plan['r']), output)
    env.update({key: str(output / 'tmp') for key in ('TMPDIR', 'TMP', 'TEMP')})
    receipt = dict(plan, status='INSTALLING', completed_packages=[], stages=[])

    def save():
        temporary = output / 'receipt.tmp'
        temporary.write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
        temporary.replace(output / 'receipt.json')

    def execute(name, command):
        stage = {'name': name, 'argv': command, 'started_epoch': time.time()}
        receipt['stages'].append(stage)
        save()
        with (output / 'logs' / (name + '.stdout')).open('xb') as stdout, \
                (output / 'logs' / (name + '.stderr')).open('xb') as stderr:
            result = run(command, env=env, cwd=output, stdout=stdout, stderr=stderr)
        stage.update(exit_code=result.returncode, finished_epoch=time.time())
        save()
        if result.returncode:
            raise RuntimeError('Stage failed; files retained: ' + name)

    try:
        for name, command in zip(plan['packages'], plan['commands'], strict=True):
            execute(name, command)
            receipt['completed_packages'].append(name)
        receipt['status'] = 'INSTALL_PASS_LOAD_CHECK_PENDING'
        save()
        execute('verify_library', plan['acceptance_command'])
        accepted = json.loads((output / 'acceptance/RESULT.json').read_text(encoding='utf-8'))
        if accepted.get('status') != '70_INSTALLED_VERSIONS_AND_EIGHT_ROOT_LOADS_PASS':
            raise RuntimeError('Installed-version and namespace load check did not pass')
        receipt['status'] = 'RESTORED_CORE_LOAD_CHECK_PASS'
        save()
    except BaseException as error:
        receipt.update(status='INCOMPLETE_FILES_RETAINED', error=str(error))
        save()
        raise
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources', type=Path, required=True)
    parser.add_argument('--r', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    plan = make_plan(args.sources, args.r, args.output)
    if args.execute:
        if plan['free_bytes'] < 10 * 1024**3:
            raise RuntimeError('Keep at least 10 GiB free before starting the core build')
        plan = install(plan)
    print(json.dumps(plan, indent=2))


if __name__ == '__main__':
    main()
