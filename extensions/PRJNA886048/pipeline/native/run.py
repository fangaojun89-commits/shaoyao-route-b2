"""Plan or execute one fresh fixed-cohort library on an inspected Linux host.

No downloader, index construction, deletion, historical-stage resume or model fit.
Planning does not open raw reads or tools. Execution is deliberately explicit.
"""
import argparse
import csv
import datetime
import hashlib
import importlib.util
import json
import os
from pathlib import Path, PurePosixPath
import platform
import shutil
import subprocess
import tarfile


def need(value, message):
    if not value:
        raise ValueError(message)


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
GATE = load_module('native_execution_gate', HERE / 'validate_PRJNA886048_native_execution_20261028.py')


def digest(path, algorithm='sha256'):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, algorithm).hexdigest()


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def posix(value):
    need(isinstance(value, str) and value.startswith('/') and
         not any(c in value for c in ('\\', '\x00', '\r', '\n', ':')),
         'Expected an absolute Linux path')
    path = PurePosixPath(value)
    need('..' not in path.parts, 'Parent traversal is not allowed')
    return path


def overlap(a, b):
    return a == b or a in b.parents or b in a.parents


def build_plan(run, input_dir, fastp, kallisto, index, tx2gene, work_dir, archive_dir):
    ledger = ROOT / 'data/raw_ledger/ENA_quantification_manifest.tsv'
    need(digest(ledger) == GATE.LEDGER_SHA, 'Unknown source ledger')
    with ledger.open(encoding='utf-8-sig') as stream:
        rows = [r for r in csv.DictReader(stream, delimiter='\t') if r['run_accession'] == run]
    need(len(rows) == 1 and run in {f'SRR217758{i}' for i in range(46, 58)}, 'Unknown fixed-cohort library')
    paths = {k: posix(v) for k, v in dict(input_dir=input_dir, fastp=fastp, kallisto=kallisto,
              index=index, tx2gene=tx2gene, work_dir=work_dir, archive_dir=archive_dir).items()}
    out = paths['work_dir'] / ('linux_full_' + run + '_01')
    for protected in ('input_dir', 'fastp', 'kallisto', 'index', 'tx2gene'):
        need(not overlap(out, paths[protected]), 'Output overlaps a protected input')
        need(not overlap(paths['archive_dir'], paths[protected]), 'Archive directory overlaps a protected input')
    need(not overlap(out, paths['archive_dir']), 'Archive directory overlaps the work output')
    row = rows[0]
    sizes, hashes = row['fastq_bytes'].split(';'), row['fastq_md5'].split(';')
    need(len(sizes) == len(hashes) == 2, 'Source ledger is not paired')
    names = [f'{run}_{m}.fastq.gz' for m in (1, 2)]
    specs = {n: {'bytes': int(s), 'md5': h} for n, s, h in zip(names, sizes, hashes)}
    sources = [str(paths['input_dir'] / n) for n in names]
    prefixes = [str(out / f'orientation_R{m}.fastq.gz') for m in (1, 2)]
    cleaned = [str(out / f'R{m}.fastq.gz') for m in (1, 2)]
    commands = []
    for label, options in [('unstranded', []), ('rf', ['--rf-stranded']), ('fr', ['--fr-stranded'])]:
        commands.append([str(paths['kallisto']), 'quant', '-i', str(paths['index']), '-o',
                         str(out / ('orientation_' + label)), '-t', '4', '--plaintext'] + options + prefixes)
    commands.append([str(paths['fastp']), '-i', sources[0], '-I', sources[1], '-o', cleaned[0], '-O', cleaned[1]] +
                    GATE.FASTP_ARGS + ['--json', str(out / 'fastp.json'), '--html', str(out / 'fastp.html')])
    commands.append([str(paths['kallisto']), 'quant', '-i', str(paths['index']), '-o', str(out / 'quant'),
                     '-t', '8', '--plaintext'] + cleaned)
    total = sum(s['bytes'] for s in specs.values())
    return {'status': 'PLAN_ONLY_NO_NATIVE_EXECUTION', 'run': run,
            'paths': {k: str(v) for k, v in paths.items()}, 'output': str(out),
            'archive': str(paths['archive_dir'] / (run + '_results_portable.tar.gz')),
            'specs': specs, 'stage_names': GATE.STAGES, 'commands': commands,
            'output_envelope_bytes': int(1.2 * total), 'disk_reserve_bytes': 1024**3,
            'memory_reserve_bytes': 18 * 1024**3,
            'raw_files_read': False, 'native_execution': False}


def output_safety(plan):
    out, archive = Path(plan['output']), Path(plan['archive'])
    need(not out.exists() and not out.is_symlink(), 'Existing output requires inspection; no restart')
    need(not archive.exists() and not archive.is_symlink(), 'Existing archive requires inspection')
    need(not archive.with_suffix(archive.suffix + '.partial').exists(), 'Partial archive requires inspection')
    for folder in ('work_dir', 'archive_dir'):
        need(Path(plan['paths'][folder]).is_dir(), 'Explicit parent directory must already exist: ' + folder)


def file_pin(path, expected_sha):
    need(path.is_file() and digest(path) == expected_sha, 'Frozen file identity mismatch: ' + str(path))


def resource_admission(plan):
    v2 = Path('/sys/fs/cgroup/memory.max')
    v1 = Path('/sys/fs/cgroup/memory/memory.limit_in_bytes')
    if v2.is_file():
        limit_text = v2.read_text().strip()
        need(limit_text != 'max', 'A finite container memory limit is required')
        limit, current = int(limit_text), int(Path('/sys/fs/cgroup/memory.current').read_text())
    else:
        need(v1.is_file(), 'Container memory accounting is unavailable')
        limit = int(v1.read_text())
        current = int(Path('/sys/fs/cgroup/memory/memory.usage_in_bytes').read_text())
    need(0 <= current < limit < 2**60, 'A finite container memory limit is required')
    need(limit - current > plan['output_envelope_bytes'] + plan['memory_reserve_bytes'],
         'Insufficient container memory headroom')
    free = shutil.disk_usage(plan['paths']['work_dir']).free
    need(free > plan['output_envelope_bytes'] + plan['disk_reserve_bytes'], 'Insufficient output space')
    need(shutil.disk_usage(plan['paths']['archive_dir']).free > 64 * 1024**2, 'Insufficient archive space')
    cores = len(os.sched_getaffinity(0))
    cpu_v2 = Path('/sys/fs/cgroup/cpu.max')
    if cpu_v2.is_file():
        quota, period = cpu_v2.read_text().split()
        if quota != 'max':
            cores = min(cores, int(quota) / int(period))
    else:
        for controller in ('cpu', 'cpu,cpuacct'):
            base = Path('/sys/fs/cgroup') / controller
            if (base / 'cpu.cfs_quota_us').is_file():
                quota = int((base / 'cpu.cfs_quota_us').read_text())
                if quota > 0:
                    cores = min(cores, quota / int((base / 'cpu.cfs_period_us').read_text()))
                break
    need(cores >= 8, 'Frozen native commands require at least eight admitted CPU cores')
    return {'limit_bytes': limit, 'current_bytes': current, 'work_free_bytes': free, 'admitted_cpu_cores': cores}


def preflight(plan):
    need(platform.system() == 'Linux', 'Native execution/preflight requires Linux')
    output_safety(plan)
    for writable in (Path(plan['output']), Path(plan['paths']['archive_dir'])):
        need(not overlap(writable, ROOT.resolve()), 'Output overlaps the package')
    for key in ('input_dir', 'fastp', 'kallisto', 'index', 'tx2gene', 'work_dir', 'archive_dir'):
        path = Path(plan['paths'][key])
        need(str(path.resolve()) == str(path), 'Use resolved paths without symlink aliases: ' + key)
    for key in ('fastp', 'kallisto'):
        path = Path(plan['paths'][key])
        file_pin(path, GATE.TOOLS[key])
        need(os.access(path, os.X_OK), 'Tool is not executable: ' + key)
    file_pin(Path(plan['paths']['index']), GATE.INDEX_SHA)
    file_pin(Path(plan['paths']['tx2gene']), GATE._content.TX2GENE_SHA256)
    for name, spec in plan['specs'].items():
        path = Path(plan['paths']['input_dir']) / name
        need(path.is_file() and not path.is_symlink() and path.stat().st_size == spec['bytes'],
             'Raw mate size or file type mismatch: ' + name)
    return resource_admission(plan)


def bundle(out, target, receipt):
    members = ['receipt.json', 'fastp.json', 'fastp.html', 'quant/run_info.json', 'quant/abundance.tsv']
    partial = target.with_suffix(target.suffix + '.partial')
    need(not target.exists() and not partial.exists(), 'Archive exists; no overwrite')
    need(shutil.disk_usage(target.parent).free > sum((out / n).stat().st_size for n in members) + 8 * 1024**2,
         'Insufficient space for verified result archive')
    with partial.open('xb') as stream, tarfile.open(fileobj=stream, mode='w:gz') as tar:
        for name in members:
            tar.add(out / name, arcname=name, recursive=False)
    with tarfile.open(partial, 'r:gz') as tar:
        need(tar.getnames() == members, 'Wrong five-member bundle')
        for name in members:
            member = tar.getmember(name)
            need(member.isfile(), 'Nonregular result member')
            with tar.extractfile(member) as stream:
                need(hashlib.file_digest(stream, 'sha256').hexdigest() == digest(out / name),
                     'Archive readback differs: ' + name)
    # Same-directory hard-link publication is atomic and refuses an existing
    # destination on Linux; Path.rename would overwrite a racing destination.
    os.link(partial, target)
    partial.unlink()
    return {'path': str(target), 'bytes': target.stat().st_size, 'sha256': digest(target), 'members': members}


def execute_plan(plan):
    if not __debug__:
        raise ValueError('Python optimization disables required checks; run without -O')
    manifest_path = ROOT / 'MANIFEST.json'
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    for name, pin in manifest['files'].items():
        path = ROOT / name
        need(path.is_file() and path.stat().st_size == pin['bytes'] and digest(path) == pin['sha256'],
             'Package component mismatch: ' + name)
    admission = preflight(plan)
    lock = Path(plan['paths']['work_dir']) / '.PRJNA886048_native.lock'
    owner = json.dumps({'pid': os.getpid(), 'started_utc': now(), 'run': plan['run'], 'argv': __import__('sys').argv})
    with lock.open('x', encoding='utf-8') as stream:
        stream.write(owner)
    out = Path(plan['output'])
    receipt = {'status': 'VERIFYING_COMPLETE_INPUTS', 'pid': os.getpid(), 'library': plan['run'],
               'started_utc': now(), 'biological_models_run': 0, 'index_rebuilt': False,
               'manifest_sha256': digest(manifest_path), 'commands': [],
               'hostname': platform.node(), 'platform': platform.platform(),
               'runner_sha256': digest(Path(__file__)),
               'scientific_core_sha256': digest(HERE / 'scientific_core.py'),
               'source_manifest_sha256': GATE.LEDGER_SHA,
               'executable_sha256': GATE.TOOLS, 'index_sha256': GATE.INDEX_SHA,
               'cgroup_admission': admission, 'input_checks': {}}
    child = None
    created_output = False

    def save():
        receipt['updated_utc'] = now()
        temporary = out / 'receipt.tmp'
        temporary.write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
        temporary.replace(out / 'receipt.json')

    def execute(label, command):
        nonlocal child
        number = len(receipt['commands'])
        need(number < 5 and label == plan['stage_names'][number] and command == plan['commands'][number],
             'Frozen scientific command differs from the validated plan')
        record = {'stage': label, 'command': command, 'started_utc': now()}
        receipt['status'] = label
        receipt['commands'].append(record)
        save()
        with (out / (label + '.stdout')).open('w') as stdout, (out / (label + '.stderr')).open('w') as stderr:
            child = subprocess.Popen(command, stdout=stdout, stderr=stderr)
            record['pid'] = child.pid
            save()
            record['exit_code'] = child.wait()
        record['finished_utc'] = now()
        save()
        need(record['exit_code'] == 0, 'Native stage stopped: ' + label)

    try:
        output_safety(plan)
        out.mkdir(exist_ok=False)
        created_output = True
        save()
        names = list(plan['specs'])
        sources = {name: Path(plan['paths']['input_dir']) / name for name in names}
        for name in names:
            actual = digest(sources[name], 'md5')
            need(actual == plan['specs'][name]['md5'], 'Raw input MD5 mismatch: ' + name)
            receipt['input_checks'][name] = dict(plan['specs'][name], md5=actual)
        save()
        resource_admission(plan)
        core = load_module('fixed_scientific_stages', HERE / 'scientific_core.py')
        core.run_stages(run=plan['run'], names=names, sources=sources, out=out, receipt=receipt,
                        execute=execute, save=save, digest=digest,
                        fastp=Path(plan['paths']['fastp']), kallisto=Path(plan['paths']['kallisto']),
                        index=Path(plan['paths']['index']))
        receipt['finished_utc'] = now()
        save()
        checked = GATE.validate(out, Path(plan['paths']['tx2gene']),
                                ROOT / 'data/raw_ledger/ENA_quantification_manifest.tsv', plan['run'])
        (out / 'record_content_validation.json').write_text(json.dumps(checked, indent=2) + '\n', encoding='utf-8')
        archived = bundle(out, Path(plan['archive']), receipt)
        print(json.dumps({'status': 'NATIVE_COMPLETE_BUNDLE_READBACK_PASSED', 'archive': archived,
                          'local_R_acceptance': False, 'biological_models_run': 0}))
    except BaseException as error:
        if created_output:
            receipt.update(status='STOPPED_REQUIRES_INSPECTION', error_type=type(error).__name__, error=str(error))
            receipt['finished_utc'] = now()
            save()
        raise
    finally:
        # Interrupted observations do not prove that the original child stopped.
        # Leave ownership in place while that exact child is still live.
        if ((child is None or child.poll() is not None) and lock.exists()
                and lock.read_text(encoding='utf-8') == owner):
            lock.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', required=True)
    for key in ('input-dir', 'fastp', 'kallisto', 'index', 'tx2gene', 'work-dir', 'archive-dir'):
        parser.add_argument('--' + key, required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--plan', action='store_true')
    mode.add_argument('--preflight', action='store_true')
    mode.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    plan = build_plan(**{k: getattr(args, k) for k in ('run', 'input_dir', 'fastp', 'kallisto',
                                                      'index', 'tx2gene', 'work_dir', 'archive_dir')})
    if args.plan:
        print(json.dumps(plan, indent=2))
    elif args.preflight:
        print(json.dumps({'status': 'PREFLIGHT_ONLY_NO_NATIVE_EXECUTION', 'resources': preflight(plan),
                          'raw_files_hashed': False}))
    else:
        execute_plan(plan)


if __name__ == '__main__':
    main()
