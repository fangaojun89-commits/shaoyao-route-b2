"""Authenticate fixed inputs, then run unchanged modules or reuse verified checkpoints.

The package root is located from this script, never from the working directory.
--check-only reads inputs and installed package metadata; it creates no run directory.
"""
from pathlib import Path
import argparse
import csv
import datetime
import hashlib
import importlib.metadata
import importlib.util
import json
import os
import shutil
import subprocess
import sys


def sha(p):
    with Path(p).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def read(p):
    return json.loads(Path(p).read_text(encoding='utf-8-sig'))


def dump(p, value):
    p = Path(p)
    temp = p.with_name(p.name + '.partial')
    if temp.exists():
        raise ValueError('Unresolved partial metadata: ' + str(temp))
    temp.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    temp.replace(p)


def pin(p, expected):
    if not p.is_file() or p.stat().st_size != expected['bytes'] or sha(p) != expected['sha256']:
        raise ValueError('Input identity mismatch: ' + str(p))


def inside(a, b):
    return a == b or a.is_relative_to(b)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', type=Path, required=True)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--check-only', action='store_true')
    mode.add_argument('--run', action='store_true')
    mode.add_argument('--resume', action='store_true')
    ap.add_argument('--output', type=Path)
    ap.add_argument('--rscript', default='Rscript')
    ap.add_argument('--r-library', action='append', default=[], help='Explicit R library; repeat in precedence order if needed')
    ap.add_argument('--adopt-completed-checkpoints', type=Path)
    ap.add_argument('--native-record-checks', action='store_true', help='Bind actual native files after canonical/historical record checks; forbids historical checkpoint adoption')
    a = ap.parse_args()
    if a.native_record_checks and a.adopt_completed_checkpoints:
        raise ValueError('Native-record mode forbids historical checkpoint adoption')
    root = Path(__file__).resolve().parents[1]
    manifest = read(root / 'MANIFEST.json')
    for name, expected in manifest['files'].items():
        p = (root / name).resolve()
        if not inside(p, root):
            raise ValueError('Manifest path escapes the package')
        pin(p, expected)
    config_path = a.config.resolve(strict=True)
    config = read(config_path)
    roles = {'native_root', 'tx2gene', 'mapped_edges', 'reference_terms'}
    if set(config) != roles or not all(isinstance(v, str) and v for v in config.values()):
        raise ValueError('Config accepts exactly four input paths; statistical overrides are not supported')
    inputs = {k: (config_path.parent / v).resolve(strict=True) for k, v in config.items()}
    input_contract = read(root / 'provenance/EXTERNAL_INPUTS.json')
    native_identity = None
    if a.native_record_checks:
        spec = importlib.util.spec_from_file_location('cohort_identity', root / 'pipeline/native/build_PRJNA886048_native_input_identity_20261028.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        native_identity = module.build(root, inputs)
    active_files = native_identity['payload']['files'] if native_identity else input_contract['files']
    for item in active_files:
        path = inputs[item['role']]
        if item.get('relative'):
            path = path / item['relative']
        pin(path, item)
    py_contract = read(root / 'environment/Python_dependency_contract.json')
    if sys.version_info[:2] != tuple(py_contract['python_major_minor']):
        raise ValueError('Python version differs from the locked major/minor runtime')
    for pkg, version in py_contract['packages'].items():
        if importlib.metadata.version(pkg) != version:
            raise ValueError('Python dependency version mismatch: ' + pkg)
    env_cmd = [a.rscript, '--vanilla', str(root / 'pipeline/engine.R'), str(root),
               '--environment-only', '-', os.pathsep.join(a.r_library)]
    env = dict(os.environ, LC_ALL='C', LANG='C', TZ='UTC', PYTHONDONTWRITEBYTECODE='1',
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    er = subprocess.run(env_cmd, capture_output=True, text=True, env=env)
    if er.returncode:
        sys.stderr.write(er.stdout + er.stderr)
        raise ValueError('R dependency closure check failed')
    adopt = {}
    if a.adopt_completed_checkpoints:
        source = a.adopt_completed_checkpoints.resolve(strict=True)
        old = read(root / 'provenance/ACCEPTED_CHECKPOINTS.json')
        if old['external_input_manifest_sha256'] != sha(root / 'provenance/EXTERNAL_INPUTS.json'):
            raise ValueError('Checkpoint contract is not bound to the authenticated input set')
        for key, item in old['files'].items():
            p = source / item['relative']
            pin(p, item)
            adopt[key] = {'path': str(p), 'sha256': item['sha256']}
    binding_names = ['provenance/EXTERNAL_INPUTS.json', 'environment/R_dependency_contract.json',
                     'environment/Python_dependency_contract.json', 'environment/renv.lock', 'uv.lock',
                     'pipeline/run.py', 'pipeline/engine.R', 'pipeline/verify_tables.R',
                     'pipeline/compare_outputs.R', 'data/derived/samples.tsv',
                     'data/derived/retained_PROGENy_comparator_basis.tsv',
                     'analysis/original/compare_PRJNA886048_retained_programs_20261025.py']
    binding_names += ['analysis/original/analyze_PRJNA886048_' + s + '_20261025.R'
                      for s in ('gene_contrast', 'programs', 'signed_activity')]
    bindings = {n: sha(root / n) for n in binding_names}
    if native_identity:
        bindings['validated_actual_native_input_identity'] = native_identity['binding']
        bindings['pipeline/native/build_PRJNA886048_native_input_identity_20261028.py'] = sha(root / 'pipeline/native/build_PRJNA886048_native_input_identity_20261028.py')
    binding = hashlib.sha256(json.dumps(bindings, sort_keys=True).encode()).hexdigest()
    if a.check_only:
        if a.output:
            raise ValueError('--check-only does not accept an output directory')
        print(json.dumps({'status': 'PASS_NATIVE_RECORD_INPUTS_AND_INSTALLED_ENVIRONMENT' if native_identity else 'PASS_PINNED_INPUTS_AND_INSTALLED_ENVIRONMENT',
                          'external_files': len(input_contract['files']), 'binding': binding,
                          'adoption_candidates': len(adopt), 'models_fitted': 0,
                          'raw_inputs_read_or_modified': 0, 'cold_environment_restored': False}))
        return
    if not a.output:
        raise ValueError('--run/--resume requires --output')
    out = a.output.resolve()
    protected = [root, config_path, *inputs.values()]
    if a.adopt_completed_checkpoints:
        protected.append(a.adopt_completed_checkpoints.resolve())
    if any(inside(out, p) or inside(p, out) for p in protected):
        raise ValueError('Output overlaps a protected source, package or input location')
    if a.run:
        out.mkdir(parents=True, exist_ok=False)
        for name in ('checkpoints', 'logs', 'provenance', 'data/intermediate', 'data/derived', 'data/raw_ledger', 'data/output'):
            (out / name).mkdir(parents=True, exist_ok=True)
        dump(out / 'provenance/identity.json', {'binding': binding, 'components': bindings})
        if native_identity:
            dump(out / 'provenance/native_input_identity.json', native_identity)
    elif not out.is_dir() or read(out / 'provenance/identity.json')['binding'] != binding:
        raise ValueError('Resume identity mismatch; the output will not be changed')
    done = out / 'provenance/COMPLETED.json'
    if done.exists():
        old = read(done)
        if old['binding'] != binding:
            raise ValueError('Completion binding mismatch')
        for name, expected in old['files'].items():
            pin(out / name, expected)
        print(json.dumps({'status': 'REUSED_VERIFIED_COMPLETED_PORTABLE_RUN', 'models_fitted': 0,
                          'binding': binding, 'output_files_verified': len(old['files'])}))
        return
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    log = out / 'logs' / (stamp + '.json')
    dump(log, {'status': 'STARTED', 'argv': sys.argv, 'python': sys.version,
               'binding': binding, 'package_manifest_sha256': sha(root / 'MANIFEST.json')})
    job = {'mode': 'resume' if a.resume else 'run', 'binding': binding,
           'inputs': {k: str(v) for k, v in inputs.items()}, 'adopt': adopt}
    jobpath = out / 'provenance/job.json'
    dump(jobpath, job)
    for name in ('data/derived/samples.tsv', 'data/derived/retained_PROGENy_comparator_basis.tsv',
                 'data/raw_ledger/ENA_quantification_manifest.tsv'):
        dst = out / name
        if dst.exists() and sha(dst) != sha(root / name):
            raise ValueError('Staged fixed input changed: ' + name)
        if not dst.exists():
            shutil.copyfile(root / name, dst)
    # Reconstruct QC directly from the independently accepted native outputs.
    samples = list(csv.DictReader((root / 'data/derived/samples.tsv').open(encoding='utf-8'), delimiter='\t'))
    qc = []
    for row in samples:
        native = inputs['native_root'] / ('linux_full_' + row['run'] + '_01')
        f = read(native / 'fastp.json')['summary']
        k = read(native / 'quant/run_info.json')
        if k['n_targets'] != 146190 or f['after_filtering']['total_reads'] != 2 * k['n_processed']:
            raise ValueError('Native read/target accounting mismatch')
        qc.append({'run': row['run'], 'input_pairs': f['before_filtering']['total_reads'] // 2,
                   'retained_pairs': k['n_processed'], 'pseudoaligned_pairs': k['n_pseudoaligned'],
                   'pseudoaligned_fraction': k['n_pseudoaligned'] / k['n_processed'],
                   'post_q30_fraction': f['after_filtering']['q30_rate'],
                   'post_gc_fraction': f['after_filtering']['gc_content']})
    with (out / 'data/output/complete_library_qc.tsv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(qc[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(qc)
    commands = []

    def execute(command, name):
        proc = subprocess.run(command, capture_output=True, text=True, env=env)
        (out / 'logs' / (stamp + '_' + name + '.stdout.txt')).write_text(proc.stdout, encoding='utf-8')
        (out / 'logs' / (stamp + '_' + name + '.stderr.txt')).write_text(proc.stderr, encoding='utf-8')
        commands.append({'argv': command, 'returncode': proc.returncode})
        dump(log, {'status': 'RUNNING' if proc.returncode == 0 else 'STOPPED_WITH_ERROR',
                   'binding': binding, 'commands': commands, 'python': sys.version})
        if proc.returncode:
            raise RuntimeError('Stage stopped; inspect retained logs before any retry: ' + name)

    execute([a.rscript, '--vanilla', str(root / 'pipeline/engine.R'), str(root), str(jobpath),
             str(out), os.pathsep.join(a.r_library)], 'analysis')
    basisdir = out / 'data/derived'
    if not (basisdir / 'receipt.json').exists():
        dump(basisdir / 'receipt.json', {'output_sha256': sha(basisdir / 'retained_PROGENy_comparator_basis.tsv'),
                                       'source': 'Frozen basis byte-verified against public package manifest'})
    comparison = out / 'data/output/comparison'
    cp = out / 'checkpoints/comparison.json'
    if cp.exists():
        z = read(cp)
        if z['binding'] != binding:
            raise ValueError('Comparison checkpoint binding mismatch')
        for name, expected in z['files'].items():
            pin(comparison / name, expected)
    elif comparison.exists():
        raise ValueError('Partial comparison directory requires inspection; no automatic repetition')
    else:
        execute([sys.executable, '-I', '-B', str(root / 'analysis/original/compare_PRJNA886048_retained_programs_20261025.py'),
                 '--basis-dir', str(basisdir), '--contrasts', str(out / 'data/output/signed_activity_contrasts.tsv'),
                 '--coverage', str(out / 'data/output/signed_activity_coverage.tsv'), '--out', str(comparison)], 'comparison')
        dump(cp, {'binding': binding, 'files': {p.name: {'sha256': sha(p), 'bytes': p.stat().st_size}
                                               for p in comparison.iterdir() if p.is_file()}})
    execute([a.rscript, '--vanilla', str(root / 'pipeline/verify_tables.R'), str(out)], 'content_readback')
    execute([a.rscript, '--vanilla', str(root / 'pipeline/compare_outputs.R'), str(root), str(out)], 'frozen_comparison')
    result = read(out / 'provenance/R_execution.json')
    for item in active_files:
        path = inputs[item['role']]
        if item.get('relative'):
            path = path / item['relative']
        pin(path, item)
    files = {p.relative_to(out).as_posix(): {'sha256': sha(p), 'bytes': p.stat().st_size}
             for folder in ('data', 'checkpoints') for p in sorted((out / folder).rglob('*')) if p.is_file()}
    dump(done, {'status': 'COMPLETE_PORTABLE_RUN_MATCHES_FROZEN_TABLES', 'binding': binding,
                'files': files, 'R_execution_counts': result['counts'],
                'cold_environment_restore_tested_by_this_command': False,
                'claim_scope': 'Acute SCI context only; no protected PSS validation'})
    dump(log, {'status': 'COMPLETE', 'binding': binding, 'commands': commands,
               'R_execution_counts': result['counts'], 'python': sys.version})
    print(json.dumps({'status': 'COMPLETE_PORTABLE_RUN_MATCHES_FROZEN_TABLES', 'binding': binding,
                      'R_execution_counts': result['counts'], 'output_files': len(files)}))


if __name__ == '__main__':
    main()
