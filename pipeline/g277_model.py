"""G277 original scientific sequence and non-fitting reuse of accepted outputs.

check stages exact inputs and checks the native environment without fitting.
run explicitly executes all five original stages in a new isolated project.
report-saved copies pinned accepted outputs and executes only reporting/validation.
"""
import argparse
import csv
import datetime
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess

REPO = Path(__file__).resolve().parents[1]
BASE = 'analysis/g277_execution/'


def identity(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''): h.update(b)
    return {'bytes': Path(path).stat().st_size, 'sha256': h.hexdigest()}


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def write(path, obj):
    Path(path).write_bytes((json.dumps(obj, indent=2, ensure_ascii=False)+'\n').encode())


def need(ok, message):
    if not ok: raise ValueError(message)


def copy_checked(source, dest, pin):
    need(identity(source) == pin, 'Changed input: '+str(source))
    dest.parent.mkdir(parents=True, exist_ok=True)
    need(not dest.exists(), 'Destination already exists: '+str(dest))
    shutil.copyfile(source, dest)
    need(identity(dest) == pin, 'Copy verification failed: '+str(dest))


def projections(project, output, accepted, codec_path, expected):
    spec = importlib.util.spec_from_file_location('original_g277_reporting_codec', codec_path)
    codec = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(codec)
    records = []
    for rec in accepted['projections']:
        binding = dict(rec['binding'])
        src = project/binding['source_path']
        blob = src.read_bytes()
        if expected:
            need(identity(src) == rec['source_identity'], 'Accepted source changed')
        else:
            # Future reruns receive a new binding, never a forged old receipt.
            table = codec.read_tsv(codec.gzip.decompress(blob))
            need(table[0] == codec.HEADER, 'Unexpected original source schema')
            binding.update(source_sha256=codec.digest(blob), rows=len(table)-1,
                source_cell_stream_sha256=codec.stream_digest(table[0], table[1:]),
                retained_cell_stream_sha256=codec.stream_digest(codec.VIEW_HEADER,
                    (tuple(row[i] for i in codec.KEEP) for row in table[1:])))
        artifacts = codec.project(blob, binding)
        # Independent comparison also proves every retained lexical value/order.
        original = list(csv.reader(io.StringIO(codec.gzip.decompress(blob).decode('utf-8-sig')), delimiter='\t'))
        kept = list(csv.reader(io.StringIO(artifacts['view'].decode('utf-8')), delimiter='\t'))
        selected = [i for i, name in enumerate(original[0]) if name not in {'ql_se','ci_low','ci_high'}]
        need(kept == [[row[i] for i in selected] for row in original], 'Retained values/order differ')
        rel = binding['source_path'].removeprefix('results/GSE277064/').removesuffix('.tsv.gz')
        products = {}
        for key, blob_out in artifacts.items():
            suffix = '.json' if key == 'schema' else '.tsv'
            dest = output/key/(rel+suffix)
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(blob_out)
            products[key] = identity(dest)
            if expected:
                need(products[key] == rec['expected_projection'][key], 'Accepted projection differs: '+key)
        records.append({'source': binding['source_path'], 'rows': binding['rows'], 'outputs': products,
                        'all_retained_cells_equal': True, 'accepted_projection_exact': expected})
    need(len(records) == 24, 'Expected 24 fixed gene families')
    write(output/'PROJECTION.json', {'status':'PASS_ORIGINAL_CODEC_AND_ALL_RETAINED_CELLS',
        'records':records, 'rows':sum(r['rows'] for r in records),
        'omitted_fields':['ql_se','ci_low','ci_high'], 'models_run_by_projection':False})
    (output/'README.md').write_bytes(b'# Current gene reporting view\n\n'
        b'Use view/ for gene reporting. Only ql_se, ci_low and ci_high are omitted. '
        b'Every retained value and ordering is unchanged. signed_sqrt_F is a ranking statistic, '
        b'not a coefficient standard error. quarantine/ retains excluded legacy values for traceability only. '
        b'Activity-score intervals are separate method-specific quantities and are unchanged. '
        b'Historical narrative and validation in project/results/GSE277064 do not reinstate excluded fields. '
        b'This reproduces the approved TSV projections, not the S14/S16 workbook formatting or whole manuscript.\n')
    return records


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=['check','run','report-saved'])
    p.add_argument('--resources', type=Path, help='exact rebuilt frozen resource RDS; check/run only')
    p.add_argument('--saved-results', type=Path, help='accepted GSE277064 result directory; report-saved only')
    p.add_argument('--rscript', required=True, type=Path)
    p.add_argument('--r-library', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    a = p.parse_args()
    need(not a.output.exists(), 'Use a fresh directory; existing results are preserved')
    out = a.output.resolve()
    need(not out.is_relative_to(REPO), 'Outputs belong outside the public source checkout')
    if a.mode == 'report-saved':
        need(a.saved_results is not None and a.resources is None, 'report-saved requires only --saved-results')
    else:
        need(a.resources is not None and a.saved_results is None, 'check/run require only --resources')
    rs = a.rscript.resolve(strict=True)
    lib = a.r_library.resolve(strict=True)
    need(rs.is_file() and lib.is_dir(), 'Explicit Rscript and R library required')
    manifest = read(REPO/'MANIFEST.json')
    contract = read(REPO/BASE/'CONTRACT.json')
    accepted = read(REPO/BASE/'ACCEPTED_OUTPUTS.json')
    inputs = read(REPO/'data/g277_inputs/PROVENANCE.json')['files']
    dependencies = ['pipeline/g277_model.py', BASE+'CONTRACT.json',BASE+'ACCEPTED_OUTPUTS.json',
                    BASE+'preflight.R', 'data/g277_inputs/PROVENANCE.json',
                    'analysis/GSE16710/expected/S1_GATE.tsv']
    dependencies += [r['public_path'] for r in contract['sources']]
    dependencies += [r['public_path'] for r in inputs]
    dependencies += ['analysis/study_sources/config/'+name for name in ('contrasts.yaml','gene_sets.yaml','thresholds.yaml')]
    for rel in dependencies:
        need(identity(REPO/rel) == manifest['files'][rel], 'Release dependency changed: '+rel)
    input_pins = {str(REPO/rel): manifest['files'][rel] for rel in dependencies}
    gate = list(csv.DictReader((REPO/'analysis/GSE16710/expected/S1_GATE.tsv').open(encoding='utf-8'), delimiter='\t'))
    need(len(gate) == 1771 and all(r['S1_pass'] == 'FALSE' for r in gate), 'This sequence requires the retained zero-admission S1')
    if a.mode == 'report-saved':
        saved = a.saved_results.resolve(strict=True)
        need(not out.is_relative_to(saved), 'Do not place outputs inside the accepted input directory')
        actual_names = {f.relative_to(saved).as_posix() for f in saved.rglob('*') if f.is_file()}
        need(actual_names == set(accepted['files']), 'Accepted output set must contain exactly the pinned 169 files')
        for rel, pin in accepted['files'].items():
            need(identity(saved/rel) == pin, 'Accepted output identity mismatch: '+rel)
            input_pins[str(saved/rel)] = pin
    else:
        res = a.resources.resolve(strict=True)
        need(identity(res) == contract['resource_identity'], 'Frozen resource identity mismatch')
        input_pins[str(res)] = contract['resource_identity']
    out.mkdir(parents=True)
    receipts = out/'receipts'; receipts.mkdir()
    project = out/'project'; project.mkdir()
    staged = {}
    for rec in contract['sources']:
        copy_checked(REPO/rec['public_path'], project/rec['runtime_path'], rec['identity'])
        staged[rec['runtime_path']] = rec['identity']
    for rec in inputs:
        copy_checked(REPO/rec['public_path'], project/rec['runtime_path'], rec['identity'])
        staged[rec['runtime_path']] = rec['identity']
    for name in ('contrasts.yaml','gene_sets.yaml','thresholds.yaml'):
        rel = 'analysis/study_sources/config/'+name
        copy_checked(REPO/rel, project/'config'/name, manifest['files'][rel])
        staged['config/'+name] = manifest['files'][rel]
    if a.mode == 'report-saved':
        for rel, pin in accepted['files'].items():
            # The original successor summary intentionally refuses this pre-existing output.
            # Keep it at its source; the fresh project creates it from the exact template.
            if rel == 'manuscript_ready_methods_results.md': continue
            copy_checked(saved/rel, project/'results/GSE277064'/rel, pin)
    else:
        copy_checked(res, project/contract['resource_path'], contract['resource_identity'])
        staged[contract['resource_path']] = contract['resource_identity']
    env = os.environ.copy()
    env.update(LANG='en_US.UTF-8', LC_ALL='en_US.UTF-8', LC_CTYPE='en_US.UTF-8', TZ='Asia/Shanghai',
        R_LIBS_USER=str(lib), R_LIBS=str(lib), RENV_CONFIG_AUTOLOADER_ENABLED='FALSE',
        OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', R_DATATABLE_NUM_THREADS='1',
        GSE277064_ROAST_NROT='9999', GSE277064_RUN_FGSEA='TRUE')
    spec = {'repo':REPO.as_posix(),'project':project.as_posix(),'receipts':receipts.as_posix(),
            'contract':(REPO/BASE/'CONTRACT.json').as_posix(),'staged':staged,'mode':a.mode}
    write(receipts/'STAGED.json', spec)
    commands = [[str(rs),'--vanilla',str(REPO/BASE/'preflight.R'),str(receipts/'STAGED.json')]]
    if a.mode != 'check':
        steps = contract['sequence'][2:] if a.mode == 'report-saved' else contract['sequence']
        commands += [[str(rs),'--vanilla',str(project/script)] for script in steps]
    write(receipts/'COMMANDS.json', commands)
    result = {'mode':a.mode,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'python':platform.python_version(), 'platform':platform.platform(),'steps':[],
              'model_execution_requested':a.mode == 'run', 'whole_study_cold_restoration':False}
    try:
        for n, cmd in enumerate(commands):
            with (receipts/f'{n:02d}_stdout.txt').open('wb') as stdout, (receipts/f'{n:02d}_stderr.txt').open('wb') as stderr:
                proc = subprocess.run(cmd,cwd=project,env=env,stdout=stdout,stderr=stderr,check=False)
            result['steps'].append({'command':cmd,'exit_code':proc.returncode})
            write(receipts/'EXECUTION.json',result)
            need(proc.returncode == 0, 'Native step failed; retained logs at '+str(receipts))
        if a.mode == 'check':
            need(not (project/'results').exists(), 'No scientific outputs may be created by check')
            result['status'] = 'PASS_INPUT_STAGING_AND_NATIVE_ENVIRONMENT_NOT_A_MODEL_RUN'
        else:
            view = out/'reporting'; view.mkdir()
            recs = projections(project, view, accepted, project/'code/python/gse277064_reporting_core_v1.py', a.mode == 'report-saved')
            result.update(status='PASS_POSTPROCESSING_AND_REPORTING_PROJECTION',projection_families=len(recs),
                projection_rows=sum(r['rows'] for r in recs),
                actual_validation=list(csv.DictReader((project/'results/GSE277064/diagnostics/output_validation.tsv').open(encoding='utf-8'),delimiter='\t')))
            need(all(r['status'] == 'PASS' for r in result['actual_validation']), 'Native validation did not pass')
    except Exception as exc:
        result['status'] = 'FAILED'; result['error'] = str(exc)
        raise
    finally:
        result['inputs_unchanged'] = all(identity(Path(path)) == pin for path,pin in input_pins.items())
        result['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        write(receipts/'EXECUTION.json',result)
    need(result['inputs_unchanged'], 'Inputs changed during execution')
    print(json.dumps({'status':result['status'],'mode':a.mode,'output':out.as_posix()}))


if __name__ == '__main__': main()
