"""Connect frozen JGW03 inputs/resources to the original reader without fitting."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

REPO = Path(__file__).resolve().parents[1]
RESOURCES = (
    'progeny_human_top500.tsv', 'progeny_mouse_top500.tsv',
    'dorothea_human_AB.tsv', 'dorothea_mouse_AB.tsv',
    'dorothea_rat_AB_one_to_one.tsv', 'reactome97_human_long.tsv',
    'quickgo_syngo_mammalian_direct_human.gmt',
    'ncbi_human_mouse_rat_one_to_one.tsv',
)


def identity(path):
    value = hashlib.sha256()
    size = 0
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            value.update(block)
            size += len(block)
    return {'bytes': size, 'sha256': value.hexdigest()}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--resource-root', type=Path, required=True,
                        help='Completed build_resources.py output root')
    parser.add_argument('--rscript', required=True)
    parser.add_argument('--r-library', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output must be a new directory')
    library = args.r_library.resolve(strict=True)
    if not library.is_dir():
        parser.error('R library must be a directory')
    manifest = json.loads((REPO / 'MANIFEST.json').read_text(encoding='utf-8'))
    admission = 'analysis/reference_sources/code/R/JGW03_input_admission_v2_20261015.R'
    core = 'analysis/reference_sources/code/R/JGW03_label_reference_core_20261015.R'
    dependencies = ['pipeline/prepare_jgw_inputs.py', 'analysis/jgw_preparation/prepare.R',
                    admission, core, 'data/jgw03_inputs/PROVENANCE.json',
                    'analysis/resource_build/CONTRACT.json']
    source_pins = {}
    for name in dependencies:
        actual = identity(REPO / name)
        if actual != manifest['files'][name]:
            raise ValueError('Release source identity mismatch: ' + name)
        source_pins[name] = actual
    provenance = json.loads((REPO / 'data/jgw03_inputs/PROVENANCE.json').read_text(encoding='utf-8'))
    resources = json.loads((REPO / 'analysis/resource_build/CONTRACT.json').read_text(encoding='utf-8'))
    expected_resources = {r['path']: r['identity'] for r in resources['outputs']}
    staging = []
    for record in provenance['files']:
        source = REPO / record['public_path']
        expected = {k: record[k] for k in ('bytes', 'sha256')}
        relative = record['accepted_runtime_namespace']
        # Archive namespaces contain one project prefix for the nine table inputs.
        destination = relative if relative.startswith('project/') else 'project/' + relative
        if identity(source) != expected:
            raise ValueError('Frozen JGW input mismatch: ' + record['public_path'])
        staging.append((source, destination, expected))
    for name in RESOURCES:
        relative = 'data/processed/gene_sets/' + name
        source = (args.resource_root / relative).resolve(strict=True)
        expected = expected_resources[relative]
        if identity(source) != expected:
            raise ValueError('Frozen JGW resource mismatch: ' + name)
        staging.append((source, 'project/' + relative, expected))
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output / 'receipts').mkdir()
    (output / 'prepared').mkdir()
    input_pins = {}
    for source, relative, expected in staging:
        destination = output / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        if identity(destination) != expected:
            raise ValueError('Staged input mismatch: ' + relative)
        input_pins[str(source)] = expected
    spec = dict(repo=REPO.as_posix(), output=output.as_posix(),
                project=(output / 'project').as_posix(), admission=admission,
                core=core, source_pins=source_pins,
                datasets=['GSE289929', 'GSE227494', 'GSE299063'])
    spec_path = output / 'receipts/RUN.json'
    write_json(spec_path, spec)
    command = [args.rscript, '--vanilla', (REPO / 'analysis/jgw_preparation/prepare.R').as_posix(),
               spec_path.as_posix()]
    env = os.environ.copy()
    locale = 'English_United States.utf8' if os.name == 'nt' else 'C.UTF-8'
    env.update(R_LIBS_USER=library.as_posix(), RENV_CONFIG_AUTOLOADER_ENABLED='FALSE',
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               LANG=locale, LC_ALL=locale, LC_CTYPE=locale)
    start = datetime.datetime.now(datetime.timezone.utc).isoformat()
    with (output / 'receipts/stdout.txt').open('wb') as stdout, (output / 'receipts/stderr.txt').open('wb') as stderr:
        result = subprocess.run(command, env=env, cwd=output, stdout=stdout, stderr=stderr, check=False)
    sources_same = all(identity(REPO / name) == pin for name, pin in source_pins.items())
    inputs_same = all(identity(Path(name)) == pin for name, pin in input_pins.items())
    qa_path = output / 'receipts/PREPARATION_QA.json'
    qa = json.loads(qa_path.read_text(encoding='utf-8')) if qa_path.exists() else {}
    success = result.returncode == 0 and sources_same and inputs_same and qa.get('status') == 'PASS_ORIGINAL_JGW_INPUT_READER_WITHOUT_SCORING_OR_FITS'
    outputs = {p.relative_to(output).as_posix(): identity(p) for p in output.rglob('*') if p.is_file()}
    write_json(output / 'receipts/PARENT.json', dict(status='PASS' if success else 'FAIL',
        started_utc=start, finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        command=command, returncode=result.returncode, input_pins=input_pins,
        source_pins=source_pins, original_sources_unchanged=sources_same,
        original_inputs_unchanged=inputs_same, outputs=outputs,
        real_data_models_run=0, new_allocations_generated=0, scores_recomputed=0))
    print(json.dumps({'status': 'PASS' if success else 'FAIL', 'output': str(output)}))
    if not success:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
