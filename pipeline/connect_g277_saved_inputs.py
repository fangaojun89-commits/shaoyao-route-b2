"""Connect saved G277 preprocessing to original annotation/resource/input consumers."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess

REPO = Path(__file__).resolve().parents[1]


def identity(path):
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'bytes': path.stat().st_size, 'sha256': digest}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared-root', required=True, type=Path,
                        help='Existing successful prepare_g277_inputs.py output')
    parser.add_argument('--resource-root', required=True, type=Path,
                        help='Existing successful build_resources.py output')
    parser.add_argument('--rscript', required=True)
    parser.add_argument('--r-library', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Use a new output directory')
    library = args.r_library.resolve(strict=True)
    if not library.is_dir():
        parser.error('R library must be a directory')
    contract_name = 'data/g277_saved_inputs/PROVENANCE.json'
    contract = json.loads((REPO / contract_name).read_text(encoding='utf-8'))
    manifest = json.loads((REPO / 'MANIFEST.json').read_text(encoding='utf-8'))
    names = ['pipeline/connect_g277_saved_inputs.py', 'analysis/g277_saved_inputs/connect.R',
             contract_name, contract['original_source'], contract['gene_method_core']]
    names.extend(r['public_path'] for r in contract['files'].values())
    pins = {name: identity(REPO / name) for name in names}
    for name, actual in pins.items():
        if actual != manifest['files'][name]:
            raise ValueError('Release source/input identity mismatch: ' + name)
    prepared = (args.prepared_root / 'prepared/GSE277064_unfitted_inputs.rds').resolve(strict=True)
    resource = (args.resource_root / contract['resource_path']).resolve(strict=True)
    if identity(prepared) != contract['prepared_input_identity']:
        raise ValueError('Prepared count/design object differs from tested original-reader output')
    if identity(resource) != contract['resource_identity']:
        raise ValueError('Frozen resource differs; reuse the verified original resource build')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    for name in ('prepared', 'receipts'):
        (output / name).mkdir()
    spec = dict(repo=str(REPO), output=str(output), prepared_input=str(prepared),
                resource_file=str(resource), source_pins=pins)
    spec_path = output / 'receipts/INVOCATION.json'
    spec_path.write_text(json.dumps(spec, indent=2) + '\n', encoding='utf-8')
    env = os.environ.copy()
    env['R_LIBS_USER'] = str(library)
    env['R_LIBS'] = str(library)
    env['LANG'] = 'en_US.UTF-8'
    env['LC_ALL'] = 'en_US.UTF-8'
    command = [args.rscript, '--vanilla', str(REPO / 'analysis/g277_saved_inputs/connect.R'), str(spec_path)]
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    result = subprocess.run(command, cwd=output, env=env, capture_output=True)
    (output / 'receipts/stdout.txt').write_bytes(result.stdout)
    (output / 'receipts/stderr.txt').write_bytes(result.stderr)
    receipt = dict(command=command, started_utc=started, native_exit=result.returncode,
                   source_pins=pins, prepared_input=identity(prepared), resource=identity(resource),
                   scientific_model_fits=0, scores_recomputed=False)
    if result.returncode == 0:
        products = list((output / 'prepared').glob('*.rds')) + [output / 'annotation/audit/ensembl_to_mouse_annotation.tsv']
        receipt['outputs'] = {str(p.relative_to(output)).replace('\\', '/'): identity(p) for p in products}
    (output / 'receipts/EXECUTION.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    if result.returncode:
        raise RuntimeError('Original annotation/resource/input connection failed; inspect retained receipts')
    print(json.dumps(dict(status='PASS_G277_SAVED_INPUT_CONNECTION', outputs=len(receipt['outputs']),
                          scientific_model_fits=0, scores_recomputed=False)))


if __name__ == '__main__':
    main()
