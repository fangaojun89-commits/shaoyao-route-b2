"""Check one canonical native record/content bundle without launching any job."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', type=Path, required=True)
    parser.add_argument('--tx2gene', type=Path, required=True)
    parser.add_argument('--run', required=True)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    root = here.parents[1]
    manifest = json.loads((root / 'MANIFEST.json').read_text(encoding='utf-8-sig'))
    names = ['pipeline/native/verify.py',
             'pipeline/native/validate_PRJNA886048_native_content_20261028.py',
             'pipeline/native/validate_PRJNA886048_native_execution_20261028.py',
             'data/raw_ledger/ENA_quantification_manifest.tsv']
    for name in names:
        path, expected = root / name, manifest['files'][name]
        if path.stat().st_size != expected['bytes'] or hashlib.sha256(path.read_bytes()).hexdigest() != expected['sha256']:
            raise ValueError('Package component identity differs: ' + name)
    spec = importlib.util.spec_from_file_location('native_record_gate', here / 'validate_PRJNA886048_native_execution_20261028.py')
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)
    result = gate.validate(args.native, args.tx2gene, root / 'data/raw_ledger/ENA_quantification_manifest.tsv', args.run)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
