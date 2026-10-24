"""Verify or restore exact accepted count-input bytes; never run scientific code."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def relative_path(value):
    path = Path(value)
    if path.is_absolute() or '..' in path.parts or '\\' in value or ':' in value:
        raise ValueError('Unsafe input destination')
    return path


def verified_payloads():
    record = json.loads((REPO / 'data/count_inputs/PROVENANCE.json').read_text(encoding='utf-8'))
    payloads = []
    destinations = set()
    for item in record['files']:
        source = relative_path(item['public_path'])
        destination = relative_path(item['runtime_relative_path'])
        if str(destination) in destinations:
            raise ValueError('Repeated destination')
        destinations.add(str(destination))
        packed = (REPO / source).read_bytes()
        if (len(packed), hashlib.sha256(packed).hexdigest()) != (item['bytes'], item['sha256']):
            raise ValueError('Compressed input identity mismatch: ' + str(source))
        raw = gzip.decompress(packed)
        if (len(raw), hashlib.sha256(raw).hexdigest()) != (item['original_bytes'], item['original_sha256']):
            raise ValueError('Original input identity mismatch: ' + str(source))
        payloads.append((destination, raw))
    if len(payloads) != 11:
        raise ValueError('Expected eleven accepted input files')
    return payloads


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('verify', 'materialize'))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.command == 'materialize' and args.output is None:
        parser.error('materialize requires --output')
    if args.output is not None and args.output.exists():
        parser.error('output must be a new directory')
    payloads = verified_payloads()
    if args.command == 'materialize':
        args.output.mkdir(parents=True, exist_ok=False)
        for relative, raw in payloads:
            target = args.output / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            with target.open('xb') as handle:
                handle.write(raw)
    print(json.dumps({'status': 'PASS_ACCEPTED_INPUT_BYTES', 'files': len(payloads),
                      'original_bytes': sum(len(raw) for _, raw in payloads),
                      'materialized': args.command == 'materialize',
                      'scientific_model_runs': 0,
                      'scope': 'Exact accepted input identities; not whole-study restoration or scientific revalidation'}))


if __name__ == '__main__':
    main()
