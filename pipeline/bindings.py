"""Verify runtime source aliases or stage exact public sources without execution."""
import argparse
import hashlib
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def identity(path):
    raw = path.read_bytes()
    return len(raw), hashlib.sha256(raw).hexdigest()


def relative_path(name):
    p = Path(name)
    if p.is_absolute() or '..' in p.parts or ':' in name or '\\' in name:
        raise ValueError('Unsafe relative namespace')
    return p


def bindings():
    data = json.loads((REPO / 'analysis/runtime_bindings/G205_BINDINGS.json').read_text(encoding='utf-8'))
    rows = []
    for source in data['new_unchanged_sources']:
        if source['transformation'] != 'NONE' or identity(REPO / relative_path(source['public_path'])) != (source['bytes'], source['sha256']):
            raise ValueError('Original source identity mismatch')
    for row in data['runtime_file_pins']:
        for source in row['public_matches']:
            if identity(REPO / relative_path(source)) != (row['bytes'], row['sha256']):
                raise ValueError('Runtime source pin mismatch')
            rows.append({'destination': row['asset_namespace'], 'source': source,
                         'bytes': row['bytes'], 'sha256': row['sha256'],
                         'evidence': row['source_record'] + ':' + row['field']})
    for row in data['source_assertions']:
        asserting = REPO / relative_path(row['asserting_public_source'])
        if identity(asserting)[1] != row['asserting_source_sha256']:
            raise ValueError('Asserting source changed')
        import re
        pattern = r'pin\(file.path\(run, "' + re.escape(row['alias']) + r'"\)\)\$sha256 == "' + row['sha256'] + '"'
        if not re.search(pattern, asserting.read_text(encoding='utf-8-sig')):
            raise ValueError('Literal source assertion absent')
        parents = {str(relative_path(r['asset_namespace']).parent) for r in data['runtime_file_pins']
                   if r['run'] == row['run'] and r['field'] == '/worker'}
        if len(parents) != 1:
            raise ValueError('Ambiguous asserted runtime namespace')
        for source in row['public_matches']:
            if identity(REPO / relative_path(source)) != (row['bytes'], row['sha256']):
                raise ValueError('Asserted source identity mismatch')
            rows.append({'destination': next(iter(parents)).replace('\\', '/') + '/' + row['alias'],
                         'source': source, 'bytes': row['bytes'], 'sha256': row['sha256'],
                         'evidence': row['binding_kind']})
    selected = {}
    for row in rows:
        relative_path(row['destination'])
        previous = selected.get(row['destination'])
        if previous and (previous['bytes'], previous['sha256']) != (row['bytes'], row['sha256']):
            raise ValueError('Conflicting identities for one runtime namespace')
        if previous:
            previous['evidence'].append(row['evidence'])
        else:
            selected[row['destination']] = dict(row, evidence=[row['evidence']])
    return data, sorted(selected.values(), key=lambda r: r['destination'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['verify', 'stage-sources'])
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    data, rows = bindings()
    result = {'status': 'RUNTIME_SOURCE_BINDINGS_VERIFIED', 'retained_record_count': len(data['receipt_provenance']),
              'runtime_file_pins': len(data['runtime_file_pins']), 'unique_source_destinations': len(rows),
              'external_assets_read': False, 'scientific_model_runs': 0, 'complete_execution_layout': False}
    if args.command == 'stage-sources':
        if args.output is None:
            parser.error('--output required')
        output = args.output.absolute()
        if output.exists() or output.is_symlink() or not output.parent.is_dir():
            parser.error('Use a new destination under an existing parent')
        output.mkdir()
        for row in rows:
            destination = output / relative_path(row['destination'])
            destination.parent.mkdir(parents=True, exist_ok=True)
            with destination.open('xb') as handle:
                handle.write((REPO / row['source']).read_bytes())
            if identity(destination) != (row['bytes'], row['sha256']):
                raise ValueError('Staged source identity mismatch')
        result['status'] = 'SOURCE_LAYOUT_STAGED_NOT_EXECUTED'
        with (output / 'SOURCE_LAYOUT_RECEIPT.json').open('x', encoding='utf-8') as handle:
            json.dump(dict(result, mappings=rows), handle, indent=2)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
