"""Verify original branch-source identities without evaluating scientific code."""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path
import re


def identity(path):
    raw = path.read_bytes()
    return len(raw), hashlib.sha256(raw).hexdigest()


def check(root, record, key='path'):
    relative = Path(record[key])
    if relative.is_absolute() or '..' in relative.parts:
        raise ValueError('Unsafe relative path')
    if identity(root / relative) != (record['bytes'], record['sha256']):
        raise ValueError('Identity mismatch: ' + str(relative))


def verify(repo):
    manifest = json.loads((repo / 'MANIFEST.json').read_text(encoding='utf-8'))
    for name, record in manifest['files'].items():
        check(repo, dict(record, path=name))
    base = repo / 'analysis/study_sources'
    provenance = json.loads((base / 'PROVENANCE.json').read_text(encoding='utf-8'))
    records = provenance['source_files']
    if len(records) != 43 or len({x['path'] for x in records}) != 43:
        raise ValueError('Original source selection changed')
    for record in records:
        check(base, record)
        if record['transformation'] != 'NONE':
            raise ValueError('Original source transformed')
    roles = json.loads((base / 'G205_SOURCE_ROLES.json').read_text(encoding='utf-8'))['roles']
    pins = {}
    for name in ('gse205906_ensembl101_checked_consumer_v1_20261006.R',
                 'gse205906_ensembl101_checked_model_RAM_worker_v4_20261006.R'):
        source = (base / 'code/R' / name).read_text(encoding='utf-8-sig')
        for role, size, digest in re.findall(r'(\w+)\s*=\s*pin\((\d+)L,\s*"([0-9a-f]{64})"\)', source):
            if role in pins:
                raise ValueError('Duplicate source role')
            pins[role] = (int(size), digest)
    if len(roles) != 20 or {x['role'] for x in roles} != set(pins):
        raise ValueError('G205 source roles changed')
    for role in roles:
        check(base, role)
        if (role['bytes'], role['sha256']) != pins[role['role']]:
            raise ValueError('G205 embedded source pin mismatch')
    sample_path = base / 'metadata/samples.tsv'
    if identity(sample_path)[1] != provenance['sample_projection']['public_sha256']:
        raise ValueError('Sample projection changed')
    with sample_path.open(encoding='utf-8-sig', newline='') as handle:
        samples = list(csv.DictReader(handle, delimiter='\t'))
    counts = dict(collections.Counter(x['dataset_id'] for x in samples))
    if counts != {'GSE277064':63, 'GSE289929':12, 'GSE227494':12, 'GSE299063':12}:
        raise ValueError('Sample projection membership changed')
    if len({(x['dataset_id'], x['sample_id']) for x in samples}) != 99:
        raise ValueError('Duplicate submitted unit')
    return {'status':'PASS_RELEASE_SOURCE_IDENTITIES_AND_ROLE_PINS',
            'release_payloads':len(manifest['files']), 'unchanged_source_and_config_files':43,
            'R_sources':38, 'G205_source_roles':20, 'submitted_units':99,
            'scientific_model_runs':0,
            'scope':'Integrity and role binding only; no evaluation, environment restore or model reproduction'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['verify'])
    parser.parse_args()
    print(json.dumps(verify(Path(__file__).resolve().parents[1]), indent=2))
