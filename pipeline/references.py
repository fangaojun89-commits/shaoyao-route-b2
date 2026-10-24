"""Read-only published reference-source integrity and selected binding checks."""
import argparse
import hashlib
import json
from pathlib import Path


def identity(path):
    raw=path.read_bytes()
    return len(raw),hashlib.sha256(raw).hexdigest()


def safe(root, name):
    rel=Path(name)
    if rel.is_absolute() or '..' in rel.parts:
        raise ValueError('Unsafe relative source path')
    return root/rel


def verify(repo):
    manifest=json.loads((repo/'MANIFEST.json').read_text(encoding='utf-8'))
    for name,info in manifest['files'].items():
        if identity(safe(repo,name))!=(info['bytes'],info['sha256']):
            raise ValueError('Release identity mismatch: '+name)
    graph=json.loads((repo/'analysis/reference_sources/SOURCE_GRAPH.json').read_text(encoding='utf-8'))
    records=graph['sources']
    if len(records)!=151 or len({r['original_path'] for r in records})!=151:
        raise ValueError('Source inventory changed')
    texts={}
    derivatives=[]
    supplements=[]
    for record in records:
        path=safe(repo,record['public_path'])
        if identity(path)!=(record['public_bytes'],record['public_sha256']):
            raise ValueError('Public source identity mismatch')
        texts[record['original_path']]=path.read_text(encoding='utf-8-sig')
        if record['transformation']=='NONE':
            if (record['original_bytes'],record['original_sha256'])!=(record['public_bytes'],record['public_sha256']):
                raise ValueError('Unchanged-source identity mismatch')
        elif record['transformation']=='ONE_PERSONAL_PROJECT_PATH_REPLACED_BY_PROJECT_ROOT_PLACEHOLDER':
            derivatives.append(record['original_path'])
            if texts[record['original_path']].count('<PROJECT_ROOT>')!=1:
                raise ValueError('Explicit path substitution changed')
        else:
            raise ValueError('Unrecognized source transformation')
        if record['origin']=='current_explicit_reference_missing_from_r02':
            supplements.append(record)
    if set(derivatives)!={'code/R/verify_JGW03_complete_reference_v2_20261015.R',
                          'code/python/verify_gse16710_anchor_baselines_independent_20261014.py'}:
        raise ValueError('Path derivative inventory changed')
    if len(supplements)!=3:
        raise ValueError('Additional source inventory changed')
    bindings=[]
    for record in supplements:
        callers=[name for name,text in texts.items() if name!=record['original_path'] and record['original_sha256'] in text]
        if not callers:
            raise ValueError('Additional source lacks literal digest binding')
        bindings.append({'source':record['original_path'],'callers':callers})
    for edge in graph['references']:
        if edge['from'] not in texts or edge['to'] not in texts:
            raise ValueError('Unmapped static reference edge')
    return {'status':'PASS_PUBLIC_REFERENCE_SOURCE_IDENTITIES_AND_SELECTED_BINDINGS',
            'selected_sources':151,'unchanged_sources':149,'explicit_path_derivatives':2,
            'new_source_files':sum(not r['reuses_prior_public_file'] for r in records),
            'additional_source_bindings':bindings,'scientific_model_runs':0,
            'scope':'Public identity and static source bindings; not complete dynamic or scientific reconstruction'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['verify'])
    parser.parse_args()
    print(json.dumps(verify(Path(__file__).resolve().parents[1]),indent=2))
