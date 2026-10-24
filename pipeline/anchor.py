"""Read-only integrity and external-input checks for the published anchor source."""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def rows(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def verify(repo):
    manifest = json.loads((repo/'MANIFEST.json').read_text(encoding='utf-8'))
    for name, info in manifest['files'].items():
        p = Path(name)
        if p.is_absolute() or '..' in p.parts:
            raise ValueError('Invalid manifest path')
        p = repo/p
        if p.stat().st_size != info['bytes'] or sha(p) != info['sha256']:
            raise ValueError('Release file mismatch: '+name)
    anchor = repo/'analysis/GSE16710'
    provenance = json.loads((anchor/'PROVENANCE.json').read_text(encoding='utf-8'))
    for record in provenance['source_files']:
        if sha(anchor/record['public_path']) != record['public_sha256']:
            raise ValueError('Source provenance mismatch')
    samples = rows(anchor/'metadata/samples.tsv')
    qc = rows(anchor/'data/processed/GSE16710/qc_inclusion_decision.tsv')
    ids = [x['sample_id'] for x in samples]
    if len(ids) != 27 or len(set(ids)) != 27 or {x['dataset_id'] for x in samples}!={'GSE16710'}:
        raise ValueError('Submitted-unit ledger mismatch')
    if len(qc)!=27 or {x['sample_id'] for x in qc}!=set(ids):
        raise ValueError('QC registry mismatch')
    excluded = [x['sample_id'] for x in qc if x['primary_analysis_decision']!='retain']
    if excluded!=['GSM418682']:
        raise ValueError('Frozen QC selection changed')
    gate = rows(anchor/'expected/S1_GATE.tsv')
    families = collections.Counter(x['collection'] for x in gate)
    if dict(families)!={'DoRothEA':101,'PROGENy':14,'Reactome':1656}:
        raise ValueError('Retained S1 family counts changed')
    if any(x['S1_pass']!='FALSE' for x in gate):
        raise ValueError('Retained S1 admission changed')
    return {'status':'PASS_RELEASE_AND_RETAINED_METADATA',
            'manifest_payloads_checked':len(manifest['files']),
            'submitted_arrays':27,'retained_primary_arrays':26,
            'retained_S1_rows':len(gate),'retained_S1_pass':0,
            'scope':'Hash and metadata verification only; no model or raw-to-result validation'}


def inputs(repo, root):
    verify(repo)
    records=json.loads((repo/'analysis/GSE16710/EXTERNAL_INPUTS.json').read_text())['files']
    findings=[]
    for record in records:
        p=root/record['relative_path']
        status='MISSING'
        if p.is_file():
            status='MATCH' if p.stat().st_size==record['bytes'] and sha(p)==record['sha256'] else 'MISMATCH'
        findings.append({'relative_path':record['relative_path'],'status':status})
    ok=all(x['status']=='MATCH' for x in findings)
    return {'status':'PASS_ENUMERATED_INPUT_IDENTITIES' if ok else 'INPUTS_INCOMPLETE_OR_DIFFERENT',
            'files':findings,'scope':'Enumerated raw/processed identities only; software installation and complete resource rebuilding are not checked'},0 if ok else 3


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['verify','inspect-inputs'])
    parser.add_argument('--input-root',type=Path)
    args=parser.parse_args()
    repo=Path(__file__).resolve().parents[1]
    if args.command=='verify':
        result,code=verify(repo),0
    else:
        if args.input_root is None:
            parser.error('inspect-inputs requires --input-root')
        result,code=inputs(repo,args.input_root.resolve())
    print(json.dumps(result,indent=2))
    raise SystemExit(code)


if __name__=='__main__':
    main()
