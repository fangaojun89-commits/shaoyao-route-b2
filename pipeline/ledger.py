"""Publish and check the dated submitted-unit ledger, without refitting models.

The export operation requires the preserved private source archive. Public
verification needs only the distributed metadata. It does not regenerate
sample decisions from raw observations or prove biological independence.
"""
import argparse
import csv
import hashlib
import io
import json
import re
import zipfile
from collections import Counter
from pathlib import Path

SOURCE_ZIP_SHA256 = '89f65a8ceb6c720797da9fb553293627d3b40283938f81132610017a7b143109'
SOURCE_PINS = {
    'pipeline/expected_membership/CURRENT_ANALYSIS_MEMBERSHIP_20261015.tsv': '52d7caf471ee96a9f3098327d6d6a349fd31b47ae12d570b67c51c993711262d',
    'pipeline/expected_membership/CURRENT_REPORTING_DATASET_LEDGER_20261015.tsv': 'c2fdcbd446ed3ecafd50687ee18435db69934b2a5498603d175285de76f35a51',
}
MEMBERS = 'ANALYSIS_MEMBERSHIP_20261015.tsv'
CONTRACTS = 'DATASET_REPORTING_20261015_public.tsv'
COUNTS = {'GSE128623':30,'GSE16710':27,'GSE205906':9,'GSE227494':12,
          'GSE277064':63,'GSE289929':12,'GSE299063':12}

def digest(raw):
    return hashlib.sha256(raw).hexdigest()

def read_table(raw):
    return list(csv.DictReader(io.StringIO(raw.decode('utf-8-sig')),delimiter='\t'))

def serialize(rows):
    f=io.StringIO(newline='')
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
    w.writeheader();w.writerows(rows)
    return f.getvalue().encode()

def check(rows, contracts):
    assert len(rows)==165
    assert len({(r['dataset_id'],r['sample_id']) for r in rows})==165
    assert Counter(r['dataset_id'] for r in rows)==Counter(COUNTS)
    branches=Counter(b for r in rows for b in r['membership_branch_ids'].split(';'))
    expected={}
    for c,n in [('primary26',26),('all27',27)]:
        for a in ['current','original']: expected[f'GSE16710_{c}_{a}']=n
    for p in ['raw','adapter']:
        for c,n in [('all30',30),('predefined29',29)]: expected[f'GSE128623_{p}_{c}']=n
    for c in ['v2a','astrocyte','microglia']: expected[f'GSE277064_natural_{c}']=15
    for c in ['v2a','astrocyte']: expected[f'GSE277064_rotarod_{c}']=12
    for d in ['GSE289929','GSE227494','GSE299063']: expected[f'{d}_canonical12']=12
    for b in ['September_20260912','October_Ensembl101_20261006']: expected[f'GSE205906_{b}']=9
    assert dict(branches)==expected
    excluded=[r['sample_id'] for r in rows if r['dataset_id']=='GSE16710' and 'GSE16710_primary26_current' not in r['membership_branch_ids'].split(';')]
    assert excluded==['GSM418682']
    excluded=[r['sample_id'] for r in rows if r['dataset_id']=='GSE128623' and 'GSE128623_raw_predefined29' not in r['membership_branch_ids'].split(';')]
    assert excluded==['GSM3681230']
    assert Counter(r['rotarod_role'] for r in rows if r['dataset_id']=='GSE277064')==Counter({'TRAINED':18,'SHARED_UNTRAINED_D14_COMPARATOR':6,'NOT_A_TRAINING_BRANCH':39})
    assert len(contracts)==9 and len({(r['dataset_id'],r['reporting_branch']) for r in contracts})==9
    assert len({r['dataset_id'] for r in contracts})==8
    for r in contracts:
        if r['dataset_id'] in COUNTS: assert int(r['submitted_member_n_20261015'])==COUNTS[r['dataset_id']]
        else: assert r['dataset_id']=='GSE85871' and r['submitted_member_n_20261015']==''
    return {'scope':'DATED_METADATA_ONLY_NOT_FULL_ANALYSIS_REPLAY','submitted_identifiers':165,
            'accessions_with_membership':7,'reporting_contracts':9,'reporting_accessions':8,
            'branch_counts':dict(sorted(branches.items())),
            'GSE85871_selected_membership':'UNKNOWN_NOT_ZERO','models_run':0,
            'biological_independence_proven':False}

def export(source_zip, output):
    assert not output.exists(), 'Use a fresh metadata directory'
    assert digest(source_zip.read_bytes())==SOURCE_ZIP_SHA256
    with zipfile.ZipFile(source_zip) as z:
        raw={n:z.read(n) for n in SOURCE_PINS}
        assert all(digest(raw[n])==h for n,h in SOURCE_PINS.items())
        code=z.read('code/python/build_analysis_membership_ledger_20261015.py')
    member_raw=raw[next(n for n in raw if 'ANALYSIS_MEMBERSHIP' in n)]
    original=read_table(raw[next(n for n in raw if 'DATASET_LEDGER' in n)])
    # Remove only the storage-location field; preserve all scientific values.
    contracts=[{k:v for k,v in row.items() if k!='current_large_data_location'} for row in original]
    contract_raw=serialize(contracts)
    check(read_table(member_raw),contracts)
    files={MEMBERS:member_raw,CONTRACTS:contract_raw,'historical_builder.py':code}
    for data in files.values():
        assert not re.search(r'(?<![A-Za-z])[A-Za-z]:[\\/]|/root/|autodl-container-|-----BEGIN .*PRIVATE KEY-----|github_pat_|ghp_',data.decode('utf-8-sig'))
    output.mkdir(parents=True)
    for name,data in files.items(): (output/name).write_bytes(data)
    provenance={'metadata_snapshot_date':'2026-10-15','source_archive_sha256':SOURCE_ZIP_SHA256,
        'original_sources':SOURCE_PINS,'transformation':{'membership':'BYTE_IDENTICAL',
        'contracts':'OMIT_ONLY_current_large_data_location_COLUMN','rows_dropped':0,
        'scientific_fields_changed':0,'historical_builder':'BYTE_IDENTICAL_CODE_ONLY_INPUT_BUNDLE_NOT_INCLUDED'},
        'files':{n:{'bytes':len(v),'sha256':digest(v)} for n,v in files.items()}}
    (output/'PROVENANCE.json').write_text(json.dumps(provenance,indent=2)+'\n')
    return verify(output)

def verify(directory):
    p=json.loads((directory/'PROVENANCE.json').read_text())
    for name,ident in p['files'].items():
        raw=(directory/name).read_bytes()
        assert len(raw)==ident['bytes'] and digest(raw)==ident['sha256'],name
    r=check(read_table((directory/MEMBERS).read_bytes()),read_table((directory/CONTRACTS).read_bytes()))
    r['file_identities_verified']=len(p['files'])
    return r

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    sp=p.add_subparsers(dest='action',required=True)
    a=sp.add_parser('export');a.add_argument('source_zip',type=Path);a.add_argument('output',type=Path)
    a=sp.add_parser('verify');a.add_argument('directory',type=Path)
    args=p.parse_args()
    result=export(args.source_zip,args.output) if args.action=='export' else verify(args.directory)
    print(json.dumps(result,indent=2))
