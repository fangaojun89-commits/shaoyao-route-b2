"""Read-only cohort identity construction; does not authorize checkpoint adoption."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('native_execution', HERE/'validate_PRJNA886048_native_execution_20261028.py')
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)


def read(p):
    return json.loads(p.read_text(encoding='utf-8-sig'))


def identity(p):
    return {'bytes':p.stat().st_size,'sha256':gate.sha(p)}


def build(package, inputs):
    if gate.sha(package/'provenance/EXTERNAL_INPUTS.json')!='236e7a1bf6cfc8d49a9234816ebceb15b840bd5448fedbaebee9f33ad40cbac1':
        raise ValueError('Historical input contract identity differs')
    contract = read(package/'provenance/EXTERNAL_INPUTS.json')
    if set(inputs) != {'native_root','tx2gene','mapped_edges','reference_terms'}:
        raise ValueError('Exactly four explicit inputs required')
    files=[]
    for item in contract['files']:
        path=inputs[item['role']]
        if item.get('relative'):
            path=path/item['relative']
        observed=identity(path)
        if item['role']!='native_root' and observed!={k:item[k] for k in ('bytes','sha256')}:
            raise ValueError('Frozen annotation/resource identity differs')
        files.append({'role':item['role'],**({'relative':item['relative']} if item.get('relative') else {}),**observed})
    evidence={}
    for n in range(46,58):
        run=f'SRR217758{n}'
        native=inputs['native_root']/f'linux_full_{run}_01'
        receipt=read(native/'receipt.json')
        if [x['stage'] for x in receipt['commands']]==gate.STAGES:
            evidence[run]=gate.validate(native,inputs['tx2gene'],package/'data/raw_ledger/ENA_quantification_manifest.tsv',run)
        else:
            if n not in (46,47,54):
                raise ValueError('Unexpected noncanonical execution')
            prefix=f'linux_full_{run}_01/'
            expected=[x for x in contract['files'] if x['role']=='native_root' and x.get('relative','').startswith(prefix)]
            if len(expected)!=5:
                raise ValueError('Incomplete historical identity contract')
            for item in expected:
                actual=next(x for x in files if x['role']=='native_root' and x.get('relative')==item['relative'])
                if any(actual[k]!=item[k] for k in ('bytes','sha256')):
                    raise ValueError('Historical variant changed; independent review required')
            evidence[run]={'status':'EXACT_HISTORICAL_VARIANT_REUSED',
                           'content':gate.validate_content(native,inputs['tx2gene'])}
    # Reject files that changed during the record/content check.
    for item in files:
        path=inputs[item['role']]
        if item.get('relative'):
            path=path/item['relative']
        if identity(path)!={k:item[k] for k in ('bytes','sha256')}:
            raise ValueError('Input changed during cohort validation')
    payload={'schema':1,'files':files,'historical_contract_sha256':gate.sha(package/'provenance/EXTERNAL_INPUTS.json'),
             'execution_gate_sha256':gate.sha(HERE/'validate_PRJNA886048_native_execution_20261028.py'),
             'content_gate_sha256':gate.sha(HERE/'validate_PRJNA886048_native_content_20261028.py')}
    binding=hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return {'status':'PASS_COHORT_INPUT_IDENTITY_ONLY','binding':binding,'payload':payload,'libraries':evidence,
            'historical_checkpoint_adoption_authorized':False,'analysis_launched':False,
            'remote_execution_attested':False}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--package',required=True,type=Path)
    for name in ('native-root','tx2gene','mapped-edges','reference-terms'):
        p.add_argument('--'+name,required=True,type=Path)
    a=p.parse_args()
    inputs={k:getattr(a,k).resolve(strict=True) for k in ('native_root','tx2gene','mapped_edges','reference_terms')}
    print(json.dumps(build(a.package.resolve(strict=True),inputs),indent=2))
