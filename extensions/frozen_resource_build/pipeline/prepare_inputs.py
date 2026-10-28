"""Assemble eight frozen inputs from explicit locations, with no network or models."""
import argparse, hashlib, json, shutil
from pathlib import Path

CONTRACT_SHA = '369dbf9b8a1d5c2f9ee8b28f04d81a89961de80ad3ffba6067acf4658f15d5e4'

def ident(p):
    b=p.read_bytes(); return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['contract','ncbi-snapshot','quickgo','reactome-gmt','matrisome-human','matrisome-mouse','output']:
        parser.add_argument('--'+name,type=Path,required=True)
    a=parser.parse_args()
    if ident(a.contract)['sha256']!=CONTRACT_SHA: raise ValueError('Unknown resource contract')
    contract=json.loads(a.contract.read_text())
    sources={
        'download_summary.tsv':a.ncbi_snapshot/'download_summary.portable.tsv',
        'ReactomePathways.gmt':a.reactome_gmt,
        'QuickGO_SynGO_assignedBy_2026-08-17.tsv':a.quickgo,
        'Hs_Matrisome_Masterlist_rev2014.csv':a.matrisome_human,
        'Mm_Matrisome_Masterlist_rev2014.csv':a.matrisome_mouse,
    }
    for label in ['human','mouse','rat']:
        name=f'NCBI_{label}_protein_coding_gene_groups_2026-08-18.tsv'
        sources[name]=a.ncbi_snapshot/name
    out=a.output.resolve()
    if out.exists(): raise ValueError('Output must be a new directory')
    staged=[]
    for item in contract['inputs']:
        relative=Path(item['path'])
        if relative.is_absolute() or '..' in relative.parts or ':' in item['path'] or '\\' in item['path']:
            raise ValueError('Unsafe input layout')
        src=sources[relative.name].resolve(strict=True)
        if src.is_relative_to(out) or out.is_relative_to(src.parent):
            raise ValueError('Output overlaps source area')
        if ident(src)!=item['identity']: raise ValueError('Frozen input mismatch: '+relative.name)
        staged.append((src,relative,item['identity']))
    if len(staged)!=8: raise ValueError('Unexpected input count')
    out.mkdir(parents=True,exist_ok=False)
    for src,relative,pin in staged:
        target=out/relative; target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(src,target)
        if ident(target)!=pin or ident(src)!=pin: raise ValueError('Copy or source changed')
    receipt={'status':'PASS_EIGHT_PINNED_INPUTS_STAGED','contract_sha256':CONTRACT_SHA,'inputs':[{'path':p.as_posix(),'identity':i} for _,p,i in staged],'source_files_unchanged':True,'models_run':0,'network_requests':0}
    (out/'INPUT_PREPARATION.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))

if __name__=='__main__': main()
