"""Reuse 72 files and select 54 additional archived gene tables for new single-herb contrasts."""
from pathlib import Path
import csv,tarfile,shutil,hashlib,json,collections
base=Path('R:/Shaoyao_storage/data/GSE298414_preparation')
out=base/'single_herb_20261025_01'
out.mkdir(exist_ok=False)
extra=out/'additional_rsem_54';extra.mkdir()
source=base/'inputs_20261025/sample_ledger_supplement_verified.tsv'
rows=list(csv.DictReader(source.open(encoding='utf-8'),delimiter='\t'))
samples=[r for r in rows if r['recipe']]
assert len(samples)==126 and len({r['sample_id'] for r in samples})==126
assert len({r['group'] for r in samples})==42 and set(collections.Counter(r['group'] for r in samples).values())=={3}
assert all(r['batch']=='A' for r in samples)
members=[];reuse=0
with tarfile.open(base/'GSE298414_RAW.tar') as tar:
 for row in samples:
  member=row['member'];assert Path(member).name==member
  previous=base/'selected_rsem_72'/member
  if previous.exists():p=previous;reuse+=1
  else:
   p=extra/member
   with tar.extractfile(member) as f,p.open('wb') as g:shutil.copyfileobj(f,g)
  row['rsem_path']=p.as_posix()
  members.append({'sample_id':row['sample_id'],'member':member,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'reused':previous.exists()})
assert reuse==72
def write(name,rows):
 with (out/name).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
write('samples.tsv',samples)
contrasts=[]
for solvent in ('W','E'):
 for mix in ('02','04','06','07','08'):
  for level in range(1,4):
   for single,label in [('01','mixture_minus_PR'),('05','mixture_minus_GR')]:
    pos=f'JG{solvent}{mix}_L{level}';neg=f'JG{solvent}{single}_L{level}'
    contrasts.append({'contrast_id':pos+'_minus_'+neg,'contrast_class':label,'positive_group':pos,'negative_group':neg,'ordinal_level':level,'solvent':solvent,'mixture_code':mix,'comparator_code':single})
old=list(csv.DictReader((base/'analysis_20261025_03/program_and_stratified_extension_01/contrast_inventory.tsv').open(),delimiter='\t'))
assert len(contrasts)==60 and not ({r['contrast_id'] for r in old}&{r['contrast_id'] for r in contrasts})
write('contrast_inventory.tsv',contrasts)
(out/'input_receipt.json').write_text(json.dumps({'status':'INPUTS_SELECTED','libraries':126,'conditions':42,'reused_rsem_files':72,'newly_unpacked_rsem_files':54,'new_contrasts':60,'prior_contrast_overlap':0,'source_ledger_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'members':members},indent=2)+'\n',encoding='utf-8')
print('Prepared 126 libraries and 60 new contrasts; reused 72 compressed inputs and unpacked 54.')
