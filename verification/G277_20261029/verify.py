"""Check distributed comparison evidence; no model fitting."""
import csv,hashlib,json
from pathlib import Path
p=Path(__file__).resolve().parent
m=json.loads((p/'EVIDENCE_MANIFEST.json').read_text(encoding='utf-8'))
for n,pin in m.items():
 b=(p/n).read_bytes()
 if {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}!=pin:raise ValueError('Evidence changed: '+n)
def rows(n):
 with (p/n).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
s=json.loads((p/'SUMMARY.json').read_text(encoding='utf-8'))
t=rows('TABLE_COMPONENT_COMPARISON.tsv');r=rows('RDS_COMPONENT_COMPARISON.tsv');v=rows('SELECTED_METADATA_AND_DESIGNS.tsv')
if (len({x['file'] for x in t}),len(t),sum(x['pass']=='FALSE' for x in t))!=(125,2789,41):raise ValueError('Table coverage')
if (len({x['file'] for x in r}),len(r),sum(x['pass']=='FALSE' for x in r))!=(5,41981,24):raise ValueError('RDS coverage')
if not all(x['same_missingness']=='TRUE' and x['exact_character_values']=='TRUE' for x in v):raise ValueError('Selected metadata/design mismatch')
g=json.loads((p/'GRAPHICS_COMPARISON.json').read_text());q=json.loads((p/'PROJECTION_COMPARISON.json').read_text())
if g['files']!=33 or not g['visual_review_complete'] or q['families']!=24:raise ValueError('Graphics/projection coverage')
if not all(x['new_source_to_view_all_cells_exact'] for x in q['records']):raise ValueError('Projection fidelity')
print(json.dumps({'status':'PASS_DISTRIBUTED_EVIDENCE_IDENTITIES_AND_COVERAGE','table_checks':len(t),'RDS_checks':len(r),'strict_failures_preserved':65,'scientific_model_runs':0}))
