"""Offline standard-library validation; does not write source data."""
import csv, hashlib, io, json
from pathlib import Path

def need(ok, why):
    if not ok: raise ValueError(why)

root = Path(__file__).resolve().parent
manifest = json.loads((root/'MANIFEST.json').read_text())
expected = {x['path'] for x in manifest['files']}
actual = {p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file()} - {'MANIFEST.json'}
need(actual == expected, 'Missing or unlisted payload')
for item in manifest['files']:
    p = root/item['path']
    need(p.resolve().is_relative_to(root.resolve()) and not p.is_symlink(), 'Invalid payload path')
    b = p.read_bytes()
    need(len(b) == item['bytes'] and hashlib.sha256(b).hexdigest() == item['sha256'], 'Payload identity mismatch: '+item['path'])
meta = json.loads((root/'PROVENANCE.json').read_text())
rows = list(csv.DictReader(io.StringIO((root/'download_summary.portable.tsv').read_text()), delimiter='\t'))
results = []
for row in rows:
    label = row['label']
    pages = sorted((root/'pages'/label).glob('page_*.tsv'))
    need([p.name for p in pages] == [f'page_{i:03d}.tsv' for i in range(1,int(row['pages'])+1)], 'Page sequence mismatch')
    lines, slots, blank, ids = [], 0, 0, []
    for p in pages:
        part = p.read_text(encoding='utf-8-sig').splitlines()
        need(bool(part), 'Empty page')
        if not lines: lines.append(part[0])
        need(part[0] == lines[0], 'Column mismatch')
        slots += len(part)-1
        for s in part[1:]:
            if not s.strip(): blank += 1; continue
            cells = s.split('\t')
            need(len(cells) == 7 and cells[4] == row['tax_id'] and cells[6] == 'PROTEIN_CODING', 'Schema or species mismatch')
            lines.append(s); ids.append(cells[0])
    need(len(ids) == len(set(ids)), 'Duplicate gene')
    need(slots == int(row['response_slots']) and blank == int(row['empty_gene_reports']) and len(ids) == int(row['parseable_gene_records']), 'Count mismatch')
    b = (root/row['output_file']).read_bytes()
    need(('\r\n'.join(lines)+'\r\n').encode() == b, 'Combined file differs from page reconstruction')
    pin = next(x for x in meta['combined_files'] if x['label'] == label)
    need(hashlib.sha256(b).hexdigest() == pin['sha256'], 'Historical combined hash mismatch')
    results.append({'taxon':row['tax_id'], 'pages':len(pages), 'records':len(ids), 'empty_reports':blank})
print(json.dumps({'status':'PASS_OFFLINE_SNAPSHOT_RECONSTRUCTION', 'taxa':results, 'payloads':len(expected)}, indent=2))
