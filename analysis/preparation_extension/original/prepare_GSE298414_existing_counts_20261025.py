"""Prepare deposited RSEM estimates and source-labelled recipe metadata; no model."""
import collections
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import tarfile

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BASE = Path('R:/Shaoyao_storage/data/GSE298414_preparation')
SOURCE = Path('R:/Shaoyao_storage/literature/GSE298414_preparation_20261025')
OUT = BASE / 'inputs_20261025'


def tsv(path, rows, fields):
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def main():
    meta = json.loads((SOURCE / 'sample_metadata.json').read_text(encoding='utf-8'))
    assert len(meta) == 141 and len({r['GSM'] for r in meta}) == 141
    ledger = []
    recipes = {1: ('3:0', 'single_herb_PR'), 2: ('2:1', 'CEM'), 4: ('1:2', 'CEM'),
               5: ('0:3', 'single_herb_GR'), 6: ('2:1', 'IEM'), 7: ('1:1', 'IEM'), 8: ('1:2', 'IEM')}
    for r in meta:
        title = r['!Sample_title'][0]
        treatment = next(x.removeprefix('treatment: ') for x in r['!Sample_characteristics_ch1'] if x.startswith('treatment: '))
        code = re.search(r'JG[WE]\d{2}', title)
        rep = int(re.search(r'_rep(\d+)$', title)[1])
        assert rep in (1, 2, 3)
        if code:
            code = code[0]
            assert code in treatment
            dose = int(re.search(r'\((\d+) ug/mL dose\)', title)[1])
            ordinal = {20: 1, 100: 2, 500: 3}[dose]
            ratio, extraction = recipes[int(code[-2:])]
            solvent = 'water' if code[2] == 'W' else '70_percent_ethanol'
            group = f'{code}_L{ordinal}'
        else:
            code = ''; dose = ''; ordinal = ''; ratio = ''; extraction = ''; solvent = ''
            group = title.rsplit('_rep', 1)[0]
        relations = r['!Sample_relation']
        biosample = next(x.rsplit('/', 1)[-1] for x in relations if x.startswith('BioSample:'))
        srx = next(x.rsplit('=', 1)[-1] for x in relations if x.startswith('SRA:'))
        file_url = r['!Sample_supplementary_file_1'][0]
        ledger.append(dict(sample_id=r['GSM'], biosample=biosample, sra_experiment=srx,
          title=title, source_treatment=treatment, group=group, recipe=code,
          nominal_GEO_dose_ug_ml=dose, ordinal_level=ordinal, replicate_label=rep,
          ratio_PR_GR=ratio, extraction=extraction, solvent=solvent,
          dose_status='ABSOLUTE_DOSE_CONFLICT' if code else 'NOT_APPLICABLE',
          batch='NOT_RESOLVED', control_status='UNTREATED_NOT_CONFIRMED_MATCHED_VEHICLE' if group=='Non-treat(N.T)' else 'NOT_A_VEHICLE',
          member=file_url.rsplit('/', 1)[-1], source_url=file_url))
    assert len({r['biosample'] for r in ledger}) == len(ledger)
    assert set(collections.Counter(r['group'] for r in ledger).values()) == {3}
    assert len({r['group'] for r in ledger}) == 47
    counts, qc, fingerprints = [], [], {}
    gene_ids = None
    with tarfile.open(BASE / 'GSE298414_RAW.tar', 'r:') as archive:
        members = archive.getmembers()
        assert len(members) == len(ledger) and all(m.isfile() for m in members)
        assert set(m.name for m in members) == {r['member'] for r in ledger}
        for r in ledger:
            raw = gzip.decompress(archive.extractfile(r['member']).read())
            data = list(csv.DictReader(io.StringIO(raw.decode()), delimiter='\t'))
            ids = [x['gene_id'] for x in data]
            assert len(set(ids)) == len(ids)
            if gene_ids is None: gene_ids = ids
            else: assert ids == gene_ids
            x = np.array([float(d['expected_count']) for d in data], dtype=np.float64)
            tpm = np.array([float(d['TPM']) for d in data])
            assert np.isfinite(x).all() and (x >= 0).all() and x.sum() > 0
            assert np.isfinite(tpm).all() and (tpm >= 0).all()
            counts.append(x)
            fingerprint = hashlib.sha256(x.astype('<f8').tobytes()).hexdigest()
            fingerprints[r['sample_id']] = fingerprint
            qc.append(dict(sample_id=r['sample_id'], expected_count_sum=float(x.sum()),
                           nonzero_genes=int((x > 0).sum()), TPM_sum=float(tpm.sum()),
                           count_vector_sha256=fingerprint))
    assert len(set(fingerprints.values())) == len(fingerprints)
    canonical = ROOT / 'data/processed/jgw03_three_cell_lines/GSE289929'
    annotation = list(csv.DictReader((canonical / 'gene_annotation.tsv').open(encoding='utf-8-sig'), delimiter='\t'))
    assert gene_ids == [r['raw_gene_id'] for r in annotation]
    old = list(csv.DictReader((canonical / 'locked_samples.tsv').open(encoding='utf-8-sig'), delimiter='\t'))
    assert not ({r['sample_id'] for r in old} & {r['sample_id'] for r in ledger})
    assert not ({r['biosample'] for r in old} & {r['biosample'] for r in ledger})
    old_fingerprints = set()
    for p in (canonical / 'selected_rsem').glob('*.gz'):
        with gzip.open(p, 'rt') as f: data = list(csv.DictReader(f, delimiter='\t'))
        assert [r['gene_id'] for r in data] == gene_ids
        x = np.array([float(r['expected_count']) for r in data], dtype='<f8')
        old_fingerprints.add(hashlib.sha256(x.tobytes()).hexdigest())
    assert len(old_fingerprints) == 12 and not (old_fingerprints & set(fingerprints.values()))
    matrix = np.column_stack(counts)
    OUT.mkdir(exist_ok=False)
    tsv(OUT / 'sample_ledger.tsv', ledger, list(ledger[0]))
    tsv(OUT / 'library_summary.tsv', qc, list(qc[0]))
    with gzip.open(OUT / 'expected_counts.tsv.gz', 'wt', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['gene_id'] + [r['sample_id'] for r in ledger])
        for gid, row in zip(gene_ids, matrix): w.writerow([gid] + [format(v, '.2f') for v in row])
    manifest = [dict(path=str(p), bytes=p.stat().st_size, sha256=hashlib.sha256(p.read_bytes()).hexdigest())
                for p in [BASE/'GSE298414_RAW.tar',SOURCE/'GSE298414_family.soft.gz',
                          SOURCE/'sample_metadata.json',canonical/'gene_annotation.tsv']]
    result = dict(status='DEPOSITED_COUNTS_PREPARED_NOT_MODELED', samples=len(ledger), genes=len(gene_ids),
      groups=47, replicates_per_source_group=3, recipe_conditions=14, recipe_libraries=126,
      untreated_libraries=3, named_comparator_drug_libraries=12,
      exact_count_duplicates_within_source=0, exact_count_duplicates_with_canonical_PC12=0,
      annotation_reused=str(canonical/'gene_annotation.tsv'), input_manifest=manifest,
      remaining=['Resolve or explicitly account for sequencing-batch metadata',
                 'No absolute-dose inference while source conflict remains',
                 'Do not substitute untreated samples for a verified DMSO vehicle',
                 'Fix admissible within-accession contrasts before new fitting'],
      biological_models_run=0, completed_analysis_reruns=0)
    (OUT/'receipt.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('input_manifest',)},indent=2))


if __name__ == '__main__':main()
