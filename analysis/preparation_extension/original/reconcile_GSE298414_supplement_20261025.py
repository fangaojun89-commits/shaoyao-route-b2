"""Reconcile the deposited metadata to the original article's sample table."""
from pathlib import Path
import csv
import hashlib
import json
import openpyxl

BASE = Path('R:/Shaoyao_storage/literature/GSE298414_preparation_20261025')
DATA = Path('R:/Shaoyao_storage/data/GSE298414_preparation/inputs_20261025')


def save_tsv(path, rows):
    assert not path.exists()
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, rows[0].keys(), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def main():
    path = BASE / '41597_2026_6759_MOESM3_ESM.xlsx'
    assert path.stat().st_size == 27080  # Size in the article's XML metadata.
    book = openpyxl.load_workbook(path, read_only=True, data_only=True)
    raw = list(book['PC12'].values)
    source = [dict(zip(raw[1], x)) for x in raw[2:]]
    indexed = {x['Sample']: x for x in source}
    assert len(indexed) == len(source) == 162
    with (DATA / 'sample_ledger.tsv').open(encoding='utf-8') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    conflicts = []
    for x in rows:
        y = indexed[x['sample_id']]
        assert y['Dataset'] == 'GSE298414' and y['Title'] == x['title'] and y['batch'] == 'A'
        assert y['treatment'] == x['source_treatment']
        x['batch'] = 'A'
        x['batch_source'] = 'Published Table S3'
        x['preparation_class'] = x['extraction']
        x['published_extraction_label'] = y['Extraction method ']
        x['extraction_interpretation'] = 'NOT_APPLICABLE'
        if x['recipe']:
            assert y['Combination ratio(PR:GR)'].removeprefix('PR:GR=') == x['ratio_PR_GR']
            assert y['Solvent'] == ('Water' if x['solvent'] == 'water' else 'Ethanol')
            expected = 'IEM' if int(x['recipe'][-2:]) >= 6 else 'CEM'
            expected_phrase = 'extraction -> mix' if expected == 'IEM' else 'mix-> extraction'
            assert expected_phrase in x['source_treatment']
            x['extraction_interpretation'] = 'GEO_ARTICLE_CODE_AND_SUPPLEMENT_AGREE'
            if y['Extraction method '] != expected:
                x['extraction_interpretation'] = 'SUPPLEMENT_LABEL_CONFLICT_GEO_AND_ARTICLE_CODE_AGREE'
                conflicts.append(dict(sample_id=x['sample_id'], recipe=x['recipe'],
                  published_extraction_label=y['Extraction method '], interpreted_extraction=expected,
                  evidence='GEO treatment operation order plus article Tables 2/3 code definition',
                  author_confirmed_correction=False))
            x['extraction'] = expected
        x['supplement_source'] = path.name
    assert {x['sample_id'] for x in conflicts} == {'GSM9013914', 'GSM9013917'}
    save_tsv(DATA / 'sample_ledger_supplement_verified.tsv', rows)
    save_tsv(DATA / 'extraction_label_conflicts.tsv', conflicts)
    with (DATA / 'cem_iem_contrast_inventory.tsv').open() as f:
        contrasts = list(csv.DictReader(f, delimiter='\t'))
    for c in contrasts:
        c['status'] = 'ORDINAL_ONLY_WITH_RECORDED_LABEL_ADJUDICATION_NOT_MODELED'
    save_tsv(DATA / 'cem_iem_contrasts_supplement_verified.tsv', contrasts)
    receipt = dict(status='BATCH_VERIFIED_TWO_EXTRACTION_LABEL_CONFLICTS_RETAINED', supplement_path=str(path),
      sha256=hashlib.sha256(path.read_bytes()).hexdigest(), table_rows=162,
      GSE298414_rows=141, GSE289929_rows=21, GSE298414_batches=['A'], title_matches=141,
      treatment_matches=141, recipe_ratio_solvent_matches=126, extraction_label_matches=124,
      extraction_conflicts=conflicts,
      absolute_dose_conflict_remains=True, controls_not_borrowed=True, new_models_run=0,
      source='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12988177/supplementaryFiles')
    (DATA / 'supplement_reconciliation.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
