"""Check the published frozen JGW matrices, annotations and sample order without fitting."""
import csv
from decimal import Decimal
import gzip
import hashlib
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def main():
    base = REPO / 'data/jgw03_inputs'
    provenance = json.loads((base / 'PROVENANCE.json').read_text(encoding='utf-8'))
    for record in provenance['files']:
        path = Path(record['public_path'])
        if path.is_absolute() or '..' in path.parts:
            raise ValueError('Invalid public input path')
        raw = (REPO / path).read_bytes()
        if (len(raw), hashlib.sha256(raw).hexdigest()) != (record['bytes'], record['sha256']):
            raise ValueError('Input identity mismatch: ' + str(path))
    results = []
    for expected in provenance['datasets']:
        directory = base / expected['dataset']
        with (directory / 'locked_samples.tsv').open(encoding='utf-8', newline='') as handle:
            samples = list(csv.DictReader(handle, delimiter='\t'))
        with (directory / 'gene_annotation.tsv').open(encoding='utf-8', newline='') as handle:
            annotations = list(csv.DictReader(handle, delimiter='\t'))
        order = [r['sample_id'] for r in samples]
        if order != expected['sample_order'] or len(set(order)) != 12:
            raise ValueError('Sample order/membership mismatch')
        groups = {}
        for row in samples:
            groups[row['dose_group']] = groups.get(row['dose_group'], 0) + 1
        if groups != expected['dose_groups'] or set(groups.values()) != {3}:
            raise ValueError('Dose-group membership mismatch')
        genes, fractional, zeros = [], 0, 0
        with gzip.open(directory / 'locked_expected_counts.tsv.gz', 'rt', encoding='utf-8', newline='') as handle:
            rows = csv.reader(handle, delimiter='\t')
            if next(rows) != ['raw_gene_id'] + order:
                raise ValueError('Count header mismatch')
            for row in rows:
                if len(row) != 13:
                    raise ValueError('Malformed count row')
                genes.append(row[0])
                for value in row[1:]:
                    number = Decimal(value)
                    if not number.is_finite() or number < 0:
                        raise ValueError('Invalid expected count')
                    fractional += int(number != number.to_integral_value())
                    zeros += int(number == 0)
        if len(genes) != len(set(genes)) or len(genes) != expected['gene_rows']:
            raise ValueError('Gene membership changed')
        if genes != [r['raw_gene_id'] for r in annotations]:
            raise ValueError('Annotation row order mismatch')
        if (fractional, zeros) != (expected['noninteger_expected_count_cells'], expected['zero_count_cells']):
            raise ValueError('Count content changed')
        results.append({'dataset': expected['dataset'], 'genes': len(genes), 'libraries': 12,
                        'noninteger_expected_count_cells': fractional})
    lineage = json.loads((REPO / 'data/derived/SAVED_OUTPUT_LINEAGE.json').read_text(encoding='utf-8'))
    for key, name in [('G128','G128_FAMILY_SUMMARY.tsv'), ('JGW','JGW_FAMILY_SUMMARY.tsv')]:
        raw = (REPO / 'data/derived' / name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != lineage[key]['published_table_sha256']:
            raise ValueError('Accepted aggregate table changed')
    print(json.dumps({'status': 'FROZEN_INPUT_BYTES_CONTENT_AND_SAVED_TABLE_BINDINGS_PASS',
                      'input_files': len(provenance['files']), 'datasets': results,
                      'RDS_label_contents_checked': False, 'RDS_label_bytes_checked': True,
                      'next_label_check': 'Rscript --vanilla data/jgw03_inputs/verify_labels.R data/jgw03_inputs',
                      'models_run': 0, 'new_random_allocations': 0}, indent=2))


if __name__ == '__main__':
    main()
