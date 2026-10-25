"""Rebuild the seven-table reading package from public retained results; standard library only."""
import argparse
import csv
import hashlib
import io
import json
import platform
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REL = Path('supplement/late_results')
NAMES = {
    'A1': 'A1_G128_all_families.tsv', 'A2': 'A2_G128_family_overview.tsv',
    'A3': 'A3_G128_gene_counts.tsv', 'A4': 'A4_G128_chronic_programs.tsv',
    'A5': 'A5_G128_day28_reference_spaces.tsv', 'A6': 'A6_JGW_Linux_all_families.tsv',
    'A7': 'A7_G205_matched_alignment.tsv',
}

def require(ok, message):
    if not ok:
        raise ValueError(message)

def digest(data):
    return hashlib.sha256(data).hexdigest()

def table(rows, columns=None):
    columns = columns or list(rows[0])
    def escape(value):
        return str(value).replace('|', '\\|').replace('\n', ' ')
    return ('| ' + ' | '.join(columns) + ' |\n|' + '|'.join(['---'] * len(columns)) + '|\n'
            + ''.join('| ' + ' | '.join(escape(row[c]) for c in columns) + ' |\n' for row in rows)).rstrip()

def render(data):
    """Descriptive aggregation and formatting only; no P/q computation."""
    groups = {}
    for row in data['A6']:
        groups.setdefault((row['accession'], row['space']), []).append(row)
    cell_names = {'GSE289929': 'PC12', 'GSE227494': 'C2C12', 'GSE299063': 'HepG2'}
    overview = []
    for (acc, space), rows in groups.items():
        overview.append({'Cell': cell_names[acc], 'Reference': space, 'Families': len(rows),
                         'Method-feature records': sum(int(r['feature_n']) for r in rows),
                         'Unavailable point P': sum(int(r['conditional_P_unavailable_n']) for r in rows),
                         'Unavailable feature-assignment comparisons': sum(int(r['null_unavailable_feature_draws']) for r in rows)})
    comparisons = []
    def rounded(value):
        return value if value == 'NA' else format(float(value), '.4f')
    for r in data['A7']:
        comparisons.append({'Method': r['table'], 'Collection': r['collection'],
                            'Contrast': {'local_disease': 'SNI−sham', 'intervention': 'PLM−SNI'}[r['contrast']],
                            'Eligible HISAT2 / STAR': r['left_eligible_n'] + ' / ' + r['right_eligible_n'],
                            'q<0.05 HISAT2 / STAR': r['left_q_strict_lt_05_n'] + ' / ' + r['right_q_strict_lt_05_n'],
                            'Shared hits': r['q_hit_intersection_n'], 'Hit Jaccard': rounded(r['q_hit_Jaccard']),
                            'Paired continuous effects': r['paired_estimable_continuous_n'],
                            'Effect or NES Spearman rho': rounded(r['effect_or_NES_spearman'])})
    index = [{'ID': k, 'File': f'[{NAMES[k]}](tables/{NAMES[k]})', 'Data rows': len(data[k])} for k in NAMES]
    return {'INDEX': table(index), 'A2': table(data['A2']), 'A3': table(data['A3']),
            'JGW': table(overview), 'A7': table(comparisons)}

def reconcile(data):
    """Compare accepted compact views with the complete accepted G128 family table."""
    checks = []
    def check(location, actual, expected):
        require(str(actual) == str(expected), 'Retained summary mismatch: ' + location)
        checks.append({'location': location, 'value': str(actual), 'status': 'MATCH'})
    g = data['A1']
    lookup = {(r['cell_id'], r['reference_space'], r['method'], r['collection'], r['contrast']): r for r in g}
    require(len(lookup) == 448, 'Duplicate or incomplete G128 family coverage')
    for row in data['A2']:
        sub = [r for r in g if r['cell_id'] == row['cell_id'] and r['reference_space'] == 'six_group']
        values = {'sampled_family_n': len(sub),
                  'nonzero_empirical_q_family_n': sum(int(r['empirical_q_le_05_n']) > 0 for r in sub),
                  'zero_empirical_q_family_n': sum(int(r['empirical_q_le_05_n']) == 0 for r in sub),
                  'feature_family_discovery_events_NOT_UNIQUE_FEATURES': sum(int(r['empirical_q_le_05_n']) for r in sub),
                  'unavailable_feature_family_events': sum(int(r['empirical_unavailable_n']) for r in sub)}
        for col, val in values.items():
            check('A2/' + row['cell_id'] + '/' + col, row[col], val)
    for row in data['A3']:
        for contrast in ('chronic_combined', 'D28_vs_control', 'D42_vs_control'):
            src = lookup[(row['cell_id'], 'six_group', row['method'], 'GENE', contrast)]
            for dest, original in ((contrast, 'empirical_q_le_05_n'), ('tested_genes', 'n'), ('reference_allocation_n', 'required_allocation_n')):
                check('A3/' + '/'.join((row['cell_id'], row['method'], contrast, dest)), row[dest], src[original])
    for row in data['A4']:
        src = lookup[(row['cell_id'], 'six_group', row['method'], row['collection'], 'chronic_combined')]
        for dest, original in (('tested_features', 'n'), ('empirical_q_le_05_n', 'empirical_q_le_05_n'), ('native_fraction_allocations_any_q_le_05', 'native_fraction_allocations_any_q_le_05'), ('reference_allocation_n', 'required_allocation_n')):
            check('A4/' + '/'.join((row['cell_id'], row['method'], row['collection'], dest)), row[dest], src[original])
    for row in data['A5']:
        for prefix, space in (('six_group', 'six_group'), ('day_control_exact', 'D28_vs_control_exact')):
            src = lookup[(row['cell_id'], space, row['method'], row['collection'], 'D28_vs_control')]
            for dest, original in ((prefix + '_empirical_q_le_05_n', 'empirical_q_le_05_n'), (prefix + '_allocation_n', 'required_allocation_n'), ('tested_features', 'n')):
                check('A5/' + '/'.join((row['cell_id'], row['method'], row['collection'], dest)), row[dest], src[original])
    require(len(checks) == 440, 'Unexpected reconciliation scope')
    return checks

def run(output):
    require(not output.exists(), 'Use a new output directory; existing outputs are preserved')
    release = json.loads((ROOT / 'MANIFEST.json').read_text(encoding='utf-8'))
    # Verify only this command's explicit inputs, not unrelated scientific stages.
    def load(relative):
        raw = (ROOT / relative).read_bytes()
        pin = release['files'][str(relative).replace('\\', '/')]
        require(len(raw) == pin['bytes'] and digest(raw) == pin['sha256'], 'Release input identity differs: ' + str(relative))
        return raw
    prov_raw = load(REL / 'PROVENANCE.json')
    prov = json.loads(prov_raw)
    template_raw = load(REL / 'reading_template.md')
    records = prov['tables']
    require([r['table_id'] for r in records] == list(NAMES), 'Table role/order mismatch')
    data, payloads, identities = {}, {}, []
    for record in records:
        k = record['table_id']
        raw = load(Path(record['public_path']))
        require(len(raw) == record['bytes'] and digest(raw) == record['sha256'], 'Source pin differs: ' + k)
        data[k] = list(csv.DictReader(io.StringIO(raw.decode('utf-8-sig')), delimiter='\t'))
        require(len(data[k]) == record['rows'], 'Row coverage differs: ' + k)
        name = 'tables/' + NAMES[k]
        payloads[name] = raw
        identities.append(dict(table_id=k, source_public_path=record['public_path'], output_path=name,
                               bytes=len(raw), sha256=digest(raw), rows=len(data[k]), transformation='none; exact byte copy'))
    require(sum(len(v) for v in data.values()) == 669, 'Table row coverage differs')
    checks = reconcile(data)
    text = template_raw.decode('utf-8')
    for key, content in render(data).items():
        token = '@@' + key + '@@'
        require(text.count(token) == 1, 'Template slot mismatch: ' + key)
        text = text.replace(token, content)
    require('@@' not in text, 'Unfilled template slot')
    payloads['SUPPLEMENTARY_RESULTS.md'] = text.encode('utf-8')
    payloads['SOURCE_TABLES.json'] = (json.dumps({'scope': 'Accepted result byte copies; no scientific recomputation', 'tables': identities}, indent=2) + '\n').encode()
    buf = io.StringIO(newline='')
    w = csv.DictWriter(buf, fieldnames=('location', 'value', 'status'), delimiter='\t', lineterminator='\n')
    w.writeheader(); w.writerows(checks)
    payloads['SUMMARY_RECONCILIATION.tsv'] = buf.getvalue().encode()
    output.mkdir(parents=True, exist_ok=False)
    for name, raw in payloads.items():
        dest = output / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        with dest.open('xb') as f:
            f.write(raw)
    receipt = {'status': 'PASS_PUBLIC_SAVED_TABLE_REPORTING', 'release_version': release['version'],
               'python': platform.python_version(), 'platform': platform.system(),
               'implementation': platform.python_implementation(), 'argv': sys.argv,
               'source_sha256': digest(Path(__file__).read_bytes()),
               'template_sha256': digest(template_raw), 'provenance_sha256': digest(prov_raw),
               'tables': 7, 'table_rows': 669, 'G128_summary_field_checks': len(checks),
               'outputs': {n: {'bytes': len(v), 'sha256': digest(v)} for n, v in payloads.items()},
               'models_run': 0, 'P_q_recalculated': False, 'random_draws': 0,
               'scope': 'Retained result copies, summary reconciliation and display generation; not upstream scientific reconstruction'}
    (output / 'RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    return receipt

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.output), indent=2))
