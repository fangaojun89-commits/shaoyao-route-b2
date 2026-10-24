"""Restore accepted counts and call original R input readers; no model fitting."""
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import re
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MISSING = '__ROUTE_B2_MISSING_GENE_NAME__'
ATTR = re.compile(r'([A-Za-z_][A-Za-z0-9_]*)[ \t]+"([^"]*)";')
GTF_COLUMNS = ['gene_id', 'gene_name', 'gene_biotype', 'source_gene_record_ordinal',
               'source_noncomment_record_ordinal', 'seqname', 'source', 'feature',
               'start', 'end', 'score', 'strand', 'frame', 'raw_gtf_attributes']


def identity(path):
    h = hashlib.sha256()
    size = 0
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            size += len(block)
            h.update(block)
    return {'bytes': size, 'sha256': h.hexdigest()}


def gene_row(fields, ordinal, noncomment):
    attributes = fields[8]
    if ATTR.sub('', attributes).strip():
        raise ValueError('Unparsed GTF attribute text')
    pairs = ATTR.findall(attributes)
    values = dict(pairs)
    if len(values) != len(pairs) or not {'gene_id', 'gene_biotype'} <= values.keys():
        raise ValueError('Repeated or missing GTF gene attributes')
    if not re.fullmatch(r'ENSRNOG[0-9]+', values['gene_id']):
        raise ValueError('Unexpected gene identifier')
    if not values['gene_biotype'] or values['gene_biotype'] != values['gene_biotype'].strip():
        raise ValueError('Invalid gene biotype')
    if MISSING in values.values():
        raise ValueError('Missing-value marker collision')
    return [values['gene_id'], values.get('gene_name', MISSING), values['gene_biotype'],
            str(ordinal), str(noncomment), *fields[:8], attributes]


def extract_attributes(source, target, expected_plain):
    digest = hashlib.sha256()
    size = 0
    noncomment = 0
    genes = set()
    with gzip.open(source, 'rb') as src, target.open('x', encoding='utf-8', newline='') as dst:
        writer = csv.writer(dst, delimiter='\t', lineterminator='\n')
        writer.writerow(GTF_COLUMNS)
        for raw in src:
            digest.update(raw)
            size += len(raw)
            line = raw.decode('utf-8').rstrip('\n').removesuffix('\r')
            if line.startswith('#'):
                continue
            noncomment += 1
            fields = line.split('\t')
            if len(fields) != 9:
                raise ValueError('GTF record must contain nine fields')
            if fields[2] != 'gene':
                continue
            row = gene_row(fields, len(genes) + 1, noncomment)
            if row[0] in genes:
                raise ValueError('Duplicate GTF gene record')
            genes.add(row[0])
            writer.writerow(row)
    observed = {'bytes': size, 'sha256': digest.hexdigest()}
    if observed != expected_plain or len(genes) != 32883:
        raise ValueError('Complete GTF identity or gene count mismatch')
    return {'plain_GTF': observed, 'gene_records': len(genes), 'noncomment_records': noncomment,
            'attribute_table': identity(target), 'rule': 'All original gene attributes and source record order preserved'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gtf-gz', type=Path, required=True)
    parser.add_argument('--rscript', required=True)
    parser.add_argument('--r-library', default='')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output must be a new directory')
    contract = json.loads((REPO / 'analysis/input_bridge/CONTRACT.json').read_text(encoding='utf-8'))
    if identity(args.gtf_gz) != contract['GTF_gzip']:
        raise ValueError('Expected the frozen Ensembl101 GTF gzip')
    for record in contract['sources']:
        if identity(REPO / record['path']) != {k: record[k] for k in ('bytes', 'sha256')}:
            raise ValueError('Original source identity mismatch: ' + record['path'])
    # Verify the versioned bridge and restored-count utility against the release manifest.
    manifest = json.loads((REPO / 'MANIFEST.json').read_text(encoding='utf-8'))
    for name in ('analysis/input_bridge/prepare.R', 'analysis/input_bridge/CONTRACT.json',
                 'pipeline/count_inputs.py', 'pipeline/prepare_count_inputs.py',
                 'data/count_inputs/PROVENANCE.json'):
        if identity(REPO / name) != manifest['files'][name]:
            raise ValueError('Bridge dependency identity mismatch: ' + name)
    spec = importlib.util.spec_from_file_location('accepted_counts', REPO / 'pipeline/count_inputs.py')
    counts = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(counts)
    payloads = counts.verified_payloads()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    for directory in ('original_inputs', 'derived', 'prepared', 'receipts'):
        (output / directory).mkdir()
    for relative, raw in payloads:
        target = output / 'original_inputs' / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('xb') as handle:
            handle.write(raw)
    extraction = extract_attributes(args.gtf_gz, output / 'derived/G205_gene_attributes.tsv', contract['GTF_plain'])
    command = [args.rscript, '--vanilla', str(REPO / 'analysis/input_bridge/prepare.R'),
               str(REPO), str(output), args.r_library, str(args.gtf_gz.resolve())]
    (output / 'receipts/COMMAND.json').write_text(json.dumps(command, indent=2) + '\n', encoding='utf-8')
    with (output / 'receipts/R.stdout.txt').open('wb') as stdout, (output / 'receipts/R.stderr.txt').open('wb') as stderr:
        result = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
    if result.returncode:
        raise RuntimeError('Input-reader process failed; inspect retained receipts/R.stderr.txt')
    summary = json.loads((output / 'receipts/INPUT_READER_QA.json').read_text(encoding='utf-8'))
    if summary['status'] != 'PASS_ORIGINAL_INPUT_READERS_WITHOUT_FITTING':
        raise ValueError('Unexpected input-reader outcome')
    # Completed outputs are new portable input objects, not historical fitted checkpoints.
    receipt = {'status': 'PASS_PORTABLE_INPUT_PREPARATION', 'R_exit_code': result.returncode,
               'original_files_restored': len(payloads), 'GTF_extraction': extraction,
               'source_contract': identity(REPO / 'analysis/input_bridge/CONTRACT.json'),
               'outputs': {p.relative_to(output).as_posix(): identity(p) for p in sorted((output / 'prepared').glob('*.rds'))},
               'scientific_model_runs': 0, 'whole_study_restoration_verified': False}
    (output / 'receipts/PREPARE_COMPLETE.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
