"""Parse the exact Ensembl107 GRCm39 GTF into read-only-source annotation dictionaries."""
from pathlib import Path
import argparse, csv, gzip, hashlib, json, re


def digest(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--gtf', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    source = a.gtf.resolve(strict=True)
    out = a.out.resolve()
    expected = '15af16b7bef88c4309a82d49579fdd967f6ac8769b7f0bf657bbf7ee73b779e3'
    if out.exists() or source.is_relative_to(out):
        raise ValueError('Output must be a new directory outside the source')
    if digest(source) != expected:
        raise ValueError('GTF differs from the frozen Ensembl107 source')
    genes, transcripts, gids, tids = [], [], set(), set()
    with gzip.open(source, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            cells = line.rstrip('\n').split('\t')
            if len(cells) != 9:
                raise ValueError('Invalid GTF row')
            if cells[2] not in ('gene', 'transcript'):
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)";', cells[8]))
            if cells[2] == 'gene':
                gid = attrs['gene_id']
                if gid in gids:
                    raise ValueError('Duplicate gene')
                gids.add(gid)
                genes.append(dict(gene_id=gid, gene_version=attrs['gene_version'],
                                  gene_name=attrs.get('gene_name', ''), gene_biotype=attrs['gene_biotype'], seqname=cells[0]))
            else:
                tid = attrs['transcript_id']
                if tid in tids:
                    raise ValueError('Duplicate transcript')
                tids.add(tid)
                transcripts.append(dict(transcript_id=tid, transcript_version=attrs['transcript_version'],
                                        gene_id=attrs['gene_id'], gene_name=attrs.get('gene_name', ''),
                                        transcript_biotype=attrs['transcript_biotype']))
    if len(genes) != 56748 or len(transcripts) != 149478 or not {r['gene_id'] for r in transcripts} <= gids:
        raise ValueError('Frozen annotation counts or gene linkage mismatch')
    out.mkdir(parents=True, exist_ok=False)
    outputs = {}
    for name, rows in [('genes.tsv', genes), ('transcript_to_gene.tsv', transcripts)]:
        p = out / name
        with p.open('w', encoding='utf-8', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
            w.writeheader(); w.writerows(rows)
        outputs[name] = {'sha256': digest(p), 'bytes': p.stat().st_size, 'rows': len(rows)}
    result = {'status': 'PARSED_FROZEN_GTF', 'source_sha256': expected, 'outputs': outputs,
              'source_modified': False, 'expression_data_read': False, 'models_fitted': 0}
    (out/'parser_receipt.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
