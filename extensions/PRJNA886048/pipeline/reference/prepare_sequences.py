"""Verify downloaded Ensembl 107 inputs and prepare their explicit target universe."""
from collections import Counter
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

import argparse
import sys
if not __debug__:
    raise SystemExit('Optimized Python is unsupported: validation assertions are required')
parser = argparse.ArgumentParser(description='Reconstruct the frozen Ensembl107 sequence universe without building an index')
parser.add_argument('--source', required=True, type=Path)
parser.add_argument('--annotation', required=True, type=Path)
parser.add_argument('--out', required=True, type=Path)
args = parser.parse_args()
SOURCE = args.source.resolve()
ANNOTATION = args.annotation.resolve()
ROOT = args.out.resolve()
if ROOT.exists() or ROOT == SOURCE or ROOT in SOURCE.parents or SOURCE in ROOT.parents or ROOT in ANNOTATION.parents:
    raise SystemExit('Output must be new and separate from all source paths')
PINS = {'Mus_musculus.GRCm39.cdna.all.fa.gz': '34499a727ff72bdddb2f32e5831b14fb5a9d20218cd79b8639dc8736b2c596c0', 'Mus_musculus.GRCm39.ncrna.fa.gz': '0472b7d912e84d838485ecf236670a68319a32fb4580062cd9768c7bce85d9ab', 'Mus_musculus.GRCm39.107.gtf.gz': '15af16b7bef88c4309a82d49579fdd967f6ac8769b7f0bf657bbf7ee73b779e3', 'cdna_CHECKSUMS.txt': 'e2bc0007cd6f880df526eee5548c6e986e4d91f898ec46f0a7967f69dd5e27b0', 'ncrna_CHECKSUMS.txt': '3d80a6c4784561b65834456322f76b9b72cc3e1376a6e8e8832cc73be6ce5329', 'gtf_CHECKSUMS.txt': '6572e8b44e10f1e43d774782d6fa1e191d0530d21693538c47af4c3b7febcb18'}
for name, expected_sha in PINS.items():
    if hashlib.sha256((SOURCE / name).read_bytes()).hexdigest() != expected_sha:
        raise SystemExit('Source SHA256 mismatch: ' + name)
if hashlib.sha256(ANNOTATION.read_bytes()).hexdigest() != '89176b6187573094bc4d3b7139290880870562c4ff44e9ff1b860fb3f77eb03d':
    raise SystemExit('Annotation SHA256 mismatch')
ROOT.mkdir(parents=True)


def write_tsv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

checks = []
for name, checksum in [
    ('Mus_musculus.GRCm39.cdna.all.fa.gz', 'cdna_CHECKSUMS.txt'),
    ('Mus_musculus.GRCm39.ncrna.fa.gz', 'ncrna_CHECKSUMS.txt'),
    ('Mus_musculus.GRCm39.107.gtf.gz', 'gtf_CHECKSUMS.txt'),
]:
    data = (SOURCE / name).read_bytes()
    total = 0
    for value in data:
        total = (((total >> 1) | ((total & 1) << 15)) + value) & 65535
    expected = next(line.split() for line in (SOURCE / checksum).read_text().splitlines() if line.split()[-1] == name)
    assert total == int(expected[0]) and math.ceil(len(data) / 1024) == int(expected[1]), name
    checks.append(dict(file=name, official_BSD_sum=total, official_blocks_1024=int(expected[1]),
                       checksum_matches=True, sha256=hashlib.sha256(data).hexdigest()))

with ANNOTATION.open(encoding='utf-8') as f:
    annotation = {row['transcript_id']: row for row in csv.DictReader(f, delimiter='\t')}
seen, source_counts = {}, {}
target = ROOT / 'transcripts_Ensembl107_GRCm39.fa'
temporary = target.with_suffix('.fa.part')
with temporary.open('w', encoding='ascii', newline='\n') as output:
    for filename in ['Mus_musculus.GRCm39.cdna.all.fa.gz', 'Mus_musculus.GRCm39.ncrna.fa.gz']:
        source_counts[filename] = 0
        with gzip.open(SOURCE / filename, 'rt', encoding='ascii') as f:
            for line in f:
                if line.startswith('>'):
                    fields = line[1:].split()
                    full_id = fields[0]
                    transcript = full_id.rsplit('.', 1)[0]
                    gene = next(x[5:] for x in fields if x.startswith('gene:')).rsplit('.', 1)[0]
                    assert transcript in annotation and transcript not in seen
                    expected = annotation[transcript]
                    assert gene == expected['gene_id']
                    assert full_id == transcript + '.' + expected['transcript_version']
                    seen[transcript] = dict(transcript_id_versioned=full_id, transcript_id=transcript,
                                            gene_id=gene, gene_name=expected['gene_name'], length=0)
                    source_counts[filename] += 1
                else:
                    seen[transcript]['length'] += len(line.strip())
                output.write(line)
temporary.replace(target)
missing = [annotation[k] for k in sorted(set(annotation) - set(seen))]
write_tsv(ROOT / 'kallisto_tx2gene.tsv', list(seen.values()))
write_tsv(ROOT / 'GTF_transcripts_not_in_published_FASTA.tsv', missing)
receipt = dict(
    status='PUBLISHED_SEQUENCE_UNIVERSE_VERIFIED_NATIVE_INDEX_PENDING',
    selected_release='Ensembl 107 / GRCm39', original_article_exact_annotation_release='not reported',
    checksums=checks, sequence_counts=source_counts, reference_transcripts=len(seen),
    represented_genes=len({r['gene_id'] for r in seen.values()}),
    all_sequence_ids_versions_and_gene_ids_match_GTF=True,
    GTF_transcripts=len(annotation), GTF_transcripts_without_published_sequence=len(missing),
    missing_transcript_biotypes=dict(Counter(r['transcript_biotype'] for r in missing)),
    all_GTF_transcripts_represented=not missing,
    transcripts_shorter_than_k31=sum(r['length'] < 31 for r in seen.values()),
    combined_fasta_bytes=target.stat().st_size,
    combined_fasta_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
    interpretation='Use the explicitly represented sequence universe. Missing FASTA targets are unavailable, not zero-expression measurements. No new sequences are imputed.',
)
(ROOT / 'reference_concordance_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
print(json.dumps(receipt, indent=2))
