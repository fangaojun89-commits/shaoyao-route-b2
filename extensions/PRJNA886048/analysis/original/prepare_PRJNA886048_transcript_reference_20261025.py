"""Verify downloaded Ensembl 107 inputs and prepare their explicit target universe."""
from collections import Counter
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

ROOT = Path('R:/Shaoyao_storage/data/PRJNA886048/reference_preparation')

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
    data = (ROOT / name).read_bytes()
    total = 0
    for value in data:
        total = (((total >> 1) | ((total & 1) << 15)) + value) & 65535
    expected = next(line.split() for line in (ROOT / checksum).read_text().splitlines() if line.split()[-1] == name)
    assert total == int(expected[0]) and math.ceil(len(data) / 1024) == int(expected[1]), name
    checks.append(dict(file=name, official_BSD_sum=total, official_blocks_1024=int(expected[1]),
                       checksum_matches=True, sha256=hashlib.sha256(data).hexdigest()))

with (ROOT / 'transcript_to_gene.tsv').open(encoding='utf-8') as f:
    annotation = {row['transcript_id']: row for row in csv.DictReader(f, delimiter='\t')}
seen, source_counts = {}, {}
target = ROOT / 'transcripts_Ensembl107_GRCm39.fa'
temporary = target.with_suffix('.fa.part')
with temporary.open('w', encoding='ascii', newline='\n') as output:
    for filename in ['Mus_musculus.GRCm39.cdna.all.fa.gz', 'Mus_musculus.GRCm39.ncrna.fa.gz']:
        source_counts[filename] = 0
        with gzip.open(ROOT / filename, 'rt', encoding='ascii') as f:
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
