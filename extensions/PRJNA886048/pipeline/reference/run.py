"""Rebuild frozen annotation, sequence and memberships without expression data."""
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--source', type=Path, required=True)
p.add_argument('--resources', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
here = Path(__file__).resolve().parent
source, resources, out = a.source.resolve(strict=True), a.resources.resolve(strict=True), a.out.resolve()
if out.exists():
    raise SystemExit('Existing output refused; inspect partial or completed output rather than restarting')
for protected in (source, resources, here):
    if out == protected or out in protected.parents or protected in out.parents:
        raise SystemExit('Output overlaps protected input/code')
out.mkdir(parents=True)
steps = [
    [sys.executable, '-I', '-B', str(here/'parse_gtf.py'), '--gtf', str(source/'Mus_musculus.GRCm39.107.gtf.gz'), '--out', str(out/'annotation')],
    [sys.executable, '-I', '-B', str(here/'prepare_sequences.py'), '--source', str(source), '--annotation', str(out/'annotation/transcript_to_gene.tsv'), '--out', str(out/'sequences')],
    [sys.executable, '-I', '-B', str(here/'map_resources.py'), '--reference', str(out/'mapping_reference'), '--resources', str(resources), '--out', str(out/'memberships')],
]
for i, cmd in enumerate(steps):
    if i == 2:
        bridge = out/'mapping_reference'
        bridge.mkdir()
        shutil.copyfile(out/'annotation/genes.tsv', bridge/'genes.tsv')
        shutil.copyfile(out/'sequences/kallisto_tx2gene.tsv', bridge/'kallisto_tx2gene.tsv')
    (out/f'step_{i+1}_argv.json').write_text(json.dumps(cmd, indent=2)+'\n')
    with (out/f'step_{i+1}.stdout.txt').open('w',encoding='utf-8') as stdout, (out/f'step_{i+1}.stderr.txt').open('w',encoding='utf-8') as stderr:
        subprocess.run(cmd, stdout=stdout, stderr=stderr, check=True)
files = {}
for sub in ('annotation', 'sequences', 'memberships'):
    for f in sorted((out/sub).glob('*')):
        if f.is_file():
            files[f.relative_to(out).as_posix()] = dict(bytes=f.stat().st_size, sha256=hashlib.file_digest(f.open('rb'), 'sha256').hexdigest())
(out/'COMPLETED.json').write_text(json.dumps(dict(status='REFERENCE_RECONSTRUCTION_COMPLETE', files=files, index_builds=0, quantifications=0, models=0),indent=2)+'\n')
print('REFERENCE_RECONSTRUCTION_COMPLETE')
