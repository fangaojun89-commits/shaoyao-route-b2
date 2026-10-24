#!/usr/bin/env python3
"""Direct, staged October Rnor6/Ensembl101 STAR -> featureCounts comparison.

Default: print commands only. No downloads, installs, old driver, R or models.
Future execution order: admit -> index -> reservoir -> largest -> remaining ->
matrix. Use --execute --quiet-window-confirmed --reserved-growth-bytes N
--disk-budget-bytes N for each heavy stage, only after actual scheduling review.
These resource arguments are a measured-window contract, not a scheduler or an
assertion that H/Atlas is idle. Existing successful library stages are reused.
Failed/partial leaves are preserved and never overwritten or auto-retried.
"""
import argparse
from collections import Counter
import csv
import datetime
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import threading
import time

RUN = Path('/root/autodl-tmp/GSE205906_matched_STAR_Ensembl101_20261015r01')
F = Path('/root/autodl-tmp/GSE205906_full_read_counts_fasp_v18e_20261005/formal')
A = Path('/root/autodl-tmp/GSE205906_full_read_reconstruction_156_20261003')
META = RUN / 'metadata_preparation_01/METADATA_AND_RESERVOIR_BINDING.json'
META_SHA = '6dd7b53d40433baaacdeedf285204dee83f2d960c090e6a25f489f96d20fc5cb'
MANIFEST_SHA = '751e43b3aba8c2245981dc0041fec273c455c6390e85ef7ec4514e56118bf3b0'
SYNTHETIC_SHA = 'dc3f3975ee00cb8856034dd01d1b0997edea12d37a81a135f661a9348dc805f1'
PLAN_SHA = 'bb5207f1a0f616f0a2c2d9c4d7bbd5bdf32a17041b631a5a0c7a8dbdc3e8d622'
AUTH_SHA = 'b16c12f01912bbebd2e4a3a781ac2f83e53d7e478785539df4fbf503d8ab8a3b'
OLD_MATRIX_SHA = 'e98296fdf308dbaa82325a3ca6b0b041c6f2569e98be75a45757141b0a1f77a2'
TOOLS = {
    'STAR': (Path('/usr/lib/rna-star/bin/STAR-avx2'), 'fcc68116ebe3b77371b27eb524a538f8bcfe67e788edeee6ee4432b0831197e1'),
    'featureCounts': (A / 'runtime/bin/featureCounts', '357dbf889d747b7ea062568321adec91824fa8f59ce27b52c45e56369d37fd75'),
    'samtools': (A / 'runtime/bin/samtools', '83221bf89c8b74820f4f5c49b1c3aabf9aedb744da97ab113e5333c2bad3b6b5'),
}
FASTA = A / 'references/rat_transferred/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa'
GTF = A / 'references/derived/Rattus_norvegicus.Rnor_6.0.101.gtf'
GIB, MEMORY = 1024 ** 3, 40 * 1024 ** 3
INDEX = RUN / 'index_01/index'
STAGES = ('admit', 'index', 'reservoir', 'largest', 'remaining', 'matrix')
FC_OPTIONS = ['-p', '--countReadPairs', '-B', '-C', '-Q', '10', '-t', 'exon', '-g', 'gene_id', '-T', '2']


def need(value, message):
    if not value:
        raise ValueError(message)


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def pin(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(b)
    return {'path': str(path), 'bytes': Path(path).stat().st_size, 'sha256': h.hexdigest()}


def stat_file(path):
    p = Path(path); s = p.lstat()
    need(p.is_file() and not p.is_symlink(), 'Regular nonsymlink file required')
    return {'path': str(p), 'bytes': s.st_size, 'mtime_ns': s.st_mtime_ns, 'inode': s.st_ino}


def put(path, value):
    with Path(path).open('x', encoding='utf-8') as f:
        json.dump(value, f, indent=2, allow_nan=False); f.write('\n'); f.flush(); os.fsync(f.fileno())


def load(path, expected=None):
    if expected:
        need(pin(path)['sha256'] == expected, 'Pinned metadata changed')
    return json.loads(Path(path).read_bytes())


def unchanged(source_sha):
    need(pin(Path(__file__))['sha256'] == source_sha and pin(META)['sha256'] == META_SHA,
         'Worker source or frozen metadata changed')


def metadata(path=META):
    x = load(path, META_SHA)
    need(x['admitted'] is True and x['fixed_reservoir_pairs'] == 1000000
         and x['uniform_QC_read_length'] == 151 and x['sjdbOverhang'] == 150, 'Frozen metadata/QC contract')
    need(x['largest_full_library_pilot'] == {'label': 'SNI3', 'SRA_run': 'SRR19631202', 'registered_pairs': 30778167}, 'Fixed largest library')
    need(len(x['libraries']) == 9 and len(x['original_sample_order']) == 9, 'Nine-library domain')
    need([r['label'] for r in x['libraries']] == x['original_sample_order'] ==
         [f'{g}{i}' for g in ('JA', 'SNI', 'PLM') for i in (1, 2, 3)]
         and [r['SRA_run'] for r in x['libraries']] == [f'SRR{i}' for i in range(19631207, 19631198, -1)]
         and [r['GEO_sample'] for r in x['libraries']] == [f'GSM{i}' for i in range(6235292, 6235301)], 'Exact inverse GEO/SRA crosswalk')
    return x


def full_order(meta):
    return ['SRR19631202'] + [x['SRA_run'] for x in meta['libraries'] if x['SRA_run'] != 'SRR19631202']


def commands(stage, folder, mates=(), compressed=False, python=None, source=None, source_sha=None):
    star, sam, fc = [str(TOOLS[k][0]) for k in ('STAR', 'samtools', 'featureCounts')]
    if stage == 'index':
        return [[star, '--runMode', 'genomeGenerate', '--runThreadN', '2', '--runRNGseed', '155', '--genomeDir', str(INDEX),
                 '--genomeFastaFiles', str(FASTA), '--sjdbGTFfile', str(GTF), '--sjdbOverhang', '150']]
    folder = Path(folder)
    alignment = [star, '--runThreadN', '2', '--genomeDir', str(INDEX), '--sjdbOverhang', '150', '--readFilesIn', *map(str, mates),
                 '--twopassMode', 'Basic', '--runRNGseed', '155', '--alignIntronMin', '20',
                 '--alignIntronMax', '500000', '--outSAMtype', 'BAM', 'Unsorted',
                 '--outFileNamePrefix', str(folder / 'alignment.')]
    if compressed:
        # No --option tokens inside STAR's readFilesCommand argument list.
        alignment += ['--readFilesCommand', python or sys.executable, '-I', '-B',
                      str(source or Path(__file__).resolve()), 'stream', str(folder / 'streams'),
                      source_sha or pin(Path(__file__))['sha256']]
    bam = folder / 'sorted.bam'
    return [alignment,
            [sam, 'sort', '-@', '1', '-m', '1G', '-T', str(folder / 'sort_tmp'), '-o', str(bam),
             str(folder / 'alignment.Aligned.out.bam')],
            [sam, 'quickcheck', str(bam)],
            *[[fc, *FC_OPTIONS, '-s', str(s), '-a', str(GTF), '-o', str(folder / f'counts_s{s}.tsv'), str(bam)] for s in (0, 1, 2)]]


def expected_md5():
    m = load(F / 'provenance/run_manifest.json', MANIFEST_SHA)
    return {str(F / 'raw' / f"{r['SRA_run']}_{mate}.fastq.gz"): md5
            for r in m['metadata_rows']
            for mate, md5 in enumerate(r['ENA_generated_FASTQ_MD5'].split(';'), 1)}


def inflate(path, expected, result):
    """Only gzip bytes go to stdout; CRC and frozen compressed MD5 are checked."""
    need(str(path) in expected_md5() and expected_md5()[str(path)] == expected
         and Path(result).resolve().is_relative_to((RUN / 'libraries').resolve()), 'Frozen internal inflater paths')
    before = stat_file(path); md5 = hashlib.md5(); sha = hashlib.sha256(); consumed = 0
    class DigestReader:
        def __init__(self, raw): self.raw = raw
        def read(self, n=-1):
            nonlocal consumed
            b = self.raw.read(n); md5.update(b); sha.update(b); consumed += len(b); return b
    receipt = {'input': before, 'full_gzip_EOF_CRC': False, 'frozen_MD5_match': False,
               'source_sha256': pin(Path(__file__))['sha256']}
    code = 1
    try:
        with Path(path).open('rb') as raw, gzip.GzipFile(fileobj=DigestReader(raw)) as gz:
            for block in iter(lambda: gz.read(1024 * 1024), b''):
                sys.stdout.buffer.write(block)
        sys.stdout.buffer.flush()
        need(before == stat_file(path) and consumed == before['bytes'] and md5.hexdigest() == expected, 'Full frozen gzip identity')
        receipt.update(full_gzip_EOF_CRC=True, frozen_MD5_match=True, compressed_bytes=consumed,
                       compressed_MD5=md5.hexdigest(), compressed_SHA256=sha.hexdigest())
        code = 0
    except Exception as exc:
        receipt['exception_class_only'] = type(exc).__name__
    put(result, receipt)
    return code


def stream(directory, source_sha, path):
    """Each Basic pass has its own atomic leaf and waited inflater exit.

    The outer process retains the FIFO stdout descriptor until after the child
    wait and EXIT receipt. Thus a normal all-input EOF cannot precede those
    receipts. STAR's closeReadsFiles kills its command-shell PID; an early
    close/missing receipt is a failed stream, never inferred successful.
    """
    directory, path = Path(directory), Path(path)
    need(directory.resolve().is_relative_to((RUN / 'libraries').resolve()) and directory.name == 'streams', 'Owned stream directory only')
    need(pin(Path(__file__))['sha256'] == source_sha, 'Inflate source changed')
    md5 = expected_md5(); need(str(path) in md5, 'Only frozen raw mates')
    mate_directory = directory / path.name; mate_directory.mkdir(exist_ok=True)
    invocation = None
    for pass_number in (1, 2):
        candidate = mate_directory / f'pass{pass_number:02d}'
        try:
            candidate.mkdir(); invocation = candidate; break
        except FileExistsError:
            pass
    need(invocation is not None, 'Exactly two Basic read opens per mate; never overwrite/retry')
    result, terminal = invocation / 'CRC.json', invocation / 'EXIT.json'
    p = subprocess.run([sys.executable, '-I', '-B', str(Path(__file__).resolve()), '_inflate', str(path), md5[str(path)], str(result)],
                       stdout=sys.stdout.buffer, stderr=subprocess.PIPE)
    child = load(result) if result.exists() else None
    accepted = p.returncode == 0 and child is not None and child['full_gzip_EOF_CRC'] is True and child['frozen_MD5_match'] is True
    put(terminal, {'input_path': str(path), 'actual_child_exit': p.returncode, 'accepted': accepted,
                   'Basic_read_open_ordinal': pass_number,
                   'child_receipt': pin(result) if result.exists() else None, 'source_sha256': source_sha,
                   'stderr_bytes': len(p.stderr), 'stderr_sha256': hashlib.sha256(p.stderr).hexdigest(),
                   'raw_condition_text_not_retained': True})
    return 0 if accepted else 1


def check_streams(folder, mates, source_sha):
    md5 = expected_md5()
    for path in mates:
        mate_directory = Path(folder) / 'streams' / Path(path).name
        need(sorted(p.name for p in mate_directory.iterdir()) == ['pass01', 'pass02'], 'Both Basic passes required per mate')
        for ordinal in (1, 2):
            invocation = mate_directory / f'pass{ordinal:02d}'; result = invocation / 'CRC.json'
            terminal, crc = load(invocation / 'EXIT.json'), load(result)
            need(terminal['Basic_read_open_ordinal'] == ordinal and terminal['input_path'] == str(path)
                 and terminal['actual_child_exit'] == 0 and terminal['accepted'] is True
                 and terminal['source_sha256'] == source_sha and terminal['child_receipt'] == pin(result)
                 and crc['source_sha256'] == source_sha and crc['full_gzip_EOF_CRC'] is True
                 and crc['frozen_MD5_match'] is True and crc['compressed_MD5'] == md5[str(path)]
                 and crc['input'] == stat_file(path), 'Any mate/pass gzip child failure invalidates STAR exit zero')


def tree_bytes(path):
    return sum(p.stat().st_size for p in Path(path).rglob('*') if p.is_file() and not p.is_symlink())


def resources(args):
    need(args.quiet_window_confirmed and args.reserved_growth_bytes is not None
         and args.reserved_growth_bytes >= 0 and args.disk_budget_bytes is not None
         and args.disk_budget_bytes > 0, 'Explicit quiet window, other-job reserve and stage disk budget required')
    limit = Path('/sys/fs/cgroup/memory.max').read_text().strip()
    current = int(Path('/sys/fs/cgroup/memory.current').read_text())
    need(limit != 'max' and int(limit) - current >= MEMORY, '40 GiB process budget lacks current memory headroom')
    free = shutil.disk_usage(RUN).free
    need(free >= args.reserved_growth_bytes + args.disk_budget_bytes, 'Reserved disk headroom unavailable')
    return {'free_bytes': free, 'cgroup_memory_current': current, 'cgroup_memory_limit': int(limit),
            'declared_other_job_growth_reserve': args.reserved_growth_bytes, 'declared_stage_disk_budget': args.disk_budget_bytes,
            'quiet_window_confirmed_by_caller': True, 'H_or_Atlas_idle_not_inferred_by_code': True}


def run_command(argv, events, folder, reserve):
    import resource
    def limits(): resource.setrlimit(resource.RLIMIT_AS, (MEMORY, MEMORY))
    before = time.monotonic(); low_free = shutil.disk_usage(RUN).free; peak = tree_bytes(folder)
    p = subprocess.Popen(argv, cwd=folder, env=dict(os.environ, LC_ALL='C', LANG='C', OMP_NUM_THREADS='1',
                         OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         preexec_fn=limits)
    while True:
        try:
            stdout, stderr = p.communicate(timeout=1); break
        except subprocess.TimeoutExpired:
            low_free = min(low_free, shutil.disk_usage(RUN).free); peak = max(peak, tree_bytes(folder))
            # Do not signal other jobs. A depleted reservation is recorded and
            # blocks the next command; no automatic retry/deletion is performed.
    low_free = min(low_free, shutil.disk_usage(RUN).free); peak = max(peak, tree_bytes(folder))
    event = {'command': list(map(str, argv)), 'actual_exit': p.returncode, 'elapsed_seconds': time.monotonic() - before,
             'cumulative_native_child_maxrss_kib': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
             'maximum_owned_stage_bytes_observed': peak, 'minimum_cloud_free_bytes_observed': low_free,
             'stdout_bytes': len(stdout), 'stdout_sha256': hashlib.sha256(stdout).hexdigest(),
             'stderr_bytes': len(stderr), 'stderr_sha256': hashlib.sha256(stderr).hexdigest()}
    events.append(event)
    need(p.returncode == 0 and low_free >= reserve, 'Native command exit/resource reservation failed')
    return stdout


def reference_dictionary(meta):
    before = [stat_file(p) for p in (FASTA, GTF)]
    dictionary = {}; h = hashlib.sha256(); current = None
    with FASTA.open('rb') as f:
        for line in f:
            h.update(line)
            if line.startswith(b'>'):
                current = line[1:].split()[0].decode('ascii')
                need(current not in dictionary, 'Duplicate FASTA name'); dictionary[current] = 0
            else:
                seq = line.strip(); need(current is not None and re.fullmatch(b'[A-Za-z]*', seq), 'FASTA sequence domain')
                dictionary[current] += len(seq)
    need(h.hexdigest() == meta['expected_reference_pins_not_freshly_hashed'][0]['sha256'], 'Complete FASTA identity')
    h = hashlib.sha256(); genes = set(); exon_n = 0
    with GTF.open('rb') as f:
        for line in f:
            h.update(line)
            if line.startswith(b'#'): continue
            x = line.decode('utf-8').rstrip('\r\n').split('\t')
            need(len(x) == 9 and x[0] in dictionary and 1 <= int(x[3]) <= int(x[4]) <= dictionary[x[0]], 'GTF coordinate/dictionary mismatch')
            if x[2] == 'exon':
                ids = re.findall(r'(?:^|;\s*)gene_id "([^";]+)";', x[8])
                need(len(ids) == 1, 'Exact exon gene_id required'); genes.add(ids[0]); exon_n += 1
    need(h.hexdigest() == meta['expected_reference_pins_not_freshly_hashed'][1]['sha256'] and len(genes) == 32883, 'Complete GTF identity/domain')
    need(before == [stat_file(p) for p in (FASTA, GTF)], 'Reference changed during full admission')
    return {'sequence_lengths': dictionary, 'exon_genes': len(genes), 'exon_rows': exon_n,
            'reference_pins': meta['expected_reference_pins_not_freshly_hashed'], 'reference_stat': before}


def tools_admission(events, folder, reserve):
    answer = {}
    for name, (path, expected) in TOOLS.items():
        p = pin(path); need(p['sha256'] == expected, 'Exact selected native tool')
        run_command([str(path), {'STAR': '--version', 'featureCounts': '-v', 'samtools': '--version-only'}[name]], events, folder, reserve)
        out = run_command(['/usr/bin/ldd', str(path)], events, folder, reserve).decode('utf-8')
        need('not found' not in out, 'Missing native shared dependency')
        dependencies = sorted(set(re.findall(r'(?:=>\s+|^\s*)(/[^\s]+)', out, re.M)))
        answer[name] = {'binary': p, 'linked_libraries': [pin(Path(x).resolve()) for x in dependencies]}
    return answer


def completed(folder, source_sha):
    r = load(Path(folder) / 'COMPLETE.json')
    need(r['status'] == 'COMPLETE' and r['source_sha256'] == source_sha and r['metadata_sha256'] == META_SHA, 'Bound completed stage required')
    for p in r.get('small_products', []): need(pin(Path(p['path'])) == p, 'Completed compact product changed')
    for p in r.get('current_bams_stat_only', []): need(stat_file(Path(p['path'])) == p, 'Completed BAM stat drift')
    return r


def fc_table(path):
    with Path(path).open() as f:
        reader = csv.DictReader((line for line in f if not line.startswith('#')), delimiter='\t')
        need(reader.fieldnames[:6] == ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'] and len(reader.fieldnames) == 7, 'One-BAM native FC table')
        answer = {}
        for row in reader:
            value = row[reader.fieldnames[6]]; need(re.fullmatch(r'[0-9]+', value) and row['Geneid'] not in answer, 'Unique ASCII nonnegative integer fragment counts')
            answer[row['Geneid']] = int(value)
    need(len(answer) == 32883, 'Complete unfiltered gene domain')
    return answer


def bam_qc(path, folder, events, reserve, dictionary):
    sam = str(TOOLS['samtools'][0])
    header = run_command([sam, 'view', '-H', str(path)], events, folder, reserve).decode()
    observed = {}
    for line in header.splitlines():
        if line.startswith('@SQ\t'):
            x = dict(s.split(':', 1) for s in line.split('\t')[1:]); observed[x['SN']] = int(x['LN'])
    need(observed == dictionary, 'BAM/index/reference dictionary equality')
    p = subprocess.Popen([sam, 'view', str(path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    err = {'bytes': 0, 'sha': hashlib.sha256()}
    def drain():
        for block in iter(lambda: p.stderr.read(65536), b''): err['bytes'] += len(block); err['sha'].update(block)
    thread = threading.Thread(target=drain); thread.start()
    counts, mapq, nh = Counter(), Counter(), Counter()
    started = time.monotonic()
    try:
        for line in p.stdout:
            x = line.decode().rstrip('\n').split('\t'); flag = int(x[1]); counts['alignment_records'] += 1
            if flag & 4: counts['unmapped_records'] += 1; continue
            tags = dict(t.split(':', 2)[::2] for t in x[11:] if t.count(':') >= 2)
            need(all(k in tags for k in ('NH', 'HI', 'AS', 'nM')), 'Native STAR tags retained')
            mapq[int(x[4])] += 1; nh[int(tags['NH'])] += 1
            if not flag & (256 | 2048): counts['primary_mapped_records'] += 1
    finally:
        p.stdout.close(); code = p.wait(); thread.join()
        import resource
        events.append({'command': [sam, 'view', str(path)], 'actual_exit': code,
                       'elapsed_seconds': time.monotonic() - started,
                       'cumulative_native_child_maxrss_kib': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                       'maximum_owned_stage_bytes_observed': tree_bytes(folder),
                       'minimum_cloud_free_bytes_observed': shutil.disk_usage(RUN).free})
    need(code == 0, 'Native BAM QC stream exited nonzero')
    return {'alignment_record_counts': dict(counts), 'MAPQ_alignment_record_histogram': dict(mapq),
            'NH_alignment_record_histogram': dict(nh), 'samtools_view_actual_exit': code,
            'samtools_view_stderr_bytes': err['bytes'], 'samtools_view_stderr_sha256': err['sha'].hexdigest(),
            'STAR_default_unmapped_BAM_output': 'None; input-pair denominators come from STAR, not the FC table total',
            'equal_MAPQ_error_calibration_claimed': False}


def library(meta, library, folder, args, reference, source_sha, reservoir=False):
    if (folder / 'COMPLETE.json').exists(): return completed(folder, source_sha)
    resources(args); folder.mkdir(parents=True, exist_ok=False); events = []; started = time.monotonic()
    receipt = {'status': 'INCOMPLETE', 'source_sha256': source_sha, 'metadata_sha256': META_SHA,
               'started_utc': utc(), 'events': events, 'technical_reservoir_only': reservoir,
               'SRA_run': library['SRA_run'], 'label': library['label'], 'resource_before': resources(args)}
    try:
        if reservoir:
            mates = [Path(p['path']) for p in meta['reservoir_fresh_file_pins'][:2]]; pairs = 1000000
            for p in meta['reservoir_fresh_file_pins']: need(pin(Path(p['path'])) == p, 'Same fixed reservoir and ordinals')
        else:
            mates = [Path(p['path']) for p in library['raw_mates']]; pairs = library['registered_pairs']; (folder / 'streams').mkdir()
            for p, m in zip(mates, library['raw_mates']): need(stat_file(p)['bytes'] == m['bytes'], 'Frozen raw mate size')
        cmds = commands('library', folder, mates, not reservoir, source_sha=source_sha)
        run_command(cmds[0], events, folder, args.reserved_growth_bytes)
        stream_products = []
        if not reservoir:
            check_streams(folder, mates, source_sha)
            stream_products = [pin(folder / 'streams' / mate.name / f'pass{ordinal:02d}' / name)
                               for mate in mates for ordinal in (1, 2) for name in ('CRC.json', 'EXIT.json')]
        native = {}
        for line in (folder / 'alignment.Log.final.out').read_text().splitlines():
            if '|' in line:
                k, v = line.split('|', 1); native[k.strip()] = v.strip()
        need(int(native['Number of input reads']) == pairs, 'STAR must consume every registered pair')
        for cmd in cmds[1:]: run_command(cmd, events, folder, args.reserved_growth_bytes)
        qc = bam_qc(folder / 'sorted.bam', folder, events, args.reserved_growth_bytes, reference['sequence_lengths'])
        qc.update(input_pairs=pairs, STAR_Log_final=native, counts={}, assignment={}, assigned_fraction={})
        products = stream_products
        for strand in (0, 1, 2):
            path = folder / f'counts_s{strand}.tsv'; table = fc_table(path)
            with Path(str(path) + '.summary').open() as f:
                rows = list(csv.reader(f, delimiter='\t'))
            assignment = {row[0]: int(row[1]) for row in rows[1:]}
            need(sum(table.values()) == assignment['Assigned'] <= pairs, 'Fragment assigned/input denominator')
            qc['counts'][str(strand)] = {'genes': len(table), 'sum': sum(table.values())}
            qc['assignment'][str(strand)] = assignment; qc['assigned_fraction'][str(strand)] = assignment['Assigned'] / pairs
            products += [pin(path), pin(Path(str(path) + '.summary'))]
        put(folder / 'QC.json', qc); products += [pin(folder / 'QC.json'), pin(folder / 'alignment.Log.final.out')]
        receipt.update(status='COMPLETE', input_pairs=pairs, small_products=products,
                       retained_stage_bytes=tree_bytes(folder), peak_stage_bytes_observed=max(e['maximum_owned_stage_bytes_observed'] for e in events),
                       stage_elapsed_seconds=time.monotonic() - started, native_counting_primary_strand=0,
                       current_bams_stat_only=[stat_file(folder / n) for n in ('alignment.Aligned.out.bam', 'sorted.bam')],
                       index_alignment_sort_not_overlapped=True, new_model_run=False)
        for p in reference['reference_stat']: need(stat_file(Path(p['path'])) == p, 'Reference changed during library processing')
        unchanged(source_sha)
        put(folder / 'COMPLETE.json', receipt)
    except Exception as exc:
        receipt.update(status='FAILED', exception_class_only=type(exc).__name__)
        put(folder / 'FAILED.json', receipt); raise
    return receipt


def matrix(meta, source_sha):
    source = F / 'output/GSE205906_fixed9_gene_fragments_s0.tsv'
    need(pin(source)['sha256'] == OLD_MATRIX_SHA, 'Preserved HISAT2 comparison identity')
    with source.open() as f:
        rows = csv.reader(f, delimiter='\t'); header = next(rows); genes = [r[0] for r in rows]
    need(header == ['gene_id', *meta['original_sample_order']] and len(genes) == len(set(genes)) == 32883
         and genes == sorted(genes) and all(re.fullmatch(r'ENSRNOG[0-9]+', g) for g in genes), 'Frozen full lexical ENSRNOG/header order')
    data, inputs = {}, []
    for lib in meta['libraries']:
        folder = RUN / 'libraries' / lib['SRA_run']; done = completed(folder, source_sha)
        need(not done['technical_reservoir_only'] and done['input_pairs'] == lib['registered_pairs'], 'Full library, not reservoir')
        inputs.append(pin(folder / 'COMPLETE.json'))
        for strand in (0, 1, 2):
            table = fc_table(folder / f'counts_s{strand}.tsv'); need(set(table) == set(genes), 'Exact full gene-ID domain')
            data[(lib['label'], strand)] = table
    output = RUN / 'matrix_01'; output.mkdir(exist_ok=False)
    products = []
    for strand in (0, 1, 2):
        path = output / f'GSE205906_STAR_fixed9_gene_fragments_s{strand}.tsv'
        with path.open('x', newline='', encoding='ascii') as f:
            writer = csv.writer(f, delimiter='\t', lineterminator='\n'); writer.writerow(['gene_id', *meta['original_sample_order']])
            writer.writerows([g, *[data[(label, strand)][g] for label in meta['original_sample_order']]] for g in genes)
        products.append(pin(path))
    unchanged(source_sha)
    put(output / 'COMPLETE.json', {'status': 'COMPLETE', 'source_sha256': source_sha, 'metadata_sha256': META_SHA,
        'dimensions': [32883, 9], 'full_libraries': inputs, 'small_products': products, 'primary_strand': 0,
        'original_matrix_preserved': pin(source), 'new_models_or_candidate_admission': False})


def execute(args, meta, source_sha):
    need(sys.platform.startswith('linux') and Path(args.metadata) == META, 'Execute only the fixed cloud inputs')
    need(RUN.is_dir() and not RUN.is_symlink(), 'Existing declared independent root required')
    if args.stage == 'matrix': return matrix(meta, source_sha)
    before = resources(args)
    if args.stage == 'admit':
        folder = RUN / 'reference_admission_01'; folder.mkdir(exist_ok=False); events = []
        try:
            synthetic = load(RUN / 'native_count_synthetic_01/NATIVE_SYNTHETIC_COUNT_RECEIPT.json', SYNTHETIC_SHA)
            need(synthetic['accepted'] is True and len(synthetic['checks']) == 18 and all(x['passed'] for x in synthetic['checks']), 'Retain accepted native tests')
            for p in meta['expected_reference_pins_not_freshly_hashed']:
                need(stat_file(Path(p['path']))['bytes'] == p['bytes'], 'Full reference, never the short sibling')
            reference = reference_dictionary(meta); tools = tools_admission(events, folder, args.reserved_growth_bytes)
            unchanged(source_sha)
            put(folder / 'REFERENCE.json', reference)
            put(folder / 'COMPLETE.json', {'status': 'COMPLETE', 'source_sha256': source_sha, 'metadata_sha256': META_SHA,
                'reference': reference, 'tools': tools, 'resource_before': before, 'events': events,
                'small_products': [pin(folder / 'REFERENCE.json')], 'plan_sha256': PLAN_SHA, 'authorization_sha256': AUTH_SHA,
                'Python_interpreter': pin(Path(sys.executable).resolve()), 'Python_version': sys.version,
                'Python_stdlib_bytes_not_individually_authenticated': True})
        except Exception as exc:
            put(folder / 'FAILED.json', {'status': 'FAILED', 'stage': 'admit', 'source_sha256': source_sha,
                'metadata_sha256': META_SHA, 'resource_before': before, 'events': events,
                'exception_class_only': type(exc).__name__, 'condition_text_not_retained': True})
            raise
        return
    admission = completed(RUN / 'reference_admission_01', source_sha)
    for item in admission['tools'].values():
        for p in [item['binary'], *item['linked_libraries']]: need(pin(Path(p['path'])) == p, 'Admitted tool/dependency drift')
    reference = admission['reference']
    for p in reference['reference_stat']: need(stat_file(Path(p['path'])) == p, 'Admitted full reference stat drift')
    if args.stage == 'index':
        folder = RUN / 'index_01'; folder.mkdir(exist_ok=False); events = []
        try:
            INDEX.mkdir()
            run_command(commands('index', folder)[0], events, folder, args.reserved_growth_bytes)
            index_dict = dict(line.rstrip().split('\t') for line in (INDEX / 'chrNameLength.txt').read_text().splitlines())
            need({k: int(v) for k, v in index_dict.items()} == reference['sequence_lengths'], 'STAR index dictionary equality')
            products = [pin(p) for p in sorted(INDEX.iterdir()) if p.is_file()]
            need(all((INDEX / x).stat().st_size > 0 for x in ('Genome', 'SA', 'SAindex', 'genomeParameters.txt')), 'Complete native STAR index')
            for p in reference['reference_stat']: need(stat_file(Path(p['path'])) == p, 'Reference changed during index generation')
            unchanged(source_sha)
            put(folder / 'COMPLETE.json', {'status': 'COMPLETE', 'source_sha256': source_sha, 'metadata_sha256': META_SHA,
                'index_files': products, 'index_stat': [stat_file(Path(p['path'])) for p in products], 'events': events,
                'reference_admission': pin(RUN / 'reference_admission_01/COMPLETE.json'), 'small_products': []})
        except Exception as exc:
            put(folder / 'FAILED.json', {'status': 'FAILED', 'stage': 'index', 'source_sha256': source_sha,
                'metadata_sha256': META_SHA, 'resource_before': before, 'events': events,
                'exception_class_only': type(exc).__name__, 'condition_text_not_retained': True})
            raise
        return
    index = completed(RUN / 'index_01', source_sha)
    for p in index['index_stat']: need(stat_file(Path(p['path'])) == p, 'Admitted index stat drift')
    if args.stage == 'reservoir':
        lib = next(x for x in meta['libraries'] if x['SRA_run'] == 'SRR19631199')
        library(meta, lib, RUN / 'reservoir_01', args, reference, source_sha, reservoir=True); return
    pilot = completed(RUN / 'reservoir_01', source_sha)
    largest = next(x for x in meta['libraries'] if x['SRA_run'] == 'SRR19631202')
    if args.stage == 'largest':
        projection = math.ceil(pilot['peak_stage_bytes_observed'] * largest['registered_pairs'] / 1000000)
        need(args.disk_budget_bytes >= projection, 'Budget must cover measured reservoir extrapolation; not a guarantee')
        library(meta, largest, RUN / 'libraries/SRR19631202', args, reference, source_sha); return
    large_done = completed(RUN / 'libraries/SRR19631202', source_sha)
    remaining = [x for x in meta['libraries'] if x['SRA_run'] != 'SRR19631202']
    pending = [x for x in remaining if not (RUN / 'libraries' / x['SRA_run'] / 'COMPLETE.json').exists()]
    projection = math.ceil(large_done['retained_stage_bytes'] * sum(x['registered_pairs'] for x in pending) / largest['registered_pairs']
                          + max(0, large_done['peak_stage_bytes_observed'] - large_done['retained_stage_bytes']))
    need(args.disk_budget_bytes >= projection, 'Remaining budget must cover measured largest-library projection')
    remaining_budget = args.disk_budget_bytes
    for lib in remaining:
        folder = RUN / 'libraries' / lib['SRA_run']; reused = (folder / 'COMPLETE.json').exists()
        child_args = argparse.Namespace(**vars(args)); child_args.disk_budget_bytes = remaining_budget
        done = library(meta, lib, folder, child_args, reference, source_sha)
        if not reused: remaining_budget = max(1, remaining_budget - done['retained_stage_bytes'])


def main():
    if len(sys.argv) > 1 and sys.argv[1] == '_inflate': return inflate(Path(sys.argv[2]), sys.argv[3], Path(sys.argv[4]))
    if len(sys.argv) > 1 and sys.argv[1] == 'stream': return stream(*sys.argv[2:])
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=STAGES)
    parser.add_argument('--metadata', type=Path, default=META, help='Exact metadata mirror allowed for dry-run only')
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--quiet-window-confirmed', action='store_true')
    parser.add_argument('--reserved-growth-bytes', type=int)
    parser.add_argument('--disk-budget-bytes', type=int)
    args = parser.parse_args(); meta = metadata(args.metadata); source_sha = pin(Path(__file__))['sha256']
    if not args.execute:
        largest = next(x for x in meta['libraries'] if x['SRA_run'] == 'SRR19631202')
        folder = RUN / ('reservoir_01' if args.stage == 'reservoir' else 'libraries/SRR19631202')
        mates = [x['path'] for x in meta['reservoir_fresh_file_pins'][:2]] if args.stage == 'reservoir' else [x['path'] for x in largest['raw_mates']]
        print(json.dumps({'status': 'DRY_RUN_NOT_RESOURCE_ADMITTED', 'stage': args.stage, 'plan_sha256': PLAN_SHA,
            'authorization_sha256': AUTH_SHA, 'full_library_order': full_order(meta), 'sample_order': meta['original_sample_order'],
            'commands': commands('index' if args.stage == 'index' else 'library', folder, mates, args.stage != 'reservoir', source_sha=source_sha)
                        if args.stage in ('index', 'reservoir', 'largest') else [],
            'remaining_serial_commands': [{'label': x['label'], 'SRA_run': x['SRA_run'],
                'commands': commands('library', RUN / 'libraries' / x['SRA_run'], [p['path'] for p in x['raw_mates']], True, source_sha=source_sha)}
                for x in meta['libraries'] if x['SRA_run'] != 'SRR19631202'] if args.stage == 'remaining' else [],
            'matrix_excludes_reservoir': True, 'source_sha256': source_sha, 'actual_heavy_execution': False}))
        return 0
    execute(args, meta, source_sha)
    print(json.dumps({'stage': args.stage, 'actual_stage_call_returned': True, 'new_model_run': False}))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
