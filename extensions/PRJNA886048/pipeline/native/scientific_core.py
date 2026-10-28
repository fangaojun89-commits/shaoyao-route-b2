"""Frozen scientific stages extracted from the actually executed R3 entry.

Only the function boundary/indentation changes. Native execution is unverified
through the portable wrapper; no completed library should be reprocessed here.
"""
import gzip
import hashlib


def run_stages(*, run, names, sources, out, receipt, execute, save, digest, fastp, kallisto, index):
    # One technical prefix per new library, before treatment effects are fitted.
    prefixes = [out / f'orientation_R{m}.fastq.gz' for m in (1, 2)]
    hashes = [hashlib.sha256(), hashlib.sha256()]
    with gzip.open(sources[names[0]], 'rb') as a, gzip.open(sources[names[1]], 'rb') as b, \
            gzip.open(prefixes[0], 'wb', compresslevel=1) as x, gzip.open(prefixes[1], 'wb', compresslevel=1) as y:
        for _ in range(250000):
            records = [[a.readline() for _ in range(4)], [b.readline() for _ in range(4)]]
            for m, lines in enumerate(records):
                assert all(lines) and lines[0].startswith(b'@') and lines[2].startswith(b'+')
                assert len(lines[1].rstrip()) == len(lines[3].rstrip())
                hashes[m].update(b''.join(lines))
            assert records[0][0].split()[0].removesuffix(b'/1') == records[1][0].split()[0].removesuffix(b'/2')
            x.writelines(records[0])
            y.writelines(records[1])
    receipt['orientation_input_record_sha256'] = [h.hexdigest() for h in hashes]
    receipt['orientation_run_info'] = {}
    for label, options in [('unstranded', []), ('rf', ['--rf-stranded']), ('fr', ['--fr-stranded'])]:
        target = out / ('orientation_' + label)
        execute('ORIENTATION_' + label, [str(kallisto), 'quant', '-i', str(index), '-o', str(target),
                                          '-t', '4', '--plaintext'] + options + [str(p) for p in prefixes])
        info = json.loads((target / 'run_info.json').read_text())
        assert info['n_processed'] == 250000 and info['n_targets'] == 146190
        receipt['orientation_run_info'][label] = info
        save()
    count = receipt['orientation_run_info']['unstranded']['n_pseudoaligned']
    assert count > 0
    assert all(0.45 < receipt['orientation_run_info'][s]['n_pseudoaligned'] / count < 0.55 for s in ['rf', 'fr']), 'Different orientation pattern requires review before full quantification'
    receipt['orientation'] = 'unstranded'
    receipt['orientation_basis'] = 'First 250000 deposited pairs from this complete MD5-verified library; neither RF nor FR dominates'
    command = [str(fastp), '-i', str(sources[names[0]]), '-I', str(sources[names[1]]),
               '-o', str(out / 'R1.fastq.gz'), '-O', str(out / 'R2.fastq.gz'), '-w', '8',
               '--detect_adapter_for_pe', '--cut_front', '--cut_tail', '--cut_front_window_size', '1',
               '--cut_tail_window_size', '1', '--cut_mean_quality', '20', '--length_required', '37',
               '--qualified_quality_phred', '15', '--unqualified_percent_limit', '40', '--n_base_limit', '5',
               '--trim_poly_g', '--poly_g_min_len', '10', '--json', str(out / 'fastp.json'), '--html', str(out / 'fastp.html')]
    execute('FASTP', command)
    qc = json.loads((out / 'fastp.json').read_text())
    receipt['qc_summary'] = qc['summary']
    receipt['filtering_result'] = qc['filtering_result']
    retained_reads = qc['summary']['after_filtering']['total_reads']
    assert retained_reads > 0 and retained_reads % 2 == 0
    execute('KALLISTO', [str(kallisto), 'quant', '-i', str(index), '-o', str(out / 'quant'),
                         '-t', '8', '--plaintext', str(out / 'R1.fastq.gz'), str(out / 'R2.fastq.gz')])
    info = json.loads((out / 'quant/run_info.json').read_text())
    assert info['n_processed'] == retained_reads // 2 and info['n_targets'] == 146190
    receipt['run_info'] = info
    receipt['outputs'] = {str(p.relative_to(out)): dict(bytes=p.stat().st_size, sha256=digest(p))
                          for p in [out / 'fastp.json', out / 'fastp.html', out / 'quant/abundance.tsv', out / 'quant/run_info.json']}
    receipt['status'] = 'COMPLETE_FULL_LIBRARY_QUANTIFICATION'
