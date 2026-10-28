"""Read-only consistency gate for canonical five-stage Linux native receipts.

This checks retained records and outputs, not remote-host attestation. Historical
46/47 orientation reuse and 54 continuation require their historical acceptance.
No transfer, native run, acceptance publication or statistical fit is performed.
"""
import argparse
import csv
import datetime as dt
import hashlib
import importlib.util
import json
import re
import shlex
from pathlib import Path, PurePosixPath

_content_path = Path(__file__).resolve().with_name('validate_PRJNA886048_native_content_20261028.py')
_spec = importlib.util.spec_from_file_location('prj_native_content', _content_path)
_content = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_content)
validate_content = _content.validate

LEDGER_SHA = '8c9dcbb9d2b83cb4442081bbdaf5c4f5c0d5fa0342a4dc2594a8c6ace54bc157'
TOOLS = {'fastp': '1144920d6f29f1ec56d36ea69bf58c50bc10cfe2c6421209689bad2c19787485',
         'kallisto': '78f60499a28d5957cffb95bb1d8d9cce6c2cafd0097f687640ed0b43f60b630e'}
INDEX_SHA = '8345052ff38d43971d50cbff8f5b767ad438cfb6e868c820d836bc056d2b1f73'
STAGES = ['ORIENTATION_unstranded', 'ORIENTATION_rf', 'ORIENTATION_fr', 'FASTP', 'KALLISTO']
OUTPUTS = {'fastp.json', 'fastp.html', 'quant/run_info.json', 'quant/abundance.tsv'}
FASTP_ARGS = ['-w', '8', '--detect_adapter_for_pe', '--cut_front', '--cut_tail',
              '--cut_front_window_size', '1', '--cut_tail_window_size', '1',
              '--cut_mean_quality', '20', '--length_required', '37',
              '--qualified_quality_phred', '15', '--unqualified_percent_limit', '40',
              '--n_base_limit', '5', '--trim_poly_g', '--poly_g_min_len', '10']


def require(value, message):
    if not value:
        raise ValueError(message)


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def instant(text):
    value = dt.datetime.fromisoformat(text)
    require(value.tzinfo is not None, 'Execution time lacks timezone')
    return value


def linux_path(text):
    require(isinstance(text, str) and text.startswith('/') and not any(c in text for c in '\x00\r\n'),
            'Expected an absolute Linux execution path')
    value = PurePosixPath(text)
    require('..' not in value.parts and str(value) == text, 'Noncanonical execution path')
    return value


def validate(native, tx2gene, ledger, run):
    require(run in {f'SRR217758{n}' for n in range(46, 58)}, 'Unknown run')
    require(sha(ledger) == LEDGER_SHA, 'ENA ledger identity differs')
    with ledger.open(encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    require(len(rows) == 12 and len({x['run_accession'] for x in rows}) == 12, 'Invalid sample ledger')
    row = next(x for x in rows if x['run_accession'] == run)
    require(row['scientific_name'] == 'Mus musculus' and row['library_layout'] == 'PAIRED',
            'Wrong source species/layout')
    expected = {}
    for url, size, md5 in zip(row['fastq_ftp'].split(';'), row['fastq_bytes'].split(';'),
                               row['fastq_md5'].split(';'), strict=True):
        expected[PurePosixPath(url).name] = {'bytes': int(size), 'md5': md5}
    require(set(expected) == {f'{run}_{m}.fastq.gz' for m in (1, 2)}, 'Wrong paired input identities')
    receipt = read(native / 'receipt.json')
    require(receipt.get('library') == run and receipt.get('status') == 'COMPLETE_FULL_LIBRARY_QUANTIFICATION',
            'Wrong or incomplete library receipt')
    require(receipt.get('input_checks') == expected, 'Original input MD5/size record differs')
    require(receipt.get('executable_sha256') == TOOLS and receipt.get('index_sha256') == INDEX_SHA,
            'Tool or index identity record differs')
    require(type(receipt.get('biological_models_run')) is int and receipt['biological_models_run'] == 0
            and receipt.get('index_rebuilt') is False, 'Unexpected native-stage scope')
    commands = receipt['commands']
    require([x['stage'] for x in commands] == STAGES,
            'Canonical five-stage record required; historical reuse/continuation needs separate acceptance')
    begin = instant(receipt['started_utc'])
    end = instant(receipt['finished_utc'])
    require(begin <= end, 'Reversed native execution times')
    prev = begin
    for stage in commands:
        require(type(stage.get('exit_code')) is int and stage['exit_code'] == 0, 'Unsuccessful native stage')
        require(type(stage.get('pid')) is int and stage['pid'] > 0, 'Missing child process identity')
        start, finish = instant(stage['started_utc']), instant(stage['finished_utc'])
        require(prev <= start <= finish <= end, 'Overlapping, reversed or out-of-window native stage')
        prev = finish
        require(isinstance(stage['command'], list) and all(isinstance(s, str) for s in stage['command']),
                'Command must be a token array')
    orientation = commands[0]['command']
    fastp = commands[3]['command']
    require(len(orientation) == 11 and len(fastp) == 13 + len(FASTP_ARGS), 'Unexpected native command shape')
    kallisto_path = str(linux_path(orientation[0]))
    fastp_path = str(linux_path(fastp[0]))
    index_path = str(linux_path(orientation[3]))
    output_root = linux_path(fastp[6]).parent
    inputs = [str(linux_path(fastp[i])) for i in (2, 4)]
    require([PurePosixPath(x).name for x in inputs] == [f'{run}_{m}.fastq.gz' for m in (1, 2)],
            'Command uses wrong raw mates or run')
    prefixes = [str(output_root / f'orientation_R{m}.fastq.gz') for m in (1, 2)]
    cleaned = [str(output_root / f'R{m}.fastq.gz') for m in (1, 2)]
    require(len(set(inputs + prefixes + cleaned)) == 6, 'Input/output path collision')
    for i, (label, option) in enumerate([('unstranded', []), ('rf', ['--rf-stranded']), ('fr', ['--fr-stranded'])]):
        expected_command = [kallisto_path, 'quant', '-i', index_path, '-o',
                            str(output_root / ('orientation_' + label)), '-t', '4', '--plaintext'] + option + prefixes
        require(commands[i]['command'] == expected_command, 'Orientation command differs: ' + label)
    expected_fastp = [fastp_path, '-i', inputs[0], '-I', inputs[1], '-o', cleaned[0], '-O', cleaned[1]]
    expected_fastp += FASTP_ARGS + ['--json', str(output_root / 'fastp.json'), '--html', str(output_root / 'fastp.html')]
    require(fastp == expected_fastp, 'FASTP scientific parameters or paths differ')
    expected_quant = [kallisto_path, 'quant', '-i', index_path, '-o', str(output_root / 'quant'),
                      '-t', '8', '--plaintext'] + cleaned
    require(commands[4]['command'] == expected_quant, 'Full quantification command differs')
    orient_info = receipt['orientation_run_info']
    require(set(orient_info) == {'unstranded', 'rf', 'fr'} and receipt.get('orientation') == 'unstranded',
            'Incomplete orientation evidence')
    for i, label in enumerate(('unstranded', 'rf', 'fr')):
        info = orient_info[label]
        require(info['n_processed'] == 250000 and type(info['n_processed']) is int
                and info['n_targets'] == 146190 and type(info['n_targets']) is int
                and info['kallisto_version'] == '0.46.1' and info['index_version'] == 10
                and type(info['n_bootstraps']) is int and info['n_bootstraps'] == 0,
                'Orientation sample size/tool/reference differs')
        require(all(type(info[k]) is int and info[k] >= 0 for k in ('n_unique', 'n_pseudoaligned'))
                and info['n_unique'] <= info['n_pseudoaligned'] <= 250000, 'Invalid orientation counts')
        require(shlex.split(info['call']) == commands[i]['command'], 'Orientation call record differs')
    denominator = orient_info['unstranded']['n_pseudoaligned']
    require(denominator > 0 and all(0.45 < orient_info[s]['n_pseudoaligned'] / denominator < 0.55 for s in ('rf', 'fr')),
            'Orientation requires scientific review')
    record_hashes = receipt.get('orientation_input_record_sha256')
    require(isinstance(record_hashes, list) and len(record_hashes) == 2
            and all(isinstance(h, str) and re.fullmatch('[0-9a-f]{64}', h) for h in record_hashes),
            'Missing orientation prefix record hashes')
    require(set(receipt['outputs']) == OUTPUTS, 'Wrong native output set')
    for name, identity in receipt['outputs'].items():
        target = native / name
        require(type(identity['bytes']) is int and target.stat().st_size == identity['bytes']
                and sha(target) == identity['sha256'], 'Output identity record differs: ' + name)
    qc, quant = read(native / 'fastp.json'), read(native / 'quant/run_info.json')
    require(receipt['qc_summary'] == qc['summary'] and receipt['filtering_result'] == qc['filtering_result']
            and receipt['run_info'] == quant, 'Receipt and actual QC/quantification differ')
    require(shlex.split(qc['command']) == fastp and shlex.split(quant['call']) == expected_quant,
            'Actual output command differs from successful stage record')
    content = validate_content(native, tx2gene)
    require(content['raw_pairs'] == int(row['read_count']), 'Raw pair count differs from source ledger')
    return {'status':'PASS_CANONICAL_NATIVE_RECORD_AND_CONTENT', 'run':run,
            'receipt_sha256':sha(native / 'receipt.json'), 'ledger_sha256':LEDGER_SHA,
            'recorded_original_inputs':expected, 'content':content,
            'stage_count':5, 'remote_execution_attested':False, 'raw_files_rehashed':False,
            'acceptance_published':False, 'historical_checkpoint_adoption_authorized':False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('native', 'tx2gene', 'ledger'):
        parser.add_argument('--' + key, type=Path, required=True)
    parser.add_argument('--run', required=True)
    args = parser.parse_args()
    print(json.dumps(validate(args.native, args.tx2gene, args.ledger, args.run), indent=2))
