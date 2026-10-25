#!/usr/bin/env python3
"""Bind retained metadata/QC and the existing reservoir, without real alignment.

Reads small successful provenance records, FastQC archives and reservoir hashes.
Does not read the full raw FASTQs or hash/rebuild the full genome during other
cloud jobs. Reference pins here remain expected historical pins until admission
in the later quiet execution window.
"""
import csv
import datetime
import hashlib
import json
from pathlib import Path
import shutil
import time
import zipfile


FORMAL = Path('/root/autodl-tmp/GSE205906_full_read_counts_fasp_v18e_20261005/formal')
ASSET = Path('/root/autodl-tmp/GSE205906_full_read_reconstruction_156_20261003')
RUN = Path('/root/autodl-tmp/GSE205906_matched_STAR_Ensembl101_20261015r01')
MANIFEST_SHA = '751e43b3aba8c2245981dc0041fec273c455c6390e85ef7ec4514e56118bf3b0'
SYNTHETIC_SHA = 'dc3f3975ee00cb8856034dd01d1b0997edea12d37a81a135f661a9348dc805f1'


def pin(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return {'path': str(path), 'bytes': path.stat().st_size, 'sha256': digest.hexdigest()}


def read_json(path, expected=None):
    observed = pin(path)
    if expected is not None and observed['sha256'] != expected:
        raise ValueError('Bound successful metadata identity differs')
    return json.loads(path.read_bytes()), observed


def main():
    started = time.monotonic()
    output = RUN / 'metadata_preparation_01'
    output.mkdir(exist_ok=False)
    receipt = {'schema': 'gse205906-matched-quantifier-metadata-preparation-v1',
               'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'source_pin': pin(Path(__file__)), 'admitted': False,
               'real_alignment_count_or_model_run': False,
               'full_raw_FASTQ_read_or_download': False,
               'full_reference_fresh_hash_or_index_admitted': False}
    try:
        manifest, manifest_pin = read_json(FORMAL / 'provenance/run_manifest.json', MANIFEST_SHA)
        synthetic, synthetic_pin = read_json(RUN / 'native_count_synthetic_01/NATIVE_SYNTHETIC_COUNT_RECEIPT.json', SYNTHETIC_SHA)
        assert synthetic['accepted'] is True and len(synthetic['checks']) == 18
        assert all(x['passed'] and x['actual'] == x['expected'] for x in synthetic['checks'])
        assert all(x['actual_exit_code'] == 0 for x in synthetic['events'])
        records = manifest['metadata_rows']
        expected_labels = [f'{group}{i}' for group in ('JA', 'SNI', 'PLM') for i in (1, 2, 3)]
        assert manifest['sample_order'] == expected_labels
        assert [row['label'] for row in records] == expected_labels
        assert [row['GEO_sample'] for row in records] == [f'GSM{n}' for n in range(6235292, 6235301)]
        assert [row['SRA_run'] for row in records] == [f'SRR{n}' for n in range(19631207, 19631198, -1)]
        qc = list(csv.DictReader((FORMAL / 'output/QC_summary.tsv').open(), delimiter='\t'))
        assert [row['label'] for row in qc] == expected_labels
        libraries, full_success_pins = [], []
        total_bytes = 0
        for row, q in zip(records, qc):
            run = row['SRA_run']
            pairs = int(row['total_spots'])
            assert q['SRA_run'] == run and int(q['input_pairs']) == pairs
            success, success_pin = read_json(FORMAL / 'validation' / (run + '_full_FASTQ_and_reservoir.json'))
            assert success['status'] == 'PASS_COMPLETE_PAIRED_FASTQ_AND_FULL_FLOW_RESERVOIR'
            assert success['SRA_run'] == run and success['full_gzip_EOF_CRC_and_pairing'] is True
            assert success['total_pairs'] == pairs
            full_success_pins.append(success_pin)
            declared_bytes = list(map(int, row['ENA_generated_FASTQ_bytes'].split(';')))
            mates, quality = [], []
            for mate in (1, 2):
                raw = FORMAL / 'raw' / f'{run}_{mate}.fastq.gz'
                assert raw.is_file() and not raw.is_symlink() and raw.stat().st_size == declared_bytes[mate - 1]
                total_bytes += raw.stat().st_size
                mates.append({'path': str(raw), 'bytes': raw.stat().st_size, 'fresh_hash': False})
                archive = FORMAL / 'QC/FastQC' / run / f'{run}_{mate}_fastqc.zip'
                with zipfile.ZipFile(archive) as z:
                    names = [name for name in z.namelist() if name.endswith('/fastqc_data.txt')]
                    assert len(names) == 1
                    basic = {}
                    for line in z.read(names[0]).decode('utf-8').splitlines():
                        fields = line.split('\t')
                        if len(fields) == 2 and fields[0] in ('Filename', 'Total Sequences', 'Sequence length'):
                            basic[fields[0]] = fields[1]
                assert basic == {'Filename': raw.name, 'Total Sequences': str(pairs), 'Sequence length': '151'}
                quality.append({'archive_pin': pin(archive), 'basic_statistics': basic})
            libraries.append({'label': row['label'], 'group': row['group'],
                              'GEO_sample': row['GEO_sample'], 'SRA_run': run,
                              'registered_pairs': pairs, 'raw_mates': mates, 'FastQC': quality})
        assert total_bytes == 33141066868
        reservoir_record, reservoir_record_pin = read_json(FORMAL / 'validation/SRR19631199_full_FASTQ_and_reservoir.json')
        assert reservoir_record['sampled_pairs'] == 1000000
        assert reservoir_record['sample_original_headers_and_order_preserved'] is True
        reservoir_files = []
        expected_reservoir = FORMAL / 'derived/samples/SRR19631199'
        for name in ('SRR19631199_reservoir_1.fastq', 'SRR19631199_reservoir_2.fastq', 'sample_pair_ordinals.tsv'):
            path = expected_reservoir / name
            historical = reservoir_record['files'][str(path)]
            assert not path.is_symlink()
            fresh = pin(path)
            assert fresh['bytes'] == historical['bytes'] and fresh['sha256'] == historical['sha256']
            reservoir_files.append(fresh)
        derived, derived_pin = read_json(ASSET / 'references/derived/DERIVED_REFERENCE_156.json',
                                         '82bf95aacbb7253b671ea40972a6c659b4b7f2b89ba79e1ebf66ce9416cf7d63')
        expected_reference = [dict(manifest['reference_files']['rat']['fasta'])]
        gtf = ASSET / 'references/derived/Rattus_norvegicus.Rnor_6.0.101.gtf'
        expected_reference.append({'path': str(gtf), **derived['files'][gtf.name]})
        for ref in expected_reference:
            path = Path(ref['path'])
            assert path.is_file() and not path.is_symlink() and path.stat().st_size == ref['bytes']
        largest = max(libraries, key=lambda item: item['registered_pairs'])
        assert largest['SRA_run'] == 'SRR19631202' and largest['registered_pairs'] == 30778167
        receipt.update(admitted=True, formal_manifest_pin=manifest_pin,
                       synthetic_receipt_pin=synthetic_pin, full_FASTQ_success_record_pins=full_success_pins,
                       libraries=libraries, total_compressed_raw_bytes=total_bytes,
                       original_sample_order=expected_labels, uniform_QC_read_length=151, sjdbOverhang=150,
                       expected_reference_pins_not_freshly_hashed=expected_reference,
                       derived_reference_record_pin=derived_pin, reservoir_record_pin=reservoir_record_pin,
                       reservoir_fresh_file_pins=reservoir_files,
                       fixed_reservoir_pairs=1000000, reservoir_RNG=reservoir_record['RNG'],
                       reservoir_seed=reservoir_record['seed'],
                       largest_full_library_pilot={'label': largest['label'], 'SRA_run': largest['SRA_run'],
                                                   'registered_pairs': largest['registered_pairs']},
                       available_cloud_bytes_snapshot=shutil.disk_usage(RUN).free,
                       next_step='QUIET_WINDOW_REFERENCE_HASH_AND_INDEX_THEN_FIXED_RESERVOIR_RESOURCE_PILOT',
                       limitations='Metadata/reservoir identity only; no STAR or full-nine result acceptance. Current capacity is not a reservation.')
    except Exception as exc:
        receipt['exception_class_only'] = type(exc).__name__
    receipt.update(finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   elapsed_seconds=time.monotonic() - started)
    path = output / 'METADATA_AND_RESERVOIR_BINDING.json'
    with path.open('x', encoding='utf-8') as stream:
        json.dump(receipt, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({'admitted': receipt['admitted'], 'scope': 'metadata_and_reservoir_only',
                      'receipt_pin': pin(path), 'full_reference_or_STAR_admitted': False}))
    return 0 if receipt['admitted'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
