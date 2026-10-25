"""Pure, filesystem-free reporting codec; NOT an authorized exporter.

Callers must independently authenticate their inventory, labels and destination
bindings. This module cannot grant authority, capture files, author workbooks,
make checkpoints or adopt manuscript text. Tests use synthetic bytes only.
"""
import csv
import gzip
import hashlib
import io
import json
import re
from pathlib import PurePosixPath

HEADER = (
    'gene_id_versioned', 'gene_id_base', 'symbol', 'entrez_id', 'gene_name',
    'logFC', 'logCPM', 'F', 'PValue', 'FDR', 'signed_sqrt_F', 'ql_se',
    'ci_low', 'ci_high', 'df_total', 'contrast', 'cell_type', 'domain',
)
KEEP = tuple(range(11)) + tuple(range(14, 18))
OMIT = (11, 12, 13)
VIEW_HEADER = tuple(HEADER[i] for i in KEEP)
STATUS = 'UNVALIDATED_LEGACY_GENE_UNCERTAINTY'
SIDE_HEADER = ('source_path', 'source_sha256', 'data_row_ordinal_1based',
               'gene_id_versioned', 'status', 'ql_se', 'ci_low', 'ci_high')
BIND_KEYS = {'source_path', 'source_sha256', 'rows', 'source_cell_stream_sha256',
             'retained_cell_stream_sha256', 'labels', 'workbook', 'worksheet',
             'worksheet_member'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def frame(row):
    return (json.dumps(list(row), ensure_ascii=False, separators=(',', ':')) + '\n').encode('utf-8')


def stream_digest(header, rows):
    h = hashlib.sha256(frame(header))
    for row in rows:
        h.update(frame(row))
    return h.hexdigest()


def relative_path(value, suffix):
    require(type(value) is str and bool(value), 'Missing relative path')
    parts = value.split('/')
    require(not PurePosixPath(value).is_absolute() and '\\' not in value
            and ':' not in value and all(p not in ('', '.', '..') for p in parts)
            and not any(ord(c) < 32 for c in value) and value.endswith(suffix),
            'Noncanonical relative path')


def validate_binding(binding):
    require(type(binding) is dict and set(binding) == BIND_KEYS, 'Binding key mismatch')
    relative_path(binding['source_path'], '.tsv.gz')
    relative_path(binding['workbook'], '.xlsx')
    relative_path(binding['worksheet_member'], '.xml')
    require(binding['worksheet_member'].startswith('xl/worksheets/'), 'Wrong worksheet member')
    require(type(binding['worksheet']) is str and 0 < len(binding['worksheet']) <= 31
            and not any(c in binding['worksheet'] for c in '[]:*?/\\')
            and not any(ord(c) < 32 for c in binding['worksheet'])
            and binding['worksheet'].casefold() != 'reporting_revision', 'Invalid target sheet')
    require(type(binding['rows']) is int and binding['rows'] > 0, 'Invalid expected row count')
    for key in ('source_sha256', 'source_cell_stream_sha256', 'retained_cell_stream_sha256'):
        require(type(binding[key]) is str and re.fullmatch('[0-9a-f]{64}', binding[key]), 'Invalid digest')
    labels = binding['labels']
    require(type(labels) is dict and set(labels) == {'contrast', 'cell_type', 'domain'}
            and all(type(v) is str and v not in ('', 'NA') for v in labels.values()),
            'Explicit fixed labels required')


def read_tsv(data):
    require(type(data) is bytes, 'Immutable bytes required')
    try:
        return tuple(tuple(row) for row in csv.reader(
            io.StringIO(data.decode('utf-8-sig'), newline=''), delimiter='\t', strict=True))
    except (UnicodeError, csv.Error) as exc:
        raise ValueError('Malformed UTF-8 TSV') from exc


def write_tsv(header, rows):
    stream = io.StringIO(newline='')
    writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
    writer.writerow(header)
    writer.writerows(rows)
    return stream.getvalue().encode('utf-8')


def source_rows(blob, binding):
    validate_binding(binding)
    require(type(blob) is bytes, 'Immutable captured source bytes required')
    require(digest(blob) == binding['source_sha256'], 'Compressed source identity drift')
    try:
        table = read_tsv(gzip.decompress(blob))
    except (OSError, EOFError) as exc:
        raise ValueError('Invalid compressed source') from exc
    require(table and table[0] == HEADER, 'Full ordered 18-column schema required')
    rows = table[1:]
    require(len(rows) == binding['rows'], 'Source row count mismatch')
    keys = set()
    for row in rows:
        require(len(row) == 18, 'Malformed source row width')
        require(row[0] and row[0] not in keys, 'Missing or duplicate source gene key')
        keys.add(row[0])
        require(row[15:] == tuple(binding['labels'][k] for k in HEADER[15:]), 'Source labels differ from binding')
    require(stream_digest(HEADER, rows) == binding['source_cell_stream_sha256'], 'Full cell-stream drift')
    require(stream_digest(VIEW_HEADER, (tuple(r[i] for i in KEEP) for r in rows))
            == binding['retained_cell_stream_sha256'], 'Retained cell-stream drift')
    return rows


def schema_for(binding):
    # This records the caller's checked binding; it is not authentication of it.
    return {
        'format': 'GSE277064_NON_V9_REPORTING_CODEC_V1',
        'binding': json.loads(json.dumps(binding)),
        'status': STATUS,
        'row_ordinal_excludes_header': True,
        'columns': [dict(source_position_1based=i + 1, name=name,
                         destination='view' if i in KEEP else 'quarantine',
                         destination_position_1based=KEEP.index(i) + 1 if i in KEEP else OMIT.index(i) + 6)
                    for i, name in enumerate(HEADER)],
    }


def no_duplicate_keys(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, 'Duplicate JSON key')
        result[key] = value
    return result


def verify_projection(blob, binding, artifacts):
    """Parse actual output bytes and reconstruct every source cell independently.

    Does not regenerate the projection for comparison. Byte identity of source
    compression is bound separately; reconstruction concerns lexical cell rows,
    not source TSV quoting/newlines or gzip headers.
    """
    rows = source_rows(blob, binding)
    require(type(artifacts) is dict and set(artifacts) == {'view', 'quarantine', 'schema'}, 'Artifact set mismatch')
    require(all(type(v) is bytes for v in artifacts.values()), 'Artifact bytes required')
    try:
        schema = json.loads(artifacts['schema'].decode('utf-8'), object_pairs_hook=no_duplicate_keys)
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError('Invalid schema JSON') from exc
    require(json.dumps(schema, sort_keys=True, ensure_ascii=False, allow_nan=False)
            == json.dumps(schema_for(binding), sort_keys=True, ensure_ascii=False, allow_nan=False),
            'Schema or destination binding changed')
    view, side = read_tsv(artifacts['view']), read_tsv(artifacts['quarantine'])
    require(view and view[0] == VIEW_HEADER, 'Wrong view schema')
    require(side and side[0] == SIDE_HEADER, 'Wrong quarantine schema')
    require(len(view) == len(side) == len(rows) + 1, 'Output row loss or duplication')
    nonmissing = [0, 0, 0]
    for ordinal, (original, kept, quarantined) in enumerate(zip(rows, view[1:], side[1:]), 1):
        require(len(kept) == 15 and len(quarantined) == 8, 'Output row width mismatch')
        require(quarantined[:5] == (binding['source_path'], binding['source_sha256'],
                                   str(ordinal), original[0], STATUS), 'Sidecar provenance or row order mismatch')
        restored = [None] * 18
        for p, value in zip(KEEP, kept):
            restored[p] = value
        for j, (p, value) in enumerate(zip(OMIT, quarantined[5:])):
            restored[p] = value
            nonmissing[j] += value not in ('', 'NA')
        require(tuple(restored) == original, 'Lossy lexical reconstruction')
    return {'status': 'PASS_IN_MEMORY_LEXICAL_ROUNDTRIP_ONLY', 'rows': len(rows),
            'retained_cells': len(rows) * 15, 'nonmissing_each': nonmissing,
            'source_sha256': digest(blob),
            'artifact_sha256': {k: digest(v) for k, v in artifacts.items()}}


def project(blob, binding):
    """Return serialized derivative bytes in memory only. No IO or authority."""
    rows = source_rows(blob, binding)
    view = [tuple(row[i] for i in KEEP) for row in rows]
    side = [(binding['source_path'], binding['source_sha256'], str(n), row[0], STATUS)
            + tuple(row[i] for i in OMIT) for n, row in enumerate(rows, 1)]
    artifacts = {'view': write_tsv(VIEW_HEADER, view), 'quarantine': write_tsv(SIDE_HEADER, side),
                 'schema': (json.dumps(schema_for(binding), ensure_ascii=False, sort_keys=True,
                                       indent=2) + '\n').encode('utf-8')}
    verify_projection(blob, binding, artifacts)
    return artifacts
