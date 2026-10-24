"""Portable descriptive projection of two accepted family tables; stdlib only.

Does not regenerate native saved objects, infer P values, adjust BH or fit models.
Input content identities are fixed; their locations are caller supplied.
"""
import argparse
import csv
import hashlib
import io
import json
import platform
from pathlib import Path

PINS = {
    'g128': ('81bd5bd9fb018e714b06fcf1e7c11cd167865dccfd8e98db44c3285270cf0512', 448),
    'jgw': ('829f3ea818cd9e3b10a5f4ee6c6804ae27fd5bc980aa752be78218c21836645a', 123),
}

def require(ok, message):
    if not ok:
        raise ValueError(message)

def load(path, role):
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    require(digest == PINS[role][0], role + ': accepted input identity differs')
    rows = list(csv.DictReader(io.StringIO(raw.decode('utf-8-sig')), delimiter='\t'))
    require(len(rows) == PINS[role][1], role + ': family coverage differs')
    return rows, {'sha256': digest, 'bytes': len(raw), 'rows': len(rows)}

def tsv(rows):
    buf = io.StringIO(newline='')
    w = csv.DictWriter(buf, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
    w.writeheader(); w.writerows(rows)
    return buf.getvalue().encode()

def run(g128, jgw, output):
    require(not output.exists(), 'Fresh output directory required')
    a, ai = load(g128, 'g128'); b, bi = load(jgw, 'jgw')
    summaries = []
    for cell, space in sorted({(r['cell_id'], r['reference_space']) for r in a}):
        group = [r for r in a if (r['cell_id'], r['reference_space']) == (cell, space)]
        require(len({r['family_id'] for r in group}) == len(group), 'Duplicate G128 family')
        summaries.append(dict(cell_id=cell, reference_space=space, family_n=len(group),
            families_with_saved_empirical_q_le_005=sum(int(r['empirical_q_le_05_n']) > 0 for r in group),
            feature_family_records=sum(int(r['n']) for r in group),
            unavailable_feature_family_records=sum(int(r['empirical_unavailable_n']) for r in group)))
    jrows = []
    for cell, space in sorted({(r['accession'], r['space']) for r in b}):
        group = [r for r in b if (r['accession'], r['space']) == (cell, space)]
        require(len({r['family_id'] for r in group}) == len(group), 'Duplicate JGW family')
        total, nonid = (10001, 10000) if space == 'four_group_MC' else (20, 19)
        require(all(int(r['reference_N']) == total and int(r['nonidentity_N']) == nonid for r in group), 'Reference denominator differs')
        jrows.append(dict(accession=cell, space=space, family_n=len(group),
            feature_family_records=sum(int(r['feature_n']) for r in group),
            conditional_P_unavailable_n=sum(int(r['conditional_P_unavailable_n']) for r in group),
            null_unavailable_feature_draws=sum(int(r['null_unavailable_feature_draws']) for r in group),
            reference_N=total, nonidentity_N=nonid))
    require(len(summaries) == 24 and len(jrows) == 6, 'Incomplete reporting groups')
    require(sum(r['feature_family_records'] for r in jrows) == 484723, 'JGW feature coverage')
    files = {'G128_GROUP_SUMMARY.tsv': tsv(summaries), 'JGW_GROUP_SUMMARY.tsv': tsv(jrows)}
    files['READING_NOTES.md'] = (
        '# Accepted saved reference tables: descriptive projection\n\n'
        'G128 counts use stored empirical q<=0.05 fields within their original families. '
        'JGW summaries retain conditional-reference missingness and fixed denominators; '
        'they are not native or empirical-BH discovery counts. Feature-family records are '
        'not unique genes or independent samples. The branches must not be pooled.\n\n'
        'This reproduces descriptive summaries from accepted table exports, not R objects, '
        'original archive arithmetic, stochastic allocations or raw-to-model computation. '
        'No new biological validation, S1 admission or antispasticity mechanism is established.\n'
    ).encode()
    # Confirm input stability before publishing generated output.
    require(load(g128, 'g128')[1] == ai and load(jgw, 'jgw')[1] == bi, 'Input changed')
    output.mkdir(parents=False)
    for name, data in files.items():
        with (output / name).open('xb') as f: f.write(data)
    receipt = {'status':'PASS_ACCEPTED_TABLE_DESCRIPTIVE_PROJECTION_ONLY', 'inputs':{'g128':ai,'jgw':bi},
        'outputs':{n:{'bytes':len(v),'sha256':hashlib.sha256(v).hexdigest()} for n,v in files.items()},
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'python':platform.python_version(), 'models_run':0,'P_q_recalculated':False,
        'native_object_reconstruction':False}
    with (output/'RECEIPT.json').open('x',encoding='utf-8') as f: json.dump(receipt,f,indent=2)
    return receipt

if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--g128',type=Path,required=True);p.add_argument('--jgw',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    print(json.dumps(run(args.g128,args.jgw,args.output),indent=2))
