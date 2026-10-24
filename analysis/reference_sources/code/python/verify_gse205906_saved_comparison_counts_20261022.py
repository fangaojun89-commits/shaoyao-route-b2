"""Independent streaming reconciliation of saved comparison exports; no fits."""
import collections
import csv
import hashlib
import json
import math
import pathlib
import sys


def main(root):
    root = pathlib.Path(root)
    source = root / 'full_source_union.tsv'
    summary = root / 'native_table_summary.tsv'
    groups = collections.defaultdict(collections.Counter)
    seen = set()
    with source.open(newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = tuple(row[x] for x in ('table', 'collection', 'contrast'))
            identity = key + (row['feature_id'],)
            assert identity not in seen, 'DUPLICATE_FAMILY_FEATURE'
            seen.add(identity)
            count = groups[key]
            count['union_domain_n'] += 1
            hits = []
            for side in ('left', 'right'):
                for field, target in [('domain_present', 'domain_n'),
                                      ('eligible', 'eligible_n'),
                                      ('native_present', 'returned_n'),
                                      ('native_not_returned', 'not_returned_n')]:
                    assert row[side+'_'+field] in ('TRUE', 'FALSE')
                    count[side+'_'+target] += row[side+'_'+field] == 'TRUE'
                q = row[side+'_native_q']
                present = q != 'NA'
                if present:
                    q = float(q)
                    assert math.isfinite(q) and 0 <= q <= 1, 'INVALID_Q'
                count[side+'_q_available_n'] += present
                hit = present and q < .05
                hits.append(hit)
                count[side+'_q_strict_lt_05_n'] += hit
                count[side+'_q_exactly_05_n'] += present and q == .05
            count['q_hit_intersection_n'] += all(hits)
            count['q_hit_union_n'] += any(hits)
            if row['paired_estimable_continuous'] == 'TRUE':
                count['paired_estimable_continuous_n'] += 1
                a, b = (float(row[s+'_effect_or_NES']) for s in ('left','right'))
                assert math.isfinite(a) and math.isfinite(b)
                nonzero = a != 0 or b != 0
                both = a != 0 and b != 0
                same = a*b > 0
                count['nonzero_union_n'] += nonzero
                count['both_nonzero_n'] += both
                count['same_sign_nonzero_union_n'] += nonzero and same
                count['same_sign_both_nonzero_n'] += both and same
                count['both_zero_n'] += not nonzero
    checked = 0
    summary_keys = set()
    with summary.open(newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = tuple(row[x] for x in ('table','collection','contrast'))
            assert key not in summary_keys, 'DUPLICATE_SUMMARY_FAMILY'
            summary_keys.add(key)
            count = groups[key]
            for field in ('union_domain_n','left_domain_n','right_domain_n',
                          'left_eligible_n','right_eligible_n','left_returned_n',
                          'right_returned_n','left_not_returned_n','right_not_returned_n',
                          'left_q_available_n','right_q_available_n',
                          'left_q_strict_lt_05_n','right_q_strict_lt_05_n',
                          'left_q_exactly_05_n','right_q_exactly_05_n',
                          'q_hit_intersection_n','q_hit_union_n',
                          'paired_estimable_continuous_n','nonzero_union_n',
                          'both_nonzero_n','same_sign_nonzero_union_n',
                          'same_sign_both_nonzero_n','both_zero_n'):
                assert int(row[field]) == count[field], (key, field)
                checked += 1
            union = count['q_hit_union_n']
            assert (row['q_hit_Jaccard'] == 'NA' if not union else
                    abs(float(row['q_hit_Jaccard'])-count['q_hit_intersection_n']/union) < 1e-12)
            checked += 1
    assert summary_keys == set(groups) and len(groups) == 22
    return dict(status='PASS_SAVED_EXPORT_COUNT_RECONCILIATION', families=len(groups),
                unique_family_features=len(seen), assertions=checked,
                iterative_fits=0, rows_removed=0,
                scope='Export counts and hit/sign partitions, not full-chain numerical acceptance',
                inputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in (source,summary)})


if __name__ == '__main__':
    print(json.dumps(main(sys.argv[1]),sort_keys=True))
