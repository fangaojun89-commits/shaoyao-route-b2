"""Exact conditional references for the accepted pooled-nine rank statistic.

Standard library only. No count/RDS, model fitting, normalization, RNG or gate
mutation. Biological exchangeability is assumed, not authenticated here.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path
import platform
import sys

GROUPS = ("sham", "SNI", "PLM")
OBSERVED = (0, 0, 0, 1, 1, 1, 2, 2, 2)
SAMPLES = tuple("GSM" + str(i) for i in range(6235292, 6235301))
CONTRASTS = {"local_disease": (1, 0), "intervention": (2, 1)}
PINS = {
    "scores.tsv": "29c095651f42498b45be273526e0ffea23eb6b41c0b66ca97a363e3d824339ac",
    "metadata.tsv": "5f69dde8c63c9a2d395af6b578e117843e899e60cf41ba94829cb321c3c130f1",
    "BASELINES.tsv": "ee191aac799f84388154bdc71b5558ca3a46a776667c5d6b8bd3801dd3290feb",
    "RECOVERY_COMPLETE.json": "d98e815bde1f5ef157e5761d523e3766795cac7af069a69b5299f7dc342ca408",
    "BASELINE_COMPLETE.json": "9d19ac0cc742df4fe5c8e03c96a882800239f88c7e6d1b198c418f8ae3aaee5e",
    "INDEPENDENT_COUNT_SCORE_BRIDGE.json": "2e60d2d30e2b7d2c3b4967bdd34be07ef30f2b455a73cba9989d4cd7387c5b9c",
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    raw = Path(path).read_bytes()
    return {"path": str(Path(path).resolve()), "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames
        require(fields is not None and len(set(fields)) == len(fields), "Unique header required")
        return fields, list(reader)


def write_tsv(path, rows, fields=None):
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def twice_ranks(values):
    require(len(values) > 0 and all(v.is_finite() for v in values), "Finite scores required")
    order = sorted(range(len(values)), key=lambda i: values[i])
    result = [None] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        # Positions start+1 ... end: twice the average is start+1+end.
        for index in order[start:end]:
            result[index] = start + 1 + end
        start = end
    return tuple(result)


def all_allocations():
    indices = tuple(range(9))
    result = []
    for sham in combinations(indices, 3):
        remaining = tuple(i for i in indices if i not in sham)
        for sni in combinations(remaining, 3):
            result.append(tuple(0 if i in sham else 1 if i in sni else 2 for i in indices))
    require(len(result) == len(set(result)) == 1680 and result.count(OBSERVED) == 1, "Complete unique orbit")
    return result


def numerator(ranks, allocation, contrast):
    pos, neg = CONTRASTS[contrast]
    return sum(r for r, g in zip(ranks, allocation) if g == pos) - sum(r for r, g in zip(ranks, allocation) if g == neg)


def tail_probabilities(statistics):
    histogram = Counter(abs(s) for s in statistics)
    count = 0
    tails = {}
    for value in sorted(histogram, reverse=True):
        count += histogram[value]
        tails[value] = Fraction(count, len(statistics))
    return [tails[abs(s)] for s in statistics]


def exact_bh(values):
    n = len(values)
    order = sorted(range(n), key=lambda i: values[i])
    result = [None] * n
    running = Fraction(1)
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        running = min(running, values[i] * n / rank)
        result[i] = running
    return result


def frac_fields(prefix, value):
    return {prefix + "_numerator": value.numerator, prefix + "_denominator": value.denominator,
            prefix + "_display": format(float(value), ".17g")}


def selftest():
    allocations = all_allocations()
    require(all(Counter(a) == Counter({0: 3, 1: 3, 2: 3}) for a in allocations), "Group-size preservation")
    require(twice_ranks([Decimal(x) for x in (1, 1, 4)]) == (3, 3, 6), "Average ties")
    require(twice_ranks([Decimal(1)] * 9) == (10,) * 9, "Constant ranks")
    constant = [numerator((10,) * 9, a, "local_disease") for a in allocations]
    require(set(tail_probabilities(constant)) == {Fraction(1)}, "Constant inclusive P=1")
    values = [Decimal(x) for x in (1, 2, 9, 3, 4, 5, 6, 7, 8)]
    ranks = twice_ranks(values)
    require(numerator(ranks, OBSERVED, "local_disease") == 0, "Exact zero")
    require(twice_ranks([v * 3 + 7 for v in values]) == ranks, "Positive affine rank invariance")
    require(twice_ranks([-v for v in values]) == tuple(20 - r for r in ranks), "Sign reversal")
    for name, (pos, neg) in CONTRASTS.items():
        third = ({0, 1, 2} - {pos, neg}).pop()
        subset = [a for a in allocations if all((a[i] == third) == (OBSERVED[i] == third) for i in range(9))]
        require(len(subset) == 20 and subset.count(OBSERVED) == 1, "Companion complete orbit")
        stats = [numerator(ranks, a, name) for a in allocations]
        require(Counter(stats) == Counter(-v for v in stats), "Sign complement symmetry")
        require(min(tail_probabilities(stats)) >= Fraction(1, 840), "Global two-sided resolution")
        require(min(tail_probabilities([numerator(ranks, a, name) for a in subset])) >= Fraction(1, 10), "Pairwise two-sided resolution")
    require(exact_bh([Fraction(1, 100), Fraction(1, 10), Fraction(1)]) ==
            [Fraction(3, 100), Fraction(3, 20), Fraction(1)], "Exact rational BH")
    require(tail_probabilities([-2, 0, 2, 2]) == [Fraction(3, 4), Fraction(1), Fraction(3, 4), Fraction(3, 4)], "Inclusive tied tails")
    rejected = False
    try:
        twice_ranks([Decimal("NaN")])
    except ValueError:
        rejected = True
    require(rejected, "Nonfinite rejection")
    return {"status": "PASS_SYNTHETIC_ONLY", "checks": ["1680_unique_identity_group_sizes", "average_ties",
        "constant_P1", "integer_exact_zero", "affine_invariance", "sign_reversal", "two_companion20",
        "sign_complement", "two_sided_resolution", "exact_BH", "inclusive_tails", "nonfinite_rejection"]}


def run(source, output, project):
    require(not output.exists() and str(output.resolve()).startswith("R:\\Shaoyao_storage\\"), "Unused output on R required")
    records = {name: pin(source / name) for name in PINS}
    require(all(records[name]["sha256"] == sha for name, sha in PINS.items()), "Accepted source bytes changed")
    protocol = project / "protocol/GSE205906_POOLED_RANK_LABEL_REFERENCE_20261015.md"
    protocol_record = pin(protocol)
    fields, score_rows = read_tsv(source / "scores.tsv")
    _, metadata = read_tsv(source / "metadata.tsv")
    _, old_rows = read_tsv(source / "BASELINES.tsv")
    require(fields == ["collection", "program_id", *SAMPLES], "Exact sample order")
    require([r["sample_id"] for r in metadata] == list(SAMPLES) and
            tuple(GROUPS.index(r["group"]) for r in metadata) == OBSERVED, "Fixed group alignment")
    keys = [(r["collection"], r["program_id"]) for r in score_rows]
    require(len(keys) == len(set(keys)) == 95 and Counter(c for c, _ in keys) == Counter({"PROGENy": 14, "DoRothEA_AB": 81}), "Complete native programs")
    old = {(r["collection"], r["program_id"], r["contrast"]): r for r in old_rows}
    require(len(old) == len(old_rows) == 190 and set(old) == {(c, p, k) for c, p in keys for k in CONTRASTS}, "Complete parent contrast keys")
    rank_rows = [twice_ranks([Decimal(r[s]) for s in SAMPLES]) for r in score_rows]
    tied_n = [sum(Counter(Decimal(r[s]) for s in SAMPLES)[Decimal(r[s])] > 1 for s in SAMPLES) for r in score_rows]
    allocations = all_allocations()
    observed_id = allocations.index(OBSERVED)
    stats = {(contrast, *key): [numerator(rank_rows[j], a, contrast) for a in allocations]
             for contrast in CONTRASTS for j, key in enumerate(keys)}
    require(all(Counter(v) == Counter(-x for x in v) for v in stats.values()), "Every sign-complement distribution")
    zero_n = 0
    for j, (collection, program) in enumerate(keys):
        for contrast in CONTRASTS:
            effect = Fraction(stats[contrast, collection, program][observed_id], 6)
            parent = old[collection, program, contrast]
            require(abs(Decimal(parent["rank_effect"]) - Decimal(effect.numerator) / Decimal(effect.denominator)) < Decimal("1e-13"), "Observed parent rank-effect export mismatch")
            require(int(parent["rank_tied_observations"]) == tied_n[j], "Parent ties mismatch")
            zero_n += effect == 0
    require(zero_n == 6, "All six exact-zero observed rows retained")
    features, families = [], []
    for contrast, (pos, neg) in CONTRASTS.items():
        third = ({0, 1, 2} - {pos, neg}).pop()
        spaces = {"global": list(range(1680)), "pairwise": [i for i, a in enumerate(allocations)
                  if all((a[s] == third) == (OBSERVED[s] == third) for s in range(9))]}
        require(len(spaces["pairwise"]) == 20, "Companion cardinality")
        for space, selected in spaces.items():
            n = len(selected)
            local_obs = selected.index(observed_id)
            for collection in ("PROGENy", "DoRothEA_AB"):
                selected_keys = [key for key in keys if key[0] == collection]
                p_rows = [tail_probabilities([stats[contrast, *key][i] for i in selected]) for key in selected_keys]
                q_columns = [exact_bh([p[i] for p in p_rows]) for i in range(n)]
                tail_counts = [sum(p[i] <= Fraction(1, 20) for p in p_rows) for i in range(n)]
                bh_counts = [sum(q < Fraction(1, 20) for q in col) for col in q_columns]
                for j, (c, program) in enumerate(selected_keys):
                    key_index = keys.index((c, program))
                    effect_num = stats[contrast, c, program][observed_id]
                    p = p_rows[j][local_obs]
                    q = q_columns[local_obs][j]
                    features.append({"space": space, "collection": c, "program_id": program,
                        "contrast": contrast, "library_n": 9, "rank_effect_numerator": effect_num,
                        "rank_effect_denominator": 6, "rank_effect_display": format(effect_num / 6, ".17g"),
                        "rank_tied_observations": tied_n[key_index], "rank_exact_zero": effect_num == 0,
                        "reference_N": n, "tail_count": int(p * n),
                        **frac_fields("exact_P", p), **frac_fields("exact_BH", q),
                        "rank_CI_status": "NOT_ESTIMATED_NO_EFFECT_MODEL_INVERSION",
                        "interpretation": "CONDITIONAL_EXCHANGEABILITY_REFERENCE_NOT_BIOLOGICAL_CALIBRATION"})
                sorted_tails = sorted(tail_counts)
                median_count = Fraction(sorted_tails[(n - 1) // 2] + sorted_tails[n // 2], 2)
                families.append({"space": space, "collection": collection, "contrast": contrast,
                    "program_n": len(selected_keys), "allocation_n": n,
                    "observed_P_le_005_n": tail_counts[local_obs],
                    "observed_BH_lt_005_n": bh_counts[local_obs],
                    "observed_min_P": str(min(p[local_obs] for p in p_rows)),
                    "observed_min_BH": str(min(q_columns[local_obs])),
                    "reference_min_P": str(min(min(p) for p in p_rows)),
                    "allocations_any_BH_lt_005": sum(x > 0 for x in bh_counts),
                    "fraction_allocations_any_BH": str(Fraction(sum(x > 0 for x in bh_counts), n)),
                    "reference_tail05_min_fraction": str(Fraction(min(tail_counts), len(selected_keys))),
                    "reference_tail05_median_fraction": str(median_count / len(selected_keys)),
                    "reference_tail05_max_fraction": str(Fraction(max(tail_counts), len(selected_keys))),
                    "tail05_count_histogram": json.dumps(dict(sorted(Counter(tail_counts).items()))),
                    "BH_count_histogram": json.dumps(dict(sorted(Counter(bh_counts).items())))})
    require(len(features) == 380 and len(families) == 8, "Complete output families")
    require(all(Fraction(r["exact_P_numerator"], r["exact_P_denominator"]) >= Fraction(1, 10)
                for r in features if r["space"] == "pairwise"), "Pairwise resolution")
    require(all(pin(source / n) == record for n, record in records.items()) and pin(protocol) == protocol_record, "Source preserved")
    output.mkdir(parents=True, exist_ok=False)
    statistic_columns = {f"{contrast}::{c}::{p}": (contrast, c, p) for contrast in CONTRASTS for c, p in keys}
    allocation_rows = ({"allocation_id": i + 1, "labels": "".join(map(str, a)),
                        **{column: stats[key][i] for column, key in statistic_columns.items()}}
                       for i, a in enumerate(allocations))
    write_tsv(output / "ALLOCATION_STATISTICS.tsv", allocation_rows, ["allocation_id", "labels", *statistic_columns])
    write_tsv(output / "FEATURE_RANK_REFERENCE.tsv", features)
    write_tsv(output / "FAMILY_RANK_REFERENCE.tsv", families)
    receipt = {"status": "ENUMERATED_PENDING_INDEPENDENT_CHECK", "utc": datetime.now(timezone.utc).isoformat(),
        "builder": pin(__file__), "protocol": protocol_record, "inputs": records,
        "outputs": {n: pin(output / n) for n in ("ALLOCATION_STATISTICS.tsv", "FEATURE_RANK_REFERENCE.tsv", "FAMILY_RANK_REFERENCE.tsv")},
        "runtime": {"python": sys.version, "executable": sys.executable, "platform": platform.platform()},
        "program_n": 95, "sample_n": 9, "score_values": 855, "global_N": 1680,
        "pairwise_N_each_contrast": 20, "observed_allocation_id": observed_id + 1,
        "feature_rows": 380, "family_rows": 8, "integer_statistic_cells": 319200,
        "observed_exact_zero_contrast_rows": zero_n, "exhaustive_including_identity": True,
        "models_run": False, "normalization_or_score_recovery_run": False, "RDS_or_counts_read": False,
        "parent_descriptive_rank_fields_modified": False, "gate_upgrade": False,
        "biological_exchangeability_authenticated": False, "whole_study_complete": False}
    with (output / "RANK_REFERENCE_COMPLETE.json").open("x", encoding="utf-8") as stream:
        json.dump(receipt, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"status": receipt["status"], "feature_rows": 380, "integer_statistic_cells": 319200,
                      "output_bytes": sum(p.stat().st_size for p in output.iterdir())}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("selftest", "run"))
    parser.add_argument("--source", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    if args.mode == "selftest":
        print(json.dumps(selftest()))
    else:
        require(args.source is not None and args.output is not None, "Source and output required")
        run(args.source.resolve(), args.output.resolve(), args.project.resolve())


if __name__ == "__main__":
    main()
