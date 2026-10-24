"""Post-outcome descriptive synthesis of named, retained GSE16710 tables.

No model, BH adjustment, raw matrix, acquisition, historical error or gate is
executed. Freeze the explicit file/metric plan before computing new summaries.
"""
import argparse
import csv
import hashlib
import itertools
import json
import math
import statistics
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


CONTRAST = "spasticity_progression_interaction"
SPECS = ("primary26_current", "primary26_original", "all27_current")
TABLES = {
    "Reactome_CAMERA": (
        "reactome_camera_primary_all_contrasts.tsv",
        "reactome_camera_original_annotation_sensitivity.tsv",
        "reactome_camera_all27_sensitivity.tsv",
    ),
    "Reactome_mean_z": (
        "reactome_mean_z_activity_all_contrasts.tsv",
        "reactome_activity_original_annotation_sensitivity.tsv",
        "reactome_activity_all27_sensitivity.tsv",
    ),
    "PROGENy_ULM": (
        "progeny_ulm_activity_all_contrasts.tsv",
        "progeny_original_annotation_sensitivity.tsv",
        "progeny_all27_sensitivity.tsv",
    ),
    "DoRothEA_ULM": (
        "dorothea_ulm_activity_all_contrasts.tsv",
        "dorothea_original_annotation_sensitivity.tsv",
        "dorothea_all27_sensitivity.tsv",
    ),
}
REFERENCES = (
    "results/GSE16710/program_level/S1_gate_all_collections.tsv",
    "results/GSE16710/program_level/S1_summary_counts.tsv",
    "results/GSE16710/diagnostics/analysis_sample_inclusion.tsv",
    "results/GSE16710/diagnostics/design_matrix_primary.tsv",
    "results/GSE16710/diagnostics/model_estimability.tsv",
    "results/GSE16710/diagnostics/contrast_matrix.tsv",
    "results/master_gates/master_gate_summary.tsv",
    "code/R/analyze_GSE16710_spasticity_anchor.R",
    "config/contrasts.yaml", "config/thresholds.yaml",
    "protocol/statistical_analysis_plan.md",
    "results/diagnostics/Route_B2_anchor_existing_specification_synthesis_20261013/METHOD_SCOPE_REVIEW.md",
)


def identity(path):
    payload = path.read_bytes()
    return {"path": path.as_posix(), "bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest()}


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def new_json(path, value):
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")


def new_tsv(path, data, columns=None):
    if columns is None:
        columns = list(data[0])
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


def sign(value):
    return 1 if value > 1e-12 else (-1 if value < -1e-12 else 0)


def average_ranks(values):
    ordered = sorted(values, key=lambda key: (values[key], key))
    ranks, start = {}, 0
    while start < len(ordered):
        end = start + 1
        while end < len(ordered) and values[ordered[end]] == values[ordered[start]]:
            end += 1
        for key in ordered[start:end]:
            ranks[key] = (start + 1 + end) / 2
        start = end
    return ranks


def tau_b(x, y):
    if len(x) != len(y):
        raise ValueError("Tau length mismatch")
    concordant = discordant = ties_x = ties_y = 0
    for i in range(len(x)):
        for j in range(i):
            sx = (x[i] > x[j]) - (x[i] < x[j])
            sy = (y[i] > y[j]) - (y[i] < y[j])
            ties_x += sx == 0
            ties_y += sy == 0
            concordant += sx * sy > 0
            discordant += sx * sy < 0
    pairs = len(x) * (len(x) - 1) // 2
    denominator = math.sqrt((pairs - ties_x) * (pairs - ties_y))
    return (concordant - discordant) / denominator if denominator else None


def top20(values):
    if len(values) < 20:
        return None
    boundary = sorted(values.values())[19]
    return {key for key, value in values.items() if value <= boundary}


def self_test():
    assert tau_b([1, 2, 3], [1, 2, 3]) == 1
    assert tau_b([1, 2, 3], [3, 2, 1]) == -1
    assert tau_b([1, 2, 2, 4], [2, 1, 1, 4]) == .2
    assert tau_b([1, 1], [1, 1]) is None
    assert average_ranks({"a": 2, "b": 1, "c": 2}) == {"b": 1, "a": 2.5, "c": 2.5}
    assert top20({str(i): i for i in range(14)}) is None
    assert len(top20({str(i): min(i, 19) for i in range(21)})) == 21
    assert sign(0) == 0 and sign(1e-13) == 0


def freeze(root, registry):
    paths = ["results/GSE16710/program_level/" + name for names in TABLES.values() for name in names]
    paths += list(REFERENCES)
    plan = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Post-outcome descriptive synthesis of 12 retained tables; no statistical refit or gate rerun",
        "root": root.as_posix(), "contrast": CONTRAST,
        "primary_specification": SPECS[0], "specifications": list(SPECS),
        "unavailable_crossed_specification": "all27_original_annotation_NOT_RECONSTRUCTED",
        "groups": TABLES, "inputs": {p: identity(root / p) for p in paths},
        "producer": identity(Path(__file__).resolve()),
        "q_rule": "stored table q_value_BH <= 0.05; not S1 selection frequency",
        "ranking": "ascending stored p_value; average ranks for ties; long table keeps full-spec ranks, summary ranks use fixed three-spec common universe",
        "pairwise_comparison": "same group/layer only, fixed three-spec common feature intersection for all pairs, explicit nonshared coverage",
        "top_k": 20, "top_k_ties": "include all at the 20th p-value boundary; fewer than 20 -> NA",
        "direction_rule": "CAMERA Up/Down; native activity sign with abs <= 1e-12 classified zero",
        "effect_reporting": "per-specification native estimates/CIs only; no pooled/median effect, no cross-spec 95% CI",
        "candidate_metrics": "NOT_IDENTIFIABLE_EMPTY_S1_FAMILY; no retrospective gate admission",
        "interpretation": "same arrays and altered score definitions; not independent validation or complete factorial multiverse",
        "original_files_modified": False, "models_run": False,
    }
    registry.parent.mkdir(parents=True, exist_ok=True)
    new_json(registry, plan)
    print(json.dumps({"status": "DESCRIPTIVE_METRIC_PLAN_FROZEN", "input_files": len(paths),
                      "registry": identity(registry)}))


def check_metadata(root):
    inclusion = rows(root / REFERENCES[2])
    ids = [r["sample_id"] for r in inclusion]
    if len(ids) != 27 or len(set(ids)) != 27:
        raise ValueError("Retained array ledger uniqueness/count mismatch")
    kept = [r["sample_id"] for r in inclusion if r["primary_analysis_decision"] == "retain"]
    excluded = [r["sample_id"] for r in inclusion if r["primary_analysis_decision"] != "retain"]
    design = rows(root / REFERENCES[3])
    if len(kept) != 26 or excluded != ["GSM418682"] or [r["sample_id"] for r in design] != kept:
        raise ValueError("Linked primary array/design membership mismatch")
    groups = Counter((r["experimental_arm"], r["time_day"]) for r in inclusion)
    if set(groups) != {("Sham", "21"), ("Sham", "60"), ("Spi", "21"), ("Spi", "60")}:
        raise ValueError("Unexpected injury-time groups")
    s1 = rows(root / REFERENCES[0])
    if len(s1) != 1771 or len({(r["collection"], r["feature_id"]) for r in s1}) != 1771:
        raise ValueError("S1 family mismatch")
    if any(r["S1_pass"].upper() != "FALSE" for r in s1):
        raise ValueError("Frozen S1 state differs from expected empty candidate family")
    master = [r for r in rows(root / REFERENCES[6]) if r["gate"] == "S1"]
    if len(master) != 1 or master[0]["evaluated_n"] != "1771" or master[0]["pass_n"] != "0":
        raise ValueError("Master gate state mismatch")
    return {"all_array_rows": 27, "primary_array_rows": 26, "primary_excluded_ids": excluded,
            "all_array_groups": {"/".join(k): v for k, v in sorted(groups.items())},
            "s1_evaluated": 1771, "s1_admitted": 0,
            "unit": "array observations; legacy independent_rat labels are not adopted as verified animal linkage",
            "qc_chronology": "retained source decision, not new verification of pre-outcome timing",
            "matrix_alignment_scope": "linked metadata and saved design only; no raw matrix read"}


def validate_table(path, group):
    data = rows(path)
    camera = group == "Reactome_CAMERA"
    required = {"feature_id", "contrast", "p_value", "q_value_BH"}
    required |= {"direction", "n_genes"} if camera else {"estimate", "ci_low", "ci_high", "standard_error", "moderated_t", "df_total"}
    if not data or not required.issubset(data[0]):
        raise ValueError("Missing result schema")
    keys = [(r["feature_id"], r["contrast"]) for r in data]
    if len(keys) != len(set(keys)) or any(not k[0] or not k[1] for k in keys):
        raise ValueError("Result duplicate/missing key")
    for row in data:
        p, q = float(row["p_value"]), float(row["q_value_BH"])
        if not (math.isfinite(p) and math.isfinite(q) and 0 <= p <= 1 and 0 <= q <= 1):
            raise ValueError("Result probability out of range")
        if camera:
            n = float(row["n_genes"])
            if row["direction"] not in {"Up", "Down"} or not n.is_integer() or n <= 0:
                raise ValueError("CAMERA direction/membership invalid")
        else:
            values = {k: float(row[k]) for k in required - {"feature_id", "contrast"}}
            if not all(math.isfinite(v) for v in values.values()):
                raise ValueError("Nonfinite native result")
            if not (values["ci_low"] <= values["estimate"] <= values["ci_high"]):
                raise ValueError("CI/estimate ordering mismatch")
            if values["standard_error"] <= 0 or values["df_total"] <= 0:
                raise ValueError("Nonpositive native SE/df")
    selected = {r["feature_id"]: r for r in data if r["contrast"] == CONTRAST}
    if not selected:
        raise ValueError("Missing retained interaction")
    return data, selected


def summarize(root, registry, output, report):
    self_test()
    plan = json.loads(registry.read_text(encoding="utf-8"))
    if plan["root"] != root.as_posix() or plan["groups"] != {k: list(v) for k, v in TABLES.items()}:
        raise ValueError("Frozen plan mismatch")
    if identity(Path(__file__).resolve()) != plan["producer"]:
        raise ValueError("Producer changed after freeze")
    for p, pin in plan["inputs"].items():
        if identity(root / p) != pin:
            raise ValueError("Source changed after freeze")
    metadata = check_metadata(root)
    tables, validation, audit, spec_summary = {}, [], [], []
    for group, names in TABLES.items():
        for spec, name in zip(SPECS, names):
            rel = "results/GSE16710/program_level/" + name
            data, selected = validate_table(root / rel, group)
            tables[group, spec] = selected
            q = [float(r["q_value_BH"]) for r in selected.values()]
            signs = [1 if r["direction"] == "Up" else -1 for r in selected.values()] if group == "Reactome_CAMERA" else [sign(float(r["estimate"])) for r in selected.values()]
            spec_summary.append({"group": group, "specification": spec, "source_rows_retained": len(data),
                                 "interaction_feature_n": len(selected), "stored_q_le_005_n": sum(v <= .05 for v in q),
                                 "min_stored_q": min(q), "positive_direction_n": signs.count(1),
                                 "negative_direction_n": signs.count(-1), "zero_direction_n": signs.count(0),
                                 "source_file": rel, "candidate_gate_status": "NOT_RERUN_EMPTY_S1_FAMILY"})
            validation.append({"file": rel, "rows_before": len(data), "rows_retained_original": len(data),
                               "rows_removed": 0, "duplicate_keys": 0, "invalid_numeric_or_ci_rows": 0,
                               "interaction_features": len(selected), "status": "PASS_RESULT_FIELD_CHECKS_ONLY"})
    curve, feature_summary, coverage, pairwise = [], [], [], []
    for group in TABLES:
        maps = {spec: tables[group, spec] for spec in SPECS}
        ids = sorted(set().union(*(set(m) for m in maps.values())))
        common_universe = sorted(set.intersection(*(set(m) for m in maps.values())))
        ranks = {spec: average_ranks({key: float(r["p_value"]) for key, r in m.items()}) for spec, m in maps.items()}
        common_ranks = {spec: average_ranks({key: float(m[key]["p_value"]) for key in common_universe}) for spec, m in maps.items()}
        for key in ids:
            available = [spec for spec in SPECS if key in maps[spec]]
            native_signs = {spec: (1 if maps[spec][key]["direction"] == "Up" else -1) if group == "Reactome_CAMERA" else sign(float(maps[spec][key]["estimate"])) for spec in available}
            primary_sign = native_signs.get(SPECS[0])
            normalized_ranks = []
            for spec in SPECS:
                if key not in maps[spec]:
                    coverage.append({"group": group, "feature_id": key, "absent_specification": spec,
                                     "status": "NOT_ESTIMATED_NO_SOURCE_ROW", "imputed": False})
                    continue
                r = maps[spec][key]
                fraction = (ranks[spec][key] - 1) / (len(maps[spec]) - 1)
                common_fraction = (common_ranks[spec][key] - 1) / (len(common_universe) - 1) if key in common_universe else None
                if common_fraction is not None:
                    normalized_ranks.append(common_fraction)
                curve.append({"group": group, "feature_id": key, "specification": spec,
                              "stored_p": r["p_value"], "stored_q": r["q_value_BH"],
                              "native_direction": native_signs[spec], "stored_estimate": r.get("estimate", "NA_NO_CAMERA_EFFECT"),
                              "stored_ci_low": r.get("ci_low", "NA_NO_CAMERA_CI"), "stored_ci_high": r.get("ci_high", "NA_NO_CAMERA_CI"),
                              "p_value_average_rank": ranks[spec][key], "full_spec_feature_n": len(maps[spec]),
                              "full_spec_rank_fraction": fraction,
                              "fixed_three_spec_common_n": len(common_universe),
                              "common_universe_rank_fraction": common_fraction if common_fraction is not None else "NA_NOT_COMMON_TO_THREE",
                              "role": "POST_OUTCOME_SOURCE_WIDE_NOT_GATE_ELIGIBLE"})
            matches = sum(v == primary_sign for v in native_signs.values()) if primary_sign is not None else None
            passes = sum(float(maps[s][key]["q_value_BH"]) <= .05 for s in available)
            feature_summary.append({"group": group, "feature_id": key, "valid_retained_specification_n": len(available),
                                    "missing_retained_specification_n": len(SPECS) - len(available), "primary_present": primary_sign is not None,
                                    "table_q_le_005_n": passes, "table_q_le_005_frequency": passes / len(available),
                                    "native_direction_matches_primary_n": matches if matches is not None else "NA_NO_PRIMARY_ROW",
                                    "native_direction_comparison_n": len(available) if matches is not None else 0,
                                    "native_direction_matches_primary_fraction": matches / len(available) if matches is not None else "NA_NO_PRIMARY_ROW",
                                    "in_fixed_three_spec_common_universe": key in common_universe,
                                    "fixed_three_spec_common_n": len(common_universe),
                                    "median_common_rank_fraction": statistics.median(normalized_ranks) if normalized_ranks else "NA_NOT_COMMON_TO_THREE",
                                    "min_common_rank_fraction": min(normalized_ranks) if normalized_ranks else "NA_NOT_COMMON_TO_THREE",
                                    "max_common_rank_fraction": max(normalized_ranks) if normalized_ranks else "NA_NOT_COMMON_TO_THREE",
                                    "pooled_effect": "NA_NATIVE_SCORE_NOT_COMMON_FIXED_SCALE",
                                    "candidate_selection_frequency": "NA_EMPTY_S1_FAMILY"})
        for left, right in itertools.combinations(SPECS, 2):
            a, b = maps[left], maps[right]
            common = common_universe
            p_a = {key: float(a[key]["p_value"]) for key in common}
            p_b = {key: float(b[key]["p_value"]) for key in common}
            tau = tau_b(list(p_a.values()), list(p_b.values()))
            sa = [1 if a[key]["direction"] == "Up" else -1 for key in common] if group == "Reactome_CAMERA" else [sign(float(a[key]["estimate"])) for key in common]
            sb = [1 if b[key]["direction"] == "Up" else -1 for key in common] if group == "Reactome_CAMERA" else [sign(float(b[key]["estimate"])) for key in common]
            ta, tb = top20(p_a), top20(p_b)
            pairwise.append({"group": group, "left_specification": left, "right_specification": right,
                             "left_feature_n": len(a), "right_feature_n": len(b), "common_feature_n": len(common),
                             "left_only_feature_n": len(set(a) - set(b)), "right_only_feature_n": len(set(b) - set(a)),
                             "left_excluded_from_fixed_common_n": len(a) - len(common),
                             "right_excluded_from_fixed_common_n": len(b) - len(common),
                             "kendall_tau_b_p_ranks": tau if tau is not None else "NA_CONSTANT_RANKS",
                             "native_direction_agreement_n": sum(x == y for x, y in zip(sa, sb)),
                             "native_direction_agreement_fraction": sum(x == y for x, y in zip(sa, sb)) / len(common),
                             "top_k": 20, "top_k_left_actual_n": len(ta) if ta is not None else "NA_LESS_THAN_20_COMMON",
                             "top_k_right_actual_n": len(tb) if tb is not None else "NA_LESS_THAN_20_COMMON",
                             "top_k_union_n": len(ta | tb) if ta is not None and tb is not None else "NA_LESS_THAN_20_COMMON",
                             "top_k_intersection_n": len(ta & tb) if ta is not None and tb is not None else "NA_LESS_THAN_20_COMMON",
                             "top_k_tie_inclusive_Jaccard": len(ta & tb) / len(ta | tb) if ta is not None and tb is not None else "NA_LESS_THAN_20_COMMON",
                             "interpretation": "shared-feature descriptive comparison, not external replication"})
    output.mkdir(parents=True, exist_ok=False)
    for name, records in (("specification_summary.tsv", spec_summary), ("specification_curve_long.tsv", curve),
                          ("feature_coverage_metrics.tsv", feature_summary), ("coverage_missing_rows.tsv", coverage),
                          ("pairwise_metrics.tsv", pairwise), ("result_field_validation.tsv", validation)):
        new_tsv(output / name, records)
    for rel, pin in plan["inputs"].items():
        audit.append({"file_role": "derived_result" if rel.startswith("results/GSE16710/program_level/") else "linked_context_not_raw_matrix",
                      "expected_path": (root / rel).as_posix(), "actual_path_exists": (root / rel).as_posix(),
                      "bytes": pin["bytes"], "sha256": pin["sha256"], "source_accession": "GSE16710 linked context",
                      "version_download_date": "retained analysis output; original download date not reauthenticated",
                      "license_data_use": "noncommercial academic study; public redistribution not cleared here",
                      "sample_metadata_match": "saved design/ledger checked; original matrix and animal linkage not reviewed",
                      "inclusion_exclusion": "26 primary/27 sensitivity; exclude GSM418682 primary only",
                      "script_reference": "code/R/analyze_GSE16710_spasticity_anchor.R named output references",
                      "duplicate_stale_check": "explicit named path, not glob/latest; frozen pre-summary identity",
                      "manifest_match": "matches new post-outcome metric registry, not proof of historical execution",
                      "status": "NEEDS REVIEW for full historical provenance; input identity checked",
                      "discrepancy": "unit linkage/historical chronology/public rights outside this scope",
                      "required_action": "retain qualification; no input substitution or automatic source modification"})
    new_tsv(report / "FILE_LEVEL_INPUT_AUDIT.tsv", audit)
    summary = {"status": "DESCRIPTIVE_EXISTING_RESULT_SYNTHESIS_COMPLETE",
               "overall_provenance_status": "NEEDS REVIEW", "validation_status": "PASS_RETAINED_RESULT_FIELDS_AND_METADATA_ONLY",
               "registry": identity(registry), "producer": identity(Path(__file__).resolve()),
               "metadata": metadata, "source_specification_rows": spec_summary,
               "long_curve_rows": len(curve), "feature_metric_rows": len(feature_summary),
               "missing_coverage_rows": len(coverage), "pairwise_rows": len(pairwise),
               "coverage_missingness_is_not_zero_response": True,
               "models_run": False, "gates_rerun": False, "input_rows_deleted": 0,
               "full_multiverse_completed": False, "whole_goal_completed": False,
               "warnings": ["post-outcome reporting rules, not pre-outcome statistical plan",
                            "score units/membership vary; no pooled effect or cross-specification confidence interval",
                            "same source arrays, unresolved animal linkage and retrospective QC chronology",
                            "missing all27/original crossed cell; no label/null/alternate-transfer/quantifier outputs supplied"],
               "outputs_on_R": {p.name: identity(p) for p in sorted(output.iterdir())},
               "all_bound_inputs_unchanged_after": all(identity(root / p) == pin for p, pin in plan["inputs"].items())}
    if not summary["all_bound_inputs_unchanged_after"]:
        raise ValueError("Source drift during descriptive synthesis")
    new_json(report / "SYNTHESIS_SUMMARY.json", summary)
    print(json.dumps({"status": summary["status"], "long_rows": len(curve), "feature_rows": len(feature_summary),
                      "coverage_missing_rows": len(coverage), "pairwise_rows": len(pairwise),
                      "sources_unchanged": summary["all_bound_inputs_unchanged_after"]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("freeze", "summarize", "self-test"))
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--registry", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    self_test()
    if args.mode == "freeze":
        freeze(args.root.resolve(), args.registry.resolve())
    elif args.mode == "summarize":
        summarize(args.root.resolve(), args.registry.resolve(), args.output.resolve(), args.report.resolve())
    else:
        print("DESCRIPTIVE_HELPER_SYNTHETIC_TESTS_PASS")
