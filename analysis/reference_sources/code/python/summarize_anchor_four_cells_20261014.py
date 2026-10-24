"""Finite four-cell reporting; reuse pinned descriptive utilities, never fit models."""
import argparse
import importlib.util
import itertools
import json
import statistics
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
UTILITY = ROOT / "code/python/summarize_anchor_existing_specifications_20261013.py"
UTILITY_SHA = "520b45c88be410c07dbf4ee5011a07b43740aecd72a380631c5d32f45ff1d88e"
import hashlib
assert hashlib.sha256(UTILITY.read_bytes()).hexdigest() == UTILITY_SHA
module_spec = importlib.util.spec_from_file_location("retained_anchor_metrics", UTILITY)
metrics = importlib.util.module_from_spec(module_spec)
module_spec.loader.exec_module(metrics)  # utilities only; its guarded CLI is not invoked
SPECS = (*metrics.SPECS, "all27_original")
NEW_NAMES = ("reactome_camera_original_all27_sensitivity.tsv", "reactome_activity_original_all27_sensitivity.tsv",
             "progeny_original_all27_sensitivity.tsv", "dorothea_original_all27_sensitivity.tsv")


def direction(row, group):
    return (1 if row["direction"] == "Up" else -1) if group == "Reactome_CAMERA" else metrics.sign(float(row["estimate"]))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("freeze", "report"))
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    native, output = args.native.resolve(), args.output.resolve()
    registry_path = native / "REPORT_REGISTRY_v2.json"
    groups = {group: [*(ROOT / "results/GSE16710/program_level" / p for p in names), native / NEW_NAMES[i]]
              for i, (group, names) in enumerate(metrics.TABLES.items())}
    old_paths = [p for paths in groups.values() for p in paths[:3]]
    old_registry = ROOT / "results/diagnostics/Route_B2_anchor_existing_specification_synthesis_20261013/METRIC_REGISTRY.json"
    historical_pins = json.loads(old_registry.read_text(encoding="utf-8"))["inputs"]
    for p in old_paths:
        assert metrics.identity(p) == historical_pins[p.relative_to(ROOT).as_posix()]
    inputs = old_paths + [UTILITY, old_registry, native / "REGISTRY.json"]
    if args.mode == "freeze":
        plan = {"status": "POST_OUTCOME_FOUR_CELL_METRIC_PLAN_FROZEN_BEFORE_NEW_CELL",
                "utc": datetime.now(timezone.utc).isoformat(), "producer": metrics.identity(Path(__file__).resolve()),
                "specifications": SPECS, "groups": {k: [p.as_posix() for p in v] for k, v in groups.items()},
                "known_inputs": [metrics.identity(p) for p in inputs],
                "ranking": "Ascending native interaction p; average exact ties; fixed four-cell common IDs for every pair",
                "top20": "Tie-inclusive boundary at rank20; fewer than20 common IDs -> NA",
                "direction": "CAMERA Up/Down; activity native coefficient tolerance1e-12; no pooled effects",
                "coverage": "Outer union separately retains missing/not estimated; no imputation",
                "interpretation": "Native q threshold frequencies are not S1 admissions; no common-scale effect or full multiverse claim"}
        metrics.new_json(registry_path, plan)
        print(json.dumps({"status": plan["status"], "registry_sha256": metrics.identity(registry_path)["sha256"]}))
        return
    plan = json.loads(registry_path.read_text(encoding="utf-8"))
    assert plan["producer"] == metrics.identity(Path(__file__).resolve())
    assert plan["groups"] == {k: [p.as_posix() for p in v] for k, v in groups.items()}
    for pin in plan["known_inputs"]:
        assert metrics.identity(Path(pin["path"])) == pin
    metrics.self_test()
    metadata = metrics.check_metadata(ROOT)
    completion = json.loads((native / "COMPLETE.json").read_text(encoding="utf-8"))
    assert completion["status"] == "COMPLETE_EXPLORATORY_ORIGINAL_ALL27_CELL"
    assert completion["inputs_unchanged"] and completion["original26_comparisons_passed"]
    for pin in completion["output_identities"]:
        observed = metrics.identity(Path(pin["path"]))
        assert observed["sha256"] == pin["sha256"] and observed["bytes"] == pin["bytes"]
    before = [metrics.identity(p) for v in groups.values() for p in v]
    families, pairwise, curve, missing, features = [], [], [], [], []
    for group, paths in groups.items():
        maps = {}
        for spec, path in zip(SPECS, paths):
            data, selected = metrics.validate_table(path, group)
            maps[spec] = selected
            q = [float(r["q_value_BH"]) for r in selected.values()]
            families.append({"layer": group, "specification": spec, "tested": len(selected),
                             "q_le_005": sum(v <= .05 for v in q), "minimum_q": min(q), "native_table": path.as_posix()})
        common = sorted(set.intersection(*(set(m) for m in maps.values())))
        union = sorted(set.union(*(set(m) for m in maps.values())))
        assert len(common) > 1
        ranks = {s: metrics.average_ranks({k: float(m[k]["p_value"]) for k in common}) for s, m in maps.items()}
        for key in union:
            available = [s for s in SPECS if key in maps[s]]
            signs = {s: direction(maps[s][key], group) for s in available}
            primary = signs.get(SPECS[0])
            common_ranks = []
            for spec in SPECS:
                if key not in maps[spec]:
                    missing.append({"layer": group, "feature_id": key, "missing_specification": spec,
                                    "status": "NOT_ESTIMATED_NO_SOURCE_ROW"})
                    continue
                r = maps[spec][key]
                rank = ranks[spec].get(key)
                if rank is not None:
                    common_ranks.append(rank)
                curve.append({"layer": group, "feature_id": key, "specification": spec,
                              "native_estimate": r.get("estimate", "NA_NO_CAMERA_COEFFICIENT"),
                              "native_ci_low": r.get("ci_low", "NA_NO_CAMERA_CI"),
                              "native_ci_high": r.get("ci_high", "NA_NO_CAMERA_CI"),
                              "native_p": r["p_value"], "native_q": r["q_value_BH"], "native_direction": signs[spec],
                              "fixed_four_cell_common_n": len(common), "average_p_rank_common": rank if rank is not None else "NA_NOT_COMMON"})
            features.append({"layer": group, "feature_id": key, "available_specs": len(available),
                             "missing_specs": 4 - len(available),
                             "native_q_threshold_frequency": sum(float(maps[s][key]["q_value_BH"]) <= .05 for s in available) / len(available),
                             "native_sign_agreement_primary": sum(v == primary for v in signs.values()) / len(available) if primary is not None else "NA_NO_PRIMARY",
                             "median_fixed_common_p_rank": statistics.median(common_ranks) if common_ranks else "NA_NOT_COMMON",
                             "minimum_fixed_common_p_rank": min(common_ranks) if common_ranks else "NA_NOT_COMMON",
                             "maximum_fixed_common_p_rank": max(common_ranks) if common_ranks else "NA_NOT_COMMON",
                             "pooled_effect": "NA_NATIVE_SCALES_DIFFER", "gate_selection_frequency": "NA_EMPTY_S1_FAMILY"})
        for a, b in itertools.combinations(SPECS, 2):
            pa = {k: float(maps[a][k]["p_value"]) for k in common}
            pb = {k: float(maps[b][k]["p_value"]) for k in common}
            ta, tb = metrics.top20(pa), metrics.top20(pb)
            pairwise.append({"layer": group, "left": a, "right": b, "four_cell_common_n": len(common),
                             "Kendall_tau_b_p_ranks": metrics.tau_b(list(pa.values()), list(pb.values())),
                             "native_direction_agreement": sum(direction(maps[a][k], group) == direction(maps[b][k], group) for k in common) / len(common),
                             "top20_left_n": len(ta) if ta is not None else "NA_FEWER_THAN20",
                             "top20_right_n": len(tb) if tb is not None else "NA_FEWER_THAN20",
                             "top20_intersection_n": len(ta & tb) if ta is not None else "NA_FEWER_THAN20",
                             "top20_union_n": len(ta | tb) if ta is not None else "NA_FEWER_THAN20",
                             "top20_Jaccard": len(ta & tb) / len(ta | tb) if ta is not None else "NA_FEWER_THAN20"})
    output.mkdir(exist_ok=False)
    for name, records in (("FAMILY_SUMMARY.tsv", families), ("PAIRWISE_METRICS.tsv", pairwise),
                          ("CURVE_LONG.tsv", curve), ("MISSING_COVERAGE.tsv", missing), ("FEATURE_METRICS.tsv", features)):
        metrics.new_tsv(output / name, records)
    report = ["# GSE16710 four-cell exploratory sensitivity", "", "14 October 2026. AI-assisted analysis and reporting, not independent biological or human scientific validation.", "",
              "## Primary specification and executed comparison", "",
              "The primary remains current annotation/26 array observations and the injury-by-time interaction. Only the missing original GPL1355/all27 cell was newly fitted. Four original26 tables were reconstructed first and matched all scientific fields under the existing absolute1e-8 + relative1e-7 tolerance. Inputs, old results and S1 gates remained unchanged. The input is retained RMA/annotation data, not a fresh CEL-to-RMA or restored-environment run.", "",
              "## Native interaction families", "", "| Layer | Specification | Tested | q≤0.05 | Minimum q |", "|---|---|---:|---:|---:|"]
    report += [f"| {r['layer']} | {r['specification']} | {r['tested']} | {r['q_le_005']} | {r['minimum_q']:.6g} |" for r in families]
    report += ["", "## Robustness interpretation", "", "| Robustness test | Rationale | Result | Consistency | Interpretation | Impact on main conclusion |", "|---|---|---|---|---|---|",
               "| Original annotation × all27 | Complete the missing annotation/sample sensitivity cell | Full native families above and complete rank/direction comparisons below; all four cells retained | Not estimable as a confirmatory equivalence class | No common-scale equivalence criterion was specified; report the observed descriptors without inventing a pass label | unchanged: no new S1 admission |", "",
               "## All fixed-common-ID pairwise descriptors", "", "| Layer | Left | Right | Common N | Kendall tau-b | Direction agreement | Tie-inclusive top20 Jaccard |", "|---|---|---|---:|---:|---:|---:|"]
    for r in pairwise:
        tau = "NA_CONSTANT_RANKS" if r["Kendall_tau_b_p_ranks"] is None else f"{r['Kendall_tau_b_p_ranks']:.6g}"
        jaccard = f"{r['top20_Jaccard']:.6g}" if isinstance(r["top20_Jaccard"], float) else r["top20_Jaccard"]
        report.append(f"| {r['layer']} | {r['left']} | {r['right']} | {r['four_cell_common_n']} | {tau} | {r['native_direction_agreement']:.6g} | {jaccard} |")
    report += ["",
               "Ranks use a fixed four-cell common-ID universe separately for each layer, with average ranks for exact P-value ties. Top20 includes boundary ties; PROGENy has fewer than20 common IDs, so its top20 fields are NA. The outer union records absence rather than zero response. Native q-threshold frequencies are descriptive and are not gate selection frequencies. Native mean-z/ULM scores are restandardized and may have different membership; effects and intervals are preserved individually, not pooled into a common-scale effect or cross-specification CI.", "",
               "## Pipeline, files and reproduction", "", "```text", "C: existing read-only processed inputs/resources + small protocol/source", "R: pre-model registry -> original26 compatibility -> original27 tables", "   -> four-cell curve/pairwise/coverage/feature reports", "```", "",
               "The native driver reuses only named functions of the pinned analyzer. The existing isolated Windows R launcher records installed versions against renv.lock; no new renv restore occurred. The report uses bundled Python -I -B and stdlib with a pinned old metric utility, no uv packages. Actual commands and stdout/exit state are retained separately. All result tables are code-generated; old files and primary figures are untouched.", "",
               "```powershell", f"& environment/run_r.ps1 -Script code/R/run_GSE16710_original_all27_sensitivity_20261014.R prepare PROJECT_ROOT FRESH_NATIVE_R_DIRECTORY", f"& environment/run_r.ps1 -Script code/R/run_GSE16710_original_all27_sensitivity_20261014.R run PROJECT_ROOT SAME_NATIVE_R_DIRECTORY", "# Freeze descriptive metrics before fitting the new cell; run report only after native completion.", "python -I -B code/python/summarize_anchor_four_cells_20261014.py freeze --native SAME_NATIVE_R_DIRECTORY --output FRESH_REPORT_R_DIRECTORY", "python -I -B code/python/summarize_anchor_four_cells_20261014.py report --native SAME_NATIVE_R_DIRECTORY --output FRESH_REPORT_R_DIRECTORY", "```", "",
               "For an actual fresh replay, prepare native registry and freeze descriptive registry before run; do not rerun into an occupied directory. Source hashes, exact inputs, current package versions, seed, locale and session are recorded. The current installed environment is not an independently restored cold runtime.", "",
               "## Unresolved risks and evidence boundary", "",
               "This is post-outcome, finite four-cell sensitivity work, not the whole ten-factor multiverse. Primary S1 remains 1,771 evaluated/0 admitted and downstream candidate families remain empty. Animal linkage, outcome-protection chronology, full raw preprocessing/dependency restoration, RNA-seq label-null/alternate quantifier coverage, true PSS holdout, exposure/mediator links and public redistribution remain unresolved. Planned unexecuted GSE277064 alternate fits are not results. No manuscript or gate was silently rewritten."]
    report_path = ROOT / "results/diagnostics/GSE16710_original_all27_sensitivity_20261014/SCIENCE_REPORT.md"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("x", encoding="utf-8") as stream:
        stream.write("\n".join(report) + "\n")
    after = [metrics.identity(p) for v in groups.values() for p in v]
    assert before == after
    result = {"status": "FOUR_CELL_DESCRIPTIVE_REPORT_COMPLETE", "registry": metrics.identity(registry_path),
              "native_completion": metrics.identity(native / "COMPLETE.json"), "metadata": metadata,
              "input_identities": before, "family_rows": len(families), "pairwise_rows": len(pairwise),
              "curve_rows": len(curve), "missing_rows": len(missing), "feature_rows": len(features),
              "output_identities": [metrics.identity(p) for p in sorted(output.iterdir())], "report": metrics.identity(report_path),
              "all_result_inputs_unchanged": True, "whole_multiverse_complete": False, "whole_goal_complete": False}
    metrics.new_json(output / "REPORT_COMPLETE.json", result)
    print(json.dumps({"status": result["status"], "families": len(families), "pairs": len(pairwise), "curve_rows": len(curve)}))


if __name__ == "__main__":
    main()
