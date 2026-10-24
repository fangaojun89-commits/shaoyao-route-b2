"""Prepare a private reporting fixture from exact reviewed input versions, not raw scientific data."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RETAINED = Path("R:/Shaoyao_storage/data/GSE16710_four_cell_synthesis_20261014_r01/REPORT_COMPLETE.json")
HISTORICAL = ROOT / "results/diagnostics/Route_B2_anchor_existing_specification_synthesis_20261013/METRIC_REGISTRY.json"
SOURCE = ROOT / "code/python/summarize_anchor_four_cells_20261014.py"
SCHEMA = "route-b2-four-cell-retained-result-report-20261014-v1"


def identity(path):
    payload = path.read_bytes()
    return {"bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}


def copy_input(source, dest, files, before):
    pin = identity(source)
    before[source.as_posix()] = pin
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        raise FileExistsError("Preserve existing staged member")
    shutil.copyfile(source, dest)
    if identity(dest) != pin:
        raise ValueError("Staged member differs")
    files[dest.relative_to(PACKAGE).as_posix()] = pin


def main():
    global PACKAGE
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    PACKAGE = args.output.resolve()
    if PACKAGE.exists() or not PACKAGE.parent.is_dir():
        raise ValueError("Fresh package directory required")
    retained = json.loads(RETAINED.read_text(encoding="utf-8"))
    if identity(RETAINED)["sha256"] != "7f23dcb9d7be3a0a399e688b1fde57c13271c9bd1f70c62ae4b8387c9668ffd1":
        raise ValueError("Reviewed completion changed")
    plan = json.loads(HISTORICAL.read_text(encoding="utf-8"))
    if identity(SOURCE)["sha256"] != "3ca6ed6733f9e6347a928c8f705931e032609a374e1465f863c3bfc088bb16ca":
        raise ValueError("Reviewed reporting source changed")
    layers = tuple(plan["groups"])
    specs = ("primary26_current", "primary26_original", "all27_current", "all27_original")
    if len(retained["input_identities"]) != 16:
        raise ValueError("Expected all 16 selected scientific result inputs")
    before, files, inputs = {}, {}, {}
    PACKAGE.mkdir()
    sources = ("summarize_anchor_existing_specifications_20261013.py", "anchor_four_cell_reporting_core_20261014.py",
               "replay_anchor_four_cell_portable_20261014.py")
    for source_name in sources:
        copy_input(ROOT / "code/python" / source_name, PACKAGE / "code/python" / source_name, files, before)
    for i, layer in enumerate(layers):
        inputs[layer] = {}
        for j, spec in enumerate(specs):
            pin = retained["input_identities"][4 * i + j]
            source = Path(pin["path"])
            if identity(source) != {"bytes": pin["bytes"], "sha256": pin["sha256"]}:
                raise ValueError("Scientific input identity changed")
            relative = f"data/derived/retained_results/{layer}/{spec}.tsv"
            copy_input(source, PACKAGE / relative, files, before)
            inputs[layer][spec] = {"package_path": relative, "historical_path_label_only": pin["path"],
                                  "role": "derived retained result, not raw matrix"}
    context_relatives = ("results/GSE16710/diagnostics/analysis_sample_inclusion.tsv",
                         "results/GSE16710/diagnostics/design_matrix_primary.tsv",
                         "results/GSE16710/program_level/S1_gate_all_collections.tsv",
                         "results/master_gates/master_gate_summary.tsv")
    for relative in context_relatives:
        pin = plan["inputs"][relative]
        source = ROOT / relative
        if identity(source) != {"bytes": pin["bytes"], "sha256": pin["sha256"]}:
            raise ValueError("Selected metadata context changed")
        copy_input(source, PACKAGE / "data/derived/retained_context" / relative, files, before)
    for name in ("README.md", "REPRODUCING.md"):
        copy_input(ROOT / "releases/anchor_four_cell_portable_20261014_v1" / name, PACKAGE / name, files, before)
    references = {}
    for pin in retained["output_identities"]:
        if identity(Path(pin["path"])) != {"bytes": pin["bytes"], "sha256": pin["sha256"]}:
            raise ValueError("Reviewed comparison table changed")
        references[Path(pin["path"]).name] = {"bytes": pin["bytes"], "sha256": pin["sha256"]}
    manifest = {"schema": SCHEMA, "files": files, "inputs": inputs,
                "configuration": {"specifications": specs, "layers": layers, "q_threshold": .05,
                                  "sign_zero_tolerance": 1e-12, "top_k": 20, "common_universe": "four_cell_intersection"},
                "context_root": "data/derived/retained_context",
                "utility": "code/python/summarize_anchor_existing_specifications_20261013.py",
                "core": "code/python/anchor_four_cell_reporting_core_20261014.py",
                "runner": "code/python/replay_anchor_four_cell_portable_20261014.py",
                "reference_outputs": references, "reference_table_payloads_in_package": False,
                "science_scope": "stored-result descriptor node only, not raw-to-model or full-study release",
                "scientific_git_commit": "NOT_ASSIGNED; exact member hashes identify this private source candidate",
                "sharing_status": "PRIVATE_NONCOMMERCIAL_RESEARCH; public redistribution and code licence pending",
                "historical_native_driver_exit_code": 1, "historical_metadata_successor_without_refit": True,
                "source_versions": {"reviewed_report_producer_sha256": identity(SOURCE)["sha256"],
                                    "retained_completion_sha256": identity(RETAINED)["sha256"],
                                    "historical_metric_registry_sha256": identity(HISTORICAL)["sha256"]}}
    with (PACKAGE / "MANIFEST.json").open("x", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2, ensure_ascii=False)
        stream.write("\n")
    if any(identity(Path(p)) != pin for p, pin in before.items()):
        raise ValueError("Original source/input mutation during preparation")
    print(json.dumps({"status": "PRIVATE_COMPLETE_FOUR_CELL_REPORT_FIXTURE_PREPARED", "files": len(files),
                      "input_tables": 16, "context_files": 4, "raw_or_models_copied": False,
                      "original_inputs_unchanged": True, "manifest": identity(PACKAGE / "MANIFEST.json")}))


if __name__ == "__main__":
    main()
