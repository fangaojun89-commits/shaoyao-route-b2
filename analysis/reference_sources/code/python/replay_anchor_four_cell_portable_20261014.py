"""Portable retained-result reporting node; stdlib only, never a scientific model launcher."""
import argparse
import hashlib
import importlib.util
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

SCHEMA = "route-b2-four-cell-retained-result-report-20261014-v1"


def identity(path):
    payload = path.read_bytes()
    return {"bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}


def utc():
    return datetime.now(timezone.utc).isoformat()


def safe_path(root, relative):
    rel = PurePosixPath(relative)
    if rel.is_absolute() or not rel.parts or any(p in {"", ".", ".."} for p in rel.parts) or ":" in relative or "\\" in relative:
        raise ValueError("Nonportable manifest path")
    path = root.joinpath(*rel.parts)
    if not path.resolve().is_relative_to(root):
        raise ValueError("Manifest path escapes the supplied package root")
    return path


def load(root):
    manifest_path = root / "MANIFEST.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest["schema"] != SCHEMA:
        raise ValueError("Unexpected package schema")
    observed = {}
    for relative, pin in manifest["files"].items():
        path = safe_path(root, relative)
        observed[relative] = identity(path)
        if observed[relative] != pin:
            raise ValueError("Changed package member: " + relative)
    actual = {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file()}
    if actual != set(manifest["files"]) | {"MANIFEST.json"}:
        raise ValueError("Package file set differs from frozen manifest")
    if Path(__file__).resolve() != safe_path(root, manifest["runner"]).resolve():
        raise ValueError("Execute the runner from the supplied isolated package")
    return manifest, observed, identity(manifest_path)


def import_source(root, relative, name):
    spec = importlib.util.spec_from_file_location(name, safe_path(root, relative))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def runtime():
    executable = Path(sys.executable).resolve()
    return {"python_version": sys.version, "executable": executable.as_posix(),
            "executable_identity": identity(executable), "platform": platform.platform(),
            "implementation": platform.python_implementation(), "isolated": bool(sys.flags.isolated),
            "no_bytecode": bool(sys.dont_write_bytecode), "randomness": "none"}


def new_json(path, value):
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


def replay(root, out):
    manifest, before, manifest_pin = load(root)
    if not sys.flags.isolated or not sys.dont_write_bytecode:
        raise ValueError("Use Python -I -B for the isolated deterministic node")
    if out.exists() or not out.parent.is_dir() or out.is_relative_to(root) or root.is_relative_to(out):
        raise ValueError("Supply a fresh output sibling with existing parent, separate from the package")
    metrics = import_source(root, manifest["utility"], "frozen_metric_utilities")
    core = import_source(root, manifest["core"], "frozen_four_cell_compute")
    config = manifest["configuration"]
    if (tuple(config["specifications"]) != core.SPECS or config["q_threshold"] != .05 or
            config["top_k"] != 20 or config["sign_zero_tolerance"] != 1e-12 or config["common_universe"] != "four_cell_intersection"):
        raise ValueError("Configuration differs from the fixed numerical node")
    metrics.self_test()
    metadata = metrics.check_metadata(safe_path(root, manifest["context_root"]))
    groups = {layer: [safe_path(root, manifest["inputs"][layer][s]["package_path"]) for s in core.SPECS]
              for layer in config["layers"]}
    labels = {layer: {s: manifest["inputs"][layer][s]["historical_path_label_only"] for s in core.SPECS}
              for layer in config["layers"]}
    out.mkdir()
    started = utc()
    environment = runtime()
    new_json(out / "RUN_STARTED.json", {"start_utc": started, "package_manifest": manifest_pin,
             "command": sys.argv, "cwd": Path.cwd().as_posix(), "runtime": environment,
             "node": "stored-result-validation-to-four-cell-descriptors", "models_run": False})
    tables = core.compute(groups, labels, metrics)
    if set(tables) != set(manifest["reference_outputs"]):
        raise ValueError("Output-table set differs from the reviewed node")
    for name, records in tables.items():
        metrics.new_tsv(out / name, records)
    actual = {name: identity(out / name) for name in tables}
    if actual != manifest["reference_outputs"]:
        raise ValueError("Regenerated complete tables differ from reviewed outputs")
    _, after, manifest_after = load(root)
    if before != after or manifest_pin != manifest_after:
        raise ValueError("Package/input mutation during replay")
    audit = []
    for layer, selections in manifest["inputs"].items():
        for spec, item in selections.items():
            p = safe_path(root, item["package_path"])
            audit.append({"File / role": "derived retained result: " + layer + "/" + spec,
                "Expected path": item["package_path"], "Actual path / exists": p.as_posix(),
                "Size / hash": str(before[item["package_path"]]["bytes"]) + " / " + before[item["package_path"]]["sha256"],
                "Version / source / accession / download date": "GSE16710 retained analysis; original download not reauthenticated",
                "License / data-use": "private noncommercial research replay; public sharing not cleared",
                "Sample ID and metadata-matrix match": "27/26 saved array/design context checked; no raw matrix or animal mapping reverified",
                "Inclusion/exclusion counts": "26 primary; GSM418682 primary-only exclusion, all27 sensitivity",
                "Notebook/script reference": manifest["core"] + ": compute; supplied relative manifest mapping",
                "Duplicate/stale/superseded check": "exact named version; bytes/hash equal selected reviewed input",
                "Manifest match": "exact", "Status": "NEEDS REVIEW for wider provenance; selected-file identity matched",
                "Discrepancy": "animal linkage, historical chronology, raw download/custody and public rights remain unresolved",
                "Required action": "retain limits; do not infer full-study or raw-to-model acceptance"})
    metrics.new_tsv(out / "FILE_LEVEL_AUDIT.tsv", audit)
    receipt = {"schema": SCHEMA, "status": "PASS_ALL_FIVE_COMPLETE_TABLES_BYTE_IDENTICAL",
               "start_utc": started, "end_utc": utc(), "command": sys.argv, "cwd": Path.cwd().as_posix(),
               "runtime": environment, "package_manifest": manifest_pin, "package_root": root.as_posix(),
               "complete_output_identities": actual, "table_rows": {k: len(v) for k, v in tables.items()},
               "metadata": metadata, "package_inputs_unchanged": True, "old_project_root_required": False,
               "historical_absolute_paths": "display-only frozen provenance labels, never input locations",
               "scientific_models_run": False, "reference_tables_available_to_generator": False,
               "raw_to_model_restoration_verified": False, "full_study_or_public_release_complete": False,
               "audit_overall_status": "NEEDS REVIEW"}
    new_json(out / "REPLAY_COMPLETE.json", receipt)
    checkpoint = {"package_manifest": manifest_pin, "runtime": environment,
                  "files": {p.name: identity(p) for p in sorted(out.iterdir()) if p.is_file()}}
    new_json(out / "CHECKPOINT.json", checkpoint)
    return {"status": receipt["status"], "table_rows": receipt["table_rows"], "models_run": False}


def check(root, out):
    manifest, _, manifest_pin = load(root)
    checkpoint = json.loads((out / "CHECKPOINT.json").read_text(encoding="utf-8"))
    if checkpoint["package_manifest"] != manifest_pin or checkpoint["runtime"] != runtime():
        raise ValueError("Source or runtime drift invalidates checkpoint reuse")
    for name, pin in checkpoint["files"].items():
        if identity(safe_path(out, name)) != pin:
            raise ValueError("Checkpoint output changed")
    actual = {name: identity(out / name) for name in manifest["reference_outputs"]}
    if actual != manifest["reference_outputs"]:
        raise ValueError("Checkpoint scientific table mismatch")
    return {"status": "PASS_READ_ONLY_CHECKPOINT_AND_ALL_FIVE_TABLES", "models_or_report_rerun": False}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("run", "check"))
    parser.add_argument("--package-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root, out = args.package_root.resolve(), args.output.resolve()
    result = replay(root, out) if args.mode == "run" else check(root, out)
    print(json.dumps(result, sort_keys=True))
