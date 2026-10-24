"""Dated metadata-only supplement to the existing corrected reporting ledger.

No count/expression/RDS, model, network, historical private error, or environment
is consumed. Legacy metadata and prior ledgers remain unchanged. Membership is
at submitted-array/library level, never a newly authenticated animal/culture map.
"""

import argparse
import csv
import hashlib
import io
import json
from collections import Counter
from pathlib import Path


PINS = {
    "metadata/samples.tsv": "367894147ee329498e6491c9397d1fe146ca8be35b696110964c1365da3bdf66",
    "metadata/datasets.tsv": "5a450674a5019f79f309a665e9e9e39bf292be42771512a6aeadb5bbccaa4ca4",
    "metadata/GSE128623_fixed30_manifest_20261009.tsv": "5e8bc210a6fb9c4c5c5682cdfcc907b20a6f6dce798d86c2afb6f10dd335212c",
    "metadata/dataset_reporting_completed_20261013/Route_B2_dataset_reporting_completed_20261013.tsv": "a6102e9e2e9bfde55495e824fc6860dde6a91561f900799955d325393ec7bf4c",
    "data/processed/GSE16710/qc_inclusion_decision.tsv": "816cdd2ec20f92b48410d7ce1ce516b30e6b7db27e7a51254ae7e4b65f9fad49",
    "results/diagnostics/GSE16710_original_all27_sensitivity_20261014/SCIENCE_REPORT_v2.md": "0b68e5035ed3f6868a97b1c9b835ee711e74401bd4d4303f1a9e464d40875b12",
    "results/diagnostics/GSE128623_temporal_dependence_20261014/SCIENTIFIC_REPORT_r02.md": "bcdf1444481f3e1f6140da754a181750ff4f492352ef12e03294d423de5009e6",
    "results/diagnostics/GSE277064_retained_temporal_20261014/SCIENCE_REPORT.md": "3ef04a9c3aa1b36799255beb1cbd723d4794b6a9d1400d584749fb626a31c3f9",
    "manuscript/supplement/Route_B2_GSE205906_full_read_lineage_addendum_20261014.md": "1536fc9c6507498f048371597d514d336e74cbcc485bcac89bb0b36c34736004",
    "results/diagnostics/JGW03_RIN_evidence_20261015/SCIENCE_REPORT.md": "986e1bc391e32442103d720b5c5720fa03a46e8e18fe532b739a4f35dc484af4",
    "results/master_gates/master_gate_summary.tsv": "faa4180715f0d536c992424a694d85be801bda935375f7699f759a0e346d6a47",
    "data/processed/jgw03_three_cell_lines/GSE289929/locked_samples.tsv": "4a1e97a06c692c5750be6895deb66673dcce22d464def2fd7c9f42c1c7212e21",
    "data/processed/jgw03_three_cell_lines/GSE227494/locked_samples.tsv": "e90434b1b1dcc38e61c958d9c8c18e5708847ef1d936db58750eee01b27bc583",
    "data/processed/jgw03_three_cell_lines/GSE299063/locked_samples.tsv": "18e414456468b3971e7ce0f9eaf2175bad0813d62892c0bc20236a4d07239bd6",
    "external_retained_metadata/SELECTED_36_PUBLISHED_RIN_CROSSWALK.tsv": "f29479b567d9fd88872303d9cd7d96579bf18c4f1ee308007dbf03cc55b3c7a6",
}
PRIOR = next(p for p in PINS if p.endswith("completed_20261013.tsv"))
RIN = "external_retained_metadata/SELECTED_36_PUBLISHED_RIN_CROSSWALK.tsv"
JGW = {"GSE289929": "PC12", "GSE227494": "C2C12", "GSE299063": "HepG2"}
EXPECTED = {"GSE128623": 30, "GSE16710": 27, "GSE205906": 9,
            "GSE227494": 12, "GSE277064": 63, "GSE289929": 12, "GSE299063": 12}
EVIDENCE = {
    d: next(p for p in PINS if token in p and p.endswith(".md"))
    for d, token in {"GSE16710": "original_all27", "GSE128623": "temporal_dependence",
                     "GSE277064": "retained_temporal", "GSE205906": "full_read_lineage"}.items()
}
EVIDENCE.update({d: "results/diagnostics/JGW03_RIN_evidence_20261015/SCIENCE_REPORT.md" for d in JGW})
EVIDENCE["GSE85871"] = "results/master_gates/master_gate_summary.tsv"
MEMBER_FILE = "CURRENT_ANALYSIS_MEMBERSHIP_20261015.tsv"
LEDGER_FILE = "CURRENT_REPORTING_DATASET_LEDGER_20261015.tsv"
SOURCE_FIELDS = ("sample_order", "source_column_or_member", "experimental_arm", "condition",
                 "time_day", "exercise_hour", "cell_type", "dose_level", "replicate", "batch")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def identity(content):
    return {"bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}


def table(content):
    return list(csv.DictReader(io.StringIO(content.decode("utf-8-sig")), delimiter="\t"))


def keyed(rows, fields):
    result = {tuple(r[f] for f in fields): r for r in rows}
    require(len(result) == len(rows), "Duplicate metadata identity: " + ",".join(fields))
    return result


def branch_definitions():
    definitions = {}
    for cohort, role in (("primary26", "PRIMARY"), ("all27", "SENSITIVITY")):
        for annotation in ("current", "original"):
            key = f"GSE16710_{cohort}_{annotation}"
            definitions[key] = {"n": 26 if cohort == "primary26" else 27,
                                "role": role if annotation == "current" else "SENSITIVITY",
                                "dependency_group": "GSE16710",
                                "note": "Same arrays; annotation/cohort specifications are not independent samples. QC timing is source-reported, not independently authenticated."}
    for process in ("raw", "adapter"):
        for cohort in ("all30", "predefined29"):
            definitions[f"GSE128623_{process}_{cohort}"] = {
                "n": 30 if cohort == "all30" else 29,
                "role": "PRIMARY" if (process, cohort) == ("raw", "all30") else "SENSITIVITY",
                "dependency_group": "GSE128623",
                "note": "Same submitted libraries; predefined29 omits GSM3681230/SRR8757507 only. Not four cohorts or a passed S2 gate."}
    for cell in ("v2a", "astrocyte", "microglia"):
        definitions[f"GSE277064_natural_{cell}"] = {
            "n": 15, "role": "SOURCE_WIDE_CONTEXT", "dependency_group": "GSE277064",
            "note": "Natural-course 5 time groups x 3 libraries, cell-enriched RiboTag; cell types are not paired subjects."}
        if cell != "microglia":
            definitions[f"GSE277064_rotarod_{cell}"] = {
                "n": 12, "role": "SEPARATE_TRAINING_CONTEXT", "dependency_group": "GSE277064",
                "note": "9 trained plus 3 shared untrained day14 comparators; not 12 new unique libraries."}
    for dataset, line in JGW.items():
        definitions[f"{dataset}_canonical12"] = {
            "n": 12, "role": "FORMULA_CONTEXT", "dependency_group": dataset,
            "note": f"{line}: 9 treated plus 3 own accession/canonical vehicle libraries; source replicate labels are not paired blocks or independent product lots. Low-RIN refit not run."}
    for branch in ("September_20260912", "October_Ensembl101_20261006"):
        definitions[f"GSE205906_{branch}"] = {
            "n": 9, "role": "CONSTITUENT_CONTEXT", "dependency_group": "GSE205906",
            "note": "Same nine submitted libraries; not independent cohorts, whole formula, PSS validation, or passed E2. October full-read lineage and September historical custody are distinct."}
    return definitions


def build(inputs):
    require(set(inputs) == set(PINS), "Only the explicit small metadata/report whitelist is allowed.")
    for name, expected in PINS.items():
        require(identity(inputs[name])["sha256"] == expected, "Pinned input changed: " + name)
    samples = table(inputs["metadata/samples.tsv"])
    sample_keys = keyed(samples, ("dataset_id", "sample_id"))
    require(Counter(r["dataset_id"] for r in samples) == Counter(EXPECTED), "165-member inventory changed.")
    inventory = table(inputs["metadata/datasets.tsv"])
    inventory_ids = {r["dataset_id"] for r in inventory}
    require(len(inventory_ids) == len(inventory) == 14, "Expected 14 accession/inventory records.")
    old = table(inputs[PRIOR])
    keyed(old, ("dataset_id", "reporting_branch"))
    require(len(old) == 9 and len({r["dataset_id"] for r in old}) == 8, "Nine contracts/eight accessions required.")
    qc = keyed(table(inputs["data/processed/GSE16710/qc_inclusion_decision.tsv"]), ("sample_id",))
    anchor_ids = {s for d, s in sample_keys if d == "GSE16710"}
    require({k[0] for k in qc} == anchor_ids, "Anchor adjudication membership mismatch.")
    exclusion_code = "exclude_primary_include_all_sample_sensitivity"
    excluded = {k[0] for k, r in qc.items() if r["primary_analysis_decision"] == exclusion_code}
    require(excluded == {"GSM418682"}, "Actual retained primary exclusion changed.")
    require(all(r["primary_analysis_decision"] in {"retain", exclusion_code} for r in qc.values()), "Unknown adjudication decision.")
    fixed30 = keyed(table(inputs["metadata/GSE128623_fixed30_manifest_20261009.tsv"]), ("sample_id",))
    require({k[0] for k in fixed30} == {s for d, s in sample_keys if d == "GSE128623"}, "Fixed30 membership mismatch.")
    require(fixed30[("GSM3681230",)]["run_accession"] == "SRR8757507", "Predefined29 identity changed.")
    rin = keyed(table(inputs[RIN]), ("dataset_id", "gsm"))
    jgw_keys = {k for k in sample_keys if k[0] in JGW}
    require(set(rin) == jgw_keys, "RIN sidecar must match precisely the 36 selected libraries.")
    locked = {}
    for d in JGW:
        rows = table(inputs[f"data/processed/jgw03_three_cell_lines/{d}/locked_samples.tsv"])
        rows_by_key = keyed(rows, ("dataset_id", "sample_id"))
        require(set(rows_by_key) == {k for k in jgw_keys if k[0] == d}, "Locked selection mismatch: " + d)
        locked.update(rows_by_key)
    gates = table(inputs["results/master_gates/master_gate_summary.tsv"])
    require(len(gates) == 11 and gates[0]["gate"] == "S1" and gates[0]["pass_n"] == "0", "Gate schema changed.")
    require(all(r["pass_n"] == "NA" for r in gates[1:]), "Downstream NA must not become zero or PASS.")
    members = []
    for key, source in sample_keys.items():
        dataset, gsm = key
        # Canonical JGW selection is taken from its locked table, not old flags.
        row = locked.get(key, source)
        require(all(row[f] == source[f] or {row[f], source[f]} <= {"", "NA"} for f in SOURCE_FIELDS), "Locked and inventory fields diverge: " + gsm)
        branches = []
        role = "NOT_A_TRAINING_BRANCH"
        group = row["experimental_arm"]
        time_note = "SOURCE_METADATA_LITERAL_NOT_A_NEW_TIMING_AUTHENTICATION"
        if dataset == "GSE16710":
            branches = [f"GSE16710_all27_{a}" for a in ("current", "original")]
            if gsm not in excluded:
                branches += [f"GSE16710_primary26_{a}" for a in ("current", "original")]
            group += "_day" + row["time_day"]
        elif dataset == "GSE128623":
            branches = [f"GSE128623_{p}_all30" for p in ("raw", "adapter")]
            if gsm != "GSM3681230":
                branches += [f"GSE128623_{p}_predefined29" for p in ("raw", "adapter")]
            group = "control" if row["experimental_arm"] == "control" else "D" + row["time_day"]
        elif dataset == "GSE277064":
            cell = row["cell_type"]
            if row["exercise_hour"]:
                require(cell in {"v2a", "astrocyte"}, "No trained microglia group exists.")
                branches = [f"GSE277064_rotarod_{cell}"]
                role = "TRAINED"
            else:
                branches = [f"GSE277064_natural_{cell}"]
                if row["time_day"] == "14" and cell in {"v2a", "astrocyte"}:
                    branches += [f"GSE277064_rotarod_{cell}"]
                    role = "SHARED_UNTRAINED_D14_COMPARATOR"
            group = cell + "_D" + row["time_day"] + ("_rotarod_" + row["exercise_hour"] + "h" if row["exercise_hour"] else "")
        elif dataset in JGW:
            branches = [f"{dataset}_canonical12"]
            group = row["dose_group"]
        elif dataset == "GSE205906":
            branches = [f"GSE205906_{b}" for b in ("September_20260912", "October_Ensembl101_20261006")]
            group = {"JA": "sham", "SNI": "SNI", "PLM": "PLM"}[row["experimental_arm"]]
            time_note = "LEGACY_BUILDER_REGIMEN_LITERAL_21_NOT_AUTHENTICATED_GEO_SAMPLING_DAY"
        member = {"dataset_id": dataset, "sample_id": gsm,
                  **{"source_" + f: source[f] for f in SOURCE_FIELDS},
                  "current_observation_unit": "submitted_array" if dataset == "GSE16710" else "submitted_library",
                  "current_group_label": group,
                  "membership_branch_ids": ";".join(sorted(branches)),
                  "rotarod_role": role,
                  "biological_unit_authentication": "NOT_INDEPENDENTLY_ESTABLISHED",
                  "reported_replication_boundary": "SOURCE_REPORTED_CULTURE_REPLICATION_NOT_PAIRING_OR_PRODUCT_LOT" if dataset in JGW else "ANIMAL_LINKAGE_POOLING_BATCH_NOT_AUTHENTICATED",
                  "time_field_qualification": time_note,
                  "RIN_evidence_status": rin[key]["mapping_status"] if dataset in JGW else "NOT_IN_THIS_RIN_AUDIT",
                  "RIN_refit_status": "NOT_RUN" if dataset in JGW else "NOT_APPLICABLE_TO_THIS_JGW_AUDIT",
                  "dated_evidence_reference": EVIDENCE[dataset]}
        if dataset in JGW:
            require(rin[key]["locked_dose_group"] == row["dose_group"], "RIN dose group mismatch.")
            require(rin[key]["locked_cell_line"] == JGW[dataset], "RIN cell line mismatch.")
            require(bool(rin[key]["RIN_raw"]) == (row["dose_group"] != "vehicle"), "Only 27 treatment RIN records are resolved.")
        members.append(member)
    members.sort(key=lambda r: (r["dataset_id"], int(r["source_sample_order"])))
    definitions = branch_definitions()
    branch_counts = Counter(b for r in members for b in r["membership_branch_ids"].split(";"))
    require(dict(branch_counts) == {b: v["n"] for b, v in definitions.items()}, "Branch counts mismatch.")
    require(Counter(r["current_group_label"] for r in members if "GSE16710_primary26_current" in r["membership_branch_ids"].split(";")) == Counter({"Spi_day21": 8, "Spi_day60": 8, "Sham_day21": 6, "Sham_day60": 4}), "Anchor group sizes changed.")
    require(Counter(r["rotarod_role"] for r in members if r["dataset_id"] == "GSE277064") == Counter({"TRAINED": 18, "SHARED_UNTRAINED_D14_COMPARATOR": 6, "NOT_A_TRAINING_BRANCH": 39}), "Training shared-control membership changed.")
    supplemented = []
    for oldrow in old:
        d = oldrow["dataset_id"]
        supplemented.append({**oldrow,
            "membership_file_20261015": MEMBER_FILE if d in EXPECTED else "NOT_PRESENT_IN_THIS_METADATA_TABLE",
            "submitted_member_n_20261015": EXPECTED.get(d),
            "dated_scope_evidence_20261015": EVIDENCE[d],
            "membership_qualification_20261015": "Same-source views share sample identities; branch n is not additional biological n. Legacy fields and flags are not new execution switches." if d in EXPECTED else "GSE85871 source inventory 212; selected membership unknown here, E3 not run, not 0 analyzed profiles."})
    core = {r["dataset_id"] for r in old}
    context_only = sorted(inventory_ids - core)
    require(len(context_only) == 6, "Six inventory-only context records required.")
    receipt = {
        "schema": "route-b2-analysis-membership-supplement-20261015-v1",
        "status": "METADATA_SCOPE_RECONCILED_NOT_SCIENTIFIC_OR_PUBLIC_RELEASE_PASS",
        "legacy_rows_and_values_preserved": True,
        "missingness_representation_note": "NA/empty equivalence only for locked-table comparison; output source fields preserve inventory literals. Not biological-value imputation.",
        "sample_rows": len(members), "sample_accessions": len(EXPECTED),
        "reporting_contract_rows": len(old), "core_accessions": len(core),
        "missing_membership": {"GSE85871": None},
        "inventory_only_not_reaudited": context_only,
        "analysis_branches": definitions,
        "branch_counts": dict(sorted(branch_counts.items())),
        "JGW_RIN": {"selected_libraries": 36, "resolved_treated": 27, "unresolved_vehicle": 9,
                    "new_exclusions": 0, "low_RIN_sensitivity_run": False,
                    "raw_numeric_RIN_not_republished_by_this_sidecar": True},
        "gate_states_retained": [{k: r[k] for k in ("gate", "eligible_input_n", "evaluated_n", "pass_n", "status")} for r in gates],
        "input_identities": {n: identity(b) for n, b in inputs.items()},
        "interpretation": "165 distinct submitted identifiers are not 165 authenticated independent biological units. Source-wide and nested model results cannot change the empty sequential gate families.",
        "scientific_vectors_read": False, "models_run": False,
        "input_membership_or_parameters_modified": False,
        "independent_animal_or_culture_linkage_proven": False,
        "licence_clearance_public_deposit_or_cold_restore_claimed": False,
    }
    return members, supplemented, receipt


def tsv_bytes(rows):
    handle = io.StringIO(newline="")
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return handle.getvalue().encode("utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project", type=Path)
    parser.add_argument("rin_crosswalk", type=Path)
    parser.add_argument("new_output", type=Path)
    args = parser.parse_args()
    project = args.project.resolve()
    output = args.new_output.resolve()
    require(not output.exists(), "Unused output directory required; never replace old ledgers.")
    inputs = {n: (args.rin_crosswalk if n == RIN else project / n).read_bytes() for n in PINS}
    members, ledger, receipt = build(inputs)
    payloads = {MEMBER_FILE: tsv_bytes(members), LEDGER_FILE: tsv_bytes(ledger)}
    receipt["output_identities"] = {n: identity(b) for n, b in payloads.items()}
    receipt["generator_identity"] = identity(Path(__file__).read_bytes())
    # Re-read only the same whitelist before writing. No discovered file walks.
    require(all((args.rin_crosswalk if n == RIN else project / n).read_bytes() == b for n, b in inputs.items()), "Input changed during metadata build.")
    output.mkdir(parents=True, exist_ok=False)
    for name, content in payloads.items():
        (output / name).write_bytes(content)
    (output / "MEMBERSHIP_RECEIPT.json").write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": receipt["status"], "members": len(members), "contracts": len(ledger), "output_bytes": sum(p.stat().st_size for p in output.iterdir())}))


if __name__ == "__main__":
    main()
