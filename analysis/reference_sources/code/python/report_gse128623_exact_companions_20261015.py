"""Report the complete saved-only exact-companion acceptance; no R or fit.

Consumes explicit small receipts/summary only. Requires the actual reader exit
and all 20 complete spaces; never substitutes this report for six-group work.
"""
import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import platform
import sys


STATUS = "PASS_ALL_20_EXACT_COMPANION_ARITHMETIC_AND_LINEAGE_ONLY"
ENTRY_SHA = "7f028ce8243fc9c3a2e7405eb2e5af98488ebd185202f8910bfedf3210baf192"
CONTRACT_SHA = "d6be674ea0e80d41c6ac60af725577391b63b0579976668a4b5e10ca78da1d02"
CELLS = ("raw_all30", "raw_predefined29", "adapter_all30", "adapter_predefined29")
DAYS = ("D4", "D7", "D14", "D28", "D42")
FAMILIES = (("edgeR", "GENE"), ("voom", "GENE"), ("DESeq2", "GENE"),
            ("activity", "DoRothEA_AB"), ("activity", "PROGENy"),
            ("activity", "Reactome"), ("activity", "SynGO_assignedBy_proxy_sensitivity"),
            ("CAMERA", "Reactome"), ("CAMERA", "SynGO_assignedBy_proxy_sensitivity"))
FIELDS = ("method", "collection", "contrast", "family_id", "n", "actual_BH_adjustment_n",
          "cell_id", "reference_space", "required_allocation_n", "empirical_unavailable_n",
          "empirical_q_le_05_n", "native_fraction_allocations_any_q_le_05")


def need(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve(strict=True)
    raw = path.read_bytes()
    return {"path": str(path), "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def integer(value):
    result = int(value)
    need(str(result) == value and result >= 0, "Expected a nonnegative integer field")
    return result


def validate(receipt, process, rows):
    need(receipt["status"] == STATUS and process["actual_R_exit_code"] == 0 and
         process["native_status"] == STATUS and process["entry_unchanged"] is True and
         process["contract_unchanged"] is True, "Actual successful reader exit required")
    for field, value in (("job_n", 20), ("archive_n", 4788), ("family_slots", 180),
                         ("scientific_model_calls", 0)):
        need(type(receipt[field]) is int and receipt[field] == value and
             process[field] == value, "Receipt/process complete-scope count differs")
    for field in ("selected_114_actual_shards_exit_zero", "all_selected_metadata_and_sources_unchanged",
                  "source_native_observed_repack_fields_match", "archive_bytes_authenticated_each_arithmetic_pass",
                  "producer_independent_complete_tails_and_BH_agree"):
        need(receipt[field] is True, "Required complete arithmetic/identity acceptance missing")
    for field in ("original_results_modified", "exact_companions_receive_outer_MC_intervals",
                  "full_24_job_finalization_invoked", "six_group_reference_accepted",
                  "full_native_fields_or_final_RNG_consumption_independently_reestimated",
                  "biological_exchangeability_authenticated", "candidate_gate_changed", "whole_study_complete"):
        need(receipt[field] is False, "Saved-only/exact-space scope was upgraded")
    need(receipt["entry"]["sha256"] == ENTRY_SHA and
         receipt["caller_contract"]["sha256"] == CONTRACT_SHA, "Reviewed entry/contract differs")
    expected_jobs = {cell + "/" + day + "_vs_control_exact" for cell in CELLS for day in DAYS}
    need(set(receipt["jobs"]) == expected_jobs and len(rows) == 180, "Require all original exact spaces")
    expected_keys = {(cell, day + "_vs_control_exact", method, collection)
                     for cell in CELLS for day in DAYS for method, collection in FAMILIES}
    seen = set()
    for row in rows:
        key = (row["cell_id"], row["reference_space"], row["method"], row["collection"])
        need(key in expected_keys and key not in seen, "Missing/extra/duplicate family key")
        seen.add(key)
        cell, space, method, collection = key
        target = space.removesuffix("_exact")
        total = 126 if cell.endswith("predefined29") and space == "D4_vs_control_exact" else 252
        job = receipt["jobs"][cell + "/" + space]
        need(job["allocations"] == total and job["families"] == 9 and
             job["producer_independent_tails_and_BH_agree"] is True and
             integer(row["required_allocation_n"]) == total and row["contrast"] == target and
             row["family_id"] == "::".join((method, collection, target)), "Family/label-space context differs")
        n, valid, missing, hits = (integer(row[k]) for k in
                                 ("n", "actual_BH_adjustment_n", "empirical_unavailable_n", "empirical_q_le_05_n"))
        need(valid + missing == n and hits <= valid, "Available/unavailable/reference counts disagree")
        fraction = Decimal(row["native_fraction_allocations_any_q_le_05"])
        need(fraction.is_finite() and 0 <= fraction <= 1, "Invalid native allocation fraction")
    need(seen == expected_keys, "Incomplete fixed family domain")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for field in ("receipt", "receipt-sha256", "process", "process-sha256", "summary", "output"):
        parser.add_argument("--" + field, required=True, type=str)
    args = parser.parse_args()
    receipt_path, process_path, summary_path = (Path(p).resolve(strict=True) for p in
                                               (args.receipt, args.process, args.summary))
    output = Path(args.output).resolve()
    companion = output.with_suffix(".json")
    need(output.suffix == ".md" and output.parent.is_dir() and
         not output.exists() and not companion.exists(), "Use a fresh explicit report/receipt output")
    need(output not in (receipt_path, process_path, summary_path) and
         companion not in (receipt_path, process_path, summary_path), "Output overlaps an input")
    inputs = [pin(p) for p in (receipt_path, process_path, summary_path, Path(__file__))]
    need(inputs[0]["sha256"] == args.receipt_sha256 and
         inputs[1]["sha256"] == args.process_sha256, "Caller receipt/process identity differs")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    process = json.loads(process_path.read_text(encoding="utf-8"))
    need(process["native_receipt"]["sha256"] == inputs[0]["sha256"] and
         process["native_receipt"]["bytes"] == inputs[0]["bytes"], "Actual exit belongs to another receipt")
    recorded = [r for r in receipt["output_records"]
                if Path(r["path"]).name == "FAMILY_REFERENCE_SUMMARY.tsv"]
    need(len(recorded) == 1 and recorded[0]["bytes"] == inputs[2]["bytes"] and
         recorded[0]["sha256"] == inputs[2]["sha256"], "Full native family summary differs")
    with summary_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        need(reader.fieldnames == list(FIELDS), "Exact source family-summary schema required")
        rows = list(reader)
    validate(receipt, process, rows)
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["method"], row["collection"]].append(row)
    lines = ["# GSE128623: complete exact day/control companion references", "",
             "Overall status: PASS for all complete exact-space archive arithmetic and lineage. "
             "These post-outcome conditional references do not authenticate animal independence "
             "or biological exchangeability; the six-group references remain separate.", "",
             "## Complete scope and interpretation", "",
             "All 20 processing/cohort/day-space tasks, 4,788 archives and 180 method/collection "
             "family slots passed producer/independent complete-tail and BH comparison after "
             "actual reader exit 0. No new model was fitted. Each exhaustive denominator includes "
             "the observed allocation once: 252, except predefined29 D4/control at 126. "
             "Unavailable observed/required-null statistics remain unavailable; no denominator "
             "was shortened. The reported BH adjustment count is the actual nonmissing P count.", "",
             "Original gene/program effects and their documented uncertainty remain in the "
             "retained observed results. This report adds label-reference counts, not replacement "
             "effect estimates, effect confidence intervals, mroast/fgsea day analyses or outer "
             "Monte Carlo intervals. No result reopens S1 or the candidate chain.", "",
             "## Descriptive comparison across the complete prespecified spaces", "",
             "Each row below covers the same 20 specification/space slots, not 20 independent "
             "samples. 'Any q' counts slots containing at least one empirical q<=0.05 within "
             "their original family; reference spaces and families were never pooled for BH. "
             "Ranges describe counts of features, not confidence intervals. Read hit counts "
             "together with unavailable counts; absence of crossings is not equivalence or "
             "evidence of biological absence.", "",
             "| Method | Collection | Slots | Slots with any empirical q<=0.05 | Hit-count min | Hit-count max | Slots with unavailable features |",
             "|---|---|---:|---:|---:|---:|---:|"]
    aggregate = []
    for method, collection in FAMILIES:
        block = grouped[method, collection]
        need(len(block) == 20, "Complete cross-specification layer required")
        counts = [integer(row["empirical_q_le_05_n"]) for row in block]
        vector = (20, sum(n > 0 for n in counts), min(counts), max(counts),
                  sum(integer(row["empirical_unavailable_n"]) > 0 for row in block))
        lines.append("| " + " | ".join(map(str, (method, collection, *vector))) + " |")
        aggregate.append({"method": method, "collection": collection, "values": list(vector)})
    lines += ["", "## Complete family-level results", "",
              f"The full [180-row native summary]({summary_path.as_posix()}) retains every "
              "processing/cohort/space/method/collection key, total features, actual BH count, "
              "unavailable count, empirical hit count and native fraction of allocations with "
              "any q<=0.05. No source row is omitted. The complete tail/native-diagnostic RDS "
              "and lineage tables remain at their cloud locations pinned in the receipt.", "",
              "## Quantitative evidence inventory", "",
              "All rows below are exploratory computational claims, not confirmatory biological claims.", "",
              "| Claim | Reported value | Source file | Source location/object | Verification status | Discrepancy | Required action |",
              "|---|---|---|---|---|---|---|"]
    for field, value in (("job_n", 20), ("archive_n", 4788), ("family_slots", 180), ("scientific_model_calls", 0)):
        lines.append(f"| Complete exact-space {field} | {value} | {receipt_path.as_posix()} | {field} | VERIFIED | None | Retain exact-space scope |")
    for item in aggregate:
        selector = f"method={item['method']}; collection={item['collection']}; all 20 fixed cell/space rows"
        values = "; ".join(f"{name}={value}" for name, value in zip(
            ("slots", "any_q_slots", "hit_min", "hit_max", "unavailable_slots"), item["values"]))
        lines.append(f"| Displayed complete-layer count vector | {values} | {summary_path.as_posix()} | {selector} | VERIFIED | None | Do not treat slots as biological N |")
    lines += ["", "## Source identities and execution", "",
              f"Native acceptance SHA-256: `{inputs[0]['sha256']}`. Actual process receipt "
              f"SHA-256: `{inputs[1]['sha256']}`. Complete family-summary SHA-256: "
              f"`{inputs[2]['sha256']}`. Reviewed saved-only reader: `{ENTRY_SHA}`. "
              f"Caller contract: `{CONTRACT_SHA}`.", "",
              "The exact source arguments/runtime and report identity are recorded in the "
              "adjacent report JSON. This stdlib-only presentation step reads no RDS, counts, "
              "raw reads or original error strings; it installs nothing and starts no R/model "
              "or cloud process. It is not environment restoration or whole-study acceptance.", "",
              "Unresolved items: sampled six-group completion and its full retained-native "
              "acceptance remain pending. Native stochastic algorithm replay and biological "
              "exchangeability are not established by these saved-field checks.", "",
              "Blocking issues: none for the complete exact-space arithmetic report. "
              "Protected phenotype validation, formulation/exposure and candidate boundaries "
              "retain their original status.", ""]
    need([pin(p) for p in (receipt_path, process_path, summary_path, Path(__file__))] == inputs,
         "Report source inputs changed")
    with output.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines))
    record = {"status": "PASS_REPORT_OF_COMPLETE_EXACT_SPACE_ARITHMETIC_ONLY", "utc": datetime.now(timezone.utc).isoformat(),
              "command": sys.argv, "runtime": {"python": sys.version, "executable": sys.executable, "platform": platform.platform()},
              "inputs": inputs, "output": pin(output), "complete_source_rows": len(rows),
              "complete_layer_vectors": aggregate, "models_run": False, "originals_modified": False,
              "six_group_reference_accepted": False, "whole_study_complete": False}
    with companion.open("x", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2); handle.write("\n")
    print(json.dumps({"status": record["status"], "source_rows": len(rows), "output": str(output)}))


if __name__ == "__main__":
    main()
