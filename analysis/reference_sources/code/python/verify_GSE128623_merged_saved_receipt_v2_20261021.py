"""Explicit-manifest small-evidence closure. No RDS, fitting, P/BH or RNG replay.

CLI: MANIFEST_JSON MANIFEST_SHA256 FRESH_RECEIPT_JSON. All manifest artifacts
are physical small-file pins; embedded producer paths remain logical identities.
The manifest is frozen by its caller only after the two pending sampled QAs
and the real JOIN_ONLY/native-exit evidence exist. No production run on import.
"""
import argparse
import csv
from decimal import Decimal, InvalidOperation
import hashlib
import json
import re
from pathlib import Path
import types

CELLS = ("raw_all30", "raw_predefined29", "adapter_all30", "adapter_predefined29")
ROOT = "/root/autodl-tmp/gse128623_label_reference_20261015r01"
JOIN_SHA = "ffff5709cb45d5b4977a5ba42e993bfef1e460da94d951560cf2e9ad65d9e2bf"
SAMPLED_SHA = "3443186f36b98ebde524d473b5b41301018a977ba7b8f4336973e86d5fa270a3"
EXACT_VALIDATOR_SHA = "cbbde071a2f00eff4eba53e3bab14a22d80088d2778d0e2b1c35dd9f88595a47"
OBS_QA_SHA = "0b9f032fa9087cab6d27df5aa9c2d45662836b4cc2d3e287364bc87d5aa40469"
SENT_QA_SHA = "5a49c7bb2ad6b0c88e359bac42ff2d392cb2b4e0d07204983ad44c2ba1ddfd62"
SENTINEL_SHAS = {"e061c95737d865df65ddba9d5cb71d9345e8353592018d1ac6ab4454bd08927d",
    "404e732ff6e91d3481700a22e9a435038830eb8c780b37a864cc5ecd3f8d4e7d",
    "aecd73b6bea44dc92731006f01ba7e39172c345d80b0eeb61c7c157a2949b766",
    "b306b639b2c7cc03d4f2603d1d17a76d74c7fbdf0c5df2e8f908b356c28ac4fc"}


def need(value, label):
    if not value: raise ValueError(label)


def count(value, expected):
    return type(value) is int and value == expected


def same_bytes(a, b):
    return a["bytes"] == b["bytes"] and a["sha256"] == b["sha256"]


def pins(value):
    if isinstance(value, dict):
        if {"path", "bytes", "sha256"} <= set(value): yield value
        else:
            for child in value.values(): yield from pins(child)
    elif isinstance(value, list):
        for child in value: yield from pins(child)


def index_pins(value):
    result = {}
    for record in pins(value):
        path = record["path"]
        need(path not in result or same_bytes(result[path], record), "ConflictingRecordedIdentity")
        result[path] = record
    return result


def add_current_observed_receipt(current, cell, item, value):
    """Resolve a current pinned small receipt, never infer identity from old QA.

    Physical mirrors may be named OBSERVED_RECEIPT.json; the logical producer
    path stays OBSERVED_COMPLETE.json. Only the already-authenticated metadata
    pin grants this receipt authority. No scientific bundle is opened/hashed.
    """
    logical=ROOT+"/observed/"+cell+"/OBSERVED_COMPLETE.json"
    need(cell in CELLS and item["logical_path"]==logical and logical in current and
         same_bytes(item["file"],current[logical]), "CurrentObservedSmallReceiptPin")
    need(value["cell_id"]==cell and value["binding"]["cell_id"]==cell and
         value["binding"]["input_manifest_sha256"]==current[ROOT+"/INPUT_MANIFEST.tsv"]["sha256"] and
         value["status"]=="OBSERVED_NATIVE_FIELD_BRIDGE_PASS_RESOURCE_NOT_ADMITTED" and
         count(value["family_n"],67) and value["actual_process_exit_confirmed_here"] is False and
         value["full_run_admitted"] is False, "CurrentObservedReceiptCellAndScope")
    bundle=value["observed_bundle"]
    need(bundle["path"]==ROOT+"/observed/"+cell+"/OBSERVED_BUNDLE.rds" and
         type(bundle["bytes"]) is int and bundle["bytes"]>0 and
         isinstance(bundle["sha256"],str) and re.fullmatch("[0-9a-f]{64}",bundle["sha256"]),
         "CurrentObservedRecordedBundlePin")
    need(bundle["path"] not in current or same_bytes(bundle,current[bundle["path"]]),
         "ConflictingCurrentObservedBundleIdentity")
    current[bundle["path"]]=bundle


def scalar_equal(a, b):
    if a == b: return True
    try: return Decimal(a) == Decimal(b)
    except InvalidOperation: return False


def row_key(row):
    return tuple(row[k] for k in ("cell_id", "reference_space", "family_id"))


def validate_join(join, merged_rows, sample_rows, exact_rows):
    need(join["schema"] == "gse128623-authenticated-saved-summary-join-v1" and join["scope"] == "JOIN_ONLY" and
         join["status"] == "JOINED_SAVED_SUMMARIES_NOT_SCIENTIFIC_ACCEPTANCE" and
         join["source"]["sha256"] == JOIN_SHA and count(join["component_n"],24) and
         count(join["joined_family_n"],448), "Actual448JoinOnlyScope")
    for key in ("original_CLI_success_asserted", "scientific_acceptance", "upstream_scientific_closure_performed_here"):
        need(join[key] is False, "JoinScopeUpgraded")
    need(join["original_CLI_success"] is None, "JoinCannotReplaceOriginalCLIExit")
    for key in ("model_calls", "RNG_calls", "P_value_or_BH_recalculations", "archive_or_scientific_RDS_reads"):
        need(count(join[key],0), "JoinPerformedExtraScience")
    expected = exact_rows + [row for cell in CELLS for row in sample_rows[cell]]
    fields = ("method", "collection", "contrast", "family_id", "n", "actual_BH_adjustment_n", "cell_id",
              "reference_space", "required_allocation_n", "empirical_unavailable_n", "empirical_q_le_05_n",
              "native_fraction_allocations_any_q_le_05")
    wanted = {row_key(r):r for r in expected}; actual = {row_key(r):r for r in merged_rows}
    need(len(expected) == len(wanted) == len(merged_rows) == len(actual) == 448 and actual.keys() == wanted.keys(),
         "Exact448UniqueFamilyUnion")
    for key in wanted:
        need(all(scalar_equal(actual[key][f], wanted[key][f]) for f in fields), "JoinedSavedValueChanged")


def validate_sampled(value, native, cell, rows):
    need(value["schema"] == "gse128623-saved-sampled-output-qa-v1" and
         value["status"] == "PASS_SELECTED_SAVED_SAMPLE_JOB_OUTPUTS_AND_PINNED_METADATA_ONLY" and
         value["cell_id"] == cell and value["job_key"] == cell+"/six_group" and
         count(value["checked_allocation_n"],10001) and count(value["checked_family_n"],67) and
         len(rows) == 67 and len({row_key(r) for r in rows}) == 67 and
         all(r["cell_id"] == cell and r["reference_space"] == "six_group" and
             r["required_allocation_n"] == "10001" for r in rows), "CompleteSampledSavedQA")
    need(count(native["actual_native_exit_code"],0) and native["source_sha256"] == SAMPLED_SHA and
         native["stop_category"] is None, "ActualSampledQANativeZero")
    need(SAMPLED_SHA in {r["sha256"] for r in pins(value["metadata_records"])}, "ExecutedSampledQASourceBinding")
    for key in ("producer_independent_saved_tails_agree", "native_source_BH_diagnostic_ID_coverage", "deterministic_saved_MC_checked"):
        need(value[key] is True, "SampledSavedComponentNotAccepted")
    for key in ("worker_fit_diagnostics_checked_here", "native_process_exit_verified_here", "original_finalizer_CLI_success_asserted",
                "full_native_acceptance_granted", "final_scientific_acceptance", "input_or_saved_results_modified"):
        need(value[key] is False, "SampledScopeUpgraded")
    for key in ("archive_payloads_opened", "model_calls", "RNG_replay_calls"):
        need(count(value[key],0), "SampledQAPerformedExtraScience")


def validate_compact(value, qa, native, rows):
    need(value["schema"] == "gse128623-compact-fit-diagnostics-v1" and
         value["status"] == "SCOPED_COMPACT_DIAGNOSTIC_READ_COMPLETE" and
         count(value["job_n"],24) and count(value["rows_checked"],44792) and count(value["reused_rows"],16) and
         count(value["failed_rows"],0) and len(value["job_records"]) == 24 and len(value["output_records"]) == 24,
         "FullCommittedCompactScope")
    need(qa["status"] == "PASS_FULL_COMMITTED_COMPACT_DIAGNOSTIC_LEDGER_AND_NATIVE_EXIT" and
         count(qa["rows"],44792) and count(qa["jobs"],24) and count(qa["review_flags"],0) and
         count(qa["failed_rows"],0) and count(qa["observed_summary_bindings"],24) and
         count(native["exit_code"],0) and native["native_exit_observed"] is True and
         native["source_sha256"] == "00ac2c361d2a89cad3ca0b0e17d397671d0430b5fb72affeb8acaec52d136aa1",
         "CompactIndependentLedgerAndRealNativeZero")
    need(count(value["original_scientific_records_modified"],0), "CompactOriginalModificationCount")
    for key in ("arithmetic_or_BH_recomputed", "full_native_model_or_stochastic_oracle", "final_scientific_acceptance"):
        need(value[key] is False, "CompactScopeUpgraded")
    need(len(rows)==24 and len({(r["cell_id"],r["reference_space"]) for r in rows})==24 and
         sum(int(r["requested_records"]) for r in rows)==44792 and
         all(r["requested_records"]==r["schema_checked_records"] and
             all(r[k]=="0" for k in ("schema_or_identity_failed_records","diagnostic_review_records","fit_type_difference_records")) for r in rows),
         "Compact24SmallSummaryRows")


def validate_adapter_producer(started, native, complete, cell, context_sha):
    actor="18b4289b31d0ed85381359ecb5ca7fd20814446d42ccd4b72ebee074f0fe16df"
    output=ROOT+"/parallel_finalization_20261021_r01/"+cell
    binding=started["binding"]
    need(cell in CELLS[2:] and started["cell"]==native["cell"]==complete["cell_id"]==cell and
         started["source_sha256"]==native["source_sha256"]==actor and
         type(binding["pid"]) is int and type(binding["birth"]) is int and binding["pid"]>0 and binding["birth"]>0 and
         binding["pid"]==native["pid"] and binding["birth"]==native["birth"] and
         started["command"][-3:]==binding["command"][-3:]==[ROOT,cell,output] and
         "--file="+ROOT+"/parallel_finalization_source_20261021_r01/finalize_GSE128623_adapter_cell_parallel_20261021.R" in binding["command"] and
         count(native["exit_code"],0) and native["native_exit_observed"] is True and
         native["actual_R_start_binding_verified"] is True and native["stop_category"] is None,
         "ActualAdapterArithmeticProducerExitAndCommand")
    need(complete["schema"]=="gse128623-original-arithmetic-cell-complete-v1" and
         complete["status"]=="ONE_ORIGINAL_SAMPLED_JOB_WRITTEN_REQUIRES_NATIVE_EXIT_AND_JOIN_QA" and
         complete["context_sha256"]==context_sha and count(complete["job_n"],1) and
         count(complete["family_n"],67) and count(complete["required_archive_n"],10001) and
         count(complete["scientific_model_calls"],0) and actor in {p["sha256"] for p in pins(complete["source_records"])} and
         all(complete[k] is False for k in ("original_results_modified","RNG_final_consumption_verified","final_scientific_acceptance")),
         "AdapterArithmeticCompleteRetainsScope")


def verify(manifest, read, table, bound):
    """read/table authenticate explicit physical pins, never follow embedded RDS paths."""
    need(manifest["schema"] == "gse128623-merged-saved-evidence-manifest-v1", "FrozenMergedManifestSchema")
    envelopes = [read(p) for p in manifest["current_contexts"]]
    contexts = [value["context"] for value in envelopes]
    need(len(contexts)==2 and {c["cell_id"] for c in contexts}==set(CELLS[2:]) and
         all(c["schema"]=="gse128623-original-arithmetic-cell-context-v1" and c["root"]==ROOT and
             c["original_sources_not_modified"] is True and count(c["scientific_model_calls"],0) for c in contexts),
         "CurrentAcceptedContextDomain")
    samples={}; sample_rows={}
    need(set(manifest["sampled"])==set(CELLS), "ExactlyFourSampledComponents")
    for cell,item in manifest["sampled"].items():
        q=read(item["qa"]); rows=table(item["summary"])
        native=read(item["native"]); started=read(item["started"])
        validate_sampled(q,native,cell,rows)
        need(type(started["pid"]) is int and started["pid"]==native["pid"] and
             started["source_sha256"]==SAMPLED_SHA and
             started["command"][-3:]==[ROOT,cell,ROOT+"/saved_sampled_output_QA_20261021_r01/"+cell],
             "SampledExitBelongsToExactCellCommand")
        need(same_bytes(q["audit_summary"],item["summary"]), "SampledAuditSummaryPin")
        samples[cell]=q; sample_rows[cell]=rows
    current=index_pins([contexts,[q["metadata_records"] for q in samples.values()]])
    need(set(manifest["current_observed_receipts"])==set(CELLS), "FourCurrentObservedSmallReceipts")
    for cell,item in manifest["current_observed_receipts"].items():
        add_current_observed_receipt(current,cell,item,read(item["file"]))
    dispatch=read(manifest["dispatch_918"])
    need(manifest["dispatch_918"]["sha256"]=="4ec5049eb50086ff1cfe14ffda2e5a2f8cf6d93e16afa5ed286786e0d3aedf11" and
         same_bytes(manifest["dispatch_918"],current[ROOT+"/dispatch/main_02/DISPATCH_EXIT.json"]) and
         dispatch["status"]=="ALL_REQUESTED_SHARDS_EXITED_ZERO_NEEDS_INDEPENDENT_CHECK" and
         all(count(dispatch[k],918) for k in ("planned_shards","started_shards","finished_shards")) and
         dispatch["whole_run_independently_accepted"] is False, "Original918RecordedSuccessfulTerminalNotScienceAcceptance")
    need(set(manifest["adapter_producers"])==set(CELLS[2:]), "BothActualAdapterProducers")
    need(manifest["adapter_launcher"]["sha256"]=="157a780c7f4815871a6ad428d05e740684f6292151e19c4eb086e240c864a1d0", "FixedAdapterObserverSource")
    bound(manifest["adapter_launcher"])
    for cell,item in manifest["adapter_producers"].items():
        complete=read(item["complete"]); native=read(item["native"])
        context_sha=next(v["context_sha256"] for v in envelopes if v["context"]["cell_id"]==cell)
        validate_adapter_producer(read(item["started"]),native,complete,cell,context_sha)
        need(native["complete_sha256"]==item["complete"]["sha256"], "AdapterNativeBindsActualComplete")
        accepted=index_pins(samples[cell]["saved_input_records"])
        need(all(p["path"] in accepted and same_bytes(p,accepted[p["path"]]) for p in pins(complete["job_records"])),
             "AdapterCompleteJobBindsSavedQAInputs")
    archive_pins={}
    for item in manifest["sentinel_manifests"]:
        need(item["logical_path"] in current and same_bytes(item["file"],current[item["logical_path"]]), "CurrentPinnedWorkerManifest")
        for row in table(item["file"]):
            archive_pins[row["file"]]={"path":row["file"],"bytes":int(row["bytes"]),"sha256":row["sha256"]}
    exact=read(manifest["exact"]["receipt"]); exact_native=read(manifest["exact"]["native"])
    exact_rows=table(manifest["exact"]["summary"])
    raw=bound(manifest["exact_validator"]).read_bytes()
    need(hashlib.sha256(raw).hexdigest()==EXACT_VALIDATOR_SHA, "ReviewedExactValidatorBeforeImport")
    module=types.ModuleType("saved_exact_validator_only")
    exec(compile(raw,str(manifest["exact_validator"]["path"]),"exec"),module.__dict__)
    need(count(exact_native["actual_R_exit_code"],0), "StrictExactNativeZero")
    module.validate(exact,exact_native,exact_rows)
    need(same_bytes(exact_native["native_receipt"],manifest["exact"]["receipt"]), "ExactExitReceiptPin")
    named=[p for p in pins(exact["output_records"]) if p["path"].endswith("/FAMILY_REFERENCE_SUMMARY.tsv")]
    need(len(named)==1 and same_bytes(named[0],manifest["exact"]["summary"]), "AcceptedExactShared180SummaryPin")
    compact=read(manifest["compact"]["receipt"]); compact_qa=read(manifest["compact"]["qa"]); compact_native=read(manifest["compact"]["native"])
    need(manifest["compact"]["receipt"]["sha256"]=="805d460aeaa336c6c5f7a0727da5d2abbe1d08dfb892fcd5c036585ceaaa0888" and
         manifest["compact"]["qa"]["sha256"]=="3fd082bf091ca9524c8d0a951e37215bd167674a936c0151ce6b02f07034bf27" and
         manifest["compact"]["native"]["sha256"]=="dd47d609f929c0be61c2eb1f68f3ecaabaead53a2abd7ad8aa8bc20fe38d32fe",
         "PreviouslyAcceptedFullCompactEvidencePins")
    validate_compact(compact,compact_qa,compact_native,table(manifest["compact"]["summary"]))
    need(compact_qa["complete_sha256"]==manifest["compact"]["receipt"]["sha256"]==compact_native["complete_sha256"] and
         compact_qa["native_exit_sha256"]==manifest["compact"]["native"]["sha256"] and
         same_bytes(compact["summary"],manifest["compact"]["summary"]), "CompactBoundReceiptNativeAndSummary")
    observed=read(manifest["observed_qa"]); sentinel=read(manifest["sentinel_qa"])
    need(manifest["observed_qa"]["sha256"]==OBS_QA_SHA and manifest["sentinel_qa"]["sha256"]==SENT_QA_SHA and
         observed["overall_status"]=="PASS_SCOPED_SMALL_EVIDENCE_CONSISTENCY_ONLY" and
         observed["reconciled_totals"]["all_observed_family_slots"]==268 and
         sentinel["cumulative_distinct_fixed_sentinels"]==20 and sentinel["remaining_fixed_sentinels"]==[] and
         sentinel["full_reference_acceptance"] is False, "FixedOldSmallQAIdentitiesAndScope")
    claims={c["cell"]:c for c in observed["claims"]}
    need(set(manifest["observed"])==set(CELLS), "FourObservedReadbacks")
    for cell,record in manifest["observed"].items():
        value=read(record)
        need(record["sha256"]==claims[cell]["receipt_sha256"] and value["protected_before"]==value["protected_after"] and
             value["full_reference_acceptance"] is False and count(value["native_model_calls"],0), "ObservedHistoricalQANotPromoted")
        for key in ("bundle","observed_receipt","input_manifest","admission"):
            p=value["protected_before"][key]
            need(p["path"] in current and same_bytes(p,current[p["path"]]), "ObservedMatchesCurrentMetadataPins")
    need({p["sha256"] for p in manifest["sentinel_aggregates"]}==SENTINEL_SHAS, "FourPreviouslyAcceptedSentinelAggregates")
    selected=set()
    for record in manifest["sentinel_aggregates"]:
        value=read(record)
        need(value["status"]=="PASS_EXPLICIT_SAVED_SENTINEL_SUBSET_NOT_FULL_REFERENCE_ACCEPTANCE" and
             value["protected_before"]==value["protected_after"] and value["full_reference_acceptance"] is False and
             value["native_stochastic_algorithms_or_posteriors_reestimated"] is False and count(value["model_calls"],0), "SentinelHistoricalScope")
        for item in value["results"].values():
            key=(item["cell_id"],item["allocation_id"]); p=item["archive"]
            need(key not in selected and item["status"]=="PASS_PINNED_SAVED_NATIVE_SENTINEL" and
                 p["path"] in archive_pins and same_bytes(p,archive_pins[p["path"]]), "SentinelCurrentManifestPin")
            selected.add(key)
    need(selected=={(c,i) for c in CELLS for i in (2,101,1001,5001,10001)}, "Exactly20FixedSentinels")
    join=read(manifest["join"]["receipt"]); join_native=read(manifest["join"]["native"])
    join_start=read(manifest["join"]["started"])
    need(count(join_native["actual_native_exit_code"],0) and join_native["source_sha256"]==JOIN_SHA and
         type(join_start["pid"]) is int and join_start["pid"]==join_native["pid"] and
         join_start["source_sha256"]==JOIN_SHA and
         join_start["command"][-3:]==[join["manifest"]["path"],join["manifest"]["sha256"],join["joined_output"]["path"].rsplit("/",1)[0]] and
         same_bytes(join["joined_output"],manifest["join"]["summary"]), "ActualJoinZeroAndOutputPin")
    validate_join(join,table(manifest["join"]["summary"]),sample_rows,exact_rows)
    need(len(join["components"])==24 and len({p["key"] for p in join["components"]})==24,
         "JoinExact24ComponentIdentity")
    for item in join["components"]:
        cell=item["cell_id"]; sampled=item["reference_space"]=="six_group"
        expected=manifest["sampled"][cell]["summary"] if sampled else manifest["exact"]["summary"]
        need(same_bytes(item["summary"],expected), "JoinComponentsBoundToAcceptedQAOutputs")
    old=read(manifest["original_serial_exit"])
    need(count(old["actual_R_exit_code"],-15) and old["aggregate"] is None and
         old["source"]["sha256"]=="2d3298699771c82a983edf4f599c0cfaeb0e210f52f824f95934d9e0220c2240" and
         old["contract"]["sha256"]=="35297c8beb13e50e571888f488f7e6a1419606e682cdbf3af00f86c5f325bbb5" and
         manifest["original_serial_exit"]["sha256"]=="061045f02458972570207060caff33e95412740b1887094010240e018ac0039a",
         "PreservedRealOriginalSerialMinus15")
    return {"status":"PASS_MERGED_SAVED_REFERENCE_EVIDENCE_ONLY", "family_slots":448,"sampled_family_slots":268,
            "exact_family_slots":180,"compact_records":44792,"observed_objects":4,"fixed_sentinels":20,
            "original_serial_actual_R_exit_code":-15,"original_CLI_success_asserted":False,
            "scientific_acceptance":False,"full_native_or_stochastic_algorithm_acceptance":False,
            "RNG_replay_or_final_consumption_verified":False,"true_FDR_or_exchangeability_verified":False,
            "models_run":0,"scientific_RDS_reads":0,"P_or_BH_recalculations":0,
            "scope":"Complete saved arithmetic/diagnostic evidence union, retaining every upstream limitation"}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest",type=Path);parser.add_argument("sha256");parser.add_argument("output",type=Path)
    args=parser.parse_args(); tracked={}
    def bound(record):
        path=Path(record["path"])
        need(path.suffix.lower() in (".json",".tsv",".py",".r") and not path.is_symlink() and
             path.stat().st_size<=8*1024**2, "SmallMetadataOrTableOnly")
        raw=path.read_bytes(); need(len(raw)==record["bytes"] and hashlib.sha256(raw).hexdigest()==record["sha256"],"ExactPhysicalEvidencePin")
        tracked[str(path)]=record
        return path
    def read(record):return json.loads(bound(record).read_text(encoding="utf-8-sig"))
    def table(record):
        with bound(record).open(encoding="utf-8-sig",newline="") as stream:return list(csv.DictReader(stream,delimiter="\t"))
    own=Path(__file__).read_bytes()
    raw=args.manifest.read_bytes();need(hashlib.sha256(raw).hexdigest()==args.sha256,"FrozenCallerManifestSHA")
    result=verify(json.loads(raw),read,table,bound)
    for record in list(tracked.values()):bound(record)
    need(args.manifest.read_bytes()==raw and Path(__file__).read_bytes()==own,"ManifestOrEntryChanged")
    result.update(manifest_sha256=args.sha256,authenticated_small_inputs=list(tracked.values()),
                  source=dict(path=str(Path(__file__).resolve()),bytes=len(own),sha256=hashlib.sha256(own).hexdigest()))
    with args.output.open("x",encoding="utf-8") as stream:json.dump(result,stream,indent=2,allow_nan=False)


if __name__=="__main__":main()
