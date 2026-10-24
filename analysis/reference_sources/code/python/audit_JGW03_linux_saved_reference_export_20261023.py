#!/usr/bin/env python3
"""Streaming descriptive-export reconciliation; no RDS decoding or model calls.

CLI: REPORT_COMPLETE_JSON ACTUAL_R_REPORT_NATIVE_EXIT_JSON FRESH_OUTPUT_DIR.
Inputs are original, canonical local paths. No relocation or mirror is inferred.
"""
import csv
import gzip
import hashlib
import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

R_SOURCE_SHA256 = "1a7069210ad804b02cfbf88bcd0bd71cbbba6499f3c768cbc912123615ce49a4"
STATUS = "PASS_SEPARATE_LINUX_SAVED_SUMMARY_DESCRIPTIVE_EXPORT_ONLY"
QA_STATUS = "PASS_INDEPENDENT_LINUX_SAVED_EXPORT_STREAM_QA_ONLY"
TOLERANCE = 2e-15
CELLS = ("GSE289929", "GSE227494", "GSE299063")
SPACES = ("four_group_MC", "high_companion_exact")
IDENTITY = ["accession", "space", "family_id", "method", "collection", "contrast"]
SAVED = ["tail_hits_nonidentity", "null_unavailable_n", "required_nonidentity_n",
         "total_denominator", "observed_available", "observed_self_hit", "conditional_P",
         "missingness_P_lower", "missingness_P_upper", "reference_complete_for_feature",
         "finite_population_P_upper", "finite_population_P_lower"]
FEATURE_SCHEMA = IDENTITY + ["native_key_field", "key_order", "feature_key", "observed_strength"] + SAVED
FAMILY_SCHEMA = IDENTITY + ["feature_n", "observed_missing_n", "null_missing_feature_n",
    "conditional_P_available_n", "conditional_P_unavailable_n", "null_unavailable_feature_draws",
    "conditional_P_min", "conditional_P_max", "conditional_P_strict_lt_005_n",
    "conditional_P_equal_005_n", "conditional_P_le_005_n", "pointwise_finite_population_lower_min",
    "pointwise_finite_population_lower_max", "pointwise_finite_population_upper_min",
    "pointwise_finite_population_upper_max", "reference_N", "nonidentity_N", "uncertainty_scope"]
UNCERTAINTY = "Pointwise conditional finite-population reference; not effect-size, inner-MC, biological, simultaneous or FDR coverage"


class AuditError(Exception):
    """A finite rule identifier, never an original error or input value."""


def need(ok, rule):
    if not ok:
        raise AuditError(rule)


def strict_json(path):
    def pairs(items):
        d = {}
        for k, v in items:
            need(k not in d, "DuplicateJSONKey")
            d[k] = v
        return d
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle, object_pairs_hook=pairs,
                         parse_constant=lambda _: (_ for _ in ()).throw(AuditError("NonfiniteJSON")))


def canonical(path):
    need(isinstance(path, str) and bool(path), "ExplicitPath")
    p = Path(path)
    need(p.is_absolute() and p.exists() and p.is_file(), "OriginalAbsoluteFile")
    need(str(p.resolve()).replace("\\", "/") == path.replace("\\", "/"), "CanonicalPath")
    need(all(not x.is_symlink() for x in (p, *p.parents)), "NoSymlinkPath")
    return p


def pin(path):
    p = canonical(str(path))
    before = p.stat()
    digest = hashlib.sha256()
    with p.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    after = p.stat()
    need((before.st_size, before.st_mtime_ns, before.st_ino) ==
         (after.st_size, after.st_mtime_ns, after.st_ino), "StableDuringHash")
    return {"path": str(p).replace("\\", "/"), "bytes": after.st_size, "sha256": digest.hexdigest()}


def bind(record, bound):
    need(isinstance(record, dict) and set(record) == {"path", "bytes", "sha256"}, "ExactPinSchema")
    need(type(record["bytes"]) is int and record["bytes"] >= 0 and
         isinstance(record["sha256"], str) and re.fullmatch("[0-9a-f]{64}", record["sha256"]), "PinTypes")
    actual = pin(record["path"])
    normalized = dict(record, path=record["path"].replace("\\", "/"))
    need(actual == normalized, "PinSizeHashMatch")
    need(actual["path"] not in bound or bound[actual["path"]] == actual, "ConsistentRepeatedPin")
    bound[actual["path"]] = actual
    return Path(actual["path"])


def integer(value):
    need(isinstance(value, str) and re.fullmatch(r"0|[1-9][0-9]*", value) is not None, "ExactNonnegativeInteger")
    return int(value)


def boolean(value):
    need(value in ("TRUE", "FALSE"), "ExplicitBoolean")
    return value == "TRUE"


def probability(value):
    if value == "NA":
        return None
    try:
        number = float(value)
    except (ValueError, TypeError):
        raise AuditError("ProbabilityType") from None
    need(math.isfinite(number) and 0 <= number <= 1, "FiniteProbabilityRange")
    return number


def close(actual, expected, rule):
    need((actual is None and expected is None) or
         (actual is not None and expected is not None and abs(actual - expected) <= TOLERANCE), rule)


def table(path, schema, compressed=False):
    opener = gzip.open if compressed else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        need(reader.fieldnames == schema, "ExactTSVSchema")
        for row in reader:
            need(None not in row and all(v is not None for v in row.values()), "ExactRowWidth")
            yield row


def family_identity(row):
    result = tuple(row[k] for k in IDENTITY)
    need(all(result) and row["accession"] in CELLS and row["space"] in SPACES, "KnownFamilyIdentity")
    need(row["family_id"] == "::".join(row[k] for k in ("method", "collection", "contrast")), "FamilyMethodCollectionContrast")
    if row["space"] == "high_companion_exact":
        need(row["contrast"] == "high_vs_vehicle", "ExactCompanionContrast")
    return result


def reconcile(family_path, feature_path):
    families = list(table(family_path, FAMILY_SCHEMA))
    features = table(feature_path, FEATURE_SCHEMA, compressed=True)
    try:
        return reconcile_tables(families, features)
    finally:
        features.close()


def reconcile_tables(families, features):
    need(len(families) == 123, "All123Families")
    seen, counts, previous_group = set(), {}, -1
    for row in families:
        key = family_identity(row)
        need(key not in seen, "UniqueFamily")
        seen.add(key)
        group = (row["accession"], row["space"])
        index = CELLS.index(group[0]) * 2 + SPACES.index(group[1])
        need(index >= previous_group, "DeclaredCellSpaceOrder")
        previous_group = index
        counts[group] = counts.get(group, 0) + 1
    for accession in CELLS:
        need(counts.get((accession, SPACES[0])) == (31 if accession == CELLS[0] else 29) and
             counts.get((accession, SPACES[1])) == (12 if accession == CELLS[0] else 11), "FixedCellFamilyCounts")
    total_rows, observed_states = 0, {"finite": 0, "positive_infinity": 0, "negative_infinity": 0, "NA": 0, "NaN": 0}
    for family in families:
        identity = family_identity(family)
        n = integer(family["feature_n"])
        need(n > 0, "NonemptyFamily")
        N = 10001 if family["space"] == SPACES[0] else 20
        need(integer(family["reference_N"]) == N and integer(family["nonidentity_N"]) == N - 1, "FamilyFixedDenominator")
        need(family["uncertainty_scope"] == UNCERTAINTY, "ExplicitUncertaintyScope")
        keys, key_field = set(), None
        totals = dict.fromkeys(["observed_missing_n", "null_missing_feature_n", "conditional_P_available_n",
            "conditional_P_unavailable_n", "null_unavailable_feature_draws", "conditional_P_strict_lt_005_n",
            "conditional_P_equal_005_n", "conditional_P_le_005_n"], 0)
        ranges = {k: [None, None] for k in ("conditional_P", "finite_population_P_lower", "finite_population_P_upper")}
        for order in range(1, n + 1):
            row = next(features, None)
            need(row is not None and family_identity(row) == identity, "FeatureFamilyDeclaredOrder")
            need(integer(row["key_order"]) == order, "SequentialFeatureKeyOrder")
            need(row["feature_key"] and row["feature_key"] not in keys, "UniqueFeatureKeyWithinFamily")
            keys.add(row["feature_key"])
            need(row["native_key_field"] in ("raw_gene_id", "feature", "term_id"), "KnownNativeKeyField")
            key_field = row["native_key_field"] if key_field is None else key_field
            need(row["native_key_field"] == key_field, "StableFamilyNativeKeyField")
            available = boolean(row["observed_available"])
            strength = row["observed_strength"]
            if strength in ("NA", "NaN"):
                observed_states[strength] += 1
                need(not available, "UnavailableObservedStrengthState")
            else:
                try:
                    obs = float(strength)
                except ValueError:
                    raise AuditError("ObservedStrengthType") from None
                need(not math.isnan(obs), "ExplicitNaNToken")
                observed_states["finite" if math.isfinite(obs) else "positive_infinity" if obs > 0 else "negative_infinity"] += 1
            hits, missing = integer(row["tail_hits_nonidentity"]), integer(row["null_unavailable_n"])
            need(integer(row["required_nonidentity_n"]) == N - 1 and integer(row["total_denominator"]) == N, "FeatureFixedDenominator")
            need(hits + missing <= N - 1 and (available or hits == 0), "BoundedHitsAndMissingCounters")
            need(row["observed_self_hit"] == ("1" if available else "NA"), "ObservedAddedOnceOrUnavailable")
            complete = available and missing == 0
            need(boolean(row["reference_complete_for_feature"]) == complete, "ReferenceAvailabilityState")
            probs = {k: probability(row[k]) for k in ("conditional_P", "missingness_P_lower", "missingness_P_upper",
                "finite_population_P_lower", "finite_population_P_upper")}
            for name in ("conditional_P", "finite_population_P_lower", "finite_population_P_upper"):
                need((probs[name] is not None) == complete, "PointProbabilityAvailability")
            for name in ("missingness_P_lower", "missingness_P_upper"):
                need((probs[name] is not None) == available, "MissingnessBoundAvailability")
            close(probs["conditional_P"], (1 + hits) / N if complete else None, "InclusiveFixedDenominatorP")
            close(probs["missingness_P_lower"], (1 + hits) / N if available else None, "MissingnessLowerArithmetic")
            close(probs["missingness_P_upper"], (1 + hits + missing) / N if available else None, "MissingnessUpperArithmetic")
            if complete:
                need(probs["finite_population_P_lower"] <= probs["finite_population_P_upper"], "OrderedFinitePopulationEndpoints")
                if N == 20:
                    close(probs["finite_population_P_lower"], probs["conditional_P"], "ExactCensusLowerEndpoint")
                    close(probs["finite_population_P_upper"], probs["conditional_P"], "ExactCensusUpperEndpoint")
            totals["observed_missing_n"] += not available
            totals["null_missing_feature_n"] += missing > 0
            totals["null_unavailable_feature_draws"] += missing
            totals["conditional_P_available_n"] += complete
            totals["conditional_P_unavailable_n"] += not complete
            point = probs["conditional_P"]
            if point is not None:
                totals["conditional_P_strict_lt_005_n"] += point < .05
                totals["conditional_P_equal_005_n"] += point == .05
                totals["conditional_P_le_005_n"] += point <= .05
            for name in ranges:
                value = probs[name]
                if value is not None:
                    low, high = ranges[name]
                    ranges[name] = [value if low is None else min(low, value), value if high is None else max(high, value)]
            total_rows += 1
        for name, value in totals.items():
            need(integer(family[name]) == value, "FamilyCountReconciliation")
        for name, limits in ranges.items():
            prefix = "conditional_P" if name == "conditional_P" else "pointwise_" + name.replace("_P_", "_")
            for suffix, expected in zip(("min", "max"), limits):
                close(probability(family[prefix + "_" + suffix]), expected, "FamilyRangeReconciliation")
    need(next(features, None) is None, "NoUndeclaredFeatureRows")
    return families, total_rows, observed_states


def audit(report_path, exit_path, output_dir, *, expected_source_sha=R_SOURCE_SHA256):
    bound = {}
    self_pin = pin(Path(__file__).resolve())
    bind(self_pin, bound)
    report_pin, exit_pin = pin(report_path), pin(exit_path)
    report = strict_json(bind(report_pin, bound))
    terminal = strict_json(bind(exit_pin, bound))
    need(report.get("status") == STATUS, "DescriptiveReportStatus")
    need(report.get("family_n") == 123 and report.get("MC_family_n") == 89 and
         report.get("companion_family_n") == 34, "ReceiptFamilyScope")
    for key in ("all_input_pins_unchanged", "saved_twelve_column_values_copied", "arithmetic_assertions_only"):
        need(report.get(key) is True, "ReceiptTrueScopeFlags")
    for key in ("P_q_BH_or_endpoints_recomputed", "unavailable_values_imputed", "denominator_reduced",
                "original_S1_or_candidate_gates_changed", "biological_or_whole_study_acceptance"):
        need(report.get(key) is False, "ReceiptFalseScopeFlags")
    for key in ("new_native_fit_calls", "RNG_calls", "archives_scanned"):
        need(type(report.get(key)) is int and report[key] == 0, "NoFitRNGOrArchiveCalls")
    need(report.get("canonical_saved_column_order") == SAVED and
         report.get("full_balanced_four_group_population") == 369600, "SavedSchemaAndPopulation")
    need(type(terminal.get("actual_native_exit_code")) is int and terminal["actual_native_exit_code"] == 0 and
         type(terminal.get("child_pid")) is int and terminal["child_pid"] > 0 and
         type(terminal.get("birth")) is int and terminal["birth"] > 0 and
         terminal.get("output") == str(Path(report_path).parent).replace("\\", "/") and
         terminal.get("raw_streams_persisted") is False and terminal.get("automatic_retry") is False and
         terminal.get("source_pins_unchanged") is True,
         "ActualReportNativeExitZero")
    need(report.get("source", {}).get("sha256") == expected_source_sha, "ExpectedReportSourceHash")
    bind(report["source"], bound)
    inputs = report.get("inputs")
    need(isinstance(inputs, dict) and bool(inputs), "ExplicitInputPinMap")
    for key, value in inputs.items():
        need(key.replace("\\", "/") == value.get("path", "").replace("\\", "/"), "InputMapPathIdentity")
        bind(value, bound)
    outputs = report.get("outputs")
    expected_names = {"FAMILY_REFERENCE_SUMMARY.tsv", "FEATURE_SAVED_REFERENCE.tsv.gz", "READING_NOTES.md"}
    need(isinstance(outputs, list) and len(outputs) == 3 and
         {Path(x["path"]).name for x in outputs} == expected_names, "ExactReportOutputs")
    paths = {}
    for record in outputs:
        p = bind(record, bound)
        need(p.parent == Path(report_path).parent, "SameReportOutputDirectory")
        paths[p.name] = p
    families, rows, states = reconcile(paths["FAMILY_REFERENCE_SUMMARY.tsv"], paths["FEATURE_SAVED_REFERENCE.tsv.gz"])
    need(type(report.get("feature_rows")) is int and report["feature_rows"] == rows, "ReceiptFeatureRowCount")
    for record in bound.values():
        need(pin(record["path"]) == record, "AllPinsUnchangedAfterStreaming")
    out = Path(output_dir)
    need(out.is_absolute() and out.name not in (".", "..") and not out.exists() and out.parent.is_dir() and
         all(not x.is_symlink() for x in (out.parent, *out.parent.parents)), "FreshCanonicalQAOutput")
    need(str(out.parent.resolve()) == str(out.parent), "CanonicalQAParent")
    out.mkdir()
    claim_path = out / "CLAIM_INVENTORY.tsv"
    with claim_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FAMILY_SCHEMA, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(families)
    receipt = {"status": QA_STATUS, "UTC": datetime.now(timezone.utc).isoformat(),
        "scope": "Independent streamed exported arithmetic/count reconciliation only; no RDS value identity or native algorithm/biological/whole-study certification",
        "claim_inventory_scope": "Verified 123-family table inventory/review matrix; not human manuscript or all-claim QA",
        "source": self_pin, "R_report_source": report["source"],
        "inputs": list(bound.values()), "family_n": 123, "MC_family_n": 89, "companion_family_n": 34,
        "feature_rows": rows, "observed_strength_states_preserved": states,
        "arithmetic_absolute_tolerance": TOLERANCE, "integer_count_tolerance": 0,
        "all_input_pins_unchanged": True, "original_rows_removed": 0,
        "new_native_fit_calls": 0, "RNG_calls": 0, "P_q_BH_or_endpoints_recomputed": False,
        "arithmetic_assertions_only": True, "denominator_reduced": False,
        "unavailable_values_imputed": False, "raw_streams_persisted": False, "automatic_retry": False,
        "original_S1_or_candidate_gates_changed": False, "biological_or_whole_study_acceptance": False,
        "outputs": [pin(claim_path)]}
    for record in bound.values():
        need(pin(record["path"]) == record, "AllPinsUnchangedAtReceipt")
    with (out / "QA_RECEIPT.json").open("w", encoding="utf-8") as handle:
        json.dump(receipt, handle, indent=2, allow_nan=False)
        handle.write("\n")
    return receipt


def main(args):
    need(len(args) == 3, "ThreeExplicitArguments")
    audit(*args)
    print(QA_STATUS)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv[1:]))
    except Exception as error:
        # No original condition, input value, or stack trace is persisted.
        print("JGW03_EXPORT_QA_STOP_" + (str(error) if isinstance(error, AuditError) else type(error).__name__), file=sys.stderr)
        sys.exit(1)
