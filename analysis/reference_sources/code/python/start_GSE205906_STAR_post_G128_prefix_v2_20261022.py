"""Post-G128 resource epoch: one original STAR admit/index call, plan-only default.

No scheduling handoff, signals, cache advice, deletion, fitting, new parameters,
or old-source edits. The current three-worker JGW controller continues. The
original 40-GiB direct-headroom gate and per-command RLIMIT_AS remain unchanged.
Only a new resource ledger retires completed G128 future-growth planning.
"""
from __future__ import annotations
import argparse
import datetime
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import socket
import subprocess
import sys
import threading
import time

GIB = 1024**3
INSTANCE = "62c8pzy536-50fb405a"
STAR = Path("/root/autodl-tmp/GSE205906_matched_STAR_Ensembl101_20261015r01")
G128 = Path("/root/autodl-tmp/gse128623_label_reference_20261015r01")
JGW = Path("/root/autodl-tmp/JGW03_full_reference_linux_20261017r03")
CLOSURE = G128 / "post_computation_saved_closure_20261022_r01"
STAR_SHA = "7fe64fcfcee554e1500f60da2cd7b4a3ce8a84fcd4e46b84c3c7e14ec1fd907d"
CHECKPOINT_SHA = "fe9f804c162a26e34c1f5fd52d5c436b8d72bd81fe19a032691a741a1323de53"
CONTROLLER_SHA = "3f1a3b9cf438ce770c0528d809e4b0c8a3cbcdb37d87e6765b64f16624774b87"
CONTROLLER, BIRTH = 2111, 2612762666
INITIAL_SHA = "175ddf0124a8e583146484930d8101d7732af8194fce19a4bd4030b667ab6d11"
CLOSURE_SHAS = dict(
    receipt="a3fc4c8bde06b5dd0ab56c91194c750db6f4b7f7992864b7a628b0acd71f461e",
    actual_exit="efe1f4a434f0485ed6fe9d6c765025d004710622a41c4aa666e9f1841547277c",
    input_manifest="570ea74f7741cbf292c5107438b9ada1103713f92cd1152fe62936fcd49fc69d",
    checker="74d947196c6294c9cec785d9becf729bce415247e5e41ec50df10dacd3136e13")
JGW_GROWTH = 261520647968  # Conservative accepted initial forecast; no later credit.
OTHER_GROWTH, DISK_FLOOR = 8*GIB, 50*GIB
RESERVED = JGW_GROWTH + OTHER_GROWTH + DISK_FLOOR
BUDGETS = {"admit": 2*GIB, "index": 64*GIB}


def need(ok, label):
    if not ok: raise ValueError(label)


def integer(value, minimum=0):
    return type(value) is int and value >= minimum


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def read(path):
    p = Path(path)
    need(p.stat().st_size <= 8*GIB//1024, "SmallMetadataOnly")
    def unique(items):
        result = {}
        for key, value in items:
            need(key not in result, "DuplicateJSONKey"); result[key] = value
        return result
    return json.loads(p.read_text(encoding="utf-8-sig"), object_pairs_hook=unique,
        parse_constant=lambda _: (_ for _ in ()).throw(ValueError("NonfiniteJSON")))


def pin(path):
    p = Path(path)
    need(p.is_absolute() and not any(q.is_symlink() for q in (p, *p.parents)) and p.is_file(),
         "RegularAbsolutePinnedFile")
    before = p.stat(); digest = hashlib.sha256()
    with p.open("rb") as stream:
        for block in iter(lambda: stream.read(1024**2), b""): digest.update(block)
    after = p.stat()
    fields = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
    need(tuple(getattr(before, x) for x in fields) == tuple(getattr(after, x) for x in fields),
         "PinnedFileChanged")
    return dict(path=str(p), bytes=after.st_size, sha256=digest.hexdigest())


def bound(record, sha=None, root=None):
    need(isinstance(record, dict) and set(record) == {"path", "bytes", "sha256"} and
         integer(record["bytes"]) and isinstance(record["sha256"], str) and
         re.fullmatch("[0-9a-f]{64}", record["sha256"]), "FullPinSchema")
    need(sha is None or record["sha256"] == sha, "FixedSourceOrEvidenceSHA")
    need(root is None or Path(record["path"]).is_relative_to(root), "FixedEvidenceNamespace")
    need(pin(record["path"]) == record, "CurrentFullPin")
    return record


def check_closure(merged, exited):
    need(merged["status"] == "PASS_MERGED_SAVED_REFERENCE_EVIDENCE_ONLY" and
         all(integer(merged[k]) and merged[k] == v for k, v in
             dict(family_slots=448, sampled_family_slots=268, exact_family_slots=180,
                  compact_records=44792, observed_objects=4, fixed_sentinels=20,
                  models_run=0, scientific_RDS_reads=0, P_or_BH_recalculations=0).items()) and
         merged["manifest_sha256"] == CLOSURE_SHAS["input_manifest"] and
         merged["source"]["sha256"] == CLOSURE_SHAS["checker"] and
         merged["original_CLI_success_asserted"] is False and
         merged["scientific_acceptance"] is False, "ActualSavedClosureScope")
    need(type(exited["actual_native_exit_code"]) is int and exited["actual_native_exit_code"] == 0 and
         exited["source_sha256"] == CLOSURE_SHAS["checker"] and
         exited["manifest_sha256"] == CLOSURE_SHAS["input_manifest"] and
         exited["platform"] == "Windows" and exited["receipt_exists"] is True and
         exited["scope"] == "ACTUAL_LOCAL_SMALL_METADATA_VERIFIER_EXIT_NOT_MODEL_EXECUTION" and
         all(type(exited[k]) is int and exited[k] == 0 for k in
             ("scientific_model_calls", "scientific_RDS_reads", "stderr_bytes")),
         "ActualLocalMetadataExitNotPretendedLinuxAcceptance")


def material(value):
    need(value["schema"] == "gse205906-post-g128-prefix-material-v1" and
         value["instance"] == INSTANCE, "FixedResourceEpoch")
    bound(value["STAR_source"], STAR_SHA, STAR)
    bound(value["checkpoint_helper"], CHECKPOINT_SHA, STAR)
    for key, sha in CLOSURE_SHAS.items(): bound(value["G128_closure"][key], sha, CLOSURE)
    check_closure(read(value["G128_closure"]["receipt"]["path"]),
                  read(value["G128_closure"]["actual_exit"]["path"]))
    initial = bound(value["JGW_initial_resources"], INITIAL_SHA, JGW)
    need(initial["path"] == str(JGW / "provenance/full_reference_02/INITIAL_RESOURCES.json") and
         initial["bytes"] == 7503, "UnchangedConservativeInitialGrowthAnchor")
    started = read(bound(value["JGW_controller_started"], root=JGW / "provenance")["path"])
    need(type(started["pid"]) is int and started["pid"] == CONTROLLER and
         type(started["start_ticks"]) is int and started["start_ticks"] == BIRTH and
         started["actual_driver"]["sha256"] == CONTROLLER_SHA, "CurrentSuccessorStart")
    need(value["JGW_controller_started"]["path"] == str(JGW /
         "provenance/incremental_reference_01/CONTROLLER_STARTED.json"), "CurrentControllerAttempt")
    bound(started["actual_driver"], CONTROLLER_SHA, JGW)
    need(started["argv"][0] == started["actual_driver"]["path"] and
         len(started["assigned_CPUs"]) == 3 and len(set(started["assigned_CPUs"])) == 3 and
         all(integer(x) for x in started["assigned_CPUs"]), "ThreeCurrentSingleCPUWorkers")
    return started


def proc(pid):
    p = Path("/proc") / str(pid)
    fields = (p / "stat").read_text().rsplit(") ", 1)[1].split()
    argv = [os.fsdecode(x) for x in (p / "cmdline").read_bytes().rstrip(b"\0").split(b"\0")]
    return dict(pid=pid, start_ticks=int(fields[19]), state=fields[0], argv=argv)


def current_controller(started, actual):
    need(actual["pid"] == CONTROLLER and actual["start_ticks"] == BIRTH and
         actual["state"] not in ("Z", "X", "T", "t") and
         actual["argv"] == [started["executable"], "-I", "-B", *started["argv"]],
         "CurrentControllerExactBirthArgvNoPause")


def no_g128_executors():
    found = []
    for directory in Path("/proc").iterdir():
        if not directory.name.isdigit(): continue
        try: actual = proc(int(directory.name))
        except (FileNotFoundError, ProcessLookupError): continue
        if actual["state"] in ("Z", "X"): continue
        if any(x.startswith(str(G128) + "/") and Path(x).suffix.lower() in (".r", ".py")
               for x in actual["argv"]): found.append(actual["pid"])
    need(not found, "NoCurrentG128ExecutionBeforeRetiringItsPlanningReserve")
    return True


def resource_plan(stage, limit, current, free, cpus, assigned, available, cpu_capacity):
    need(stage in BUDGETS and all(integer(x) for x in (limit, current, free)) and
         limit == 96*GIB and current <= limit, "Real96GiBCgroup")
    need(integer(cpu_capacity) and cpu_capacity == 20 and
         len(cpus) == len(set(cpus)) == 2 and all(integer(x) for x in cpus) and
         not set(cpus).intersection(assigned) and set(cpus) <= set(available) and
         len(available) >= 10, "TwoDisjointSTARCPUsWithinCurrentTwentyVCPUQuota")
    return dict(policy="POST_G128_COMPLETE_CONCURRENT_PREFIX_V1",
        G128_active_fit_or_reader_workers=0, G128_future_archive_growth_bytes=0,
        JGW_workers=3, JGW_total_worker_planning_bytes=9*GIB,
        STAR_native_RLIMIT_AS_bytes=40*GIB, global_memory_floor_bytes=8*GIB,
        joint_memory_planning_bytes=57*GIB, planning_is_not_aggregate_RSS_hard_limit=True,
        cgroup_limit_bytes=limit, cgroup_current_bytes=current, direct_headroom_bytes=limit-current,
        direct_memory_gate_bytes=40*GIB, memory_admitted=limit-current >= 40*GIB,
        available_cloud_bytes=free, conservative_JGW_growth_bytes=JGW_GROWTH,
        other_finite_growth_bytes=OTHER_GROWTH, one_global_disk_floor_bytes=DISK_FLOOR,
        reserved_growth_bytes=RESERVED, stage_disk_budget_bytes=BUDGETS[stage],
        required_available_cloud_bytes=RESERVED+BUDGETS[stage],
        disk_admitted=free >= RESERVED+BUDGETS[stage],
        actual_STAR_CPU_ids=cpus, actual_JGW_CPU_ids=assigned,
        cgroup_vCPU_quota=cpu_capacity, process_affinity_CPU_count=len(available),
        joint_CPU_planning_with_five_vCPU_floor=10,
        current_JGW_continues_no_idle_claim=True,
        original_quiet_window_flag_means_measured_reviewed_concurrent_window=True)


def snapshot(stage, started, cpus):
    current_controller(started, proc(CONTROLLER)); no_g128_executors()
    c = Path("/sys/fs/cgroup")
    quota, period = (c/"cpu.max").read_text().split()
    need(quota.isdigit() and period.isdigit() and int(period)>0 and
         int(quota) == 20*int(period), "ActualTwentyVCPUQuotaNotAffinityCount")
    return resource_plan(stage, int((c/"memory.max").read_text()),
        int((c/"memory.current").read_text()), shutil.disk_usage(STAR).free,
        cpus, started["assigned_CPUs"], sorted(os.sched_getaffinity(0)), 20)


def command(value, stage):
    return [sys.executable, "-I", "-B", value["STAR_source"]["path"], stage, "--execute",
            "--quiet-window-confirmed", "--reserved-growth-bytes", str(RESERVED),
            "--disk-budget-bytes", str(BUDGETS[stage])]


def put(path, value):
    with Path(path).open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, allow_nan=False); stream.write("\n")
        stream.flush(); os.fsync(stream.fileno())


def observe(call, cpus, output):
    """Direct parent waits even if /proc/STARTED fails; never signals any child."""
    child = subprocess.Popen(call, cwd=STAR, preexec_fn=lambda: os.sched_setaffinity(0, cpus),
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    streams = {}; drain_failures = {}
    def drain(name, stream):
        n = 0; h = hashlib.sha256()
        try:
            for block in iter(lambda: stream.read(65536), b""): n += len(block); h.update(block)
        except BaseException as exc:
            drain_failures[name] = type(exc).__name__
        finally:
            try: stream.close()
            except BaseException as exc: drain_failures[name] = type(exc).__name__
            streams[name] = {"bytes": n, "sha256": h.hexdigest(),
                             "complete": name not in drain_failures}
    threads = [threading.Thread(target=drain, args=(name, stream)) for name, stream in
               (("stdout", child.stdout), ("stderr", child.stderr))]
    for thread in threads: thread.start()
    identity = None; failure = None; before = time.monotonic()
    try:
        identity = proc(child.pid)
        put(output / "STARTED.json", dict(actual_pid=child.pid, identity=identity, argv=call, utc=utc()))
    except BaseException as exc: failure = type(exc).__name__
    finally:
        code = child.wait()
        for thread in threads: thread.join()
    if (drain_failures or set(streams) != {"stdout", "stderr"}) and failure is None:
        failure = "PipeDrainFailure"
    return dict(actual_native_exit_code=code, actual_pid=child.pid,
        identity=identity, argv=call, elapsed_seconds=time.monotonic()-before,
        stream_summaries=streams, raw_streams_persisted=False,
        stream_failure_classes=drain_failures, observer_failure_class_only=failure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=tuple(BUDGETS))
    parser.add_argument("--material", type=Path, required=True)
    parser.add_argument("--material-sha256", required=True)
    parser.add_argument("--cpus", type=int, nargs=2, required=True)
    parser.add_argument("--attempt", required=True)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    need(sys.platform.startswith("linux") and socket.gethostname() == "autodl-container-"+INSTANCE,
         "OnlyDeclaredNewInstance")
    need(re.fullmatch("[A-Za-z0-9][A-Za-z0-9_-]{0,100}", args.attempt), "FreshAttemptBasename")
    manifest_pin = pin(args.material)
    need(manifest_pin["sha256"] == args.material_sha256, "CallerFrozenMaterialSHA")
    value = read(args.material); started = material(value)
    spec = importlib.util.spec_from_file_location("star_checkpoint_only", value["checkpoint_helper"]["path"])
    helper = importlib.util.module_from_spec(spec); spec.loader.exec_module(helper)
    # This single source-safe helper is reused, never its old controller/scheduler.
    prior = helper.checkpoint(args.stage)
    if args.stage == "index": need(helper.checkpoint("admit") is not None, "FormalAdmissionBeforeIndex")
    plan = snapshot(args.stage, started, args.cpus)
    call = command(value, args.stage)
    report = dict(status="READ_ONLY_POST_G128_RESOURCE_PLAN_NOT_EXECUTION", stage=args.stage,
        source=pin(Path(__file__).resolve()), material=manifest_pin,
        G128_saved_closure_origin="Actual Windows metadata checker; unchanged cloud mirrors",
        scientific_acceptance=False, existing_successful_checkpoint=prior,
        resource_plan=plan, argv=call, no_old_source_or_controller_changes=True)
    if not args.execute: print(json.dumps(report)); return 0
    import fcntl
    with (STAR / "POST_G128_PREFIX.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        output = STAR / "provenance" / args.attempt
        need(not output.exists(), "FreshOperationNoOverwriteOrRetry")
        output.mkdir()
        material(value); need(pin(args.material) == manifest_pin, "MaterialChangedBeforeIntent")
        plan = snapshot(args.stage, started, args.cpus)
        report["resource_plan"] = plan
        if prior is not None:
            put(output / "INTENT.json", report)
            put(output / "RESULT.json", dict(status="REUSED_ORIGINAL_SUCCESSFUL_CHECKPOINT_NO_NATIVE_CALL",
                checkpoint=helper.checkpoint(args.stage), actual_native_exit_code=None,
                scientific_acceptance=False, new_models=0)); return 0
        need(plan["memory_admitted"] and plan["disk_admitted"], "OriginalActualResourceGatesNotMet")
        put(output / "INTENT.json", report)
        result = observe(call, args.cpus, output)
        put(output / "NATIVE_EXIT.json", result)
        need(type(result["actual_native_exit_code"]) is int and result["actual_native_exit_code"] == 0 and
             result["observer_failure_class_only"] is None, "ActualOriginalSTARStageNotSuccessful")
        completed = helper.checkpoint(args.stage)
        need(completed is not None, "OriginalNativeCompleteRequired")
        material(value); need(pin(args.material) == manifest_pin and pin(Path(__file__).resolve()) == report["source"],
                              "FixedSourcesOrMaterialChanged")
        current_controller(started, proc(CONTROLLER))
        put(output / "RESULT.json", dict(status="ORIGINAL_PREFIX_STAGE_ACTUAL_NATIVE_ZERO",
            stage=args.stage, native_exit=pin(output/"NATIVE_EXIT.json"), checkpoint=completed,
            source=report["source"], material=manifest_pin,
            scientific_acceptance=False, new_models=0, JGW_controller_not_signalled=True))
    return 0


if __name__ == "__main__":
    try: raise SystemExit(main())
    except Exception as exc:
        print(json.dumps(dict(status="NOT_STARTED_OR_STOP_RETAINED_NO_AUTOMATIC_RETRY",
                              exception_class_only=type(exc).__name__, raw_error_persisted=False)))
        raise SystemExit(1)

