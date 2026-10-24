"""One-owner, same-three-worker successor; never an edit of the old driver.

Execution needs --handoff JSON --handoff-sha256 SHA --attempt FRESH. The external
handoff authenticates the old controller's actual terminal and every completed
round. This driver never kills, resumes, clears cache, retries or replays those
rounds. MC uses the incremental R executor; exact run keeps the old CLI; finish
uses authenticated saved prefixes. Resource/native checks retain the old plan.
No installation, new scientific configuration, model/kernel/RNG change or
whole-domain scientific-acceptance claim is made here.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import types

ROOT = Path("/root/autodl-tmp/JGW03_full_reference_linux_20261017r03")
CELLS = ("GSE289929", "GSE227494", "GSE299063")
CONFIG_SHA = "eb1a72bbecd61878e61a073b025480a0683c6be76c135257ec2e7943f30bf01c"
OLD_R_SHA = "c0c610da98b6511689b46bb5e3a98e87b64b02fa1ed978718db67b412591dd2d"
OLD_DRIVER_SHA = "8f5edfa505e44c6fbf36c605e8b582dfe23c73bb29b547645644003980fec4a0"
ADMISSION_SHA = "f4cc2d610c1c6bea9d36902188230cbac803071172a961914dac243021da40da"
NEW_R_SHA = "452acfa10dbbaa29631ad41877feab120ba421dec6852bd42ab5180ef9b1b56c"
INCREMENT_SHA = "74bb8dc0d3cdb16bb9373fdf5a9e6153af7298d5c104276e58c3247fc41b972e"
FINISH_SHA = "69ba3dd31a7854462a22e03c0121d825ea8b228a0278d35d5862f1ac8a7d5613"
KEY_EXPORT_SHA = "b8bc36c3a0ab122d00b5eff5a186669709915e5fce206d9f60bae23b8246797d"
OLD_ENTRY = ROOT / "project/code/R/run_JGW03_full_reference_linux_v3_20261017.R"
RUNTIME = Path("/root/autodl-tmp/route_b2_public_analysis/runtime/bin/Rscript")
CONTROLLER_PID, CONTROLLER_BIRTH = 659622, 2576253226
MEMORY_PER_WORKER = 3 * 1024**3
AFFINITY_EXEC = "import os,sys;os.sched_setaffinity(0,{int(sys.argv[1])});os.execv(sys.argv[2],sys.argv[2:])"


def need(ok, label):
    if not ok:
        raise RuntimeError(label)


def integer(value, minimum=0):
    return type(value) is int and value >= minimum


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def pin(path):
    path = Path(path)
    need(path.is_absolute() and not any(p.is_symlink() for p in (path, *path.parents)), "ExactNonsymlinkPath")
    before = path.stat()
    need(path.is_file(), "RegularPinnedFile")
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024**2), b""):
            h.update(block)
    after = path.stat()
    fields = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
    need(all(getattr(before, k) == getattr(after, k) for k in fields), "PinChangedDuringRead")
    return {"path": str(path), "bytes": after.st_size, "sha256": h.hexdigest()}


def bound(record):
    need(isinstance(record, dict) and integer(record.get("bytes")) and
         isinstance(record.get("path"), str) and isinstance(record.get("sha256"), str) and
         re.fullmatch("[0-9a-f]{64}", record["sha256"]), "StrictPinSchema")
    actual = pin(record["path"])
    need(actual == record, "ExactFilePin")
    return actual


def read(path):
    path = Path(path)
    need(path.stat().st_size <= 16 * 1024**2, "SmallMetadataJSON")
    def unique(pairs):
        result = {}
        for k, v in pairs:
            need(k not in result, "DuplicateJSONKey")
            result[k] = v
        return result
    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=unique,
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError("NonfiniteJSON")))


def write_new(path, value):
    with Path(path).open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n"); stream.flush(); os.fsync(stream.fileno())
    return pin(path)


def load_module(record, name):
    record = bound(record)
    module = types.ModuleType(name)
    module.__file__ = record["path"]
    raw = Path(record["path"]).read_bytes()
    need(len(raw) == record["bytes"] and hashlib.sha256(raw).hexdigest() == record["sha256"], "ExactExecutedPythonSourceBuffer")
    exec(compile(raw, record["path"], "exec"), module.__dict__)
    return module


def process_identity(pid):
    raw = (Path("/proc") / str(pid) / "stat").read_text()
    fields = raw[raw.rfind(")")+2:].split()
    return {"pid": pid, "parent": int(fields[1]), "start_ticks": int(fields[19])}


def legacy_gone():
    try:
        value = process_identity(CONTROLLER_PID)
    except (FileNotFoundError, ProcessLookupError):
        return True
    need(value["start_ticks"] != CONTROLLER_BIRTH, "LegacyControllerIdentityStillLiveNoTakeover")
    return True  # A reused PID is never signalled or treated as the old actor.


def strict_child(result, legacy):
    return (type(result.get("actual_native_exit_code")) is int and
            result["actual_native_exit_code"] == 0 and
            type(result.get("actual_R_exit_code")) is int and result["actual_R_exit_code"] == 0 and
            legacy.accepted_child(result))


def same_job(value, expected):
    return (value == expected and all(type(value.get(k)) is type(expected[k]) for k in ("first", "last")))


def remaining_plan(legacy, completed):
    need(integer(completed, 1) and completed <= 201, "CompletedMCBoundaryOnly")
    full = legacy.plan()
    need(len(full) == 204 and all(len(jobs) == 3 for jobs in full), "FrozenPlanShape")
    need(all(j["mode"] == "run" and j["space"] == "four_group_MC" for group in full[:completed] for j in group),
         "NoPreviouslyCompletedExactOrFinishReexecution")
    return list(enumerate(full[completed:], completed+1))


def original_command(job, acceptance):
    call = [str(RUNTIME), "--vanilla", str(OLD_ENTRY), str(ROOT / "CONFIG.json"),
            job["mode"], job["accession"], job["space"]]
    if job["mode"] == "run":
        call += [str(job["first"]), str(job["last"]), acceptance["path"], acceptance["sha256"]]
    return call


def validate_legacy_round(record, number, jobs, h, legacy):
    """Small original command/start/individual/whole chain, no archive reads."""
    need(set(record) == {"number", "whole", "cells"} and record["number"] == number and
         type(record["number"]) is int and set(record["cells"]) == set(CELLS), "LegacyRoundDomain")
    folder = Path(h["legacy_attempt"]) / f"round_{number:04d}"
    whole_pin = bound(record["whole"])
    need(Path(whole_pin["path"]) == folder / "ROUND_EXIT.json", "LegacyWholePath")
    whole = read(whole_pin["path"])
    results, cpus = [], []
    for job in jobs:
        cell = job["accession"]
        pins = record["cells"][cell]
        need(set(pins) == {"command", "started", "terminal"}, "LegacyCellEvidence")
        files = {key: bound(value) for key, value in pins.items()}
        for key, filename in (("command", "COMMAND.json"), ("started", "STARTED.json"), ("terminal", "PROCESS_EXIT.json")):
            need(Path(files[key]["path"]) == folder / cell / filename, "LegacyCellPath")
        cmd, start, result = (read(files[k]["path"]) for k in ("command", "started", "terminal"))
        expected = original_command(job, h["legacy_acceptance"])
        need(all(same_job(v["job"], job) for v in (cmd, start, result)) and cmd["R_argv"] == expected and
             cmd["argv"] == [h["legacy_python"], "-I", "-B", "-c", AFFINITY_EXEC,
                             str(cmd["assigned_CPU"]), *expected] and
             integer(cmd["assigned_CPU"]) and cmd["driver"] == h["legacy_driver"] and
             cmd["Rscript"] == h["runtime"] and cmd["acceptance_sha256"] == h["legacy_acceptance"]["sha256"] and
             cmd["config_sha256"] == start["config_sha256"] == CONFIG_SHA and
             cmd["entry_sha256"] == start["entry_sha256"] == OLD_R_SHA,
             "LegacyExactCommandStartSource")
        need(integer(start["pid"], 1) and integer(start["start_ticks"], 1) and
             type(start["parent"]) is int and start["parent"] == CONTROLLER_PID and
             integer(result["actual_R_pid"], 1) and integer(result["start_ticks"], 1) and
             start["pid"] == result["actual_R_pid"] and start["start_ticks"] == result["start_ticks"] and
             strict_child(result, legacy), "LegacyActualThreeZero")
        cpus.append(cmd["assigned_CPU"]); results.append(result)
    need(len(set(cpus)) == 3 and whole["results"] == results and
         whole["all_actual_exits_zero_and_observed"] is True, "LegacyWholeEqualsThreeNativeTerminals")
    return whole_pin


def command(job, h, permit=None, previous=None):
    if job["mode"] == "finish":
        need(permit is not None and previous is not None, "ExplicitSavedFinishPermitAndPrefix")
        args = ["--execute", str(ROOT / "CONFIG.json"), job["accession"], job["space"],
                previous["path"], str(previous["bytes"]), previous["sha256"],
                h["legacy_acceptance"]["path"], h["legacy_acceptance"]["sha256"],
                h["new_R"]["path"], h["new_R"]["sha256"], permit["path"], permit["sha256"]]
        need(len(args) == 13, "ExactSavedFinishArgumentCount")
        return [str(RUNTIME), "--vanilla", h["finish_R"]["path"], *args]
    if job["space"] != "four_group_MC":
        return original_command(job, h["legacy_acceptance"])
    need(permit is not None and previous is not None, "NewMCPermitAndPreviousRequired")
    args = ["--execute", str(ROOT / "CONFIG.json"), job["accession"], job["space"],
            str(job["first"]), str(job["last"]), previous["path"], str(previous["bytes"]),
            previous["sha256"], h["legacy_acceptance"]["path"], h["legacy_acceptance"]["sha256"],
            permit["path"], permit["sha256"]]
    need(len(args) == 13, "ExactSuccessorRArgumentCount")
    return [str(RUNTIME), "--vanilla", h["new_R"]["path"], *args]


def job_permission(job, h, previous, boundary, resource, producer):
    if job["mode"] == "finish":
        return {"schema": "jgw03-authenticated-saved-prefix-finish-acceptance-v1",
            "status": "ACCEPTED_FOR_ONE_SAVED_PREFIX_FINISH", "resources_admitted": True,
            "independent_wiring_review_accepted": True, "previous_boundary_authenticated_by_caller": True,
            "maximum_parallel_single_core_workers": 3, "per_worker_memory_bytes": MEMORY_PER_WORKER,
            "actual_finish": h["finish_R"], "legacy_scientific_entry": h["legacy_R"],
            "successor_context_source": h["new_R"], "increment_helper": h["increment_helper"],
            "config": h["config"], "legacy_acceptance": h["legacy_acceptance"],
            "job": {"accession": job["accession"], "space": job["space"], "previous_prefix": previous},
            "output_leaf": str(producer), "evidence": {"review": h["review"], "resources": resource, "previous_boundary": boundary}}
    return {"schema": "jgw03-incremental-prefix-executor-acceptance-v1",
        "status": "ACCEPTED_FOR_ONE_BOUNDED_INCREMENTAL_RUN", "resources_admitted": True,
        "independent_wiring_review_accepted": True, "saved_prefix_boundary_comparison_passed": True,
        "previous_boundary_authenticated_by_caller": True, "maximum_parallel_single_core_workers": 3,
        "per_worker_memory_bytes": MEMORY_PER_WORKER, "actual_executor": h["new_R"],
        "legacy_scientific_entry": h["legacy_R"], "config": h["config"],
        "legacy_acceptance": h["legacy_acceptance"], "increment_helper": h["increment_helper"],
        "job": {"accession": job["accession"], "space": job["space"], "first": job["first"],
                "last": job["last"], "previous_prefix": previous}, "output_leaf": str(producer),
        "evidence": {"review": h["review"], "resources": resource, "previous_boundary": boundary}}


def run_child(job, folder, cpu, h, permit=None, previous=None):
    """Frozen run_child supervision; only command/source/actor metadata differs.

    folder is exclusively created by this driver before the permit; R creates
    folder/producer itself. Observation failure drains owned children, never
    kills a fit or retries it. The 3 GiB observation is NOT an RLIMIT promise.
    """
    r_call = command(job, h, permit, previous)
    actual_entry = h["finish_R"] if job["mode"] == "finish" else (h["new_R"] if permit else h["legacy_R"])
    call = [sys.executable, "-I", "-B", "-c", AFFINITY_EXEC, str(cpu), *r_call]
    env = dict(os.environ,
        R_LIBS_USER="/root/autodl-tmp/GSE128623_R_overlay_20260920/library:/root/autodl-tmp/route_b2_public_analysis/R-library",
        OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", LANG="C", LC_ALL="C",
        TMPDIR=str(ROOT / "tmp"))
    write_new(folder / "COMMAND.json", {"argv": call, "R_argv": r_call, "job": job,
        "assigned_CPU": cpu, "utc": utc(), "driver": h["actual_driver"], "Rscript": h["runtime"],
        "config_sha256": CONFIG_SHA, "entry_sha256": actual_entry["sha256"],
        "legacy_scientific_entry": h["legacy_R"], "executor_acceptance": permit,
        "legacy_acceptance": h["legacy_acceptance"], "fixed_thread_environment": {key: env[key]
        for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "LANG", "LC_ALL", "R_LIBS_USER", "TMPDIR")}})
    child = subprocess.Popen(call, cwd=ROOT, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    start = time.monotonic(); peak = 0; observer = None; collected = False
    stdout, stderr = b"", b""; identity = None
    try:
        identity = process_identity(child.pid)
        write_new(folder / "STARTED.json", {**identity, "utc": utc(), "job": job,
            "config_sha256": CONFIG_SHA, "entry_sha256": actual_entry["sha256"]})
        while True:
            try:
                text = (Path("/proc") / str(child.pid) / "status").read_text()
                peak = max([peak] + [int(m[1])*1024 for m in re.finditer(r"^Vm(?:HWM|RSS):\s+(\d+) kB$", text, re.M)])
            except (FileNotFoundError, ProcessLookupError):
                pass
            try:
                stdout, stderr = child.communicate(timeout=1); collected = True; break
            except subprocess.TimeoutExpired:
                pass
    except BaseException as exc:
        observer = type(exc).__name__
    finally:
        if not collected:
            stdout, stderr = child.communicate()
        child.wait()
    import resource
    children_max = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss*1024
    result = {"status": "ACTUAL_BOUNDED_EXIT_NOT_INDEPENDENT_NUMERICAL_ACCEPTANCE", "job": job,
        "actual_R_pid": child.pid, "execution_kind": "AFFINITY_BOOTSTRAP_EXECV_SAME_PID",
        "actual_native_exit_code": child.returncode, "actual_R_exit_code": 0 if child.returncode == 0 else None,
        "nonzero_exit_stage_unknown_until_diagnosed": child.returncode != 0,
        "start_ticks": identity["start_ticks"] if identity else None, "finished_utc": utc(),
        "elapsed_seconds": time.monotonic()-start, "sampled_child_peak_RSS_bytes": peak,
        "controller_children_max_RSS_bytes": children_max, "memory_budget_met": max(peak, children_max) <= MEMORY_PER_WORKER,
        "observer_failure_class_only": observer, "stdout_bytes": len(stdout), "stderr_bytes": len(stderr),
        "stdout_sha256": hashlib.sha256(stdout).hexdigest(), "stderr_sha256": hashlib.sha256(stderr).hexdigest(),
        "raw_streams_persisted": False, "automatic_native_retry": False}
    write_new(folder / "PROCESS_EXIT.json", result)
    return result


# Explicit versioned copies of the two frozen admission functions. AST tests
# require only source==actor_for(accession,key) and the credit-call rename to
# differ. Bind into an INDEPENDENT globals dictionary; never patch old module.
def resource_disk_credits(decision, plan):
    inodes = set(); original = current = 0
    pending_publications = 0
    entry = pin(ROOT / 'project' / ENTRY)
    cutoff = dt.datetime.fromisoformat(decision['recorded_utc']).timestamp()
    need(len(plan['counts']) == 24 and sum(plan['counts'].values()) == 44792 and len(plan['tasks']) == 918,
         'OriginalG128CompletePlan')
    for key, count in plan['counts'].items():
        folder = clean(G128 / 'allocations' / key)
        need(folder.is_dir(), 'OriginalArchiveParent')
        for p in folder.iterdir():
            m = re.fullmatch(r'allocation_([0-9]{5})\.rds', p.name)
            if not m:
                continue
            need(1 <= int(m[1]) <= count and p.lstat().st_mtime >= cutoff,
                 'OriginalArchiveCountOrDate')
            if g128_publication_pending(p):
                pending_publications += 1
                continue  # No credit: retain the entire original remaining-growth allowance.
            s = allocated_file(p, inodes)
            need(1 <= int(m[1]) <= count and s.st_mtime >= cutoff, 'OriginalArchiveCountOrDate')
            original += s.st_size
    for accession in CELLS:
        folder = ROOT / 'derived' / accession
        if not folder.exists():
            continue
        clean(folder); need(folder.is_dir(), 'CurrentArchiveParent')
        identity = None
        complete_n = 0
        prepared_mtime = (ROOT / 'intermediate' / accession / 'PREPARED.json').stat().st_mtime
        for leaf in folder.iterdir():
            if not leaf.is_dir():
                continue  # Prefix/summary overhead stays in its separate 8 GiB reserve.
            key = leaf.name
            need(re.fullmatch('[0-3]{12}', key) and all(key.count(c) == 3 for c in '0123'), 'CurrentArchiveLabelKey')
            start = read(leaf / 'STARTED.json'); complete = read(leaf / 'COMPLETE.json')
            source = start['execution_entry']
            need(start['status'] == 'ALLOCATION_OWNED_BEFORE_NATIVE_CALL' and
                 complete['status'] == 'IMMUTABLE_COMPLETE_COMPACT_CHECKPOINT' and
                 start['allocation_key'] == complete['allocation_key'] == key and
                 re.fullmatch('[0-9a-f]{64}', start['scientific_identity_sha256']) and
                 start['scientific_identity_sha256'] == complete['scientific_identity_sha256'] and
                 source == actor_for(accession, key),
                 'CurrentArchiveMetadataCoupling')
            if identity is None:
                identity = start['scientific_identity_sha256']
            need(identity == start['scientific_identity_sha256'], 'CurrentCellScientificIdentityChanged')
            p = leaf / 'COMPACT.rds'; s = allocated_file(p, inodes); rec = complete['archive']
            need(integer(rec['bytes']) and s.st_size == rec['bytes'] and re.fullmatch('[0-9a-f]{64}', rec['sha256']) and
                 s.st_mtime >= prepared_mtime and integer(complete['new_native_calls']) and
                 complete['new_native_calls'] in (0, 1) and type(start['new_fit_intended']) is bool and
                 int(start['new_fit_intended']) == complete['new_native_calls'], 'CurrentArchiveStatMetadata')
            current += s.st_size
            complete_n += 1
        need(complete_n <= 10020, 'CurrentCompleteSpaceArchiveCount')
    floor = math.ceil(decision['minimum_available_disk_bytes'])
    remaining = max(floor, math.ceil(decision['complete_additional_disk_scenario_and_floor_bytes']) - original)
    return {'original_publication_pending_files_zero_credit': pending_publications,
            'original_credited_logical_bytes': original, 'original_remaining_including_floor_bytes': remaining,
            'original_floor_bytes': floor, 'current_credited_logical_bytes': current,
            'credit_is_stat_and_saved_metadata_not_RDS_numeric_validation': True}


def resource_inspect(workers=3):
    """Read-only barrier inspection; resource/window waits are not identity failures."""
    calculate(workers)  # strict caller/type check before any filesystem observation
    configuration()
    release = release_evidence()
    rp, pp = pin(G128 / 'RESOURCE_DECISION.json'), pin(G128 / 'dispatch/main_02/PLAN.json')
    need(rp['sha256'] == RESOURCE_SHA and pp['sha256'] == PLAN_SHA, 'OriginalG128BudgetOrPlanIdentity')
    decision, plan = read(rp['path']), read(pp['path'])
    controller = g128_controller(plan)
    credits = resource_disk_credits(decision, plan)
    calculation = calculate(workers, credits['original_remaining_including_floor_bytes'],
                            credits['current_credited_logical_bytes'])
    waits = set()
    own = {os.getpid(), os.getppid()}
    for p in Path('/proc').iterdir():
        if not p.name.isdigit() or int(p.name) in own:
            continue
        item = proc(int(p.name))
        if item is None:
            continue
        exe = item['cmd'].split('\0')[0]
        if 'STAR' in Path(exe).name or 'atlas' in item['cmd'].lower():
            waits.add('WAIT_ATLAS_OR_STAR_WINDOW')
        if ('/R' in exe or 'Rscript' in exe) and str(G128) not in item['cmd']:
            waits.add('WAIT_OTHER_SCIENTIFIC_R_WINDOW')
    limit = int(Path('/sys/fs/cgroup/memory.max').read_text())
    used = int(Path('/sys/fs/cgroup/memory.current').read_text())
    quota, period = map(int, Path('/sys/fs/cgroup/cpu.max').read_text().split())
    need(decision['maximum_workers_accepted'] == 12, 'OriginalG128CPUReservation')
    cpu = cpu_budget(quota, period, len(os.sched_getaffinity(0)), workers)
    waits.update(cpu['wait_categories'])
    v = os.statvfs(ROOT); free = v.f_bavail*v.f_frsize
    waits.update(resource_waits(calculation, decision['complete_memory_scenario_bytes'], limit, used, free))
    return {'schema': 'jgw03-full-resource-inspection-v1', 'utc': utc(), 'resources_admitted': not waits,
        'wait_categories': sorted(waits), 'calculation': calculation, 'disk_credits': credits,
        'free_disk_bytes': free, 'memory_limit_bytes': limit, 'memory_current_bytes': used,
        'original_G128_complete_memory_scenario_bytes': decision['complete_memory_scenario_bytes'],
        'CPU': cpu, 'original_G128_controller_evidence': controller,
        'original_budget': rp, 'original_plan': pp, 'SAM_release_evidence': release,
        'window_observation_not_future_idle_guarantee': True, 'scientific_models_launched': False}


def actor_map(h, keys):
    through = min(10001, h["completed_rounds"] * 50)
    mapping = {}
    for cell in CELLS:
        mc = keys["cells"][cell]["four_group_MC"]
        exact = keys["cells"][cell]["high_companion_exact"]
        mapping[cell] = {key: h["new_R"] if i > through else h["legacy_R"]
                         for i, key in enumerate(mc, 1)}
        for key in exact:
            mapping[cell].setdefault(key, h["legacy_R"])
    def owner(cell, key):
        need(cell in mapping and key in mapping[cell], "UnknownFrozenAllocationKey")
        return mapping[cell][key]
    return owner


def resource_adapter(admission, h, keys):
    namespace = dict(vars(admission))
    namespace["actor_for"] = actor_map(h, keys)
    namespace["resource_disk_credits"] = types.FunctionType(resource_disk_credits.__code__, namespace,
                                                            "resource_disk_credits", resource_disk_credits.__defaults__)
    inspect = types.FunctionType(resource_inspect.__code__, namespace, "resource_inspect", resource_inspect.__defaults__)
    def observe():
        value = inspect(workers=3)
        return dict(value, resource_policy="UNCHANGED_THRESHOLDS_EXACT_FROZEN_KEY_ACTOR_MAPPING_V1",
                    resource_adapter_source=h["actual_driver"], legacy_admission=h["admission"],
                    frozen_keys=h["keys_sidecar"], legacy_completed_through=min(10001, h["completed_rounds"]*50))
    return observe


def no_fit_native(value):
    need(type(value.get("actual_R_exit_code")) is int and value["actual_R_exit_code"] == 0 and
         value.get("stop_reason") is None and type(value.get("new_fits")) is int and value["new_fits"] == 0 and
         integer(value.get("actual_pid"), 1), "ActualNoFitNativeZeroRequired")


def validate_keys(keys, h, config):
    need(keys["schema"] == "jgw03-frozen-allocation-keys-v1" and keys["config"] == h["config"] and
         type(keys["new_fits"]) is int and keys["new_fits"] == 0 and keys["registry_regenerated"] is False and
         set(keys["cells"]) == set(CELLS), "FrozenKeysScope")
    source = bound(keys["label_function_source"])
    relative = Path(source["path"]).relative_to(ROOT / "project").as_posix()
    need(relative == "code/R/prepare_JGW03_reference_first_nonidentity_20261015.R" and
         {k: source[k] for k in ("bytes", "sha256")} == config["source_pins"][relative], "OriginalPureLabelFunctionPin")
    export = bound(keys["export_entry"])
    need(export["sha256"] == KEY_EXPORT_SHA and Path(export["path"]).is_relative_to(ROOT) and
         not Path(export["path"]).is_relative_to(ROOT / "project"),
         "SeparateSavedKeysExportEntry")
    for cell in CELLS:
        item = keys["cells"][cell]
        registry = bound(item["registry"])
        need(registry == config["cells"][cell]["registry"], "ExactOriginalCellRegistryPin")
        mc, exact = item["four_group_MC"], item["high_companion_exact"]
        for values, total in ((mc, 10001), (exact, 20)):
            need(isinstance(values, list) and len(values) == total and len(set(values)) == total and
                 all(isinstance(k, str) and re.fullmatch("[0-3]{12}", k) and
                     all(k.count(c) == 3 for c in "0123") for k in values), "FullOrderedFrozenKeys")
        need(mc[0] == exact[0] and item["companion_MC_indices"] ==
             [mc.index(k)+1 if k in mc else None for k in exact], "ExactMCOverlapNotActorRelabelled")
    return keys


def validate_handoff(h, legacy, pre_retirement=False):
    expected_status = ("MATERIALS_PREPARED_NOT_HANDOFF_NOT_EXECUTION_ACCEPTANCE" if pre_retirement else
                       "REVIEWED_EXPLICIT_PLANNED_HANDOFF_NOT_OLD_CONTROLLER_SUCCESS")
    need(h["schema"] == "jgw03-incremental-successor-handoff-v1" and
         h["status"] == expected_status and type(pre_retirement) is bool and
         h["run_root"] == str(ROOT), "ExplicitReviewedHandoffSchema")
    remaining_plan(legacy, h["completed_rounds"])
    expected_sources = {"actual_driver": None, "legacy_driver": OLD_DRIVER_SHA,
        "admission": ADMISSION_SHA, "legacy_R": OLD_R_SHA, "new_R": NEW_R_SHA,
        "increment_helper": INCREMENT_SHA, "finish_R": FINISH_SHA, "config": CONFIG_SHA, "runtime": None}
    for name, sha in expected_sources.items():
        actual = bound(h[name])
        need(sha is None or actual["sha256"] == sha, "FixedSourceSHA:" + name)
        if name not in ("runtime",):
            need(Path(actual["path"]).is_relative_to(ROOT), "RunBoundSource:" + name)
    need(h["actual_driver"] == pin(Path(__file__).resolve()) and
         Path(h["legacy_driver"]["path"]) == ROOT / "execution_tools_02/run_JGW03_full_linux_20261017.py" and
         Path(h["admission"]["path"]) == ROOT / "execution_tools_02/admit_JGW03_full_linux_20261017.py" and
         Path(h["legacy_R"]["path"]) == OLD_ENTRY and Path(h["config"]["path"]) == ROOT / "CONFIG.json" and
         Path(h["runtime"]["path"]) == RUNTIME and
         not Path(h["new_R"]["path"]).is_relative_to(ROOT / "project"), "ExactExecutionToolPaths")
    acceptance = bound(h["legacy_acceptance"])
    old_accept = read(acceptance["path"])
    need(old_accept["schema"] == "jgw03-linux-full-execution-acceptance-v2" and
         old_accept["status"] == "ACCEPTED_FOR_FULL_LINUX_REFERENCE_V2" and
         old_accept["config_sha256"] == CONFIG_SHA and old_accept["entry_sha256"] == OLD_R_SHA and
         old_accept["execution_tools"] == {"driver": h["legacy_driver"], "admission": h["admission"]} and
         old_accept["full_MC_family_n"] == 89 and old_accept["full_companion_family_n"] == 34 and
         old_accept["MC_nonidentity_per_cell"] == 10000 and old_accept["companion_per_cell"] == 20 and
         set(old_accept["cells"]) == set(CELLS) and old_accept["resources_admitted"] is True and
         old_accept["independent_bridge_acceptance"] is True, "UnchangedLegacyScientificAcceptance")
    for name in ("review", "context_check", "context_native_exit", "context_intent", "keys_sidecar",
                 "keys_export_native_exit", "keys_export_started"):
        bound(h[name])
        need(Path(h[name]["path"]).is_relative_to(ROOT), "RunBoundNewEvidence:" + name)
    context = read(h["context_check"]["path"])
    no_fit_native(read(h["context_native_exit"]["path"]))
    need(context["scope"] == "SAVED_ONLY_SUCCESSOR_CONTEXT_CHECK_NOT_EXECUTOR_ACCEPTANCE" and
         context["successor_source"] == h["new_R"] and context["legacy_scientific_source"] == h["legacy_R"] and
         context["all_three_contexts_exact"] is True and type(context["new_fits"]) is int and context["new_fits"] == 0,
         "AuthenticatedSavedContextIdentity")
    intent = read(h["context_intent"]["path"])
    need(intent["scope"] == "SAVED_ONLY_NO_FIT_NO_EXECUTOR_CHANGE" and
         intent["sources"][Path(h["new_R"]["path"]).name] == NEW_R_SHA and
         intent["argv"][-4:] == [str(ROOT), h["new_R"]["path"], NEW_R_SHA, h["context_check"]["path"]],
         "ContextNativeInvocationBinding")
    context_cells = {v["accession"]: v for v in context["cells"]}
    need(len(context["cells"]) == 3 and set(context_cells) == set(CELLS), "ThreeContextIdentities")
    cfg = read(h["config"]["path"])
    keys = validate_keys(read(h["keys_sidecar"]["path"]), h, cfg)
    key_exit = read(h["keys_export_native_exit"]["path"]); no_fit_native(key_exit)
    key_start = read(h["keys_export_started"]["path"])
    need(key_start["actual_pid"] == key_exit["actual_pid"] and integer(key_start["start_ticks"], 1) and
         key_start["export_entry"] == keys["export_entry"] and
         key_start["argv"][-1] == h["keys_sidecar"]["path"] and
         keys["export_entry"]["path"] in key_start["argv"], "KeysExportActualProcessAssociation")
    controller = h["legacy_controller"]
    need(controller["pid"] == CONTROLLER_PID and type(controller["pid"]) is int and
         controller["start_ticks"] == CONTROLLER_BIRTH and type(controller["start_ticks"]) is int,
         "ExactLegacyControllerIdentity")
    start_pin = bound(controller["started"])
    started = read(start_pin["path"])
    need(Path(h["legacy_attempt"]).is_relative_to(ROOT / "provenance") and
         Path(start_pin["path"]) == Path(h["legacy_attempt"]) / "CONTROLLER_STARTED.json" and
         started["pid"] == CONTROLLER_PID and started["start_ticks"] == CONTROLLER_BIRTH and
         started["driver"] == h["legacy_driver"] and started["acceptance"] == acceptance and
         started["executable"] == h["legacy_python"] and
         started["argv"] == [h["legacy_driver"]["path"], "--execute", "--acceptance", acceptance["path"],
                             "--acceptance-sha256", acceptance["sha256"], "--attempt", Path(h["legacy_attempt"]).name],
         "OriginalControllerStartExactCommand")
    if pre_retirement:
        need("terminal" not in controller, "NoFabricatedFutureTerminalInPreparedMaterial")
    else:
        terminal = read(bound(controller["terminal"])["path"])
        need(terminal["schema"] == "jgw03-legacy-controller-handoff-terminal-v1" and
             terminal["status"] == "ACTUAL_PROCESS_EXIT_OBSERVED_FOR_PLANNED_HANDOFF" and
             terminal["pid"] == CONTROLLER_PID and terminal["start_ticks"] == CONTROLLER_BIRTH and
             (type(terminal["actual_native_exit_code"]) is int or terminal["actual_native_exit_code"] is None) and
             terminal["exit_observed_after_signal"] is True and
             type(terminal["termination_signal_sent"]) is int and terminal["termination_signal_sent"] == 15 and
             terminal["driver"] == h["legacy_driver"] and
             terminal["started"] == start_pin and type(terminal["completed_rounds"]) is int and
             terminal["completed_rounds"] == h["completed_rounds"],
             "RealLegacyTerminalNotPretendedSuccess")
    history = h["completed_legacy_rounds"]
    need(len(history) == h["completed_rounds"], "AllCompletedLegacyRoundsRequired")
    for number, record in enumerate(history, 1):
        validate_legacy_round(record, number, legacy.plan()[number-1], h, legacy)
    through = min(10001, 50*h["completed_rounds"])
    need(set(h["initial_previous"]) == set(CELLS), "ThreeCanonicalPreviousPrefixes")
    for cell in CELLS:
        item = h["initial_previous"][cell]
        previous = bound(item["prefix"])
        need(Path(previous["path"]) == ROOT / "derived" / cell / f"four_group_MC_TAIL_PREFIX_{through:05d}.rds" and
             item["through"] == through and type(item["through"]) is int and
             item["scientific_identity_sha256"] == context_cells[cell]["scientific_identity_sha256"] and
             context_cells[cell]["context_identity_exact"] is True and
             re.fullmatch("[0-9a-f]{64}", item["scientific_identity_sha256"]), "InitialPrefixBoundaryAndScienceIdentity")
    if pre_retirement:
        return {"scope": "MATERIALS_VALIDATED_NOT_HANDOFF_NOT_EXECUTION_ACCEPTANCE",
                "execution_authorized": False, "legacy_terminal_required_for_actual_run": True, "validated_keys": keys}
    legacy_gone()
    return keys


def validate_new_receipt(path, job, h, previous, permit, keys):
    receipt_pin = pin(path); value = read(path)
    need(value["schema"] == "jgw03-incremental-prefix-executor-result-v1" and
         value["status"] == "BOUNDED_INCREMENTAL_PREFIX_PRODUCER_COMPLETE_NOT_INDEPENDENT_ACCEPTANCE" and
         value["accession"] == job["accession"] and value["space"] == job["space"] and
         value["first"] == job["first"] and value["last"] == job["last"] and
         value["actual_executor"] == h["new_R"] and value["legacy_scientific_entry"] == h["legacy_R"] and
         value["scientific_identity_sha256"] == h["initial_previous"][job["accession"]]["scientific_identity_sha256"] and
         value["previous_prefix"] == previous and value["executor_acceptance"] == permit and
         value["legacy_acceptance"] == h["legacy_acceptance"] and
         value["old_CLI_success_asserted"] is False and value["actual_native_exit_observed_here"] is False and
         value["scientific_acceptance"] is False and value["checker_requires_exact_entry_segment_mapping"] is True and
         value["existing_scientific_source_config_or_acceptance_modified"] is False and
         type(value["maximum_parallel_single_core_workers"]) is int and value["maximum_parallel_single_core_workers"] == 3,
         "ActualNewExecutorReceiptNotScienceAcceptance")
    p = bound(value["new_prefix"])
    need(Path(p["path"]) == ROOT / "derived" / job["accession"] / f"four_group_MC_TAIL_PREFIX_{job['last']:05d}.rds",
         "NewCanonicalPrefixPin")
    records = value["allocation_records"]
    need(len(records) == job["last"]-job["first"]+1 and value["evidence"] == read(permit["path"])["evidence"],
         "CompleteRangeAndEvidenceReceipt")
    for i, row in zip(range(job["first"], job["last"]+1), records):
        need(type(row["allocation_id"]) is int and row["allocation_id"] == i and
             row["allocation_key"] == keys["cells"][job["accession"]]["four_group_MC"][i-1] and
             row["execution_entry"] == h["new_R"] and type(row["reused_completed_archive"]) is bool and
             type(row["new_native_calls_this_invocation"]) is int and row["new_native_calls_this_invocation"] in (0, 1) and
             (not row["reused_completed_archive"] or row["new_native_calls_this_invocation"] == 0),
             "ExactNewRangeActorAndNativeAccounting")
    return p, receipt_pin


def prefix_key(job):
    return job["accession"] if job["space"] == "four_group_MC" else job["accession"] + "__high_companion_exact"


def validate_finish_receipt(path, job, h, previous, permit):
    receipt = pin(path); value = read(path)
    need(value["schema"] == "jgw03-authenticated-saved-prefix-finish-result-v1" and
         value["status"] == "SAVED_PREFIX_ARITHMETIC_COMPLETE_NOT_FULL_ARCHIVE_ACCEPTANCE" and
         value["accession"] == job["accession"] and value["space"] == job["space"] and
         value["actual_finish"] == h["finish_R"] and value["legacy_scientific_entry"] == h["legacy_R"] and
         value["successor_context_source"] == h["new_R"] and
         value["scientific_identity_sha256"] == h["initial_previous"][job["accession"]]["scientific_identity_sha256"] and
         value["final_prefix"] == previous and value["executor_acceptance"] == permit and
         value["legacy_acceptance"] == h["legacy_acceptance"] and
         value["evidence"] == read(permit["path"])["evidence"] and
         value["original_CLI_success_asserted"] is False and value["actual_native_exit_observed_here"] is False and
         value["prefix_modified"] is False and value["scientific_acceptance"] is False and
         value["full_archive_independent_checker_still_required"] is True and
         all(type(value[k]) is int and value[k] == 0 for k in ("new_fits", "allocation_RNG_calls", "archives_reread")),
         "HonestAuthenticatedSavedFinishReceipt")
    summary = bound(value["summary"])
    need(Path(summary["path"]) == ROOT / "derived" / job["accession"] / (job["space"] + "_TAIL_SUMMARY.rds"),
         "SameCellSpaceNewSummaryPin")
    bound(previous)
    return receipt


def run_round(jobs, folder, cpus, h, previous, boundary, resource, keys, legacy, runner=run_child):
    need(len(jobs) == len(cpus) == 3 and len(set(cpus)) == 3 and
         {j["accession"] for j in jobs} == set(CELLS), "ExactlyThreeDistinctSingleCoreCells")
    folder.mkdir()
    permits = {}
    for job in jobs:
        cell = job["accession"]; leaf = folder / cell; leaf.mkdir()
        if job["mode"] == "finish" or job["space"] == "four_group_MC":
            permits[cell] = write_new(leaf / "EXECUTOR_ACCEPTANCE.json", job_permission(job, h,
                previous[prefix_key(job)], boundary, resource, leaf / "producer"))
    results = []
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures = [pool.submit(runner, job, folder / job["accession"], cpu, h,
                               permits.get(job["accession"]), previous.get(prefix_key(job)))
                   for job, cpu in zip(jobs, cpus)]
        for future in futures:
            try:
                results.append(future.result())
            except BaseException as exc:
                results.append({"actual_native_exit_code": None, "actual_R_exit_code": None,
                    "memory_budget_met": False, "observer_failure_class_only": type(exc).__name__})
    accepted = all(strict_child(v, legacy) for v in results)
    whole = write_new(folder / "ROUND_EXIT.json", {"utc": utc(), "results": results,
        "all_actual_exits_zero_and_observed": accepted, "actual_driver": h["actual_driver"]})
    need(accepted, "NativeRoundFailedNoAutomaticRetry")
    updated = dict(previous); bindings = {}
    for job in jobs:
        cell = job["accession"]
        if job["mode"] == "finish":
            bindings[cell] = validate_finish_receipt(folder / cell / "producer/R_COMPLETE.json", job, h,
                previous[prefix_key(job)], permits[cell])
        elif cell in permits:
            updated[cell], bindings[cell] = validate_new_receipt(folder / cell / "producer/R_COMPLETE.json",
                job, h, previous[cell], permits[cell], keys)
        else:
            need(job["mode"] == "run" and job["space"] == "high_companion_exact" and
                 job["first"] == 1 and job["last"] == 20, "OriginalFullExact20Only")
            updated[prefix_key(job)] = pin(ROOT / "derived" / cell / "high_companion_exact_TAIL_PREFIX_00020.rds")
    validation = write_new(folder / "ROUND_ACCEPTED.json", {"scope": "ACTUAL_NATIVE_ZERO_AND_BOUNDED_PRODUCER_BINDING_ONLY",
        "whole": whole, "producer_receipts": bindings, "next_previous": updated,
        "scientific_acceptance": False, "old_controller_success_asserted": False})
    return updated, whole, validation


def progress(folder, value):
    temporary = folder / "STATUS.next.json"
    write_new(temporary, dict(value, utc=utc(), owner_pid=os.getpid()))
    os.replace(temporary, folder / "STATUS.json")


def run(handoff_path, handoff_sha, attempt):
    need(sys.platform.startswith("linux") and re.fullmatch("[A-Za-z0-9_-]{1,80}", attempt), "LinuxFreshAttemptOnly")
    handoff_pin = pin(handoff_path)
    need(handoff_pin["sha256"] == handoff_sha and Path(handoff_pin["path"]).is_relative_to(ROOT / "provenance"),
         "ExternalFixedHandoffPin")
    h = read(handoff_path)
    need(bound(h["legacy_driver"])["sha256"] == OLD_DRIVER_SHA and
         bound(h["admission"])["sha256"] == ADMISSION_SHA, "FixedPythonToolsBeforeSourceExecution")
    legacy = load_module(h["legacy_driver"], "jgw_unchanged_legacy_driver")
    admission = load_module(h["admission"], "jgw_unchanged_legacy_admission")
    keys = validate_handoff(h, legacy)
    observe = resource_adapter(admission, h, keys)
    cpus = sorted(os.sched_getaffinity(0))[-3:]
    need(len(cpus) == 3, "ThreeAllowedSingleCores")
    folder = ROOT / "provenance" / attempt
    need(not folder.exists(), "FreshAttemptNoAutomaticReplay")
    import fcntl
    with (ROOT / "provenance/FULL_REFERENCE.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        legacy_gone(); folder.mkdir()
        write_new(folder / "CONTROLLER_STARTED.json", {**process_identity(os.getpid()), "utc": utc(),
            "actual_driver": h["actual_driver"], "handoff": handoff_pin,
            "argv": sys.argv, "executable": sys.executable, "python": sys.version,
            "assigned_CPUs": cpus, "old_controller_success_asserted": False})
        remaining = remaining_plan(legacy, h["completed_rounds"])
        write_new(folder / "PLAN.json", {"remaining": remaining, "maximum_workers": 3,
            "per_worker_memory_bytes": MEMORY_PER_WORKER, "legacy_completed_rounds": h["completed_rounds"],
            "exact_may_reuse_new_MC_key_actor": True, "full_domain_checker_not_yet_compatible": True})
        previous = {cell: h["initial_previous"][cell]["prefix"] for cell in CELLS}
        boundary = h["completed_legacy_rounds"][-1]["whole"]
        completed = []
        try:
            for number, jobs in remaining:
                while True:
                    bound(handoff_pin); legacy_gone()
                    for name in ("actual_driver", "legacy_driver", "admission", "legacy_R", "new_R", "finish_R",
                                 "increment_helper", "config", "runtime", "legacy_acceptance", "keys_sidecar", "review"):
                        bound(h[name])
                    state = observe()
                    if state["resources_admitted"] is True:
                        break
                    progress(folder, {"status": "WAIT_RESOURCE_WINDOW_NO_NATIVE_PROCESS_RESTART",
                        "next_round": number, "resource_state": state})
                    time.sleep(10)
                resource = write_new(folder / f"RESOURCES_{number:04d}.json", state)
                progress(folder, {"status": "RUNNING_BOUND_SUCCESSOR_ROUND", "round": number,
                                  "completed_remaining_rounds": completed})
                previous, boundary, validation = run_round(jobs, folder / f"round_{number:04d}", cpus,
                    h, previous, boundary, resource, keys, legacy)
                completed.append({"number": number, "whole": boundary, "accepted": validation})
            bound(handoff_pin)
            for name in ("actual_driver", "legacy_driver", "admission", "legacy_R", "new_R", "finish_R",
                         "increment_helper", "config", "runtime", "legacy_acceptance", "keys_sidecar", "review"):
                bound(h[name])
            result = {"status": "REMAINING_SCOPE_ACTUAL_NATIVE_ZERO_NOT_FULL_SCIENCE_ACCEPTANCE",
                "completed_remaining_rounds": completed, "legacy_completed_rounds": h["completed_rounds"],
                "old_controller_success_asserted": False, "scientific_acceptance": False,
                "actual_driver": h["actual_driver"], "handoff": handoff_pin, "utc": utc()}
            write_new(folder / "DISPATCH_EXIT.json", result); progress(folder, result)
            return 0
        except BaseException as exc:
            write_new(folder / "DISPATCH_EXIT.json", {"status": "STOP_RETAINED_NO_AUTOMATIC_NATIVE_RETRY",
                "exception_class_only": type(exc).__name__, "completed_remaining_rounds": completed,
                "old_controller_success_asserted": False, "scientific_acceptance": False,
                "raw_error_persisted": False, "utc": utc()})
            return 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--handoff", type=Path)
    parser.add_argument("--handoff-sha256")
    parser.add_argument("--attempt")
    args = parser.parse_args()
    if not args.execute:
        print(json.dumps({"status": "PLAN_ONLY_NO_HANDOFF_OR_NATIVE_ACTION", "maximum_workers": 3,
            "per_worker_memory_bytes": MEMORY_PER_WORKER, "MC_new_executor_only": True,
            "exact_run_original_CLI": True, "finish_authenticated_saved_prefix_only": True,
            "source": pin(Path(__file__).resolve())}))
        return 0
    need(args.handoff is not None and args.handoff_sha256 is not None and args.attempt is not None,
         "FixedExternalHandoffAndFreshAttemptRequired")
    return run(args.handoff, args.handoff_sha256, args.attempt)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(json.dumps({"status": "NO_NEW_DISPATCH", "exception_class_only": type(exc).__name__,
                          "raw_error_persisted": False}))
        raise SystemExit(2)
