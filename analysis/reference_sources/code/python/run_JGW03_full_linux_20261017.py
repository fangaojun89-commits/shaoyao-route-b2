"""Thin full-scope driver for the frozen JGW03 Linux R entry.

Default CLI produces a plan summary only. Execution requires an independently
assembled full acceptance with fresh resources, and rechecks resources between
bounded rounds. At most three single-core cell processes; no native retry,
installation, source mutation, Windows checkpoint import or empirical-BH change.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime as dt
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path("/root/autodl-tmp/JGW03_full_reference_linux_20261017r03")
CONFIG_SHA = "eb1a72bbecd61878e61a073b025480a0683c6be76c135257ec2e7943f30bf01c"
ENTRY_SHA = "c0c610da98b6511689b46bb5e3a98e87b64b02fa1ed978718db67b412591dd2d"
ENTRY_REL = "code/R/run_JGW03_full_reference_linux_v3_20261017.R"
RUNTIME = Path("/root/autodl-tmp/route_b2_public_analysis/runtime/bin/Rscript")
CELLS = ("GSE289929", "GSE227494", "GSE299063")
MEMORY_PER_WORKER = 3 * 1024**3
BATCH = 50
AFFINITY_EXEC = "import os,sys;os.sched_setaffinity(0,{int(sys.argv[1])});os.execv(sys.argv[2],sys.argv[2:])"


def need(ok, category):
    if not ok:
        raise RuntimeError(category)


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def pin(path):
    p = Path(path)
    need(p.is_file() and not p.is_symlink(), "RegularPinnedSourceRequired")
    h = hashlib.sha256()
    with p.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b""):
            h.update(block)
    return {"path": str(p.resolve()), "bytes": p.stat().st_size, "sha256": h.hexdigest()}


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_new(path, obj):
    with Path(path).open("x", encoding="utf-8") as stream:
        json.dump(obj, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())


def progress(folder, obj):
    temporary = folder / "STATUS.next.json"
    write_new(temporary, {"utc": utc(), "owner_pid": os.getpid(), **obj})
    os.replace(temporary, folder / "STATUS.json")


def plan():
    rounds = []
    for space, total in (("four_group_MC", 10001), ("high_companion_exact", 20)):
        for first in range(1, total + 1, BATCH):
            last = min(total, first + BATCH - 1)
            rounds.append([{"accession": accession, "mode": "run", "space": space,
                            "first": first, "last": last} for accession in CELLS])
    for space in ("four_group_MC", "high_companion_exact"):
        rounds.append([{"accession": accession, "mode": "finish", "space": space,
                        "first": None, "last": None} for accession in CELLS])
    return rounds


def command(job, acceptance_path, acceptance_sha):
    call = [str(RUNTIME), "--vanilla", str(ROOT/"project"/ENTRY_REL), str(ROOT/"CONFIG.json"),
            job["mode"], job["accession"], job["space"]]
    if job["mode"] == "run":
        call += [str(job["first"]), str(job["last"]), str(acceptance_path), acceptance_sha]
    return call


def load_acceptance(path, expected_sha):
    actual = pin(path)
    need(actual["sha256"] == expected_sha and Path(actual["path"]).is_relative_to(ROOT/"provenance"),
         "ExplicitFullAcceptanceIdentity")
    value = read(path)
    need(value["schema"] == "jgw03-linux-full-execution-acceptance-v2"
         and value["status"] == "ACCEPTED_FOR_FULL_LINUX_REFERENCE_V2"
         and value["config_sha256"] == CONFIG_SHA and value["entry_sha256"] == ENTRY_SHA
         and value["resources_admitted"] is True and value["independent_bridge_acceptance"] is True
         and value["full_MC_family_n"] == 89 and value["full_companion_family_n"] == 34
         and value["MC_nonidentity_per_cell"] == 10000 and value["companion_per_cell"] == 20
         and set(value["cells"]) == set(CELLS), "FrozenCompleteScopeRequired")
    need(pin(ROOT/"CONFIG.json")["sha256"] == CONFIG_SHA
         and pin(ROOT/"project"/ENTRY_REL)["sha256"] == ENTRY_SHA, "FrozenEntryOrConfigChanged")
    driver = value["execution_tools"]["driver"]
    helper = value["execution_tools"]["admission"]
    need(pin(__file__) == driver and Path(driver["path"]).parent == ROOT/"execution_tools_01",
         "ReviewedDriverIdentity")
    need(pin(helper["path"]) == helper and Path(helper["path"]).parent == ROOT/"execution_tools_01"
         and Path(helper["path"]).name == "admit_JGW03_full_linux_20261017.py", "ReviewedAdmissionIdentity")
    spec = importlib.util.spec_from_file_location("jgw_full_admission_frozen", helper["path"])
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return actual, value, module


def process_identity(pid):
    text = (Path("/proc")/str(pid)/"stat").read_text()
    fields = text[text.rfind(")")+2:].split()
    return {"pid": pid, "parent": int(fields[1]), "start_ticks": int(fields[19])}


def run_child(job, folder, cpu, acceptance, acceptance_sha):
    folder.mkdir()
    r_call = command(job, acceptance, acceptance_sha)
    # Avoid Python preexec_fn after fork from a multi-threaded controller.
    # A fresh isolated standard-library interpreter pins its own CPU and execs
    # the exact original R command in the SAME PID before any R code runs.
    call = [sys.executable, "-I", "-B", "-c", AFFINITY_EXEC, str(cpu), *r_call]
    env = dict(os.environ,
        R_LIBS_USER="/root/autodl-tmp/GSE128623_R_overlay_20260920/library:/root/autodl-tmp/route_b2_public_analysis/R-library",
        OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", LANG="C", LC_ALL="C",
        TMPDIR=str(ROOT/"tmp"))
    write_new(folder/"COMMAND.json", {"argv": call, "R_argv": r_call, "job": job, "assigned_CPU": cpu,
        "utc": utc(), "driver": pin(__file__), "Rscript": pin(RUNTIME),
        "config_sha256": CONFIG_SHA, "entry_sha256": ENTRY_SHA,
        "acceptance_sha256": acceptance_sha, "fixed_thread_environment":
        {key: env[key] for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "LANG", "LC_ALL", "R_LIBS_USER", "TMPDIR")}})
    child = subprocess.Popen(call, cwd=ROOT, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    start = time.monotonic()
    peak = 0
    observer = None
    collected = False
    stdout, stderr = b"", b""
    identity = None
    try:
        identity = process_identity(child.pid)
        write_new(folder/"STARTED.json", {**identity, "utc": utc(), "job": job,
            "config_sha256": CONFIG_SHA, "entry_sha256": ENTRY_SHA})
        while True:
            try:
                text = (Path("/proc")/str(child.pid)/"status").read_text()
                peak = max([peak] + [int(m[1])*1024 for m in re.finditer(r"^Vm(?:HWM|RSS):\s+(\d+) kB$", text, re.M)])
            except (FileNotFoundError, ProcessLookupError):
                pass
            try:
                stdout, stderr = child.communicate(timeout=1)
                collected = True
                break
            except subprocess.TimeoutExpired:
                pass
    except BaseException as exc:
        observer = type(exc).__name__
    finally:
        # Drain the already-owned native process; never kill a model or refit it
        # because observation failed or another cell returned nonzero.
        if not collected:
            stdout, stderr = child.communicate()
        child.wait()
    import resource
    # RUSAGE_CHILDREN is a conservative process-wide maximum, not a per-cell
    # attribution. /proc HWM is collected separately above.
    children_max = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss*1024
    result = {"status": "ACTUAL_BOUNDED_EXIT_NOT_INDEPENDENT_NUMERICAL_ACCEPTANCE",
        "job": job, "actual_R_pid": child.pid,
        "execution_kind": "AFFINITY_BOOTSTRAP_EXECV_SAME_PID",
        "actual_native_exit_code": child.returncode,
        "actual_R_exit_code": 0 if child.returncode == 0 else None,
        "nonzero_exit_stage_unknown_until_diagnosed": child.returncode != 0,
        "start_ticks": identity["start_ticks"] if identity else None,
        "finished_utc": utc(), "elapsed_seconds": time.monotonic()-start,
        "sampled_child_peak_RSS_bytes": peak, "controller_children_max_RSS_bytes": children_max,
        "memory_budget_met": max(peak, children_max) <= MEMORY_PER_WORKER,
        "observer_failure_class_only": observer, "stdout_bytes": len(stdout), "stderr_bytes": len(stderr),
        "stdout_sha256": hashlib.sha256(stdout).hexdigest(), "stderr_sha256": hashlib.sha256(stderr).hexdigest(),
        "raw_streams_persisted": False, "automatic_native_retry": False}
    write_new(folder/"PROCESS_EXIT.json", result)
    return result


def accepted_child(result):
    return (result["actual_R_exit_code"] == 0 and result["memory_budget_met"] is True
            and result["observer_failure_class_only"] is None)


def run_round(jobs, folder, cpus, acceptance, acceptance_sha, runner=run_child):
    need(len(jobs) == len(cpus) == 3 and len(set(cpus)) == 3
         and {job["accession"] for job in jobs} == set(CELLS), "AllThreeCellsAndDistinctCoresRequired")
    folder.mkdir()
    results = []
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures = [pool.submit(runner, job, folder/job["accession"], cpu, acceptance, acceptance_sha)
                   for job, cpu in zip(jobs, cpus)]
        for future in futures:
            try:
                results.append(future.result())
            except BaseException as exc:
                results.append({"actual_R_exit_code": None, "memory_budget_met": False,
                    "observer_failure_class_only": type(exc).__name__, "remote_state_requires_inspection": True})
    # The pool has joined ALL already-launched cell jobs before caller decides
    # whether a subsequent round is allowed. Never rerun a failed child here.
    write_new(folder/"ROUND_EXIT.json", {"utc": utc(), "results": results,
        "all_actual_exits_zero_and_observed": all(accepted_child(x) for x in results)})
    return results


def run(acceptance, acceptance_sha, attempt):
    need(sys.platform.startswith("linux"), "LinuxExecutionOnly")
    need(re.fullmatch(r"[A-Za-z0-9_-]+", attempt) is not None, "ExplicitFreshAttempt")
    import fcntl
    folder = ROOT/"provenance"/attempt
    need(not folder.exists(), "AttemptAlreadyExistsInspectDoNotReplay")
    actual, value, helper = load_acceptance(acceptance, acceptance_sha)
    cpus = sorted(os.sched_getaffinity(0))
    need(len(cpus) >= 3, "ThreePermittedSingleCoresRequired")
    cpus = cpus[-3:]
    with (ROOT/"provenance/FULL_REFERENCE.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        resource_state = helper.inspect_resources(workers=3)
        need(resource_state["resources_admitted"] is True, "InitialResourceWindowUnavailable")
        folder.mkdir()
        write_new(folder/"CONTROLLER_STARTED.json", {**process_identity(os.getpid()), "utc": utc(),
            "driver": pin(__file__), "acceptance": actual, "python": sys.version,
            "executable": sys.executable, "argv": sys.argv, "cwd": str(Path.cwd()), "assigned_CPUs": cpus})
        write_new(folder/"PLAN.json", {"rounds": plan(), "batch": BATCH, "workers": 3,
            "full_MC_family_n": 89, "full_companion_family_n": 34,
            "per_cell_MC_including_observed": 10001, "per_cell_companion": 20,
            "native_checkpoint_reuse": "unchanged_R_entry_checks_same_Linux_identity",
            "provisional_P_reported": False, "new_outer_empirical_BH": False})
        write_new(folder/"INITIAL_RESOURCES.json", resource_state)
        completed = 0
        try:
            for number, jobs in enumerate(plan(), 1):
                while True:
                    load_acceptance(acceptance, acceptance_sha)  # same immutable binding; no fits
                    state = helper.inspect_resources(workers=3)
                    if state["resources_admitted"] is True:
                        break
                    progress(folder, {"status": "WAIT_RESOURCE_WINDOW_NO_NATIVE_PROCESS_RESTART",
                        "next_round": number, "resource_state": state})
                    time.sleep(10)
                target = folder/f"round_{number:04d}"
                progress(folder, {"status": "RUNNING_BOUND_CELL_ROUND", "round": number,
                                  "completed_rounds": completed})
                results = run_round(jobs, target, cpus, acceptance, acceptance_sha)
                need(all(accepted_child(x) for x in results), "NativeRoundFailedNoAutomaticRetry")
                completed = number
            result = {"status": "COMPLETE_FULL_SCOPE_EXECUTION_NOT_INDEPENDENT_NUMERICAL_ACCEPTANCE",
                      "actual_controller_exit_code": 0, "completed_rounds": completed, "utc": utc(),
                      "full_MC_family_n": 89, "full_companion_family_n": 34,
                      "every_cell_MC_including_observed": 10001, "every_cell_companion": 20,
                      "independent_numerical_acceptance": False, "new_outer_empirical_BH": False}
            write_new(folder/"DISPATCH_EXIT.json", result)
            progress(folder, result)
            return 0
        except BaseException as exc:
            result = {"status": "STOP_RETAINED_NO_AUTOMATIC_NATIVE_RETRY", "actual_controller_exit_code": 1,
                      "completed_rounds": completed, "utc": utc(), "exception_class_only": type(exc).__name__,
                      "all_prior_checkpoints_retained": True, "raw_error_persisted": False}
            write_new(folder/"DISPATCH_EXIT.json", result)
            return 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--acceptance", type=Path)
    parser.add_argument("--acceptance-sha256")
    parser.add_argument("--attempt")
    args = parser.parse_args()
    if not args.execute:
        jobs = plan()
        print(json.dumps({"status": "PLAN_ONLY_NO_RESOURCE_ACCEPTANCE_OR_NATIVE_CALL", "rounds": len(jobs),
            "maximum_parallel_single_core_cells": 3, "batch": BATCH,
            "MC_allocations_per_cell": 10001, "exact_companion_per_cell": 20,
            "full_MC_family_n": 89, "full_companion_family_n": 34,
            "saved_only_finish_commands": 6, "source": pin(__file__)}))
        return 0
    need(args.acceptance is not None and args.acceptance_sha256 is not None and args.attempt is not None,
         "ActualFullAcceptanceAndFreshAttemptRequired")
    return run(args.acceptance, args.acceptance_sha256, args.attempt)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(json.dumps({"status": "NO_NEW_DISPATCH", "exception_class_only": type(exc).__name__,
                          "raw_error_persisted": False}))
        raise SystemExit(2)
