"""One finite STAR admit/index reservation at a real, closed JGW round boundary.

Plan-only by default. Standard library only; no fitting, downloads, deletion,
old-source edits, or old native-process signals. The detached lease is armed
before STOP and owns one NEW session containing only cache/STAR prefix calls.
Both the caller and lease resume the exact old Python controller on exit. If
the caller dies, the lease remains and restores scheduling by its deadline.
Existing partial/failed STAR stage directories are never retried or replaced.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
import threading
import time

sys.dont_write_bytecode = True
GIB = 1024**3
INSTANCE = "62c8pzy536-50fb405a"
HOST = "autodl-container-" + INSTANCE
JGW = Path("/root/autodl-tmp/JGW03_full_reference_linux_20261017r03")
G128 = Path("/root/autodl-tmp/gse128623_label_reference_20261015r01")
STAR = Path("/root/autodl-tmp/GSE205906_matched_STAR_Ensembl101_20261015r01")
DRIVER = JGW / "execution_tools_02/run_JGW03_full_linux_20261017.py"
DRIVER_SHA = "8f5edfa505e44c6fbf36c605e8b582dfe23c73bb29b547645644003980fec4a0"
STAR_SHA = "7fe64fcfcee554e1500f60da2cd7b4a3ce8a84fcd4e46b84c3c7e14ec1fd907d"
META_SHA = "6dd7b53d40433baaacdeedf285204dee83f2d960c090e6a25f489f96d20fc5cb"
CONTROLLER, BIRTH = 659622, 2576253226
BASE = JGW / "provenance/full_reference_02"
CELLS = ("GSE289929", "GSE227494", "GSE299063")
ABORT = threading.Event()


def need(ok, label):
    if not ok:
        raise RuntimeError(label)


def integer(value, minimum=0):
    return type(value) is int and value >= minimum


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def read(path):
    return json.loads(Path(path).read_bytes())


def published(path):
    # Handshake records are exclusive and tiny, but existence can precede the
    # final write/fsync. Do not mistake a writer's open file for corrupt JSON.
    for _ in range(50):
        try:
            return read(path)
        except json.JSONDecodeError:
            time.sleep(.02)
    raise RuntimeError("SmallHandshakeNotCompletelyWritten")


def put(path, value):
    with Path(path).open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n"); handle.flush(); os.fsync(handle.fileno())


def pin(path):
    p = Path(path)
    need(p.is_file() and not p.is_symlink() and p.resolve() == p, "RegularAbsolutePinnedFile")
    digest = hashlib.sha256()
    with p.open("rb") as handle:
        for block in iter(lambda: handle.read(1024**2), b""):
            digest.update(block)
    return {"path": str(p), "bytes": p.stat().st_size, "sha256": digest.hexdigest()}


def stat_process(pid):
    p = Path("/proc") / str(pid)
    fields = (p / "stat").read_text().rsplit(") ", 1)[1].split()
    return {"pid": pid, "birth": int(fields[19]), "state": fields[0],
            "parent": int(fields[1]), "pgrp": int(fields[2]), "session": int(fields[3])}


def proc(pid):
    value = stat_process(pid)
    value["argv"] = [os.fsdecode(x) for x in (Path("/proc") / str(pid) / "cmdline").read_bytes().rstrip(b"\0").split(b"\0")]
    return value


def children(pid):
    result = set()
    tasks = list((Path("/proc") / str(pid) / "task").iterdir())
    need(tasks, "VisibleControllerTaskList")
    for task in tasks:
        result.update(int(x) for x in (task / "children").read_text().split())
    return sorted(result)


def match_controller(current, binding):
    need(current["pid"] == CONTROLLER and integer(current["birth"], 1)
         and current["birth"] == BIRTH == binding["birth"]
         and current["argv"] == binding["argv"]
         and current["state"] not in ("Z", "X"), "ExactOriginalControllerIdentity")
    return current


def boundary(status, number, individual, whole):
    need(integer(number, 56) and number <= 200, "FrozenMCBoundary56Through200")
    need(type(status["owner_pid"]) is int and status["owner_pid"] == CONTROLLER
         and status["status"] == "WAIT_RESOURCE_WINDOW_NO_NATIVE_PROCESS_RESTART"
         and type(status.get("next_round")) is int and status["next_round"] == number,
         "RealExactNextRoundWait")
    previous = number - 1
    terminals = []
    for cell in CELLS:
        ended, started, command = individual[cell]
        job = dict(accession=cell, first=(previous - 1)*50 + 1, last=previous*50,
                   mode="run", space="four_group_MC")
        need(all(type(ended[k]) is int and ended[k] == 0 for k in
                 ("actual_R_exit_code", "actual_native_exit_code"))
             and ended["memory_budget_met"] is True
             and ended["observer_failure_class_only"] is None
             and integer(started["pid"], 1) and integer(started["start_ticks"], 1)
             and type(started["parent"]) is int and started["parent"] == CONTROLLER
             and type(ended["actual_R_pid"]) is int and ended["actual_R_pid"] == started["pid"]
             and type(ended["start_ticks"]) is int and ended["start_ticks"] == started["start_ticks"]
             and ended["job"] == started["job"] == command["job"] == job,
             "PreviousThreeIndividualStrictNativeSuccess")
        need(command["driver"]["path"] == str(DRIVER)
             and command["driver"]["sha256"] == DRIVER_SHA, "PreviousCommandDriverBinding")
        terminals.append(ended)
    need(whole["all_actual_exits_zero_and_observed"] is True
         and len(whole["results"]) == 3
         and sorted(whole["results"], key=lambda x: x["job"]["accession"]) ==
             sorted(terminals, key=lambda x: x["job"]["accession"]),
         "WholeRoundEqualsAllThreeActualReceipts")


class Linux:
    now = staticmethod(time.monotonic)
    sleep = staticmethod(time.sleep)
    proc = staticmethod(proc)
    children = staticmethod(children)

    def signal_controller(self, binding, sig):
        # A pidfd prevents the check/signal PID-reuse race; no process-group signal.
        fd = os.pidfd_open(CONTROLLER)
        try:
            match_controller(self.proc(CONTROLLER), binding)
            signal.pidfd_send_signal(fd, sig)
        finally:
            os.close(fd)


def resume(binding, ops):
    current = match_controller(ops.proc(CONTROLLER), binding)
    sent = current["state"] in ("T", "t")
    if sent:
        ops.signal_controller(binding, signal.SIGCONT)
        end = ops.now() + 2
        while ops.now() < end:
            current = match_controller(ops.proc(CONTROLLER), binding)
            if current["state"] not in ("T", "t"):
                break
            ops.sleep(.02)
        need(current["state"] not in ("T", "t"), "OriginalControllerResumeNotYetObserved")
    return {"utc": utc(), "SIGCONT_sent": sent, "pid": CONTROLLER, "birth": BIRTH,
            "observed_state_after": current["state"]}


def pause(binding, ops, recheck):
    current = match_controller(ops.proc(CONTROLLER), binding)
    need(current["state"] not in ("T", "t"), "DoNotTakeOverSomebodyElsesPause")
    need(not ops.children(CONTROLLER), "NoOldControllerChildrenBeforeStop")
    recheck()
    # Caller has already armed the independent lease and durable ownership flag.
    ops.signal_controller(binding, signal.SIGSTOP)
    try:
        end = ops.now() + 5
        while ops.now() < end:
            state = match_controller(ops.proc(CONTROLLER), binding)
            if state["state"] in ("T", "t"):
                need(not ops.children(CONTROLLER), "StopRaceNewChildAbandonWindow")
                recheck()
                return state
            ops.sleep(.05)
        raise RuntimeError("ControllerDidNotEnterStoppedState")
    except BaseException:
        resume(binding, ops)
        raise


def paused(binding, number, ops):
    current = match_controller(ops.proc(CONTROLLER), binding)
    need(current["state"] in ("T", "t") and not ops.children(CONTROLLER),
         "OnlyStoppedChildlessOriginalScheduler")
    check_boundary(number)


def check_boundary(number):
    individual = {}
    for cell in CELLS:
        folder = BASE / ("round_%04d" % (number - 1)) / cell
        individual[cell] = [read(folder / name) for name in
                            ("PROCESS_EXIT.json", "STARTED.json", "COMMAND.json")]
    boundary(read(BASE / "STATUS.json"), number, individual,
             read(BASE / ("round_%04d" % (number - 1)) / "ROUND_EXIT.json"))


def controller_binding(ops):
    need(pin(DRIVER)["sha256"] == DRIVER_SHA, "ReviewedOriginalDriverSource")
    started = read(BASE / "CONTROLLER_STARTED.json")
    current = ops.proc(CONTROLLER)
    need(type(started["pid"]) is int and started["pid"] == CONTROLLER
         and type(started["start_ticks"]) is int and started["start_ticks"] == BIRTH
         and started["driver"] == pin(DRIVER)
         and started["argv"][0] == str(DRIVER)
         and "--execute" in started["argv"], "OriginalControllerLaunchBinding")
    binding = dict(pid=CONTROLLER, birth=BIRTH,
                   argv=[started["executable"], "-I", "-B", *started["argv"]],
                   started=pin(BASE / "CONTROLLER_STARTED.json"), driver=pin(DRIVER))
    match_controller(current, binding)
    return binding


def tool(path, expected, root):
    p = Path(path)
    need(p.is_absolute() and p.is_relative_to(root), "FixedProjectToolDirectory")
    value = pin(p)
    need(value["sha256"] == expected, "ReviewedToolSourceSHA")
    return value


def require_sources(intent):
    for key in ("scheduler_source", "STAR_source", "cache_source"):
        need(pin(intent[key]["path"]) == intent[key], "WindowSourceUnchanged")
    need(pin(DRIVER) == intent["controller"]["driver"]
         and pin(BASE / "CONTROLLER_STARTED.json") == intent["controller"]["started"],
         "OriginalSourceAndLaunchReceiptUnchanged")


def headroom():
    c = Path("/sys/fs/cgroup")
    limit, used = [int((c / name).read_text()) for name in ("memory.max", "memory.current")]
    need(limit == 96*GIB, "FixedNewInstance96GiBQuota")
    return {"limit": limit, "used": used, "direct_headroom": limit-used}


def checkpoint(stage, root=STAR):
    """Small receipt/stat reuse, not a second hash pass over genome products."""
    folder = root / {"admit": "reference_admission_01", "index": "index_01"}[stage]
    if not folder.exists():
        return None
    need(folder.is_dir() and not folder.is_symlink() and not (folder / "FAILED.json").exists(),
         "FailedStageRetainedNoRetry")
    r = read(folder / "COMPLETE.json")  # Missing COMPLETE also refuses a partial leaf.
    need(r["status"] == "COMPLETE" and r["source_sha256"] == STAR_SHA
         and r["metadata_sha256"] == META_SHA, "BoundOriginalStageComplete")
    need(r.get("events") and all(type(e.get("actual_exit")) is int and e["actual_exit"] == 0
                                for e in r["events"]), "RecordedOriginalNativeCommandsAllStrictZero")
    for item in r.get("small_products", []):
        need(pin(item["path"]) == item, "SavedSmallStageProductUnchanged")
    if stage == "admit":
        for record in r["reference"]["reference_stat"]:
            p = Path(record["path"]); s = p.lstat()
            need(p.is_file() and not p.is_symlink()
                 and {"path": str(p), "bytes": s.st_size, "mtime_ns": s.st_mtime_ns,
                      "inode": s.st_ino} == record, "SavedAdmittedReferenceStatUnchanged")
    if stage == "index":
        need(r["reference_admission"] == pin(root / "reference_admission_01/COMPLETE.json"),
             "SavedIndexAdmissionBinding")
        paths = set()
        products = {x["path"]: x for x in r["index_files"]}
        need(len(products) == len(r["index_files"]) == len(r["index_stat"]),
             "SavedIndexProductAndStatDomain")
        for record in r["index_stat"]:
            p = Path(record["path"]); s = p.lstat()
            need(p.is_relative_to(folder / "index") and p.is_file() and not p.is_symlink()
                 and {"path": str(p), "bytes": s.st_size, "mtime_ns": s.st_mtime_ns,
                      "inode": s.st_ino} == record and record["path"] in products
                 and products[record["path"]]["bytes"] == s.st_size > 0
                 and isinstance(products[record["path"]]["sha256"], str)
                 and len(products[record["path"]]["sha256"]) == 64,
                 "ExistingIndexProductStatBinding")
            paths.add(p.name)
        need({"Genome", "SA", "SAindex", "genomeParameters.txt"} <= paths,
             "CompleteExistingIndexDomain")
    return pin(folder / "COMPLETE.json")


def commands(intent):
    cache = [sys.executable, "-I", "-B", intent["cache_source"]["path"],
             "--execute", "--expected-round", str(intent["expected_round"])]
    stages = {stage: [sys.executable, "-I", "-B", intent["STAR_source"]["path"], stage,
              "--execute", "--quiet-window-confirmed", "--reserved-growth-bytes",
              str(intent["reserved_growth_bytes"]), "--disk-budget-bytes",
              str(intent[stage + "_disk_budget_bytes"])] for stage in ("admit", "index")}
    return cache, stages


def lease_ok(intent, ops):
    need(not ABORT.is_set() and ops.now() < intent["deadline_monotonic"], "FiniteWindowDeadlineOrSignal")
    parent = ops.proc(intent["caller"]["pid"])
    need(parent["birth"] == intent["caller"]["birth"] and parent["state"] not in ("Z", "X"),
         "CallerGoneLeaseRestoresScheduling")


def group_members(session):
    members = []
    for p in Path("/proc").iterdir():
        if p.name.isdigit():
            try:
                item = stat_process(int(p.name))
                if item["session"] == session:
                    need(item["pgrp"] == session, "OwnedSessionHasOnlyExplicitGroup")
                    members.append(item)
            except (FileNotFoundError, ProcessLookupError):
                pass
    return members


def owned_group(identity, current):
    need(integer(identity["pid"], 1) and integer(identity["birth"], 1)
         and current["pid"] == identity["pid"] and current["birth"] == identity["birth"]
         and current["pgrp"] == current["session"] == current["pid"], "ExactNewOwnedSessionOnly")
    return current["pid"]


def stop_group(identity, binding, ops):
    """Lease owns this session; never uses the old controller's process group."""
    sid = owned_group(identity, ops.proc(identity["pid"]))
    others = [x for x in group_members(sid) if x["pid"] != os.getpid() and x["state"] not in ("Z", "X")]
    if not others:
        return
    os.killpg(sid, signal.SIGTERM)
    end = ops.now() + 10
    while ops.now() < end:
        live = [x for x in group_members(sid) if x["pid"] != os.getpid() and x["state"] not in ("Z", "X")]
        if not live:
            return
        ops.sleep(.1)
    # If the caller died, the lease must restore scheduling BEFORE killing its
    # own session. This final SIGKILL includes only our new lease/cache/STAR tree.
    resume(binding, ops)
    owned_group(identity, ops.proc(identity["pid"]))
    os.killpg(sid, signal.SIGKILL)


def run_owned(command, intent, folder, name, ops):
    lease_ok(intent, ops); paused(intent["controller"], intent["expected_round"], ops)
    require_sources(intent)
    # The detached lease is already its OWN session/group. All new commands
    # inherit it; no other project process can enter this group by accident.
    child = subprocess.Popen(command, cwd=STAR, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
    streams = {}
    def drain(label, handle):
        digest = hashlib.sha256(); size = 0
        try:
            for block in iter(lambda: handle.read(1024**2), b""):
                digest.update(block); size += len(block)
            streams[label] = {"bytes": size, "sha256": digest.hexdigest()}
        except BaseException as exc:
            streams[label] = {"observer_failure_class_only": type(exc).__name__}
    threads = [threading.Thread(target=drain, args=(k, h), daemon=True)
               for k, h in (("stdout", child.stdout), ("stderr", child.stderr))]
    for t in threads:
        t.start()
    identity, failure = None, None
    try:
        identity = proc(child.pid)
        put(folder / (name + "_STARTED.json"), {"utc": utc(), "command": command,
            "pid": child.pid, "birth": identity["birth"], "pgrp": identity["pgrp"],
            "session": identity["session"], "lease_group": os.getpgrp()})
        while True:
            lease_ok(intent, ops)
            try:
                child.wait(timeout=.25)
                break
            except subprocess.TimeoutExpired:
                pass
    except BaseException as exc:
        failure = exc
    finally:
        if failure is not None:
            stop_group(ops.proc(os.getpid()), intent["controller"], ops)
        actual = child.wait(timeout=5)
        for t in threads:
            t.join(5)
        result = {"utc": utc(), "command": command, "pid": child.pid,
                  "birth": identity["birth"] if identity else None,
                  "actual_exit_code": actual, "streams": streams,
                  "observer_failure_class_only": type(failure).__name__ if failure else None,
                  "raw_streams_persisted": False}
        put(folder / (name + "_EXIT.json"), result)
    if failure is not None:
        raise failure
    need(all(not t.is_alive() for t in threads) and
         all("observer_failure_class_only" not in x for x in streams.values())
         and len(streams) == 2, "OwnedCommandStreamObservationComplete")
    need(type(actual) is int and actual == 0, "OwnedCommandNonzeroNoAutomaticRetry")
    require_sources(intent)
    return result


def prefix(intent, folder, ops, run=run_owned, saved=checkpoint, memory=headroom):
    cache, stages = commands(intent)
    reused = {stage: saved(stage) for stage in stages}
    if "existing_checkpoints" in intent:
        need(reused == intent["existing_checkpoints"], "CheckpointBindingUnchangedSinceReservation")
    if all(reused.values()):
        return {"stages": reused, "cache_called": False, "all_existing_prefix_reused": True}
    before = memory()
    if before["direct_headroom"] < 40*GIB:
        run(cache, intent, folder, "CACHE", ops)
    need(memory()["direct_headroom"] >= 40*GIB, "ActualFortyGiBDirectHeadroomRequired")
    results = {}
    for stage, command in stages.items():
        lease_ok(intent, ops); paused(intent["controller"], intent["expected_round"], ops)
        if reused[stage]:
            results[stage] = {"reused": True, "COMPLETE": reused[stage]}
        else:
            need(memory()["direct_headroom"] >= 40*GIB, "FreshStageFortyGiBDirectHeadroom")
            actual = run(command, intent, folder, stage.upper(), ops)
            results[stage] = {"reused": False, "native": actual, "COMPLETE": saved(stage)}
            need(results[stage]["COMPLETE"] is not None, "ActualPrefixStageReceiptRequired")
    return {"stages": results, "cache_called": before["direct_headroom"] < 40*GIB,
            "new_models": 0, "libraries_run": 0}


def lease(folder, expected_sha):
    ops = Linux(); intent_path = folder / "INTENT.json"
    need(pin(intent_path)["sha256"] == expected_sha, "ExactWindowIntent")
    intent = read(intent_path); require_sources(intent)
    need(intent["schema"] == "gse205906-star-prefix-scheduler-v1" and intent["instance"] == INSTANCE
         and integer(intent["expected_round"], 56) and intent["expected_round"] <= 200
         and intent["allowed_stages"] == ["admit", "index"]
         and folder == STAR / ("prefix_window_round_%04d_20261021_r01" % intent["expected_round"])
         and isinstance(intent["deadline_monotonic"], (int, float))
         and 0 < intent["deadline_monotonic"] - ops.now() <= 4*60*60,
         "FixedFiniteLeaseNamespaceAndScope")
    need(Path(__file__).resolve() == Path(intent["scheduler_source"]["path"]), "SameLeaseEntry")
    need(socket.gethostname() == HOST and sys.platform.startswith("linux"), "OnlyNewInstanceLinux")
    identity = ops.proc(os.getpid()); owned_group(identity, identity)
    signal.signal(signal.SIGTERM, lambda *_: ABORT.set())
    signal.signal(signal.SIGINT, lambda *_: ABORT.set())
    result, code = {}, 1
    try:
        lease_ok(intent, ops)
        put(folder / "LEASE_READY.json", {"utc": utc(), **identity,
            "deadline_monotonic": intent["deadline_monotonic"]})
        while not (folder / "GO.json").exists():
            lease_ok(intent, ops)
            need(not (folder / "ABANDON.json").exists(), "CallerAbandonedWindow")
            ops.sleep(.1)
        paused(intent["controller"], intent["expected_round"], ops)
        result = prefix(intent, folder, ops)
        result["status"] = "COMPLETE_SAVED_OR_NATIVE_STAR_ADMIT_INDEX_ONLY"
        code = 0
    except BaseException as exc:
        result = {"status": "STOP_RETAINED_NO_AUTOMATIC_RETRY",
                  "exception_class_only": type(exc).__name__, "raw_error_persisted": False}
    finally:
        # Attempt restoration even if cleanup, evidence writing, or a command
        # failed. The caller provides a second independent restoration path.
        try:
            stop_group(identity, intent["controller"], ops)
        finally:
            if (folder / "PAUSE_OWNED.json").exists():
                result["controller_resume"] = resume(intent["controller"], ops)
        result.update(actual_lease_function_exit=code, utc=utc(), scope="PREFIX_ONLY_NO_MODEL_OR_LIBRARY")
        put(folder / "LEASE_RESULT.json", result)
    return code


def schedule(args, ops=None):
    ops = ops or Linux()
    need(sys.platform.startswith("linux") and socket.gethostname() == HOST, "OnlyDeclaredNewInstance")
    need(integer(args.expected_round, 56) and args.expected_round <= 200
         and integer(args.max_stop_seconds, 60) and args.max_stop_seconds <= 4*60*60,
         "FiniteClosedMCWindowAtMostFourHours")
    need(integer(args.reserved_growth_bytes) and integer(args.admit_disk_budget_bytes, 1)
         and integer(args.index_disk_budget_bytes, 1), "ExplicitCurrentStageResourceLedger")
    binding = controller_binding(ops); check_boundary(args.expected_round)
    need(ops.proc(CONTROLLER)["state"] not in ("T", "t"), "NoPreexistingPauseOwnership")
    need(not ops.children(CONTROLLER), "NoControllerChildrenAtReservation")
    source = pin(Path(__file__).resolve())
    intent = {"schema": "gse205906-star-prefix-scheduler-v1", "utc": utc(), "instance": INSTANCE,
        "expected_round": args.expected_round, "controller": binding,
        "caller": ops.proc(os.getpid()), "scheduler_source": source,
        "STAR_source": tool(args.star_tool, STAR_SHA, STAR),
        "cache_source": tool(args.cache_tool, args.cache_sha, G128),
        "reserved_growth_bytes": args.reserved_growth_bytes,
        "admit_disk_budget_bytes": args.admit_disk_budget_bytes,
        "index_disk_budget_bytes": args.index_disk_budget_bytes,
        # Reserve 30 seconds INSIDE the requested maximum for cleanup/resume.
        "deadline_monotonic": ops.now() + args.max_stop_seconds - 30,
        "maximum_controller_stop_seconds": args.max_stop_seconds, "new_models": 0,
        "allowed_stages": ["admit", "index"], "JGW_source_or_thresholds_changed": False}
    # Reject retained partial/failed stages BEFORE pause and bind reuse receipts.
    intent["existing_checkpoints"] = {stage: checkpoint(stage) for stage in ("admit", "index")}
    folder = STAR / ("prefix_window_round_%04d_20261021_r01" % args.expected_round)
    folder.mkdir(exist_ok=False)
    put(folder / "INTENT.json", intent)
    guard = None; guard_identity = None; result = {"status": "STOP_RETAINED_NO_RETRY"}; code = 1
    def interrupted(*_):
        raise InterruptedError("WindowCallerSignal")
    old_handlers = {sig: signal.signal(sig, interrupted) for sig in (signal.SIGINT, signal.SIGTERM)}
    try:
        guard = subprocess.Popen([sys.executable, "-I", "-B", source["path"], "_lease",
            str(folder), pin(folder / "INTENT.json")["sha256"]], start_new_session=True,
            stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
        guard_identity = ops.proc(guard.pid)
        end = min(ops.now() + 20, intent["deadline_monotonic"])
        while not (folder / "LEASE_READY.json").exists():
            need(ops.now() < end and ops.proc(guard.pid)["state"] not in ("Z", "X"),
                 "IndependentLeaseMustBeArmedBeforePause")
            ops.sleep(.05)
        ready = published(folder / "LEASE_READY.json")
        owned_group(ready, ops.proc(guard.pid))
        need(ready["birth"] == guard_identity["birth"], "ActualIndependentLeaseIdentity")
        put(folder / "PAUSE_OWNED.json", {"utc": utc(), "controller": binding,
             "independent_lease": ready, "restoration_owned_before_SIGSTOP": True})
        pause(binding, ops, lambda: check_boundary(args.expected_round))
        put(folder / "GO.json", {"utc": utc(), "paused_controller": ops.proc(CONTROLLER)})
        while not (folder / "LEASE_RESULT.json").exists():
            need(ops.now() < intent["deadline_monotonic"] + 15
                 and ops.proc(guard.pid)["state"] not in ("Z", "X"), "LeaseExitedOrDeadline")
            ops.sleep(.2)
        result = published(folder / "LEASE_RESULT.json")
        code = 0 if type(result.get("actual_lease_function_exit")) is int and result["actual_lease_function_exit"] == 0 else 1
    except BaseException as exc:
        result = {"status": "STOP_RETAINED_NO_AUTOMATIC_RETRY",
                  "exception_class_only": type(exc).__name__, "raw_error_persisted": False}
        if not (folder / "GO.json").exists():
            put(folder / "ABANDON.json", {"utc": utc(), "new_science_started_by_caller": False})
    finally:
        try:
            if guard is not None:
                # Do NOT reap a failed lease leader until its own group is
                # cleaned. Even a zombie leader preserves PID/session identity.
                current = ops.proc(guard.pid)
                sid = owned_group(guard_identity, current)
                if (folder / "LEASE_RESULT.json").exists():
                    # RESULT publication precedes the interpreter's native
                    # exit by a few instructions. Allow that normal exit;
                    # don't turn a successful lease into a signalled one.
                    end = ops.now() + 2
                    while current["state"] not in ("Z", "X") and ops.now() < end:
                        ops.sleep(.05)
                        current = ops.proc(guard.pid)
                        owned_group(guard_identity, current)
                if not (folder / "LEASE_RESULT.json").exists() or current["state"] not in ("Z", "X"):
                    os.killpg(sid, signal.SIGTERM)
                    end = ops.now() + 12
                    while ops.now() < end and any(x["state"] not in ("Z", "X") for x in group_members(sid)):
                        ops.sleep(.1)
                    if any(x["state"] not in ("Z", "X") for x in group_members(sid)):
                        owned_group(guard_identity, ops.proc(guard.pid))
                        os.killpg(sid, signal.SIGKILL)
                actual = guard.wait(timeout=5)
                result["actual_independent_lease_exit"] = actual
                if type(actual) is not int or actual != 0:
                    code = 1
        finally:
            try:
                if (folder / "PAUSE_OWNED.json").exists():
                    result["caller_controller_resume"] = resume(binding, ops)
            finally:
                for sig, handler in old_handlers.items():
                    signal.signal(sig, handler)
        result.update(utc=utc(), actual_scheduler_exit=code, source=source,
                      scientific_acceptance=False, original_JGW_CLI_exit_asserted=False)
        put(folder / "SCHEDULER_RESULT.json", result)
    return code


def main():
    if len(sys.argv) == 4 and sys.argv[1] == "_lease":
        return lease(Path(sys.argv[2]), sys.argv[3])
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--expected-round", required=True, type=int)
    parser.add_argument("--star-tool", required=True, type=Path)
    parser.add_argument("--cache-tool", required=True, type=Path)
    parser.add_argument("--cache-sha", required=True)
    parser.add_argument("--reserved-growth-bytes", required=True, type=int)
    parser.add_argument("--admit-disk-budget-bytes", required=True, type=int)
    parser.add_argument("--index-disk-budget-bytes", required=True, type=int)
    parser.add_argument("--max-stop-seconds", type=int, default=4*60*60)
    args = parser.parse_args()
    if not args.execute:
        print(json.dumps({"status": "PLAN_ONLY_NO_PROCESS_OR_CACHE_ACTION", "instance": INSTANCE,
            "allowed_stages": ["admit", "index"], "source": pin(Path(__file__).resolve()),
            "resource_acceptance": False, "maximum_stop_seconds": args.max_stop_seconds}))
        return 0
    return schedule(args)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(json.dumps({"status": "NO_NEW_SCIENCE_OR_RETRY", "exception_class_only": type(exc).__name__,
                          "raw_error_persisted": False}))
        raise SystemExit(2)
