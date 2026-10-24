"""Finite natural-boundary handoff, then a lightweight direct-parent observer.

Default is plan-only. No model is stopped, replayed or retried. A pinned STAR
pidfd/STOP/CONT implementation is reused, with an independent finite recovery
lease armed before pausing ONLY the childless original Python scheduler. The
old controller is retired only after complete material and fresh resources pass.
Its unavailable parent returncode stays null. After retirement the lease ends;
the successor observer waits as long as required, without a STAR time limit.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import re
import select
import signal
import socket
import subprocess
import sys
import threading
import time
import types

ROOT = Path("/root/autodl-tmp/JGW03_full_reference_linux_20261017r03")
HOST = "autodl-container-62c8pzy536-50fb405a"
SAFETY_SHA = "fe9f804c162a26e34c1f5fd52d5c436b8d72bd81fe19a032691a741a1323de53"
SUCCESSOR_SHA = "3f1a3b9cf438ce770c0528d809e4b0c8a3cbcdb37d87e6765b64f16624774b87"
SAFETY = Path("/root/autodl-tmp/GSE205906_matched_STAR_Ensembl101_20261015r01/prefix_scheduler_source_20261021_r01/schedule_GSE205906_STAR_prefix_window_20261021.py")
CACHE = Path("/root/autodl-tmp/gse128623_label_reference_20261015r01/cache_window_after_handoff_source_20261021_r03/advise_GSE128623_closed_cache_after_handoff_20261021.py")
CACHE_SHA = "31f9c503b098741bfaa579c5af0c05385ebd3a3cb540e42637962496d11dbea3"


def need(ok, label):
    if not ok:
        raise RuntimeError(label)


def load(path, sha, name):
    raw = Path(path).read_bytes()
    need(hashlib.sha256(raw).hexdigest() == sha, "ExactReviewedExecutedSource")
    value = types.ModuleType(name); value.__file__ = str(path)
    exec(compile(raw, str(path), "exec"), value.__dict__)
    return value


def build_material(template, number, m):
    """Append only authenticated, already completed old rounds/prefix pins."""
    need(type(number) is int and 67 <= number <= 200 and
         type(template["completed_rounds"]) is int and template["completed_rounds"] <= number-1 and
         template["status"] == "MATERIALS_PREPARED_NOT_HANDOFF_NOT_EXECUTION_ACCEPTANCE" and
         "terminal" not in template["legacy_controller"], "PreparedTemplateAndExplicitLaterBoundary")
    h = copy.deepcopy(template); h["completed_rounds"] = number-1
    base = Path(h["legacy_attempt"])
    h["completed_legacy_rounds"] = []
    for n in range(1, number):
        folder = base / f"round_{n:04d}"
        h["completed_legacy_rounds"].append({"number": n, "whole": m.pin(folder / "ROUND_EXIT.json"),
            "cells": {cell: {key: m.pin(folder / cell / filename) for key, filename in
                (("command", "COMMAND.json"), ("started", "STARTED.json"), ("terminal", "PROCESS_EXIT.json"))}
                for cell in m.CELLS}})
    through = (number-1)*50
    for cell in m.CELLS:
        h["initial_previous"][cell]["through"] = through
        h["initial_previous"][cell]["prefix"] = m.pin(ROOT / "derived" / cell / f"four_group_MC_TAIL_PREFIX_{through:05d}.rds")
    return h


def exact_wait(status, number, controller):
    need(type(status["owner_pid"]) is int and status["owner_pid"] == controller and
         status["status"] == "WAIT_RESOURCE_WINDOW_NO_NATIVE_PROCESS_RESTART" and
         type(status.get("next_round")) is int and status["next_round"] == number, "TrueRequestedNaturalWAITOnly")


def wait_boundary(number, seconds, s, ops):
    need(type(seconds) is int and 0 <= seconds <= 3600, "BoundedNaturalWaitAtMostOneHour")
    binding = s.controller_binding(ops)
    deadline = ops.now() + seconds
    while True:
        s.match_controller(ops.proc(s.CONTROLLER), binding)
        status = s.read(s.BASE / "STATUS.json")
        if status["status"] == "WAIT_RESOURCE_WINDOW_NO_NATIVE_PROCESS_RESTART" and status.get("next_round") == number:
            exact_wait(status, number, s.CONTROLLER); s.check_boundary(number)
            need(not ops.children(s.CONTROLLER), "NoFitChildrenAtNaturalBoundary")
            return binding
        current = status.get("round") if status["status"] == "RUNNING_BOUND_CELL_ROUND" else status.get("next_round")
        need(type(current) is int and current <= number and
             not (status["status"] == "RUNNING_BOUND_CELL_ROUND" and current == number), "RequestedBoundaryAlreadyPassedNoReplay")
        need(ops.now() < deadline, "NaturalBoundaryWaitExpiredNoSignals")
        ops.sleep(min(2, max(.01, deadline-ops.now())))


def safe_resume(binding, s, ops):
    try:
        current = ops.proc(s.CONTROLLER)
    except (FileNotFoundError, ProcessLookupError):
        return {"SIGCONT_sent": False, "reason": "ORIGINAL_IDENTITY_GONE_NO_SIGNAL"}
    if current["birth"] != s.BIRTH or current["state"] in ("Z", "X"):
        return {"SIGCONT_sent": False, "reason": "ORIGINAL_EXITED_OR_PID_REUSED_NO_SIGNAL"}
    return s.resume(binding, ops)


def old_gone(s, ops):
    try:
        item = ops.proc(s.CONTROLLER)
    except (FileNotFoundError, ProcessLookupError):
        return True
    return item["birth"] != s.BIRTH or item["state"] in ("Z", "X")


def source_pins(intent, m):
    for key in ("coordinator_source", "successor_driver", "safety_source", "prepared_template"):
        m.bound(intent[key])
    if "closed_cache_source" in intent:
        m.bound(intent["closed_cache_source"])


def finite_ok(intent, folder, m, s, ops):
    need(ops.now() < intent["deadline_monotonic"], "FinitePauseDeadline")
    source_pins(intent, m)
    ready = s.published(folder / "LEASE_READY.json")
    guard = ops.proc(ready["pid"])
    s.owned_group(ready, guard)
    need(guard["state"] not in ("Z", "X") and not (folder / "LEASE_RESULT.json").exists(), "LiveIndependentRecoveryLease")


def retire(binding, number, s, ops, timeout=20, note=lambda event: None):
    """One pidfd covers TERM, pending-TERM CONT and actual exit observation."""
    fd = ops.pidfd_open(s.CONTROLLER)
    try:
        s.paused(binding, number, ops)
        ops.pidfd_signal(fd, signal.SIGTERM)
        note("TERM_SENT")
        # Linux defers ordinary termination while the task is STOPped. Fixed
        # old driver has the default TERM action; CONT delivers that pending TERM.
        ops.pidfd_signal(fd, signal.SIGCONT)
        note("PENDING_TERM_CONT_SENT")
        need(ops.pidfd_exited(fd, timeout), "OriginalPidfdExitNotObserved")
        note("PIDFD_EXIT_OBSERVED")
    finally:
        ops.close(fd)
    return {"pid": s.CONTROLLER, "start_ticks": s.BIRTH, "actual_native_exit_code": None,
        "exit_observed_after_signal": True, "pidfd_exit_observed": True,
        "termination_signal_sent": int(signal.SIGTERM), "pending_TERM_delivered_via_SIGCONT": True}


class Linux:
    def __init__(self, safety):
        self.original = safety.Linux()
    now = staticmethod(time.monotonic)
    sleep = staticmethod(time.sleep)
    pidfd_open = staticmethod(os.pidfd_open) if hasattr(os, "pidfd_open") else None
    pidfd_signal = staticmethod(signal.pidfd_send_signal) if hasattr(signal, "pidfd_send_signal") else None
    close = staticmethod(os.close)
    @staticmethod
    def pidfd_exited(fd, timeout):
        return bool(select.select([fd], [], [], timeout)[0])
    def proc(self, pid):
        return self.original.proc(pid)
    def children(self, pid):
        return self.original.children(pid)
    def signal_controller(self, binding, sig):
        return self.original.signal_controller(binding, sig)


def stream_readers(child):
    streams = {}; failures = []
    def drain(label, handle):
        h = hashlib.sha256(); n = 0
        try:
            for block in iter(lambda: handle.read(1024**2), b""):
                h.update(block); n += len(block)
        except BaseException as exc:
            failures.append(type(exc).__name__)
        streams[label] = {"bytes": n, "sha256": h.hexdigest()}
    threads = [threading.Thread(target=drain, args=(name, handle), daemon=True)
               for name, handle in (("stdout", child.stdout), ("stderr", child.stderr))]
    for thread in threads:
        thread.start()
    return streams, failures, threads


def prepare_cache_once(number, folder, intent, binding, m, s, ops, guarded):
    """One fixed bounded maintenance child; never a scientific process group."""
    source = m.pin(CACHE)
    need(source["sha256"] == CACHE_SHA and source == intent["closed_cache_source"], "FixedPreviouslyReviewedClosedCacheTool")
    guarded(); s.paused(binding, number, ops)
    argv = [sys.executable, "-I", "-B", source["path"], "--execute", "--expected-round", str(number)]
    m.write_new(folder / "CLOSED_CACHE_INTENT.json", {"argv": argv, "source": source,
        "no_retry_or_old_cap_reset": True, "old_controller_paused_and_material_passed": True})
    child = None; streams = {}; failures = []; threads = []
    identity = None; failed = None; actual = None
    try:
        child = subprocess.Popen(argv, cwd=ROOT, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
        streams, failures, threads = stream_readers(child)
        identity = ops.proc(child.pid)
        need(identity["argv"] == argv and identity["parent"] == os.getpid(), "OwnedCacheChildOnly")
        while True:
            guarded(); s.paused(binding, number, ops)
            try:
                actual = child.wait(timeout=.25); break
            except subprocess.TimeoutExpired:
                pass
    except BaseException as exc:
        failed = type(exc).__name__
    finally:
        if child is not None and actual is None:
            # This Popen is the fixed cache utility ONLY. Never signal a group,
            # old model/reader, the controller's R descendants or the successor.
            if child.poll() is None:
                child.terminate()
            try:
                actual = child.wait(timeout=3)
            except subprocess.TimeoutExpired:
                child.kill(); actual = child.wait()
        for thread in threads:
            thread.join()
        m.write_new(folder / "CLOSED_CACHE_NATIVE_EXIT.json", {"actual_pid": child.pid if child else None,
            "start_ticks": identity["birth"] if identity else None, "argv": argv, "source": source,
            "actual_native_exit_code": actual, "observer_failure_class_only": failed or (failures[0] if failures else None),
            "streams": streams, "raw_streams_persisted": False, "automatic_retry": False})
    need(type(actual) is int and actual == 0 and not failed and not failures, "CacheNativeFailureNoRetirementOrRetry")
    m.bound(source)


def transition(h, number, folder, intent, binding, m, s, ops, guarded):
    """All checks precede retirement; exceptions before commit are resumable."""
    guarded(); s.pause(binding, ops, lambda: s.check_boundary(number))
    guarded(); s.paused(binding, number, ops)
    material = m.validate_handoff(h, m.load_module(h["legacy_driver"], "handoff_legacy"), pre_retirement=True)
    need(material["execution_authorized"] is False, "MaterialsAreNotExecutionPermission")
    admission = m.load_module(h["admission"], "handoff_admission")
    resource = m.resource_adapter(admission, h, material["validated_keys"])()
    resource_pin = m.write_new(folder / "RESOURCES_BEFORE_RETIREMENT.json", resource)
    if resource["resources_admitted"] is not True and intent.get("prepare_closed_cache") is True:
        prepare_cache_once(number, folder, intent, binding, m, s, ops, guarded)
        resource = m.resource_adapter(admission, h, material["validated_keys"])()
        resource_pin = m.write_new(folder / "RESOURCES_AFTER_CLOSED_CACHE.json", resource)
    need(resource["resources_admitted"] is True, "ResourcesUnavailableResumeOriginalNoTERM")
    guarded(); s.paused(binding, number, ops)
    m.write_new(folder / "RETIRE_INTENT.json", {"scope": "VERIFIED_CHILDLESS_SCHEDULER_ONLY_NOT_FIT",
        "expected_round": number, "prepared_material": intent["fresh_prepared"], "resources": resource_pin,
        "controller": binding, "actual_controller_success_not_asserted": True})
    terminal = retire(binding, number, s, ops, note=lambda event: m.write_new(folder / (event + ".json"),
        {"event": event, "controller": binding, "utc": m.utc(), "coordinator_source": intent["coordinator_source"]}))
    terminal.update(schema="jgw03-legacy-controller-handoff-terminal-v1",
        status="ACTUAL_PROCESS_EXIT_OBSERVED_FOR_PLANNED_HANDOFF", driver=h["legacy_driver"],
        started=h["legacy_controller"]["started"], completed_rounds=h["completed_rounds"],
        coordinator_source=intent["coordinator_source"], resource_admission=resource_pin, utc=m.utc())
    terminal_pin = m.write_new(folder / "LEGACY_NATIVE_TERMINAL.json", terminal)
    # Permit a brief orphan/zombie reaping delay; never infer an unavailable rc.
    deadline = ops.now()+10
    while True:
        try:
            m.legacy_gone(); break
        except RuntimeError:
            need(old_gone(s, ops) and ops.now() < deadline, "RetiredIdentityNotYetReapedNoLaunch")
            ops.sleep(.05)
    handed = copy.deepcopy(h)
    handed["status"] = "REVIEWED_EXPLICIT_PLANNED_HANDOFF_NOT_OLD_CONTROLLER_SUCCESS"
    handed["legacy_controller"]["terminal"] = terminal_pin
    m.validate_handoff(handed, m.load_module(h["legacy_driver"], "handoff_post_terminal_legacy"))
    handoff_pin = m.write_new(folder / "HANDOFF.json", handed)
    m.write_new(folder / "DISARM.json", {"reason": "ORIGINAL_PIDFD_EXIT_OBSERVED", "legacy_terminal": terminal_pin})
    return handed, handoff_pin


def lease_loop(intent, folder, m, s, ops):
    """Independent CONT-only recovery, with no scientific child or long limit."""
    identity = ops.proc(os.getpid()); s.owned_group(identity, identity)
    m.write_new(folder / "LEASE_READY.json", dict(identity, deadline_monotonic=intent["deadline_monotonic"]))
    result = {}; code = 1
    try:
        while True:
            source_pins(intent, m)
            if (folder / "DISARM.json").exists():
                disarm = s.published(folder / "DISARM.json")
                terminal = m.read(m.bound(disarm["legacy_terminal"])["path"])
                need(disarm["reason"] == "ORIGINAL_PIDFD_EXIT_OBSERVED" and
                     terminal["pid"] == s.CONTROLLER and terminal["start_ticks"] == s.BIRTH and
                     terminal["pidfd_exit_observed"] is True and old_gone(s, ops), "OnlyActuallyRetiredIdentityDisarmsLease")
                result = {"status": "DISARMED_OLD_IDENTITY_EXITED_NO_FURTHER_SIGNALS"}; code = 0; break
            try:
                caller = ops.proc(intent["caller"]["pid"])
                alive = caller["birth"] == intent["caller"]["birth"] and caller["state"] not in ("Z", "X")
            except (FileNotFoundError, ProcessLookupError):
                alive = False
            need(alive and ops.now() < intent["deadline_monotonic"] and
                 not (folder / "ABANDON.json").exists(), "FiniteLeaseCallerDeadlineOrAbandon")
            ops.sleep(.1)
    except BaseException as exc:
        result = {"status": "LEASE_RESTORATION_PATH", "exception_class_only": type(exc).__name__}
    finally:
        if (folder / "PAUSE_OWNED.json").exists() and code != 0:
            try:
                result["resume"] = safe_resume(intent["controller"], s, ops)
            except BaseException as exc:
                result["resume_failure_class_only"] = type(exc).__name__
        result.update(actual_lease_function_exit=code, utc=m.utc(), no_scientific_process_signals=True)
        m.write_new(folder / "LEASE_RESULT.json", result)
    return code


def observe_successor(child, argv, handed, handoff_pin, folder, successor_attempt, intent, m, s, ops):
    """Direct wait parent, bounded stream memory, deliberately NO wall timeout."""
    streams, failures, threads = stream_readers(child)
    identity = None; coordinator = None; started_pin = None
    try:
        identity = ops.proc(child.pid)
        need(identity["argv"] == argv and identity["parent"] == os.getpid(), "ActualSuccessorDirectChildCommand")
        own = ops.proc(os.getpid())
        coordinator = {"pid": own["pid"], "parent": own["parent"], "start_ticks": own["birth"], "birth": own["birth"]}
        started = {"schema": "jgw03-successor-parent-started-v1", "actual_pid": child.pid,
            "start_ticks": identity["birth"], "argv": argv, "driver": handed["actual_driver"],
            "handoff": handoff_pin, "coordinator_source": intent["coordinator_source"],
            "coordinator_identity": coordinator, "actual_native_exit_not_yet_observed": True, "utc": m.utc()}
        started_pin = m.write_new(folder / "SUCCESSOR_STARTED.json", started)
    except BaseException as exc:
        failures.append(type(exc).__name__)
    # Even failure to observe/write STARTED cannot abandon the already-created
    # child or turn a later exit into a guessed code. Its real parent still waits.
    actual = None
    try:
        actual = child.wait()  # May take days. Does NOT retain the STOP lease.
    except BaseException as exc:
        failures.append(type(exc).__name__)
    finally:
        for thread in threads:
            thread.join()
    if actual is None:
        try:
            actual = child.wait()
        except BaseException as exc:
            failures.append(type(exc).__name__)
    dispatch = ROOT / "provenance" / successor_attempt / "DISPATCH_EXIT.json"
    dispatch_pin = None
    try:
        if dispatch.exists():
            dispatch_pin = m.pin(dispatch)
    except BaseException as exc:
        failures.append(type(exc).__name__)
    result = {"schema": "jgw03-successor-parent-native-exit-v1", "actual_pid": child.pid,
        "start_ticks": identity["birth"] if identity else None, "argv": argv, "driver": handed["actual_driver"],
        "handoff": handoff_pin, "coordinator_source": intent["coordinator_source"], "started": started_pin,
        "coordinator_identity": coordinator,
        "actual_native_exit_code": actual, "dispatch_exit": dispatch_pin,
        "stdout_bytes": streams["stdout"]["bytes"], "stdout_sha256": streams["stdout"]["sha256"],
        "stderr_bytes": streams["stderr"]["bytes"], "stderr_sha256": streams["stderr"]["sha256"],
        "raw_streams_persisted": False, "observer_failure_class_only": failures[0] if failures else None,
        "old_controller_success_asserted": False, "scientific_acceptance": False, "utc": m.utc()}
    m.write_new(folder / "SUCCESSOR_NATIVE_EXIT.json", result)
    return 0 if type(actual) is int and actual == 0 and not failures and dispatch_pin else 1


def coordinate(args):
    need(sys.platform.startswith("linux") and socket.gethostname() == HOST, "OnlyDeclaredNewLinuxInstance")
    need(type(args.max_stop_seconds) is int and 30 <= args.max_stop_seconds <= 600 and
         re.fullmatch("[A-Za-z0-9_-]{1,80}", args.attempt) and
         re.fullmatch("[A-Za-z0-9_-]{1,80}", args.successor_attempt), "ExplicitFiniteFreshAttempts")
    raw = Path(args.prepared).read_bytes()
    need(hashlib.sha256(raw).hexdigest() == args.prepared_sha256, "FixedPreparedTemplateSHA")
    template = json.loads(raw)
    need(template["actual_driver"]["sha256"] == SUCCESSOR_SHA, "FixedReviewedSuccessorBeforeImport")
    m = load(template["actual_driver"]["path"], template["actual_driver"]["sha256"], "reviewed_successor")
    s = load(SAFETY, SAFETY_SHA, "reviewed_pidfd_safety"); ops = Linux(s)
    m.validate_handoff(template, m.load_module(template["legacy_driver"], "prepared_legacy"), pre_retirement=True)
    need(Path(template["legacy_attempt"]) == s.BASE, "OneOriginalControllerAttempt")
    folder = ROOT / "provenance" / args.attempt
    need(not folder.exists() and not (ROOT / "provenance" / args.successor_attempt).exists(), "NoCoordOrScientificAttemptReplay")
    import fcntl
    # A short operational mutex prevents two coordinators from claiming the
    # same pause. Scientific ownership still uses the unchanged FULL_REFERENCE.lock.
    with (ROOT / "provenance/INCREMENTAL_HANDOFF_COORDINATOR.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        binding = wait_boundary(args.expected_round, args.max_wait_seconds, s, ops)
        h = build_material(template, args.expected_round, m)
        m.validate_handoff(h, m.load_module(h["legacy_driver"], "boundary_legacy"), pre_retirement=True)
        folder.mkdir()
        material_pin = m.write_new(folder / "PREPARED_MATERIALS.json", h)
        intent = {"schema": "jgw03-incremental-natural-boundary-coordinator-v1", "expected_round": args.expected_round,
            "controller": binding, "caller": ops.proc(os.getpid()),
            "coordinator_source": m.pin(Path(__file__).resolve()), "successor_driver": h["actual_driver"],
            "safety_source": m.pin(SAFETY), "prepared_template": m.pin(Path(args.prepared)),
            "fresh_prepared": material_pin, "deadline_monotonic": ops.now()+args.max_stop_seconds-10,
            "prepare_closed_cache": args.prepare_closed_cache,
            "successor_attempt": args.successor_attempt, "utc": m.utc()}
        if args.prepare_closed_cache:
            intent["closed_cache_source"] = m.pin(CACHE)
            need(intent["closed_cache_source"]["sha256"] == CACHE_SHA, "FixedOptionalCacheSourceBeforePause")
        intent_pin = m.write_new(folder / "INTENT.json", intent)
        guard = None; handed = handoff_pin = None; phase = {"observing": False}; observed_signals = []
        def interrupted(sig, frame):
            if phase["observing"]:
                observed_signals.append(int(sig)); return
            raise InterruptedError("BoundaryCoordinatorSignal")
        handlers = {sig: signal.signal(sig, interrupted) for sig in (signal.SIGTERM, signal.SIGINT)}
        try:
            guard = subprocess.Popen([sys.executable, "-I", "-B", intent["coordinator_source"]["path"], "_lease",
                str(folder), intent_pin["sha256"]], start_new_session=True, stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
            end = min(ops.now()+15, intent["deadline_monotonic"])
            while not (folder / "LEASE_READY.json").exists():
                need(ops.now() < end and guard.poll() is None, "IndependentLeaseREADYBeforeAnySTOP")
                ops.sleep(.05)
            need(s.published(folder / "LEASE_READY.json")["pid"] == guard.pid, "ReadyFromActualIndependentLeaseChild")
            finite_ok(intent, folder, m, s, ops)
            m.write_new(folder / "PAUSE_OWNED.json", {"controller": binding,
                "independent_lease": s.published(folder / "LEASE_READY.json"), "restoration_armed_before_STOP": True})
            handed, handoff_pin = transition(h, args.expected_round, folder, intent, binding, m, s, ops,
                lambda: finite_ok(intent, folder, m, s, ops))
            # The old process has actually exited. The finite pause lease ends
            # before launching any long-running successor scientific worker.
            actual_lease = guard.wait(timeout=5)
            ready = s.published(folder / "LEASE_READY.json")
            m.write_new(folder / "LEASE_NATIVE_EXIT.json", {"actual_pid": guard.pid, "start_ticks": ready["birth"],
                "actual_native_exit_code": actual_lease, "coordinator_source": intent["coordinator_source"]})
            need(type(actual_lease) is int and actual_lease == 0, "RetiredIdentityLeaseDisarmExit")
        except BaseException as exc:
            if not (folder / "DISARM.json").exists():
                m.write_new(folder / "ABANDON.json", {"exception_class_only": type(exc).__name__})
            restore = safe_resume(binding, s, ops)
            m.write_new(folder / "COORDINATOR_STOP.json", {"exception_class_only": type(exc).__name__,
                "resume": restore, "retired_handoff_may_need_inspection": (folder / "RETIRE_INTENT.json").exists(),
                "no_automatic_retry": True, "raw_error_persisted": False})
            for sig, handler in handlers.items():
                signal.signal(sig, handler)
            return 1
        finally:
            if handed is None:
                safe_resume(binding, s, ops)
        argv = [sys.executable, "-I", "-B", handed["actual_driver"]["path"], "--execute",
                "--handoff", handoff_pin["path"], "--handoff-sha256", handoff_pin["sha256"], "--attempt", args.successor_attempt]
        m.bound(handed["actual_driver"]); m.bound(handoff_pin)
        # Direct parent remains solely to observe actual native exit. No raw
        # stream accumulation, resource loop, re-launch or model signal occurs.
        phase["observing"] = True  # No signal gap after a native child is created.
        child = subprocess.Popen(argv, cwd=ROOT, start_new_session=True, stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
        try:
            return observe_successor(child, argv, handed, handoff_pin, folder, args.successor_attempt, intent, m, s, ops)
        finally:
            for sig, handler in handlers.items():
                signal.signal(sig, handler)


def lease(folder, expected_sha):
    raw = (folder / "INTENT.json").read_bytes()
    need(hashlib.sha256(raw).hexdigest() == expected_sha, "LeaseExactIntentSHA")
    intent = json.loads(raw)
    need(intent["successor_driver"]["sha256"] == SUCCESSOR_SHA, "LeaseFixedReviewedSuccessor")
    m = load(intent["successor_driver"]["path"], intent["successor_driver"]["sha256"], "lease_successor_tools")
    s = load(SAFETY, SAFETY_SHA, "lease_pidfd_safety"); ops = Linux(s)
    need(sys.platform.startswith("linux") and socket.gethostname() == HOST and
         intent["schema"] == "jgw03-incremental-natural-boundary-coordinator-v1" and
         Path(intent["coordinator_source"]["path"]) == Path(__file__).resolve() and
         0 < intent["deadline_monotonic"]-ops.now() <= 600 and folder.parent == ROOT / "provenance", "FiniteSameSourceIndependentLease")
    return lease_loop(intent, folder, m, s, ops)


def main():
    if len(sys.argv) == 4 and sys.argv[1] == "_lease":
        return lease(Path(sys.argv[2]), sys.argv[3])
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--prepared", type=Path)
    parser.add_argument("--prepared-sha256")
    parser.add_argument("--expected-round", type=int)
    parser.add_argument("--attempt")
    parser.add_argument("--successor-attempt")
    parser.add_argument("--max-wait-seconds", type=int, default=1800)
    parser.add_argument("--max-stop-seconds", type=int, default=120)
    parser.add_argument("--prepare-closed-cache", action="store_true")
    args = parser.parse_args()
    if not args.execute:
        print(json.dumps({"status": "PLAN_ONLY_NO_PROCESS_OR_SCIENTIFIC_ACTION", "finite_pause_max_seconds": 600,
                          "post_handoff_observer_wall_timeout": None, "no_cache_or_model_signals": True}))
        return 0
    need(all(getattr(args, name) is not None for name in
             ("prepared", "prepared_sha256", "expected_round", "attempt", "successor_attempt")), "ExplicitReviewedMaterialAndFreshAttempts")
    return coordinate(args)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(json.dumps({"status": "COORDINATOR_STOP_NO_AUTOMATIC_RETRY", "exception_class_only": type(exc).__name__,
                          "raw_error_persisted": False}))
        raise SystemExit(2)
