"""Explicit one-command Linux preparation/bridge supervisor, never full dispatch.

Raw stdout/stderr remain in RAM. Preserve existing workloads, original G128
growth budget and floor; run at most one one-core bridge, no retry or kill.
"""
import argparse
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import stat
import subprocess
import time

ROOT = Path("/root/autodl-tmp/JGW03_full_reference_linux_20261017r01")
G128 = Path("/root/autodl-tmp/gse128623_label_reference_20261015r01")
GIB = 1024 ** 3
BRIDGE_MEMORY = 3 * GIB
BRIDGE_DISK = 8 * GIB
OTHER_MEMORY_HEADROOM = 8 * GIB
RUNTIME = Path("/root/autodl-tmp/route_b2_public_analysis/runtime/bin/Rscript")


def need(ok, message):
    if not ok:
        raise ValueError(message)


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def pin(path):
    p = Path(path)
    return {"bytes": p.stat().st_size, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}


def read(path):
    return json.loads(Path(path).read_text())


def write_new(path, obj):
    with Path(path).open("x", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=2)
        handle.write("\n")


def proc(pid):
    try:
        p = Path("/proc") / str(pid)
        text = (p / "stat").read_text()
        fields = text[text.rfind(")") + 2:].split()
        cmd = (p / "cmdline").read_bytes().decode("utf-8", "replace")
        return {"parent": int(fields[1]), "start_ticks": int(fields[19]), "cmd": cmd}
    except (OSError, ValueError):
        return None


def admission():
    rp = G128 / "RESOURCE_DECISION.json"
    pp = G128 / "dispatch/main_02/PLAN.json"
    need(pin(rp)["sha256"] == "c27f7db86584f0c8ea2c344b675bc6849a77d531e5678eca69059366308d2ec7", "Original budget differs")
    need(pin(pp)["sha256"] == "35c5cb4d790529854756a778db46704cc6013bac5ed7eab5997ea73bb43ae647", "Original plan differs")
    d, plan = read(rp), read(pp)
    controller = proc(522800)
    need(controller and controller["start_ticks"] == 2556696464 and str(G128) in controller["cmd"], "Original controller not bound")
    own_ancestors = {os.getpid(), os.getppid()}
    for p in Path("/proc").iterdir():
        if not p.name.isdigit() or int(p.name) in own_ancestors:
            continue
        item = proc(int(p.name))
        if not item:
            continue
        cmd = item["cmd"]
        executable = cmd.split("\0")[0]
        need(not (("/R" in executable or "Rscript" in executable) and str(G128) not in cmd), "Another scientific R owns window")
        need(not ("STAR" in Path(executable).name or "atlas" in cmd.lower()), "Atlas or STAR owns window")
    memory_limit = int(Path("/sys/fs/cgroup/memory.max").read_text())
    memory_current = int(Path("/sys/fs/cgroup/memory.current").read_text())
    need(memory_limit - memory_current >= BRIDGE_MEMORY + OTHER_MEMORY_HEADROOM, "Direct memory headroom unavailable")
    need(memory_limit >= d["complete_memory_scenario_bytes"] + BRIDGE_MEMORY + OTHER_MEMORY_HEADROOM, "Full original growth allowance unavailable")
    occupied = 0
    inodes = set()
    for key, n in plan["counts"].items():
        folder = G128 / "allocations" / key
        need(folder.is_dir() and not folder.is_symlink() and not folder.parent.is_symlink(), "Archive parent differs")
        for p in folder.iterdir():
            match = re.fullmatch(r"allocation_([0-9]{5})\.rds", p.name)
            if not match:
                continue
            s = p.lstat()
            need(stat.S_ISREG(s.st_mode) and s.st_nlink == 1 and s.st_blocks * 512 >= s.st_size
                 and (s.st_dev, s.st_ino) not in inodes and 1 <= int(match[1]) <= n, "Disk accounting not admissible")
            need(s.st_mtime >= datetime.datetime.fromisoformat(d["recorded_utc"]).timestamp(), "Pre-admission output not discounted")
            inodes.add((s.st_dev, s.st_ino))
            occupied += s.st_size
    remaining = max(0, int(d["complete_additional_disk_scenario_and_floor_bytes"]) - occupied)
    v = os.statvfs(ROOT)
    free = v.f_bavail * v.f_frsize
    need(free >= remaining + BRIDGE_DISK, "Original disk growth plus bridge allowance unavailable")
    return {"utc": utc(), "memory_limit_bytes": memory_limit, "memory_current_bytes": memory_current,
            "original_G128_memory_scenario_bytes": d["complete_memory_scenario_bytes"],
            "bridge_memory_bytes": BRIDGE_MEMORY, "additional_memory_headroom_bytes": OTHER_MEMORY_HEADROOM,
            "free_disk_bytes": free, "remaining_original_disk_scenario_including_floor_bytes": remaining,
            "bridge_disk_allowance_bytes": BRIDGE_DISK, "single_CPU_only": True,
            "Atlas_STAR_process_classification_absent_not_future_idle_guarantee": True,
            "full_Linux_dispatch_admitted": False}


def run(args):
    need(ROOT.is_dir() and not ROOT.is_symlink(), "Staged independent root required")
    config_path = ROOT / "CONFIG.json"
    need(pin(config_path)["sha256"] == args.config_sha256, "Explicit config differs")
    config = read(config_path)
    entry = Path(config["project_root"]) / "code/R/run_JGW03_full_reference_linux_20261017.R"
    need(pin(entry) == config["source_pins"]["code/R/run_JGW03_full_reference_linux_20261017.R"], "Linux entry differs")
    need(re.fullmatch(r"[A-Za-z0-9_-]+", args.attempt), "Explicit attempt required")
    call = [str(RUNTIME), "--vanilla", str(entry), str(config_path), args.mode, args.accession]
    if args.mode == "bridge":
        need(args.allocation in (1, 2), "Only six pre-fixed fits admitted")
        call.append(str(args.allocation))
    folder = ROOT / "provenance" / args.attempt
    need(not folder.exists(), "Fresh process record required; do not relaunch completed fits")
    lock_path = ROOT / "provenance/FINITE_BRIDGE.lock"
    lock_path.parent.mkdir(exist_ok=True)
    with lock_path.open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        before = admission()
        folder.mkdir()
        write_new(folder / "ADMISSION.json", before)
        cpus = sorted(os.sched_getaffinity(0))
        need(cpus, "One allowed core required")
        assigned_cpu = cpus[-1]
        env = dict(os.environ, R_LIBS_USER="/root/autodl-tmp/GSE128623_R_overlay_20260920/library:/root/autodl-tmp/route_b2_public_analysis/R-library",
                   OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", LANG="C", LC_ALL="C",
                   TMPDIR=str(ROOT / "tmp"))
        (ROOT / "tmp").mkdir(exist_ok=True)
        write_new(folder / "COMMAND.json", {"argv": call, "config": pin(config_path), "entry": pin(entry),
                                            "supervisor": pin(__file__), "assigned_CPU": assigned_cpu,
                                            "mode": args.mode, "accession": args.accession, "allocation": args.allocation,
                                            "Rscript": pin(RUNTIME), "started_utc": utc()})
        start = time.monotonic()
        child = subprocess.Popen(call, env=env, cwd=ROOT,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 preexec_fn=lambda: os.sched_setaffinity(0, {assigned_cpu}))
        min_free = before["memory_limit_bytes"] - before["memory_current_bytes"]
        peak_rss = 0
        out, err = b"", b""
        streams_collected = False
        observation_failure = None
        try:
            identity = proc(child.pid)
            write_new(folder / "STARTED.json", {"actual_R_pid": child.pid, "start_ticks": identity["start_ticks"] if identity else None,
                                               "utc": utc(), "config_sha256": args.config_sha256})
            while True:
                try:
                    status = (Path("/proc") / str(child.pid) / "status").read_text()
                    values = [int(m[1]) * 1024 for m in re.finditer(r"^Vm(?:HWM|RSS):\s+(\d+) kB$", status, re.M)]
                    peak_rss = max([peak_rss] + values)
                    min_free = min(min_free, before["memory_limit_bytes"] - int(Path("/sys/fs/cgroup/memory.current").read_text()))
                except OSError:
                    pass
                try:
                    out, err = child.communicate(timeout=1)
                    streams_collected = True
                    break
                except subprocess.TimeoutExpired:
                    pass
        except BaseException as condition:
            # Keep only the class. Drain the already-owned child while holding
            # the exclusive lock even when a metadata/observer operation fails.
            observation_failure = type(condition).__name__
        finally:
            if not streams_collected:
                out, err = child.communicate()
            child.wait()
        peak_rss = max(peak_rss, resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * 1024)
        result = {"status": "FINITE_BRIDGE_ACTUAL_EXIT_NOT_INDEPENDENT_ACCEPTANCE", "actual_R_pid": child.pid,
                  "actual_R_exit_code": child.returncode, "finished_utc": utc(), "elapsed_seconds": time.monotonic() - start,
                  "peak_child_RSS_bytes": peak_rss, "minimum_direct_memory_headroom_bytes": min_free,
                  "bridge_memory_budget_met": peak_rss <= BRIDGE_MEMORY,
                  "observer_failure_class_only": observation_failure,
                  "stdout_bytes": len(out), "stdout_sha256": hashlib.sha256(out).hexdigest(),
                  "stderr_bytes": len(err), "stderr_sha256": hashlib.sha256(err).hexdigest(),
                  "raw_streams_persisted": False, "no_automatic_retry_or_full_dispatch": True}
        write_new(folder / "PROCESS_EXIT.json", result)
        print(json.dumps(result))
        return 0 if child.returncode == 0 and peak_rss <= BRIDGE_MEMORY and observation_failure is None else 2


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config-sha256", required=True)
    p.add_argument("--mode", choices=("prepare", "bridge"), required=True)
    p.add_argument("--accession", choices=("GSE289929", "GSE227494", "GSE299063"), required=True)
    p.add_argument("--allocation", type=int)
    p.add_argument("--attempt", required=True)
    a = p.parse_args()
    try:
        code = run(a)
    except Exception as condition:
        # Raw condition is not printed or copied; retained owned artifacts stay.
        print(json.dumps({"status": "STOP_RETAINED_NO_RETRY", "exception_class": type(condition).__name__,
                          "raw_error_persisted": False}))
        code = 2
    raise SystemExit(code)
