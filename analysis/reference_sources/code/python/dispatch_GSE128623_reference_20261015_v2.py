"""Bounded dispatcher for one admitted run; never fits, repairs or deletes itself.

CLI: RUN_ROOT CONFIG_JSON CONFIG_SHA256 [--plan-only]
Each R child consumes a disjoint range of the frozen, full allocation matrix.
"""
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_new(path, value):
    with Path(path).open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2)
        handle.write("\n")


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def main():
    assert len(sys.argv) in (4, 5)
    root = Path(sys.argv[1]).resolve(strict=True)
    assert str(root) == "/root/autodl-tmp/gse128623_label_reference_20261015r01"
    config_path = Path(sys.argv[2]).resolve(strict=True)
    assert root in config_path.parents and sha(config_path) == sys.argv[3]
    config = json.loads(config_path.read_text())
    assert config["run_id"] == root.name and config["chunk_size"] == 50
    maximum = config["maximum_workers"]
    assert isinstance(maximum, int) and 1 <= maximum <= 12
    assert config["minimum_available_disk_bytes"] >= 50 * 1024**3
    jobs_path = root / "JOBS_PREPARED.json"
    assert sha(jobs_path) == config["jobs_prepared_sha256"]
    receipt = json.loads(jobs_path.read_text())
    assert receipt["status"] == "JOBS_PREPARED_NOT_LAUNCHED"
    admission = receipt["admission"]
    assert sha(admission["path"]) == admission["sha256"]
    admitted = json.loads(Path(admission["path"]).read_text())
    assert admitted["status"] == "ALL_FOUR_OBSERVED_AND_RESOURCE_PASS"
    assert admitted["resource_acceptance"] is True and admitted["run_id"] == root.name
    resource_evidence = admitted["resource_evidence"]
    assert sha(resource_evidence["path"]) == resource_evidence["sha256"]
    resource_decision = json.loads(Path(resource_evidence["path"]).read_text())
    assert maximum <= resource_decision["maximum_workers_accepted"]
    assert len(receipt["jobs"]) == 24
    for relative, expected in config["source_pins"].items():
        assert sha(root / "source" / relative) == expected
    runtime = "/root/autodl-tmp/route_b2_public_analysis/runtime/bin/Rscript"
    # This separately pinned entry normalizes JSON metadata and binds the
    # lossless storage adapter. Numerical kernels/canonical job pins are intact.
    entry = root / "source/run_GSE128623_reference_shard_transport_v2_20261015.R"
    assert entry.name in config["source_pins"]
    tasks = []
    counts = {}
    for key, job in receipt["jobs"].items():
        cell, space = key.split("/")
        path = Path(job["path"]).resolve(strict=True)
        assert path.parent == root / "jobs" and sha(path) == job["sha256"]
        total = 10001 if space == "six_group" else (
            126 if cell.endswith("predefined29") and space == "D4_vs_control_exact" else 252)
        counts[key] = total
        for first in range(1, total + 1, 50):
            last = min(first + 49, total)
            worker = f"{config['dispatch_id']}_{cell}_{space}_{first:05d}_{last:05d}"
            tasks.append(dict(key=key, first=first, last=last, worker=worker,
                              job=str(path), job_sha256=job["sha256"]))
    # First full-space range of every cell, then exact companions, then the
    # remaining full-space ranges round-robin; never select by result values.
    tasks.sort(key=lambda task: (0 if task["first"] == 1 and task["key"].endswith("/six_group")
                                else 1 if not task["key"].endswith("/six_group") else 2,
                                task["first"], task["key"]))
    cells = ("raw_all30", "raw_predefined29", "adapter_all30", "adapter_predefined29")
    spaces = ("six_group", "D4_vs_control_exact", "D7_vs_control_exact",
              "D14_vs_control_exact", "D28_vs_control_exact", "D42_vs_control_exact")
    assert set(counts) == {cell + "/" + space for cell in cells for space in spaces}
    assert sum(counts.values()) == 44792 and len(tasks) == 918
    assert len({task["worker"] for task in tasks}) == len(tasks)
    plan_only = len(sys.argv) == 5
    if plan_only:
        assert sys.argv[4] == "--plan-only"
    output = root / "dispatch" / (config["dispatch_id"] + ("_plan" if plan_only else ""))
    assert output.parent == root / "dispatch" and output.name.replace("_", "").isalnum()
    output.mkdir(parents=True, exist_ok=False)
    write_new(output / "PLAN.json", dict(created_utc=utc(), tasks=tasks,
        required_archives=sum(counts.values()), counts=counts, config_sha256=sys.argv[3],
        controller_sha256=sha(__file__), maximum_workers=maximum))
    if plan_only:
        print("PLAN_ONLY_NO_WORKER_LAUNCHED", len(tasks), flush=True)
        return
    # A preceding dispatcher may still be draining its eight valid day ranges.
    # Count those children in the same 12-worker budget and never overlap IDs.
    prior = root / "dispatch" / config["prior_dispatch_id"]
    assert config["prior_dispatch_id"] == "main_01" and prior.is_dir()
    assert sha(root / "DISPATCH_CONFIG.json") == config["prior_dispatch_config_sha256"]

    def prior_active_ranges():
        ranges = []
        for started_path in prior.glob("*/STARTED.json"):
            if not (started_path.parent / "PROCESS_EXIT.json").exists():
                ranges.append(json.loads(started_path.read_text())["task"])
        return ranges

    def overlaps(task, previous):
        return task["key"] == previous["key"] and task["first"] <= previous["last"] and previous["first"] <= task["last"]
    env = os.environ.copy()
    env.update(R_LIBS_USER="/root/autodl-tmp/GSE128623_R_overlay_20260920/library:/root/autodl-tmp/route_b2_public_analysis/R-library",
        R_PROFILE="/dev/null", R_PROFILE_USER="/dev/null", R_ENVIRON="/dev/null", R_ENVIRON_USER="/dev/null",
        RENV_CONFIG_AUTOLOADER_ENABLED="FALSE", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1", LC_ALL="C.UTF-8", LANG="C.UTF-8")
    active = {}
    pending = list(tasks)
    next_task = 0
    failure = None
    finished = 0
    write_new(output / "STARTED.json", dict(started_utc=utc(), pid=os.getpid()))
    while active or (pending and failure is None):
        previous_active = prior_active_ranges()
        while failure is None and pending and len(active) + len(previous_active) < maximum:
            if shutil.disk_usage(root).free < config["minimum_available_disk_bytes"]:
                failure = "RESOURCE_DISK_FLOOR_NO_NEW_WORKERS"
                break
            eligible = next((index for index, candidate in enumerate(pending)
                             if not any(overlaps(candidate, previous) for previous in previous_active)), None)
            if eligible is None:
                break
            task = pending.pop(eligible)
            destination = output / task["worker"]
            args = [runtime, "--vanilla", str(entry), str(root / "source"), task["job"],
                    task["job_sha256"], str(task["first"]), str(task["last"]), task["worker"]]
            stdout = stderr = None
            process = None
            try:
                destination.mkdir(exist_ok=False)
                stdout = (destination / "stdout").open("xb")
                stderr = (destination / "stderr").open("xb")
                process = subprocess.Popen(args, env=env, stdin=subprocess.DEVNULL,
                                           stdout=stdout, stderr=stderr, close_fds=True)
                started = time.monotonic()
                active[process.pid] = (process, destination, started, stdout, stderr)
                next_task += 1
                write_new(destination / "STARTED.json", dict(task=task, command=args,
                    pid=process.pid, started_utc=utc()))
            except Exception as error:
                # Any already-spawned process remains registered for wait4.
                # Never retry a scientific allocation after an I/O failure.
                failure = "LAUNCH_OR_METADATA_EXCEPTION_NO_NEW_WORKERS"
                if process is None:
                    for handle in (stdout, stderr):
                        if handle is not None:
                            handle.close()
                print(type(error).__name__, failure, flush=True)
                break
        for pid, (process, destination, started, stdout, stderr) in list(active.items()):
            waited, status, usage = os.wait4(pid, os.WNOHANG)
            if waited == 0:
                continue
            code = os.waitstatus_to_exitcode(status)
            process.returncode = code
            stdout.close()
            stderr.close()
            finished += 1
            del active[pid]
            try:
                write_new(destination / "PROCESS_EXIT.json", dict(exit_code=code, finished_utc=utc(),
                    elapsed_seconds=time.monotonic() - started, user_seconds=usage.ru_utime,
                    system_seconds=usage.ru_stime, maximum_resident_set_kib_linux=usage.ru_maxrss))
            except Exception as error:
                failure = "TERMINAL_METADATA_EXCEPTION_NO_NEW_WORKERS"
                print(type(error).__name__, failure, "pid", pid, "exit", code, flush=True)
            if code != 0:
                failure = "NONZERO_WORKER_EXIT_NO_NEW_WORKERS"
        if active or (pending and failure is None):
            time.sleep(5)
    write_new(output / "DISPATCH_EXIT.json", dict(status=failure or "ALL_REQUESTED_SHARDS_EXITED_ZERO_NEEDS_INDEPENDENT_CHECK",
        finished_utc=utc(), started_shards=next_task, finished_shards=finished,
        planned_shards=len(tasks), whole_run_independently_accepted=False))
    if failure:
        sys.exit(1)


if __name__ == "__main__":
    main()
