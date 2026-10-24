"""One fixed downstream R invocation and its genuine terminal; no retries."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path('/root/autodl-tmp/GSE205906_matched_STAR_Ensembl101_20261015r01')
PROJECT = ROOT/'custody_source_05'
CALLER = PROJECT/'code/R/run_gse205906_matched_downstream_custody_v1_20261022.R'
CALLER_SHA = '7cdafad1f838cd01ca8ae730762c059a78b2b3f515d47785fd09fe87f9e2275e'
R = Path('/root/autodl-tmp/route_b2_public_analysis/runtime/bin/Rscript')
R_SHA = '5982fd0f4e6a94b7eec43f25cd7b83d389d37a62be750110b51edff6b290bc8d'
PRODUCER = ROOT/'source_staging_20261015_01/run_gse205906_matched_STAR_20261015.py'
OUTPUT = ROOT/'downstream_custody_01'
PARENT = ROOT/'provenance/matched_downstream_custody_parent_20261022_r01'


def pin(path):
    path = Path(path)
    return dict(path=str(path), bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def put(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True)


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def execute():
    assert pin(CALLER)['sha256'] == CALLER_SHA and pin(R)['sha256'] == R_SHA
    assert pin(ROOT/'matrix_01/COMPLETE.json')['sha256'] == '1067c83043a3c7e06dc075774251318f336124477d789828cbf35ffeb950f074'
    assert not OUTPUT.exists() and not PARENT.exists()
    PARENT.mkdir()
    temporary = ROOT/'matched_downstream_custody_tmp_01'
    temporary.mkdir()
    assert 192 in os.sched_getaffinity(0)
    os.sched_setaffinity(0, {192})  # Separate from the approved fifteen-worker IDs193–207.
    argv = [str(R), '--vanilla', str(CALLER), str(PROJECT), str(PRODUCER), str(OUTPUT),
            '--execute', '--python', '/usr/bin/python3', '--quiet-window-confirmed',
            '--reserved-growth-bytes', str(240*1024**3), '--disk-budget-bytes', str(16*1024**3)]
    env = dict(os.environ, R_LIBS_USER='/root/autodl-tmp/GSE128623_R_overlay_20260920/library:/root/autodl-tmp/route_b2_public_analysis/R-library',
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               R_DATATABLE_NUM_THREADS='1', LC_ALL='C', LANG='C', TMPDIR=str(temporary))
    child = subprocess.Popen(argv, cwd=ROOT, env=env, stdin=subprocess.DEVNULL,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, start_new_session=True)
    started = PARENT/'STARTED.json'
    ticks = None
    observer_error = None
    try:
        stat = Path('/proc', str(child.pid), 'stat').read_text().rsplit(') ', 1)[1].split()
        ticks = int(stat[19])
        put(started, dict(schema='GSE205906_MATCHED_ACTUAL_R_STARTED/v1', actual_R_pid=child.pid,
            actual_parent_pid=os.getpid(), start_ticks=ticks, argv=argv,
            assigned_CPU=192,
            source=pin(CALLER), observer_source=pin(__file__), utc=now(),
            resource_window='Known bounded Shaoyao co-host; not a claim that all processes are idle'))
    except Exception as error:
        observer_error = type(error).__name__
    stdout, stderr = child.communicate()  # Always drain and wait even after a metadata failure.
    complete = OUTPUT/'COMPLETE.json'
    put(PARENT/'PARENT_NATIVE_EXIT.json', dict(schema='GSE205906_MATCHED_PARENT_PROCESS_EXIT/v1',
        run_root=str(OUTPUT), actual_R_exit=child.returncode, actual_R_pid=child.pid,
        actual_parent_pid=os.getpid(), start_ticks=ticks, started=pin(started) if started.exists() else None,
        caller_source_sha256=CALLER_SHA, observer_source=pin(__file__), argv=argv,
        complete_pin=pin(complete) if complete.exists() else None,
        stdout_bytes=len(stdout), stdout_sha256=hashlib.sha256(stdout).hexdigest(),
        stderr_bytes=len(stderr), stderr_sha256=hashlib.sha256(stderr).hexdigest(),
        raw_streams_persisted=False, automatic_refit=False,
        observer_failure_class_only=observer_error, utc=now()))
    return child.returncode


if __name__ == '__main__':
    if sys.argv[1:] == ['--execute']:
        raise SystemExit(execute())
    assert not sys.argv[1:]
    print('PLAN_ONLY_ONE_FIXED_R_CALL_NO_RETRY_NO_IO')




