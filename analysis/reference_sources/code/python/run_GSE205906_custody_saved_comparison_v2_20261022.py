"""One-use saved comparison/component checks; no scientific fitting or retries."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path('/root/autodl-tmp/GSE205906_matched_STAR_Ensembl101_20261015r01')
PROJECT = ROOT/'custody_source_05'
R = Path('/root/autodl-tmp/route_b2_public_analysis/runtime/bin/Rscript')
NEW_PARENT = ROOT/'provenance/matched_downstream_custody_parent_20261022_r01/PARENT_NATIVE_EXIT.json'
OLD_PARENT = ROOT/'downstream_source_20261022_r01/HISAT_PREPARATION_ACTUAL_PARENT_EXIT.json'
COMPARE = PROJECT/'code/R/run_gse205906_matched_saved_comparison_custody_v2_20261022.R'
AUDIT = PROJECT/'code/R/audit_GSE205906_saved_custody_components_20261022.R'


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
    assert pin(R)['sha256'] == '5982fd0f4e6a94b7eec43f25cd7b83d389d37a62be750110b51edff6b290bc8d'
    assert pin(COMPARE)['sha256'] == 'e457bf417f18c31a2a5d9ed825fcb2e975711ad06a9fe65534ba98219590f16e'
    assert pin(AUDIT)['sha256'] == '4bd41749b8a17b167c6aa9b419d43eda733ddc56a46acc595332509187305a3c'
    terminal = json.loads(NEW_PARENT.read_text())
    assert type(terminal['actual_R_exit']) is int and terminal['actual_R_exit'] == 0
    assert terminal['observer_failure_class_only'] is None
    assert terminal['run_root'] == str(ROOT/'downstream_custody_01')
    assert terminal['caller_source_sha256'] == '7cdafad1f838cd01ca8ae730762c059a78b2b3f515d47785fd09fe87f9e2275e'
    assert terminal['complete_pin'] == pin(ROOT/'downstream_custody_01/COMPLETE.json')
    assert pin(OLD_PARENT)['sha256'] == '713e404db55975c767920e5e05020cc0be3cbc9595e9752c72cccf00f83b570d'
    old_pin, new_pin = pin(OLD_PARENT), pin(NEW_PARENT)
    jobs = [
        ('comparison', COMPARE, ROOT/'saved_custody_comparison_02',
         ['--compare', '--hisat-prepared', str(ROOT/'saved_hisat_inputs_01'),
          '--hisat-exit-receipt', str(OLD_PARENT), '--hisat-exit-sha', old_pin['sha256'],
          '--star-downstream', str(ROOT/'downstream_custody_01'),
          '--star-exit-receipt', str(NEW_PARENT), '--star-exit-sha', new_pin['sha256']]),
        ('components', AUDIT, ROOT/'saved_custody_components_01', [])]
    jobs = jobs[:1]  # Only the corrected saved comparison; no component replay.
    assert all(not output.exists() and not (ROOT/f'provenance/custody_saved_{name}_parent_02').exists()
               for name, _, output, _ in jobs)
    env = dict(os.environ, R_LIBS_USER='/root/autodl-tmp/GSE128623_R_overlay_20260920/library:/root/autodl-tmp/route_b2_public_analysis/R-library',
        OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
        R_DATATABLE_NUM_THREADS='1', LC_ALL='C', LANG='C',
        TMPDIR=str(ROOT/'matched_downstream_custody_tmp_01'))
    assert 192 in os.sched_getaffinity(0)
    os.sched_setaffinity(0, {192})
    for name, source, output, extra in jobs:
        folder = ROOT/f'provenance/custody_saved_{name}_parent_02'
        folder.mkdir()
        args = ([str(PROJECT), str(output), *extra] if name == 'comparison' else
                [str(PROJECT), str(ROOT/'numeric_readback_source_20261022_r01/code/R'),
                 new_pin['sha256'], str(output)])
        argv = [str(R), '--vanilla', str(source), *args]
        child = subprocess.Popen(argv, cwd=ROOT, env=env, stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        observer_error = None
        try:
            stat = Path('/proc', str(child.pid), 'stat').read_text().rsplit(') ', 1)[1].split()
            started = dict(actual_R_pid=child.pid, birth=int(stat[19]), actual_parent_pid=os.getpid(),
                argv=argv, source=pin(source), parent_source=pin(__file__), utc=now())
            put(folder/'STARTED.json', started)
        except Exception as condition:
            observer_error = type(condition).__name__
        stdout, stderr = child.communicate()
        put(folder/'NATIVE_EXIT.json', dict(actual_R_exit=child.returncode,
            started=pin(folder/'STARTED.json') if (folder/'STARTED.json').exists() else None,
            observer_failure_class_only=observer_error,
            old_parent=old_pin, new_parent=new_pin, source=pin(source), parent_source=pin(__file__),
            stdout_bytes=len(stdout), stdout_sha256=hashlib.sha256(stdout).hexdigest(),
            stderr_bytes=len(stderr), stderr_sha256=hashlib.sha256(stderr).hexdigest(),
            raw_streams_persisted=False, automatic_retry=False, scientific_fitter_invocations=0,
            output=str(output), utc=now()))
        if child.returncode != 0 or observer_error is not None:
            return child.returncode if child.returncode != 0 else 1
    return 0


if __name__ == '__main__':
    if sys.argv[1:] == ['--execute']:
        raise SystemExit(execute())
    assert not sys.argv[1:]
    print('PLAN_ONLY_SAVED_OBJECTS_TWO_FINITE_CHECKS_NO_FIT_NO_IO')


