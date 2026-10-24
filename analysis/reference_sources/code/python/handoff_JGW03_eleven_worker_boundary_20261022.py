"""Plan-only by default; finite childless-eight-worker handoff, then actual-child wait.

No scientific inputs, seeds or completed fits are changed. Kernel and driver
paths are explicit SHA-pinned inputs. The reviewed old source is never edited.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import socket
import subprocess
import sys
import threading
import types

ROOT = Path('/root/autodl-tmp/JGW03_full_reference_linux_20261017r03')
HOST = 'autodl-container-62c8pzy536-50fb405a'
KERNEL_SHA = '8c0640997f1ad2a917b5f57fdafdcf6a65cb7d3574445f5dae0e6416673520d7'
SAFETY_SHA = 'fe9f804c162a26e34c1f5fd52d5c436b8d72bd81fe19a032691a741a1323de53'
OLD_SHA = '6f64ba466e7d969b1feb709533159ed757353fe3795a5aec68c1e69bd994b546'
OLD_DRIVER = ROOT / 'eight_worker_source_20261022_r01/code/python/run_JGW03_eight_worker_shard_successor_20261022.py'
OLD_HANDOFF = ROOT / 'provenance/eight_worker_handoff_0068_r01/HANDOFF.json'
OLD_HANDOFF_SHA = '6ea09f352fe9d4b8ff5c9f15ba42bca0d986c8db01b3950e746efa41673a9daa'
PARENT_EXIT = ROOT / 'provenance/eight_worker_handoff_0068_r01/SUCCESSOR_NATIVE_EXIT.json'


def need(ok, label):
    if not ok:
        raise RuntimeError(label)

def load(path, sha, name):
    raw = Path(path).read_bytes()
    need(hashlib.sha256(raw).hexdigest() == sha, 'ReviewedSourceSHA')
    value = types.ModuleType(name); value.__file__ = str(path)
    exec(compile(raw, str(path), 'exec'), value.__dict__)
    return value

def safety(path):
    """Explicit isolated operational adapter; scientific source is never patched."""
    s = load(path, SAFETY_SHA, 'eleven_worker_pidfd_safety')
    s.CONTROLLER, s.BIRTH = 11848, 2613447611
    s.BASE = ROOT / 'provenance/eight_worker_reference_01'
    s.DRIVER, s.DRIVER_SHA = OLD_DRIVER, OLD_SHA
    old = load(OLD_DRIVER, OLD_SHA, 'authenticated_old_eight_worker')
    need(old.pin(OLD_HANDOFF)['sha256'] == OLD_HANDOFF_SHA, 'ActualEightWorkerHandoff')
    oldh = old.read(OLD_HANDOFF)
    validated = old.validate_handoff(oldh)
    def controller_binding(ops):
        driver = s.pin(OLD_DRIVER)
        need(driver['sha256'] == OLD_SHA, 'ReviewedActualEightPredecessor')
        start_pin = s.pin(s.BASE / 'CONTROLLER_STARTED.json')
        started = s.read(start_pin['path'])
        need(type(started['pid']) is int and started['pid'] == s.CONTROLLER and
             type(started['start_ticks']) is int and started['start_ticks'] == s.BIRTH and
             started['parent'] == 11843 and started['actual_driver'] == driver and
             started['handoff'] == s.pin(OLD_HANDOFF) and
             started['argv'][0] == str(OLD_DRIVER) and '--execute' in started['argv'],
             'ActualEightWorkerLaunchSchema')
        binding = dict(pid=s.CONTROLLER, birth=s.BIRTH,
            argv=[started['executable'], '-I', '-B', *started['argv']],
            started=start_pin, driver=driver)
        s.match_controller(ops.proc(s.CONTROLLER), binding)
        return binding
    def check_boundary(number):
        status = s.read(s.BASE / 'STATUS.json')
        closed_boundary_status(status, number, s.CONTROLLER)
        n = number - 1
        need(n >= 68 and n <= 200, 'ClosedEightWorkerMCBoundary')
        if n == 68:
            previous = {cell: oldh['initial_previous'][cell]['prefix'] for cell in old.CELLS}
            boundary = oldh['completed_predecessor_rounds'][-1]['whole']
        else:
            prior = old.read(s.BASE / ('round_%04d' % (n-1)) / 'ROUND_ACCEPTED.json')
            previous, boundary = prior['next_previous'], prior['whole']
            old.bound(boundary)
            for value in previous.values(): old.bound(value)
        old.validate_saved_round(s.BASE / ('round_%04d' % n), n, oldh, previous,
            boundary, validated, controller_pid=s.CONTROLLER)
    s.controller_binding = controller_binding
    s.check_boundary = check_boundary
    return s


def exact_wait(status, number, controller):
    need(type(number) is int and 69 <= number <= 201 and
         type(status['owner_pid']) is int and status['owner_pid'] == controller and
         status['status'] == 'WAIT_RESOURCE_WINDOW_NO_NATIVE_RESTART' and
         type(status.get('next_round')) is int and status['next_round'] == number,
         'ExactEightWorkerChildlessWait')


def closed_boundary_status(status, number, controller):
    """Accept a completed prior-round scheduling gap, never the new round."""
    if status['status'] == 'RUNNING_BOUND_EIGHT_WORKER_ROUND':
        need(type(number) is int and 69 <= number <= 201 and
             type(status['owner_pid']) is int and status['owner_pid'] == controller and
             type(status.get('round')) is int and status['round'] == number-1,
             'OnlyClosedPriorRoundSchedulingGap')
    else:
        exact_wait(status, number, controller)


def wait_boundary(number, seconds, s, ops):
    need(type(seconds) is int and 0 <= seconds <= 3600, 'FiniteBoundaryObservation')
    binding = s.controller_binding(ops)
    deadline = ops.now() + seconds
    while True:
        s.match_controller(ops.proc(s.CONTROLLER), binding)
        status = s.read(s.BASE / 'STATUS.json')
        waiting = status['status'] == 'WAIT_RESOURCE_WINDOW_NO_NATIVE_RESTART' and status.get('next_round') == number
        closed_gap = (status['status'] == 'RUNNING_BOUND_EIGHT_WORKER_ROUND' and
            status.get('round') == number-1 and
            (s.BASE / ('round_%04d' % (number-1)) / 'ROUND_ACCEPTED.json').exists())
        if waiting or (closed_gap and not ops.children(s.CONTROLLER)):
            s.check_boundary(number)
            need(not ops.children(s.CONTROLLER), 'NoActiveFitsAtRequestedBoundary')
            return binding
        current = status.get('round') if status['status'] == 'RUNNING_BOUND_EIGHT_WORKER_ROUND' else status.get('next_round')
        need(type(current) is int and current <= number and
             not (status['status'] == 'RUNNING_BOUND_EIGHT_WORKER_ROUND' and current == number),
             'RequestedBoundaryPassedNoReplay')
        need(ops.now() < deadline, 'ObservationExpiredNoSignalOrRestart')
        ops.sleep(min(.05, max(.01, deadline-ops.now())))


def source_pins(intent, m):
    for name in ('coordinator_source', 'successor_driver', 'safety_source',
                 'kernel_source', 'prepared_template', 'fresh_prepared'):
        m.bound(intent[name])

def retirement_prerequisite(h, m):
    """Original observer must still own the actual eight-worker child."""
    record = h['predecessor_controller']['parent_started']
    m.bound(record)
    started = m.read(record['path'])
    identity = started['coordinator_identity']
    need(identity['pid'] == 11843 and identity['birth'] == 2613447384 and
         started['actual_pid'] == 11848 and started['start_ticks'] == 2613447611,
         'ActualEightWorkerObservingParent')
    current = m.proc(11843)
    need(current['start_ticks'] == 2613447384 and current['state'] not in ('Z','X'),
         'ActualOriginalObserverStillLiveBeforeRetirement')
    return record


def transition(h, number, folder, intent, binding, m, s, ops, k, guarded):
    retirement_prerequisite(h, m)
    guarded(); s.pause(binding, ops, lambda: s.check_boundary(number))
    guarded(); s.paused(binding, number, ops)
    validated = m.validate_handoff(h, pre_retirement=True)
    resource = m.inspect_resources(h, validated)
    resource_pin = m.write_new(folder / 'RESOURCES_BEFORE_RETIREMENT.json', resource)
    need(resource['resources_admitted'] is True, 'ResourcesUnavailableResumeWithoutTERM')
    retirement_prerequisite(h, m)
    guarded(); s.paused(binding, number, ops)
    m.write_new(folder / 'RETIRE_INTENT.json', {'scope': 'CHILDLESS_SCHEDULER_ONLY',
        'controller': binding, 'expected_round': number, 'resources': resource_pin,
        'prepared_material': intent['fresh_prepared'], 'old_success_not_asserted': True})
    terminal = k.retire(binding, number, s, ops, note=lambda event: m.write_new(
        folder / (event + '.json'), {'event': event, 'controller': binding, 'utc': m.utc()}))
    terminal.update(schema='jgw03-eleven-worker-predecessor-terminal-v1',
        status='ACTUAL_PROCESS_EXIT_OBSERVED_FOR_PLANNED_HANDOFF',
        driver=h['predecessor_driver'], started=h['predecessor_controller']['started'],
        completed_rounds=h['completed_rounds'], coordinator_source=intent['coordinator_source'],
        resource_admission=resource_pin, utc=m.utc())
    terminal_pin = m.write_new(folder / 'PREDECESSOR_NATIVE_TERMINAL.json', terminal)
    # Disarm uses the real pidfd record, not the different parent's returncode.
    m.write_new(folder / 'DISARM.json', {'reason': 'ORIGINAL_PIDFD_EXIT_OBSERVED',
        'legacy_terminal': terminal_pin})
    handed = copy.deepcopy(h); handed['status'] = 'READY_AFTER_CLOSED_PREFIX_HANDOFF'
    handed['predecessor_controller']['terminal'] = terminal_pin
    end = min(ops.now()+20, intent['deadline_monotonic'])
    while not PARENT_EXIT.exists():
        need(ops.now() < end, 'RealPredecessorParentReceiptNotYetPublished')
        ops.sleep(.05)
    s.published(PARENT_EXIT)
    parent_pin = m.pin(PARENT_EXIT)
    m.validate_predecessor_parent_exit(handed, parent_pin)
    handed['predecessor_parent_exit'] = parent_pin
    m.validate_handoff(handed)
    handoff_pin = m.write_new(folder / 'HANDOFF.json', handed)
    return handed, handoff_pin

def observe(child, argv, h, hp, folder, attempt, intent, m, ops):
    """Drain immediately, always wait even metadata/drain failures; no timeout."""
    streams = {}; failures = []; threads = []; identity = None; started = None; coordinator = None
    def drain(name, handle):
        digest = hashlib.sha256(); size = 0
        try:
            for block in iter(lambda: handle.read(1024**2), b''):
                digest.update(block); size += len(block)
        except BaseException as exc:
            failures.append(type(exc).__name__)
        finally:
            streams[name] = {'bytes': size, 'sha256': digest.hexdigest()}
            try:
                handle.close()
            except BaseException as exc:
                failures.append(type(exc).__name__)
    for name in ('stdout', 'stderr'):
        try:
            thread = threading.Thread(target=drain, args=(name, getattr(child, name)), daemon=True)
            thread.start(); threads.append(thread)
        except BaseException as exc:
            failures.append(type(exc).__name__)
            # Failure to start a drain must not retain a pipe that can fill.
            try:
                getattr(child, name).close()
            except BaseException as close_exc:
                failures.append(type(close_exc).__name__)
    try:
        identity = ops.proc(child.pid)
        need(identity['argv'] == argv and identity['parent'] == os.getpid(), 'DirectChildIdentity')
        own = ops.proc(os.getpid())
        coordinator = {'pid': own['pid'], 'parent': own['parent'], 'birth': own['birth'],
                       'start_ticks': own['birth']}
        started = m.write_new(folder / 'SUCCESSOR_STARTED.json', {
            'schema': 'jgw03-successor-parent-started-v1', 'actual_pid': child.pid,
            'start_ticks': identity['birth'], 'argv': argv, 'driver': h['actual_driver'],
            'handoff': hp, 'coordinator_source': intent['coordinator_source'],
            'coordinator_identity': coordinator,
            'actual_native_exit_not_yet_observed': True, 'utc': m.utc()})
    except BaseException as exc:
        failures.append(type(exc).__name__)
    actual = None
    # A wait interruption is metadata failure, never permission to abandon child.
    while actual is None:
        try:
            actual = child.wait()
        except BaseException as exc:
            failures.append(type(exc).__name__)
    for thread in threads:
        thread.join()
    dispatch = None
    try:
        path = ROOT / 'provenance' / attempt / 'DISPATCH_EXIT.json'
        if path.exists():
            dispatch = m.pin(path)
    except BaseException as exc:
        failures.append(type(exc).__name__)
    result = {'schema': 'jgw03-successor-parent-native-exit-v1', 'actual_pid': child.pid,
        'start_ticks': identity['birth'] if identity else None, 'argv': argv,
        'driver': h['actual_driver'], 'handoff': hp, 'started': started,
        'coordinator_identity': coordinator,
        'coordinator_source': intent['coordinator_source'], 'actual_native_exit_code': actual,
        'dispatch_exit': dispatch, 'streams': streams, 'raw_streams_persisted': False,
        'observer_failure_class_only': failures[0] if failures else None,
        'old_controller_success_asserted': False, 'scientific_acceptance': False, 'utc': m.utc()}
    for name in ('stdout', 'stderr'):
        result[name+'_bytes'] = streams.get(name, {}).get('bytes')
        result[name+'_sha256'] = streams.get(name, {}).get('sha256')
    m.write_new(folder / 'SUCCESSOR_NATIVE_EXIT.json', result)
    return 0 if type(actual) is int and actual == 0 and dispatch and not failures else 1

def tools(intent):
    m = load(intent['successor_driver']['path'], intent['successor_driver']['sha256'], 'eleven_driver')
    k = load(intent['kernel_source']['path'], KERNEL_SHA, 'reviewed_boundary_kernels')
    s = safety(intent['safety_source']['path']); ops = k.Linux(s)
    k.source_pins = source_pins
    return m, k, s, ops

def lease(folder, expected):
    raw = (folder / 'INTENT.json').read_bytes()
    need(hashlib.sha256(raw).hexdigest() == expected, 'ExactLeaseIntent')
    intent = json.loads(raw); m, k, s, ops = tools(intent)
    need(sys.platform.startswith('linux') and socket.gethostname() == HOST and
         intent['schema'] == 'jgw03-eleven-worker-boundary-coordinator-v1' and
         Path(intent['coordinator_source']['path']) == Path(__file__).resolve() and
         0 < intent['deadline_monotonic']-ops.now() <= 600 and
         folder.parent == ROOT / 'provenance', 'FiniteOwnedIndependentLease')
    return k.lease_loop(intent, folder, m, s, ops)

def coordinate(args):
    need(sys.platform.startswith('linux') and socket.gethostname() == HOST, 'DeclaredHostOnly')
    need(30 <= args.max_stop_seconds <= 600 and 69 <= args.expected_round <= 200 and
         all(re.fullmatch('[A-Za-z0-9_-]{1,80}', x) for x in
             (args.attempt, args.successor_attempt)), 'FiniteFreshDeclaredBoundary')
    raw = args.prepared.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == args.prepared_sha256, 'PreparedSHA')
    template = json.loads(raw)
    need(template['actual_driver']['sha256'] == args.driver_sha256, 'ExplicitReviewedDriverSHA')
    intent = {'successor_driver': template['actual_driver'],
        'kernel_source': {'path': str(args.kernel)}, 'safety_source': {'path': str(args.safety)}}
    m, k, s, ops = tools(intent)
    m.validate_handoff(template, pre_retirement=True)
    retirement_prerequisite(template, m)
    folder = ROOT / 'provenance' / args.attempt
    need(not folder.exists() and not (ROOT / 'provenance' / args.successor_attempt).exists(), 'NoAttemptReplay')
    import fcntl
    with (ROOT / 'provenance/ELEVEN_WORKER_HANDOFF_COORDINATOR.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        binding = wait_boundary(args.expected_round, args.max_wait_seconds, s, ops)
        h = m.build_boundary_material(template, args.expected_round)
        m.validate_handoff(h, pre_retirement=True)
        folder.mkdir(); prepared = m.write_new(folder / 'PREPARED_MATERIALS.json', h)
        intent.update(schema='jgw03-eleven-worker-boundary-coordinator-v1', controller=binding,
            caller=ops.proc(os.getpid()), expected_round=args.expected_round,
            coordinator_source=m.pin(Path(__file__).resolve()), kernel_source=m.pin(args.kernel),
            safety_source=m.pin(args.safety), prepared_template=m.pin(args.prepared),
            fresh_prepared=prepared, deadline_monotonic=ops.now()+args.max_stop_seconds-10,
            successor_attempt=args.successor_attempt, utc=m.utc())
        ip = m.write_new(folder / 'INTENT.json', intent)
        guard = None; handed = None; observing = False
        def interrupted(sig, frame):
            if not observing:
                raise InterruptedError('CoordinatorSignal')
        handlers = {sig: signal.signal(sig, interrupted) for sig in (signal.SIGTERM, signal.SIGINT)}
        try:
            guard = subprocess.Popen([sys.executable, '-I', '-B', str(Path(__file__).resolve()),
                '_lease', str(folder), ip['sha256']], start_new_session=True,
                stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            end = min(ops.now()+15, intent['deadline_monotonic'])
            while not (folder / 'LEASE_READY.json').exists():
                need(ops.now() < end and guard.poll() is None, 'LeaseREADYBeforeSTOP'); ops.sleep(.05)
            need(s.published(folder / 'LEASE_READY.json')['pid'] == guard.pid, 'ActualLeaseChild')
            guarded = lambda: k.finite_ok(intent, folder, m, s, ops)
            guarded()
            m.write_new(folder / 'PAUSE_OWNED.json', {'controller': binding,
                'restoration_armed_before_STOP': True})
            handed, hp = transition(h, args.expected_round, folder, intent, binding, m, s, ops, k, guarded)
            actual = guard.wait(timeout=5)
            m.write_new(folder / 'LEASE_NATIVE_EXIT.json', {'actual_pid': guard.pid,
                'actual_native_exit_code': actual})
            need(type(actual) is int and actual == 0, 'LeaseDisarmedActualExit')
            argv = [sys.executable, '-I', '-B', handed['actual_driver']['path'], '--execute',
                '--handoff', hp['path'], '--handoff-sha256', hp['sha256'], '--attempt', args.successor_attempt]
            m.bound(handed['actual_driver']); m.bound(hp); observing = True
            child = subprocess.Popen(argv, cwd=ROOT, start_new_session=True, stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1'))
            return observe(child, argv, handed, hp, folder, args.successor_attempt, intent, m, ops)
        except BaseException as exc:
            if not (folder / 'DISARM.json').exists():
                m.write_new(folder / 'ABANDON.json', {'exception_class_only': type(exc).__name__})
            restored = k.safe_resume(binding, s, ops)
            m.write_new(folder / 'COORDINATOR_STOP.json', {'exception_class_only': type(exc).__name__,
                'resume': restored, 'no_automatic_retry': True, 'raw_error_persisted': False})
            return 1
        finally:
            if handed is None:
                k.safe_resume(binding, s, ops)
            for sig, handler in handlers.items():
                signal.signal(sig, handler)

def main():
    if len(sys.argv) == 4 and sys.argv[1] == '_lease':
        return lease(Path(sys.argv[2]), sys.argv[3])
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--execute', action='store_true')
    for name in ('prepared', 'kernel', 'safety'):
        p.add_argument('--'+name, type=Path)
    for name in ('prepared-sha256', 'driver-sha256', 'attempt', 'successor-attempt'):
        p.add_argument('--'+name)
    p.add_argument('--expected-round', type=int); p.add_argument('--max-wait-seconds', type=int, default=1800)
    p.add_argument('--max-stop-seconds', type=int, default=120); args = p.parse_args()
    if not args.execute:
        print(json.dumps({'status': 'PLAN_ONLY_NO_SIGNALS_OR_SCIENTIFIC_ACTION', 'observer_timeout': None})); return 0
    need(all(getattr(args, name) is not None for name in ('prepared', 'kernel', 'safety',
        'prepared_sha256', 'driver_sha256', 'attempt', 'successor_attempt', 'expected_round')), 'ExplicitInputs')
    return coordinate(args)

if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(json.dumps({'status': 'COORDINATOR_STOP', 'exception_class_only': type(exc).__name__})); raise SystemExit(2)
