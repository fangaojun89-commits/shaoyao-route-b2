"""Plan-only by default; guarded childless resumed-eleven to fifteen handoff.

The original live sources are never edited. STAR must have a genuine successful
queue closure before this coordinator is armed. Only the empty scheduler may be
retired; no R fit is signalled. A pinned independent finite lease restores an
abandoned pre-retirement pause. The successor's direct parent always waits.
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
import types

ROOT = Path('/root/autodl-tmp/JGW03_full_reference_linux_20261017r03')
HOST = 'autodl-container-62c8pzy536-50fb405a'
CONTROLLER, BIRTH, PARENT, PARENT_BIRTH = 57056, 2616685761, 57055, 2616685675
BASE = ROOT/'provenance/eleven_checkpoint_resume93_01'
OLD_DRIVER = ROOT/'eleven_worker_source_20261022_r01/code/python/run_JGW03_eleven_worker_shard_successor_20261022.py'
OLD_SHA = '718d394a97b6052c427855169fd24dab68402ed96ef59cf5338eb2259151e22c'
RUNNER = ROOT/'checkpoint_resume_source_20261022_r01/resume_JGW03_eleven_worker_checkpoint_20261022.py'
RUNNER_SHA = 'cea59028f4e94d6b366d32c1006c6cdaf0a33b8fa29683e71fa4c35b91d4c998'
GENERIC = ROOT/'eleven_worker_source_20261022_r01/code/python/handoff_JGW03_eleven_worker_boundary_20261022.py'
GENERIC_SHA = '9ae73f398f899fa21d72bf8c41a59b902e9451e1abae5482a67f83242cff57b0'
KERNEL_SHA = '8c0640997f1ad2a917b5f57fdafdcf6a65cb7d3574445f5dae0e6416673520d7'
SAFETY_SHA = 'fe9f804c162a26e34c1f5fd52d5c436b8d72bd81fe19a032691a741a1323de53'
PARENT_EXIT = BASE/'ACTUAL_PARENT_NATIVE_EXIT.json'


def need(ok, label):
    if not ok:
        raise RuntimeError(label)


def load(path, sha, name):
    raw = Path(path).read_bytes()
    need(hashlib.sha256(raw).hexdigest() == sha, 'ReviewedExecutedSourceSHA')
    value = types.ModuleType(name); value.__file__ = str(path)
    exec(compile(raw, str(path), 'exec'), value.__dict__)
    return value


def exact_wait(status, number, controller):
    need(type(number) is int and 94 <= number <= 201 and
         type(status.get('owner_pid')) is int and status['owner_pid'] == controller and
         status.get('status') == 'WAIT_RESOURCE_WINDOW_NO_NATIVE_RESTART' and
         type(status.get('next_round')) is int and status['next_round'] == number,
         'ExactResumedElevenChildlessWait')


def closed_boundary_status(status, number, controller):
    if status.get('status') == 'RUNNING_BOUND_ELEVEN_WORKER_ROUND':
        need(type(number) is int and 94 <= number <= 201 and
             type(status.get('owner_pid')) is int and status['owner_pid'] == controller and
             type(status.get('round')) is int and status['round'] == number-1,
             'OnlyClosedPriorElevenRoundSchedulingGap')
    else:
        exact_wait(status, number, controller)


def safety(path, template, m):
    """Operational adapter only. Full old scientific chain uses its own validator."""
    s = load(path, SAFETY_SHA, 'fifteen_pidfd_safety')
    s.CONTROLLER, s.BIRTH, s.BASE = CONTROLLER, BIRTH, BASE
    s.DRIVER, s.DRIVER_SHA = OLD_DRIVER, OLD_SHA

    def controller_binding(ops):
        pin = m.bound(template['old_resume_controller_started'])
        started = m.read(pin['path'])
        need(started['schema'] == 'jgw03-same-eleven-epoch-checkpoint-controller-v1' and
             type(started['pid']) is int and started['pid'] == CONTROLLER and
             type(started['start_ticks']) is int and started['start_ticks'] == BIRTH and
             type(started['parent']) is int and started['parent'] == PARENT and
             started['actual_driver'] == m.pin(OLD_DRIVER) and
             started['continuation_runner'] == m.pin(RUNNER) and
             started['continuation_runner']['sha256'] == RUNNER_SHA,
             'ActualResumedElevenControllerSourceAndParent')
        binding = dict(pid=CONTROLLER, birth=BIRTH, argv=started['argv'],
                       started=pin, driver=started['actual_driver'])
        s.match_controller(ops.proc(CONTROLLER), binding)
        return binding

    def check_boundary(number):
        closed_boundary_status(m.read(BASE/'STATUS.json'), number, CONTROLLER)
        # Cheap pause eligibility only: the complete old scientific chain is
        # validated after STOP, before TERM. This avoids spending seconds in a
        # short natural scheduling gap while the next fit could be launched.
        accepted=m.read(BASE/f'round_{number-1:04d}/ROUND_ACCEPTED.json')
        whole=m.read(m.bound(accepted['whole'])['path'])
        need(accepted['scope']=='ACTUAL_NATIVE_ZERO_AND_ORDERED_PREFIX_BINDING_ONLY' and
             accepted['scientific_acceptance'] is False and
             accepted['old_controller_success_asserted'] is False and
             whole['all_actual_exits_zero_and_observed'] is True and
             whole['actual_driver']==m.pin(OLD_DRIVER) and len(whole['results'])==3 and
             all(type(x['actual_native_exit_code'])is int and x['actual_native_exit_code']==0
                 and type(x['actual_R_exit_code'])is int and x['actual_R_exit_code']==0
                 and x['observer_failure_class_only']is None for x in whole['results']),
             'KnownClosedPriorNativeBatchPauseEligibilityOnly')

    s.controller_binding, s.check_boundary = controller_binding, check_boundary
    return s


def wait_boundary(number, seconds, s, ops):
    need(type(seconds) is int and 0 <= seconds <= 3600, 'FiniteBoundaryObservation')
    binding = s.controller_binding(ops); deadline = ops.now()+seconds
    while True:
        s.match_controller(ops.proc(CONTROLLER), binding)
        status = s.read(BASE/'STATUS.json')
        waiting = status.get('status') == 'WAIT_RESOURCE_WINDOW_NO_NATIVE_RESTART' and status.get('next_round') == number
        gap = (status.get('status') == 'RUNNING_BOUND_ELEVEN_WORKER_ROUND' and
               status.get('round') == number-1 and (BASE/f'round_{number-1:04d}/ROUND_ACCEPTED.json').exists())
        if waiting or (gap and not ops.children(CONTROLLER)):
            s.check_boundary(number)
            need(not ops.children(CONTROLLER), 'NoActiveFitsAtRequestedBoundary')
            return binding
        current = status.get('round') if status.get('status') == 'RUNNING_BOUND_ELEVEN_WORKER_ROUND' else status.get('next_round')
        need(type(current) is int and current <= number and not
             (status.get('status') == 'RUNNING_BOUND_ELEVEN_WORKER_ROUND' and current == number),
             'RequestedBoundaryPassedNoReplay')
        need(ops.now() < deadline, 'ObservationExpiredNoSignalOrRestart')
        ops.sleep(min(.05, max(.01, deadline-ops.now())))


def source_pins(intent, m):
    for key in ('coordinator_source', 'successor_driver', 'safety_source', 'kernel_source',
                'generic_coordinator_source', 'prepared_template', 'fresh_prepared', 'STAR_closure'):
        m.bound(intent[key])


def retirement_prerequisite(h, m):
    record = m.bound(h['old_resume_parent_started']); started = m.read(record['path'])
    own = started['observer']
    need(started['schema'] == 'jgw03-same-eleven-epoch-parent-started-v1' and
         type(started['actual_pid']) is int and started['actual_pid'] == CONTROLLER and
         type(started['start_ticks']) is int and started['start_ticks'] == BIRTH and
         type(own['pid']) is int and own['pid'] == PARENT and
         type(own['start_ticks']) is int and own['start_ticks'] == PARENT_BIRTH and
         started['runner'] == h['old_resume_runner'] and
         started['original_driver'] == h['old_driver'], 'ActualOriginalResumedObserver')
    current = m.proc(PARENT)
    need(current['start_ticks'] == PARENT_BIRTH and current['state'] not in ('Z','X') and
         current['argv'] == own['argv'], 'OriginalObserverStillOwnsActualChild')
    return record


def _transition(h, number, folder, intent, binding, m, s, ops, k, guarded):
    # STAR closure/full stable base were checked before entering wait_boundary.
    # The armed source guard rechecks that exact closure record before STOP.
    # Full dynamic history validation follows STOP, never after TERM.
    retirement_prerequisite(h, m)
    guarded(); s.pause(binding, ops, lambda: s.check_boundary(number))
    guarded(); s.paused(binding, number, ops)
    h=m.build_boundary_material(h,number)
    validated = m.validate_handoff(h, pre_retirement=True)
    prepared=m.write_new(folder/'BOUNDARY_MATERIALS.json',h)
    state = m.inspect_resources(h, validated)
    resource = m.write_new(folder/'RESOURCES_BEFORE_RETIREMENT.json', state)
    need(state['resources_admitted'] is True, 'ResourcesUnavailableResumeWithoutTERM')
    retirement_prerequisite(h, m)
    guarded(); s.paused(binding, number, ops)
    m.write_new(folder/'RETIRE_INTENT.json', dict(scope='CHILDLESS_SCHEDULER_ONLY',
        controller=binding, expected_round=number, resources=resource,
        prepared_material=prepared, old_success_not_asserted=True))
    terminal = k.retire(binding, number, s, ops, note=lambda event: m.write_new(
        folder/(event+'.json'), dict(event=event, controller=binding, utc=m.utc())))
    terminal.update(schema='jgw03-fifteen-worker-predecessor-terminal-v1',
        status='ACTUAL_PROCESS_EXIT_OBSERVED_FOR_PLANNED_HANDOFF', driver=h['old_driver'],
        runner=h['old_resume_runner'], started=h['old_resume_controller_started'],
        completed_rounds=number-1, next_round=number, no_active_fit_children=True,
        previous_boundary=h['previous_boundary'], initial_previous=h['initial_previous'],
        coordinator_source=intent['coordinator_source'], handoff_resource_admission=resource,
        prepared_material=prepared,utc=m.utc())
    terminal_pin = m.write_new(folder/'PREDECESSOR_NATIVE_TERMINAL.json', terminal)
    m.write_new(folder/'DISARM.json', dict(reason='ORIGINAL_PIDFD_EXIT_OBSERVED', legacy_terminal=terminal_pin))
    handed = copy.deepcopy(h); handed['status'] = 'ACCEPTED_CLOSED_ELEVEN_BOUNDARY_STAR_COMPLETE'
    handed['predecessor_retirement'] = terminal_pin
    end = min(ops.now()+20, intent['deadline_monotonic'])
    while not PARENT_EXIT.exists():
        need(ops.now() < end, 'RealPredecessorParentReceiptNotYetPublished'); ops.sleep(.05)
    s.published(PARENT_EXIT)
    handed['old_resume_parent_exit'] = m.pin(PARENT_EXIT)
    # The original genuine parent must report -15, never an invented zero.
    m.validate_handoff(handed)
    return handed, m.write_new(folder/'HANDOFF.json', handed)


def transition(h, number, folder, intent, binding, m, s, ops, k, guarded):
    try:
        return _transition(h,number,folder,intent,binding,m,s,ops,k,guarded)
    except BaseException:
        # Also protect direct callers/tests, not just coordinate()'s outer
        # finally. The pinned helper never signals a gone/reused process.
        k.safe_resume(binding,s,ops)
        raise


def tools(intent):
    m = load(intent['successor_driver']['path'], intent['successor_driver']['sha256'], 'fifteen_driver')
    k = load(intent['kernel_source']['path'], KERNEL_SHA, 'fifteen_boundary_kernel')
    generic = load(GENERIC, GENERIC_SHA, 'unchanged_generic_observer')
    template = m.read(intent['prepared_template']['path'])
    s = safety(intent['safety_source']['path'], template, m); ops = k.Linux(s)
    k.source_pins = source_pins
    return m, k, s, ops, generic


def lease(folder, expected):
    raw = (folder/'INTENT.json').read_bytes()
    need(hashlib.sha256(raw).hexdigest() == expected, 'ExactLeaseIntent')
    intent = json.loads(raw); m,k,s,ops,generic = tools(intent)
    need(sys.platform.startswith('linux') and socket.gethostname() == HOST and
         intent['schema'] == 'jgw03-fifteen-worker-boundary-coordinator-v1' and
         Path(intent['coordinator_source']['path']) == Path(__file__).resolve() and
         0 < intent['deadline_monotonic']-ops.now() <= 600 and folder.parent == ROOT/'provenance',
         'FiniteOwnedIndependentLease')
    return k.lease_loop(intent, folder, m, s, ops)


def coordinate(args):
    need(sys.platform.startswith('linux') and socket.gethostname() == HOST, 'DeclaredHostOnly')
    need(type(args.expected_round) is int and 94 <= args.expected_round <= 201 and
         30 <= args.max_stop_seconds <= 600 and all(re.fullmatch('[A-Za-z0-9_-]{1,80}', x)
         for x in (args.attempt, args.successor_attempt)), 'FiniteFreshDeclaredBoundary')
    raw = args.prepared.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == args.prepared_sha256, 'PreparedSHA')
    template = json.loads(raw)
    need(template['actual_driver']['sha256'] == args.driver_sha256, 'ExplicitReviewedDriverSHA')
    intent = dict(successor_driver=template['actual_driver'], kernel_source=dict(path=str(args.kernel)),
                  safety_source=dict(path=str(args.safety)), prepared_template=dict(path=str(args.prepared)))
    m,k,s,ops,generic = tools(intent)
    # Source/stable-base preflight may coexist with the old active N+1. This
    # relaxation is never used after STOP or by executable handoff acceptance.
    m.validate_handoff(template, pre_retirement=True,allow_active_future=True)
    retirement_prerequisite(template,m)
    folder = ROOT/'provenance'/args.attempt
    need(not folder.exists() and not (ROOT/'provenance'/args.successor_attempt).exists(), 'NoAttemptReplay')
    import fcntl
    with (ROOT/'provenance/FIFTEEN_WORKER_HANDOFF_COORDINATOR.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        binding = wait_boundary(args.expected_round,args.max_wait_seconds,s,ops)
        h=copy.deepcopy(template)
        folder.mkdir(); prepared = m.write_new(folder/'PREPARED_MATERIALS.json',h)
        intent.update(schema='jgw03-fifteen-worker-boundary-coordinator-v1',controller=binding,
            caller=ops.proc(os.getpid()),expected_round=args.expected_round,
            coordinator_source=m.pin(Path(__file__).resolve()),kernel_source=m.pin(args.kernel),
            safety_source=m.pin(args.safety),generic_coordinator_source=m.pin(GENERIC),
            prepared_template=m.pin(args.prepared),fresh_prepared=prepared,STAR_closure=h['STAR_closure'],
            deadline_monotonic=ops.now()+args.max_stop_seconds-10,
            successor_attempt=args.successor_attempt,utc=m.utc())
        ip = m.write_new(folder/'INTENT.json',intent)
        guard=None;handed=None;observing=False
        def interrupted(sig,frame):
            if not observing: raise InterruptedError('CoordinatorSignal')
        handlers={sig:signal.signal(sig,interrupted) for sig in (signal.SIGTERM,signal.SIGINT)}
        try:
            guard=subprocess.Popen([sys.executable,'-I','-B',str(Path(__file__).resolve()),'_lease',
                str(folder),ip['sha256']],start_new_session=True,stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
            end=min(ops.now()+15,intent['deadline_monotonic'])
            while not (folder/'LEASE_READY.json').exists():
                need(ops.now()<end and guard.poll() is None,'LeaseREADYBeforeSTOP');ops.sleep(.05)
            need(s.published(folder/'LEASE_READY.json')['pid']==guard.pid,'ActualLeaseChild')
            guarded=lambda:k.finite_ok(intent,folder,m,s,ops)
            guarded();m.write_new(folder/'PAUSE_OWNED.json',dict(controller=binding,restoration_armed_before_STOP=True))
            handed,hp=transition(h,args.expected_round,folder,intent,binding,m,s,ops,k,guarded)
            actual=guard.wait(timeout=5)
            m.write_new(folder/'LEASE_NATIVE_EXIT.json',dict(actual_pid=guard.pid,actual_native_exit_code=actual))
            need(type(actual)is int and actual==0,'LeaseDisarmedActualExit')
            argv=[sys.executable,'-I','-B',handed['actual_driver']['path'],'--execute','--handoff',
                  hp['path'],'--handoff-sha256',hp['sha256'],'--attempt',args.successor_attempt]
            m.bound(handed['actual_driver']);m.bound(hp);observing=True
            child=subprocess.Popen(argv,cwd=ROOT,start_new_session=True,stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
            return generic.observe(child,argv,handed,hp,folder,args.successor_attempt,intent,m,ops)
        except BaseException as exc:
            if not (folder/'DISARM.json').exists():
                m.write_new(folder/'ABANDON.json',dict(exception_class_only=type(exc).__name__))
            restored=k.safe_resume(binding,s,ops)
            m.write_new(folder/'COORDINATOR_STOP.json',dict(exception_class_only=type(exc).__name__,
                resume=restored,no_automatic_retry=True,raw_error_persisted=False))
            return 1
        finally:
            if handed is None:k.safe_resume(binding,s,ops)
            for sig,handler in handlers.items():signal.signal(sig,handler)


def main():
    if len(sys.argv)==4 and sys.argv[1]=='_lease':return lease(Path(sys.argv[2]),sys.argv[3])
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--execute',action='store_true')
    for name in ('prepared','kernel','safety'):p.add_argument('--'+name,type=Path)
    for name in ('prepared-sha256','driver-sha256','attempt','successor-attempt'):p.add_argument('--'+name)
    p.add_argument('--expected-round',type=int);p.add_argument('--max-wait-seconds',type=int,default=1800)
    p.add_argument('--max-stop-seconds',type=int,default=120);args=p.parse_args()
    if not args.execute:
        print(json.dumps(dict(status='PLAN_ONLY_NO_SIGNALS_OR_SCIENTIFIC_ACTION',observer_timeout=None)));return 0
    need(all(getattr(args,k)is not None for k in ('prepared','kernel','safety','prepared_sha256',
        'driver_sha256','attempt','successor_attempt','expected_round')),'ExplicitInputs')
    return coordinate(args)


if __name__=='__main__':
    try:raise SystemExit(main())
    except Exception as exc:
        print(json.dumps(dict(status='COORDINATOR_STOP',exception_class_only=type(exc).__name__)));raise SystemExit(2)
