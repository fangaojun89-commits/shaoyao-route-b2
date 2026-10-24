"""Same eleven-worker scientific epoch, fresh observer/attempt, no completed fits.

The original immutable driver still owns scientific commands and resource policy.
This runner authenticates one stopped dispatch and resumes its first unused round.
Plan-only by default. No signals, installations, payload deletion or automatic retry.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import types

DRIVER_SHA = '718d394a97b6052c427855169fd24dab68402ed96ef59cf5338eb2259151e22c'
COORDINATOR_SHA = '9ae73f398f899fa21d72bf8c41a59b902e9451e1abae5482a67f83242cff57b0'
ROOT = Path('/root/autodl-tmp/JGW03_full_reference_linux_20261017r03')


def need(ok, label):
    if not ok:
        raise RuntimeError(label)


def load_driver(path):
    raw = Path(path).read_bytes()
    need(hashlib.sha256(raw).hexdigest() == DRIVER_SHA, 'FrozenOriginalElevenDriver')
    d = types.ModuleType('unchanged_eleven_checkpoint_driver')
    d.__file__ = str(path)
    exec(compile(raw, str(path), 'exec'), d.__dict__)
    return d


def remaining_after(plan, completed, first):
    """Pure exact-prefix selector; rejects gaps, duplicates and all-work-complete."""
    need(type(first) is int and first >= 1 and completed, 'NonemptyClosedPrefix')
    need(all(set(x) == {'number', 'whole', 'accepted'} and type(x['number']) is int
             for x in completed), 'ExactClosedRoundRecordShape')
    numbers = [x['number'] for x in completed]
    need(numbers == list(range(first, first+len(numbers))), 'ContiguousClosedPrefixOnly')
    need([n for n, _ in plan[:len(numbers)]] == numbers, 'PrefixOfFrozenPlan')
    remaining = plan[len(numbers):]
    need(remaining and remaining[0][0] == numbers[-1]+1, 'OnlyFirstUnusedRound')
    return remaining


def exact_qc_commands(source, metadata):
    """The two original bam_qc calls missing from commands(), exact paths only."""
    sam=str(source.TOOLS['samtools'][0])
    folders=[source.RUN/'reservoir_01']+[source.RUN/'libraries'/x['SRA_run'] for x in metadata['libraries']]
    return [command for folder in folders for command in
        ([sam,'view','-H',str(folder/'sorted.bam')],[sam,'view',str(folder/'sorted.bam')])]


def reviewed_star_concurrency(d, epoch):
    """Original guard plus exactly the frozen producer's two BAM QC commands.

    Tool digest, one native process, AS limit and disjoint CPU pair stay fixed.
    This is an explicit control guard correction, not a new scientific actor.
    """
    source=d.load_module(epoch['STAR_source'],'same_epoch_STAR_QC_guard')
    metadata=source.metadata();expected=source.commands('index',source.RUN/'index_01')
    expected+=source.commands('library',source.RUN/'reservoir_01',
        [p['path'] for p in metadata['reservoir_fresh_file_pins'][:2]],False,source_sha=d.STAR_SHA)
    for lib in metadata['libraries']:
        expected+=source.commands('library',source.RUN/'libraries'/lib['SRA_run'],
            [p['path'] for p in lib['raw_mates']],True,python=sys.executable,
            source=epoch['STAR_source']['path'],source_sha=d.STAR_SHA)
    expected+=exact_qc_commands(source,metadata)
    toolpaths={str(p[0]):p[1] for p in source.TOOLS.values()};found=[]
    for directory in Path('/proc').iterdir():
        if not directory.name.isdigit():continue
        try:value=d.proc(int(directory.name))
        except (FileNotFoundError,ProcessLookupError):continue
        if value['state'] in ('Z','X') or not value['argv'] or value['argv'][0] not in toolpaths:continue
        need(value['argv'] in expected,'UnknownConcurrentSTAROrCountCommand')
        need(d.pin(value['argv'][0])['sha256']==toolpaths[value['argv'][0]],'FrozenConcurrentNativeTool')
        limits=(directory/'limits').read_text()
        match=re.search(r'^Max address space\s+(\d+)\s+(\d+)\s+bytes\s*$',limits,re.M)
        need(match and int(match[1])==int(match[2])==40*d.GIB and
             set(os.sched_getaffinity(value['pid']))<=set(epoch['STAR_CPU_ids']),'ActualConcurrentASAndAffinity')
        found.append(value)
    need(len(found)<=1,'AtMostOneBoundSTARCountNativeStage')
    return found


def stopped_parent(d, controller, ps, pe, started_pin, dispatch_pin, handoff_pin):
    """Real original coordinator's actual nonzero terminal, never a pidfd claim."""
    driver = controller['actual_driver']
    need(ps['schema'] == 'jgw03-successor-parent-started-v1' and
         pe['schema'] == 'jgw03-successor-parent-native-exit-v1' and
         ps['driver'] == pe['driver'] == driver and
         ps['handoff'] == pe['handoff'] == handoff_pin and
         ps['actual_pid'] == pe['actual_pid'] == controller['pid'] and
         type(controller['pid']) is int and controller['pid'] > 0 and
         ps['start_ticks'] == pe['start_ticks'] == controller['start_ticks'] and
         type(controller['start_ticks']) is int and controller['start_ticks'] > 0 and
         ps['argv'] == pe['argv'] == [controller['executable'], '-I', '-B', *controller['argv']] and
         controller['argv'] == [driver['path'], '--execute', '--handoff', handoff_pin['path'],
             '--handoff-sha256', handoff_pin['sha256'], '--attempt', Path(dispatch_pin['path']).parent.name] and
         pe['started'] == started_pin and pe['dispatch_exit'] == dispatch_pin and
         type(pe['actual_native_exit_code']) is int and pe['actual_native_exit_code'] == 1 and
         pe['observer_failure_class_only'] is None and pe['raw_streams_persisted'] is False and
         pe['old_controller_success_asserted'] is False and pe['scientific_acceptance'] is False and
         ps['actual_native_exit_not_yet_observed'] is True and
         ps['coordinator_source'] == pe['coordinator_source'] and
         d.bound(ps['coordinator_source'])['sha256'] == COORDINATOR_SHA and
         ps['coordinator_identity'] == pe['coordinator_identity'], 'ActualOriginalParentNonzero')
    owner = ps['coordinator_identity']
    need(all(type(owner[k]) is int and owner[k] > 0 for k in ('pid','parent','birth','start_ticks')) and
         owner['birth'] == owner['start_ticks'] and owner['pid'] == controller['parent'], 'OriginalParentIdentity')
    for stream in ('stdout','stderr'):
        need(type(pe[stream+'_bytes']) is int and pe[stream+'_bytes'] >= 0 and
             re.fullmatch('[a-f0-9]{64}', pe[stream+'_sha256']), 'DigestOnlyOriginalStreams')


def authenticate(d, args):
    pins = {k: d.pin(getattr(args, k)) for k in ('dispatch','parent_started','parent_exit')}
    need(all(Path(p['path']).is_relative_to(ROOT/'provenance') for p in pins.values()), 'RunBoundPriorReceipts')
    dispatch = d.read(args.dispatch); prior = args.dispatch.parent
    pins['controller_started'] = d.pin(prior/'CONTROLLER_STARTED.json')
    controller = d.read(pins['controller_started']['path'])
    handoff_pin = d.bound(controller['handoff']); h = d.read(handoff_pin['path'])
    need(controller['actual_driver'] == h['actual_driver'] == d.pin(d.__file__) and
         h['actual_driver']['sha256'] == DRIVER_SHA, 'SameOriginalScientificEpoch')
    ps, pe = d.read(args.parent_started), d.read(args.parent_exit)
    stopped_parent(d, controller, ps, pe, pins['parent_started'], pins['dispatch'], handoff_pin)
    try:
        live = d.proc(controller['pid'])
    except (FileNotFoundError, ProcessLookupError):
        live = None
    need(live is None or live['start_ticks'] != controller['start_ticks'] or live['state'] in ('Z','X'), 'PriorControllerActuallyTerminal')
    need(dispatch['status'] == 'STOP_RETAINED_NO_AUTOMATIC_NATIVE_RETRY' and
         dispatch['old_controller_success_asserted'] is False and dispatch['scientific_acceptance'] is False and
         dispatch['raw_error_persisted'] is False, 'StoppedDispatchNotCompletedOrForgedZero')
    v = d.validate_handoff(h)
    plan = v['predecessor'].remaining_plan(v['legacy'], h['completed_rounds'])
    records = dispatch['completed_remaining_rounds']
    remaining = remaining_after(plan, records, h['completed_rounds']+1)
    need(remaining[0][0] == args.next_round, 'ExplicitExpectedFirstUnusedRound')
    previous = {c: h['initial_previous'][c]['prefix'] for c in d.CELLS}
    boundary = h['completed_predecessor_rounds'][-1]['whole']
    for record in records:
        n = record['number']
        previous, whole, accepted = d.validate_saved_round(prior/f'round_{n:04d}', n,
            h, previous, boundary, v, controller['pid'])
        need(record['whole'] == whole and record['accepted'] == accepted, 'ExactRealClosedRoundPins')
        boundary = whole
    # A consumed future directory must never be replayed, even if unaccepted.
    need(not any((prior/f'round_{n:04d}').exists() for n, _ in remaining), 'NoPartialPriorFutureRound')
    for p in pins.values():
        d.bound(p)
    return h, handoff_pin, v, previous, boundary, remaining, pins


def child(args):
    need(sys.platform.startswith('linux') and re.fullmatch('[A-Za-z0-9_-]{1,80}', args.attempt), 'LinuxFreshAttempt')
    d = load_driver(args.driver)
    folder = ROOT/'provenance'/args.attempt
    need(folder.is_dir() and not (folder/'CONTROLLER_STARTED.json').exists(), 'FreshOwnedObserverFolder')
    import fcntl
    with (ROOT/'provenance/FULL_REFERENCE.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        h,hp,v,previous,boundary,remaining,pins = authenticate(d,args)
        identity = d.proc(os.getpid())
        runner = d.pin(Path(__file__).resolve())
        # Deliberate, source-pinned control-only replacement; original R/functions
        # and original resource arithmetic are unchanged, never monkey-patch science.
        d.star_concurrency=lambda epoch:reviewed_star_concurrency(d,epoch)
        started = d.write_new(folder/'CONTROLLER_STARTED.json', dict(
            schema='jgw03-same-eleven-epoch-checkpoint-controller-v1',
            **{k: identity[k] for k in ('pid','parent','start_ticks')}, argv=identity['argv'], parsed_argv=sys.argv,
            executable=sys.executable, actual_driver=h['actual_driver'], continuation_runner=runner,
            handoff=hp, prior_receipts=pins, assigned_CPUs=h['worker_cpu_ids'],
            resource_epoch=h['resource_epoch'], first_unused_round=args.next_round,
            resource_guard_source=runner,resource_guard_policy='EXACT_ORIGINAL_STAR_COMMANDS_AND_BAM_QC_V1',
            old_controller_success_asserted=False, scientific_acceptance=False, utc=d.utc()))
        d.write_new(folder/'PLAN.json', dict(remaining=remaining, first_unused_round=args.next_round,
            same_scientific_epoch=True, maximum_workers=11, original_driver=h['actual_driver'],
            inherited_closed_rounds=d.read(args.dispatch)['completed_remaining_rounds']))
        completed=[]; phase='RESOURCE_ADMISSION'
        try:
            for number,jobs in remaining:
                phase='RESOURCE_ADMISSION'
                while True:
                    for p in (hp,runner,*pins.values()): d.bound(p)
                    for key in ('actual_driver','fit_R','merge_R','new_R','finish_R','increment_helper','review','resource_epoch'):
                        d.bound(h[key])
                    for key in ('legacy_driver','admission','legacy_R','config','runtime','legacy_acceptance','keys_sidecar'):
                        d.bound(v['context'][key])
                    try:
                        state=d.inspect_resources(h,v)
                    except (FileNotFoundError,ProcessLookupError):
                        # /proc and transient STAR files can disappear between reads.
                        # No native job exists here; wait, then remeasure everything.
                        v['predecessor'].progress(folder,dict(status='WAIT_TRANSIENT_RESOURCE_SNAPSHOT',next_round=number))
                        time.sleep(10);continue
                    if state['resources_admitted'] is True: break
                    v['predecessor'].progress(folder,dict(status='WAIT_RESOURCE_WINDOW_NO_NATIVE_RESTART',next_round=number,resource_state=state))
                    time.sleep(10)
                state.update(resource_guard_source=runner,resource_guard_policy='EXACT_ORIGINAL_STAR_COMMANDS_AND_BAM_QC_V1')
                resource=d.write_new(folder/f'RESOURCES_{number:04d}.json',state)
                v['predecessor'].progress(folder,dict(status='RUNNING_BOUND_ELEVEN_WORKER_ROUND',round=number,completed_remaining_rounds=completed))
                phase='EXECUTE_UNUSED_ROUND'
                previous,boundary,accepted=d.execute_round(number,jobs,folder/f'round_{number:04d}',h,previous,boundary,resource,v)
                completed.append(dict(number=number,whole=boundary,accepted=accepted))
            value=dict(status='REMAINING_SAME_EPOCH_SCOPE_ACTUAL_NATIVE_ZERO_NOT_FULL_SCIENCE_ACCEPTANCE',
                completed_remaining_rounds=completed,original_driver=h['actual_driver'],continuation_runner=runner,
                handoff=hp,prior_receipts=pins,controller_started=started,
                scientific_acceptance=False,old_controller_success_asserted=False,utc=d.utc())
            d.write_new(folder/'DISPATCH_EXIT.json',value);v['predecessor'].progress(folder,value);return 0
        except BaseException as exc:
            d.write_new(folder/'DISPATCH_EXIT.json',dict(status='STOP_RETAINED_NO_AUTOMATIC_NATIVE_RETRY',
                phase=phase,exception_class_only=type(exc).__name__,
                error_text_sha256=hashlib.sha256(str(exc).encode('utf-8')).hexdigest(),
                completed_remaining_rounds=completed,
                prior_receipts=pins,raw_error_persisted=False,scientific_acceptance=False,
                old_controller_success_asserted=False,utc=d.utc()))
            return 1


def observe(args):
    """Actual Popen parent waits without a terminal timeout; no native signals."""
    need(sys.platform.startswith('linux'), 'LinuxObserver')
    d=load_driver(args.driver)
    # Authenticate before consuming any new output namespace, then again under lock.
    authenticate(d,args)
    folder=ROOT/'provenance'/args.attempt
    folder.mkdir(exist_ok=False)
    argv=[sys.executable,'-I','-B',str(Path(__file__).resolve()),'_child',*sys.argv[1:]]
    p=subprocess.Popen(argv,cwd=ROOT,stdin=subprocess.DEVNULL,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
        start_new_session=True,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
    identity=None
    for _ in range(100):
        try: identity=d.proc(p.pid);break
        except (FileNotFoundError,ProcessLookupError):
            if p.poll() is not None:break
            time.sleep(.01)
    runner=d.pin(Path(__file__).resolve());original_driver=d.pin(args.driver);observer=d.proc(os.getpid())
    start=d.write_new(folder/'ACTUAL_PARENT_STARTED.json',dict(
        schema='jgw03-same-eleven-epoch-parent-started-v1',actual_pid=p.pid,
        start_ticks=identity['start_ticks'] if identity else None,argv=argv,
        runner=runner,original_driver=original_driver,
        observer=observer,actual_native_exit_not_yet_observed=True,utc=d.utc()))
    out,err=p.communicate()
    dispatch=folder/'DISPATCH_EXIT.json'
    d.write_new(folder/'ACTUAL_PARENT_NATIVE_EXIT.json',dict(
        schema='jgw03-same-eleven-epoch-parent-native-exit-v1',actual_pid=p.pid,
        start_ticks=identity['start_ticks'] if identity else None,argv=argv,started=start,
        runner=runner,original_driver=original_driver,observer=observer,
        actual_native_exit_code=p.returncode,observer_failure_class_only=None,
        stdout_bytes=len(out),stderr_bytes=len(err),stdout_sha256=hashlib.sha256(out).hexdigest(),
        stderr_sha256=hashlib.sha256(err).hexdigest(),raw_streams_persisted=False,
        dispatch_exit=d.pin(dispatch) if dispatch.exists() else None,
        scientific_acceptance=False,old_controller_success_asserted=False,utc=d.utc()))
    return p.returncode


def main():
    is_child=len(sys.argv)>1 and sys.argv[1]=='_child'
    if is_child:del sys.argv[1]
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--execute',action='store_true')
    for name in ('driver','dispatch','parent-started','parent-exit'):p.add_argument('--'+name,type=Path)
    p.add_argument('--attempt');p.add_argument('--next-round',type=int)
    args=p.parse_args()
    if not args.execute:
        print(json.dumps(dict(status='PLAN_ONLY_NO_FITS_OR_WRITES',same_scientific_epoch=True)));return 0
    need(all(getattr(args,k) is not None for k in ('driver','dispatch','parent_started','parent_exit','attempt','next_round')),
        'ExplicitFrozenPriorReceipts')
    need(re.fullmatch('[A-Za-z0-9_-]{1,80}',args.attempt) is not None,'FreshRunBoundAttemptName')
    return child(args) if is_child else observe(args)


if __name__=='__main__':
    try:raise SystemExit(main())
    except Exception as exc:
        print(json.dumps(dict(status='CHECKPOINT_CONTINUATION_STOP',exception_class_only=type(exc).__name__,
            error_text_sha256=hashlib.sha256(str(exc).encode('utf-8')).hexdigest(),raw_error_persisted=False)))
        raise SystemExit(2)
