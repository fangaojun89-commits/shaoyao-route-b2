"""Saved small-metadata gate for the exact four-epoch JGW03 execution union.

No model, archive decode, allocation, resource remeasurement or signal.
Reviewed driver validation is reused; a prepared handoff is not completion.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import types

DRIVER_SHA = '718d394a97b6052c427855169fd24dab68402ed96ef59cf5338eb2259151e22c'
COORDINATOR_SHA = '9ae73f398f899fa21d72bf8c41a59b902e9451e1abae5482a67f83242cff57b0'


def need(ok, label):
    if not ok:
        raise RuntimeError(label)


def load_driver(path):
    raw = Path(path).read_bytes()
    need(hashlib.sha256(raw).hexdigest() == DRIVER_SHA, 'FixedReviewedElevenWorkerDriver')
    value = types.ModuleType('eleven_worker_saved_gate_driver')
    value.__file__ = str(path)
    exec(compile(raw, str(path), 'exec'), value.__dict__)
    return value


def check_parent(d, started, terminal, controller, driver, handoff, dispatch, started_pin):
    need(started['schema'] == 'jgw03-successor-parent-started-v1' and
         terminal['schema'] == 'jgw03-successor-parent-native-exit-v1' and
         started['driver'] == terminal['driver'] == driver and
         started['handoff'] == terminal['handoff'] == handoff and
         started['actual_pid'] == terminal['actual_pid'] == controller['pid'] and
         started['start_ticks'] == terminal['start_ticks'] == controller['start_ticks'] and
         d.integer(controller['pid'], 1) and d.integer(controller['start_ticks'], 1) and
         started['argv'] == terminal['argv'] == [controller['executable'], '-I', '-B', *controller['argv']] and
         controller['argv'] == [driver['path'], '--execute', '--handoff', handoff['path'],
            '--handoff-sha256', handoff['sha256'], '--attempt', Path(dispatch['path']).parent.name] and
         terminal['started'] == started_pin and terminal['dispatch_exit'] == dispatch and
         type(terminal['actual_native_exit_code']) is int and terminal['actual_native_exit_code'] == 0 and
         terminal['observer_failure_class_only'] is None and terminal['raw_streams_persisted'] is False and
         terminal['old_controller_success_asserted'] is False and terminal['scientific_acceptance'] is False and
         started['actual_native_exit_not_yet_observed'] is True and
         started['coordinator_source'] == terminal['coordinator_source'] and
         d.bound(started['coordinator_source'])['sha256'] == COORDINATOR_SHA and
         started['coordinator_identity'] == terminal['coordinator_identity'], 'ActualElevenWorkerDirectParentZero')
    owner = started['coordinator_identity']
    need(all(d.integer(owner[k], 1) for k in ('pid', 'parent', 'birth', 'start_ticks')) and
         owner['birth'] == owner['start_ticks'] and owner['pid'] == controller['parent'],
         'SameActualCoordinatorPIDAndBirth')
    for stream in ('stdout', 'stderr'):
        need(d.integer(terminal[stream+'_bytes']) and
             isinstance(terminal[stream+'_sha256'], str) and
             re.fullmatch('[0-9a-f]{64}', terminal[stream+'_sha256']), 'ParentDigestOnlyStreams')


def process_count(legacy, predecessor_closed, eight_closed):
    need(type(predecessor_closed) is int and 67 <= predecessor_closed < eight_closed and
         type(eight_closed) is int and eight_closed <= 200, 'ExactFourEpochBoundaries')
    count = 3*predecessor_closed + 11*(eight_closed-predecessor_closed)
    for number, jobs in enumerate(legacy.plan(), 1):
        if number <= eight_closed: continue
        if number <= 201:
            n=jobs[0]['last']-jobs[0]['first']+1
            count += 3+min(4,n)+min(3,n)+min(4,n)
        else: count += len(jobs)
    need(count == 1679-8*(predecessor_closed-67)+3*(200-eight_closed), 'ExactFinitePlanProcessCount')
    return count


def metadata_closure(d, records):
    """Re-pin only referenced small JSON and source documents, never RDS payload."""
    tracked = {}; pending=list(records)
    def discover(value):
        if isinstance(value,dict):
            if set(value)=={'path','bytes','sha256'}:
                if Path(value['path']).suffix.lower() in ('.json','.py','.r','.md'): pending.append(value)
            else:
                for item in value.values(): discover(item)
        elif isinstance(value,list):
            for item in value: discover(item)
    while pending:
        item=pending.pop()
        if item['path'] in tracked:
            need(tracked[item['path']]==item,'MetadataReferenceIdentity');continue
        d.bound(item);tracked[item['path']]=item
        if Path(item['path']).suffix.lower()=='.json':discover(d.read(item['path']))
    return tracked


def verify(d, dispatch_path, parent_start_path, parent_exit_path):
    dispatch_pin = d.pin(dispatch_path); dispatch = d.read(dispatch_path); attempt = Path(dispatch_path).parent
    controller_pin = d.pin(attempt/'CONTROLLER_STARTED.json')
    controller = d.read(controller_pin['path'])
    handoff_pin = d.bound(controller['handoff']); h = d.read(handoff_pin['path'])
    need(h['actual_driver'] == d.pin(d.__file__) and h['actual_driver']['sha256'] == DRIVER_SHA and
         controller['actual_driver'] == h['actual_driver'] and controller['old_controller_success_asserted'] is False,
         'ActualReviewedControllerNotPreparedMaterial')
    tracked = metadata_closure(d, [dispatch_pin,controller_pin,handoff_pin,d.pin(parent_start_path),d.pin(parent_exit_path)])
    v = d.validate_handoff(h, pre_retirement=False)
    need(v['execution_authorized'] is True and v['old_handoff']['completed_rounds'] == 66,
         'AuthenticatedOld66AndClosedIncrementalEpoch')
    predecessor_terminal = d.read(d.bound(h['predecessor_controller']['terminal'])['path'])
    need(predecessor_terminal['coordinator_source']['sha256'] == COORDINATOR_SHA,
         'FixedReviewedRealRetirementObserver')
    d.bound(predecessor_terminal['coordinator_source'])
    ps_pin, pe_pin = d.pin(parent_start_path), d.pin(parent_exit_path)
    parent_start, parent_end = d.read(ps_pin['path']), d.read(pe_pin['path'])
    check_parent(d, parent_start, parent_end, controller, h['actual_driver'], handoff_pin, dispatch_pin, ps_pin)
    need(predecessor_terminal['coordinator_source'] == parent_start['coordinator_source'],
         'SameReviewedRetirementAndSuccessorObserver')
    remaining = v['predecessor'].remaining_plan(v['legacy'], h['completed_rounds'])
    plan_pin = d.pin(attempt/'PLAN.json')
    plan = d.read(plan_pin['path'])
    tracked.update(metadata_closure(d,[plan_pin]))
    need(plan == dict(remaining=[[number, jobs] for number, jobs in remaining], maximum_workers=11,
         per_worker_memory_bytes=d.MEMORY_PER_WORKER, ordered_merge_workers=3,
         exact_and_finish_workers=3, completed_before_this_epoch=h['completed_rounds'], scientific_acceptance=False),
         'ExactRemainingFrozenFullPlan')
    records = dispatch['completed_remaining_rounds']
    need(dispatch['status'] == 'REMAINING_SCOPE_ACTUAL_NATIVE_ZERO_NOT_FULL_SCIENCE_ACCEPTANCE' and
         dispatch['actual_driver'] == h['actual_driver'] and dispatch['handoff'] == handoff_pin and
         dispatch['completed_before_this_epoch'] == h['completed_rounds'] and
         dispatch['old_controller_success_asserted'] is False and dispatch['scientific_acceptance'] is False and
         len(records) == len(remaining), 'ActualEntireRemainingScopeNotPreparedOrPartial')
    previous = {cell: h['initial_previous'][cell]['prefix'] for cell in d.CELLS}
    boundary = h['completed_predecessor_rounds'][-1]['whole']
    for (number, _), item in zip(remaining, records):
        need(set(item) == {'number', 'whole', 'accepted'} and type(item['number']) is int and
             item['number'] == number, 'ExactRemainingRoundDomainAndOrder')
        previous, whole, accepted = d.validate_saved_round(attempt/f'round_{number:04d}', number, h,
            previous, boundary, v, controller['pid'])
        need(item['whole'] == whole and item['accepted'] == accepted, 'DispatchMatchesActualSavedRound')
        boundary = whole
    count = process_count(v['legacy'], v['eight_handoff']['completed_rounds'], h['completed_rounds'])
    # Re-pin tracked small metadata/sources. No allocation scientific RDS is decoded.
    for item in list(tracked.values()):
        need(d.pin(item['path']) == item, 'AllAuthenticatedMetadataUnchanged')
    return dict(status='PASS_EXACT_FOUR_EPOCH_EXECUTION_UNION_NOT_NUMERICAL_ACCEPTANCE',
        actual_rounds=204, legacy_rounds=66, successor_rounds=138,
        predecessor_closed_rounds=v['eight_handoff']['completed_rounds'],
        eight_worker_closed_rounds=h['completed_rounds'],
        incremental_rounds=v['eight_handoff']['completed_rounds']-66,
        eight_worker_rounds=h['completed_rounds']-v['eight_handoff']['completed_rounds'],
        eleven_worker_rounds=204-h['completed_rounds'], legacy_processes=198, successor_processes=count-198,
        legacy_entry=v['context']['legacy_R'], incremental_entry=h['new_R'], fit_entry=h['fit_R'],
        merge_entry=h['merge_R'], eight_fit_entry=v['eight_handoff']['fit_R'],
        eight_merge_entry=v['eight_handoff']['merge_R'], finish_entry=h['finish_R'], handoff=handoff_pin,
        pins=list(tracked.values()), scientific_acceptance=False, old_controller_success_asserted=False,
        resource_scope='Authenticated recorded admission, not current resource remeasurement',
        predecessor_actual_native_exit_code=predecessor_terminal['actual_native_exit_code'],
        predecessor_actual_parent_exit_code=-15)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('driver', 'dispatch', 'parent_started', 'parent_exit'):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    try:
        print(json.dumps(verify(load_driver(args.driver), args.dispatch, args.parent_started, args.parent_exit),
                         sort_keys=True, allow_nan=False))
    except Exception as error:
        raise SystemExit('SAVED_ELEVEN_WORKER_GATE_STOP_CLASS_ONLY:'+type(error).__name__)

