"""Saved-only four-science-epoch gate with two honest eleven-worker controllers.

Authenticates stopped original88..92/parent1 and resumed93..204/parent0.
No model, resource measurement or scientific archive decode. The frozen
continuation authenticator only checks the historical old PID is terminal.
"""
import argparse
import hashlib
from pathlib import Path
import re
import types

DRIVER_SHA='718d394a97b6052c427855169fd24dab68402ed96ef59cf5338eb2259151e22c'
RUNNER_SHA='cea59028f4e94d6b366d32c1006c6cdaf0a33b8fa29683e71fa4c35b91d4c998'
OLD_GATE_SHA='4e3a1bdc8a681540926350e7c9d1d5eaaa00e1f044d320f52219bc8f70ef9fd1'
FIRST_UNUSED=93
GUARD='EXACT_ORIGINAL_STAR_COMMANDS_AND_BAM_QC_V1'


def need(ok,label):
    if not ok:raise RuntimeError(label)


def load_fixed(path,sha,name):
    raw=Path(path).read_bytes()
    need(hashlib.sha256(raw).hexdigest()==sha,'FrozenSavedGateDependency')
    m=types.ModuleType(name);m.__file__=str(path)
    exec(compile(raw,str(path),'exec'),m.__dict__);return m


def exact_flags(values,expected):
    """Reject duplicate/unknown flags rather than silently accepting last value."""
    result={};i=0
    while i<len(values):
        key=values[i]
        need(key in expected and key not in result,'UniqueKnownContinuationFlag')
        if key=='--execute':result[key]=True;i+=1
        else:
            need(i+1<len(values),'ContinuationFlagValue')
            result[key]=values[i+1];i+=2
    need(result==expected,'ExactActualContinuationArguments')


def check_new_parent(d,controller,ps,pe,start_pin,dispatch_pin,runner,original_driver,prior):
    need(controller['schema']=='jgw03-same-eleven-epoch-checkpoint-controller-v1' and
         controller['continuation_runner']==runner and controller['actual_driver']==original_driver and
         controller['first_unused_round']==FIRST_UNUSED and type(controller['first_unused_round']) is int and
         controller['resource_guard_source']==runner and controller['resource_guard_policy']==GUARD and
         controller['prior_receipts']==prior and controller['old_controller_success_asserted'] is False and
         controller['scientific_acceptance'] is False,'HonestSameEpochNewController')
    need(ps['schema']=='jgw03-same-eleven-epoch-parent-started-v1' and
         pe['schema']=='jgw03-same-eleven-epoch-parent-native-exit-v1' and
         ps['actual_pid']==pe['actual_pid']==controller['pid'] and d.integer(controller['pid'],1) and
         ps['start_ticks']==pe['start_ticks']==controller['start_ticks'] and d.integer(controller['start_ticks'],1) and
         ps['runner']==pe['runner']==runner and ps['original_driver']==pe['original_driver']==original_driver and
         ps['observer']==pe['observer'] and
         ps['actual_native_exit_not_yet_observed'] is True and pe['started']==start_pin and
         pe['dispatch_exit']==dispatch_pin and type(pe['actual_native_exit_code']) is int and
         pe['actual_native_exit_code']==0 and pe['observer_failure_class_only'] is None and
         pe['raw_streams_persisted'] is False and pe['old_controller_success_asserted'] is False and
         pe['scientific_acceptance'] is False,'ActualNewControllerParentZero')
    owner=ps['observer']
    need(d.integer(owner['pid'],1) and d.integer(owner['start_ticks'],1) and
         owner['pid']==controller['parent'] and owner['start_ticks']<=controller['start_ticks'],
         'ActualContinuationDirectParentIdentity')
    parsed=controller['parsed_argv'];actual=controller['argv']
    need(isinstance(parsed,list) and parsed and parsed[0]==runner['path'] and
         actual==ps['argv']==pe['argv']==[controller['executable'],'-I','-B',runner['path'],'_child',*parsed[1:]],
         'ActualProcArgvNotPretendedOriginalDriver')
    exact_flags(parsed[1:],{'--execute':True,'--driver':original_driver['path'],
        '--dispatch':prior['dispatch']['path'],'--parent-started':prior['parent_started']['path'],
        '--parent-exit':prior['parent_exit']['path'],'--attempt':Path(dispatch_pin['path']).parent.name,
        '--next-round':str(FIRST_UNUSED)})
    for stream in ('stdout','stderr'):
        need(d.integer(pe[stream+'_bytes']) and isinstance(pe[stream+'_sha256'],str) and
             re.fullmatch('[a-f0-9]{64}',pe[stream+'_sha256']),'DigestOnlyActualParentStreams')


def exact_remaining_records(records,remaining):
    need(len(records)==len(remaining) and all(
        set(record)=={'number','whole','accepted'} and type(record['number']) is int and record['number']==n
        for record,(n,_) in zip(records,remaining)),'AllAndOnlyUnfinishedRoundsInOrder')


def verify(d,dispatch_path,parent_start_path,parent_exit_path):
    dispatch_pin=d.pin(dispatch_path);dispatch=d.read(dispatch_path);attempt=Path(dispatch_path).parent
    need(attempt.is_relative_to(d.ROOT/'provenance') and dispatch_path.name=='DISPATCH_EXIT.json',
         'RunBoundActualContinuationDispatch')
    controller_pin=d.pin(attempt/'CONTROLLER_STARTED.json');controller=d.read(controller_pin['path'])
    runner_pin=d.bound(controller['continuation_runner'],RUNNER_SHA)
    runner=load_fixed(runner_pin['path'],RUNNER_SHA,'saved_same_eleven_epoch_runner')
    need(dispatch['status']=='REMAINING_SAME_EPOCH_SCOPE_ACTUAL_NATIVE_ZERO_NOT_FULL_SCIENCE_ACCEPTANCE' and
         dispatch['original_driver']==controller['actual_driver']==d.pin(d.__file__) and
         dispatch['continuation_runner']==runner_pin and dispatch['controller_started']==controller_pin and
         dispatch['prior_receipts']==controller['prior_receipts'] and
         dispatch['old_controller_success_asserted'] is False and dispatch['scientific_acceptance'] is False,
         'ActualFullContinuationNotPartialOrSyntheticUnion')
    prior=dispatch['prior_receipts']
    need(set(prior)=={'dispatch','parent_started','parent_exit','controller_started'},'ExactPriorClosurePins')
    for p in prior.values():d.bound(p)
    args=argparse.Namespace(dispatch=Path(prior['dispatch']['path']),parent_started=Path(prior['parent_started']['path']),
        parent_exit=Path(prior['parent_exit']['path']),next_round=FIRST_UNUSED)
    h,hp,v,previous,boundary,remaining,authenticated=runner.authenticate(d,args)
    need(authenticated==prior and controller['handoff']==dispatch['handoff']==hp and
         controller['resource_epoch']==h['resource_epoch'] and controller['assigned_CPUs']==h['worker_cpu_ids'] and
         v['execution_authorized'] is True and v['old_handoff']['completed_rounds']==66,
         'SameFourScientificEpochsAndPriorClosure')
    original_records=d.read(args.dispatch)['completed_remaining_rounds']
    need([x['number'] for x in original_records]==list(range(h['completed_rounds']+1,FIRST_UNUSED)) and
         remaining[0][0]==FIRST_UNUSED and remaining[-1][0]==204,'ExactTwoActualElevenAttempts')
    ps_pin,pe_pin=d.pin(parent_start_path),d.pin(parent_exit_path)
    need(Path(ps_pin['path'])==attempt/'ACTUAL_PARENT_STARTED.json' and
         Path(pe_pin['path'])==attempt/'ACTUAL_PARENT_NATIVE_EXIT.json','ActualOwnedContinuationParentPaths')
    check_new_parent(d,controller,d.read(ps_pin['path']),d.read(pe_pin['path']),ps_pin,dispatch_pin,
        runner_pin,h['actual_driver'],prior)
    old_gate_path=Path(h['actual_driver']['path']).parent/'verify_JGW03_eleven_worker_execution_gate_20261022.py'
    old_gate_pin=d.pin(old_gate_path);old_gate=load_fixed(old_gate_path,OLD_GATE_SHA,'saved_four_science_epoch_helpers')
    plan_pin=d.pin(attempt/'PLAN.json');plan=d.read(plan_pin['path'])
    need(plan==dict(remaining=[[n,j] for n,j in remaining],first_unused_round=FIRST_UNUSED,
        same_scientific_epoch=True,maximum_workers=11,original_driver=h['actual_driver'],
        inherited_closed_rounds=original_records),'ExactContinuationPlanOnlyRemainingIDs')
    records=dispatch['completed_remaining_rounds'];exact_remaining_records(records,remaining);resource_pins=[]
    for (n,_),record in zip(remaining,records):
        resource_pin=d.pin(attempt/f'RESOURCES_{n:04d}.json');resource_pins.append(resource_pin)
        resource=d.read(resource_pin['path'])
        need(resource['resource_guard_source']==runner_pin and resource['resource_guard_policy']==GUARD and
             resource['direct_memory_gate_bytes']==41*d.GIB and resource['resources_admitted'] is True,
             'RecordedUnchangedAdmissionWithExplicitQCGuard')
        previous,whole,accepted=d.validate_saved_round(attempt/f'round_{n:04d}',n,h,previous,boundary,v,controller['pid'])
        need(record['whole']==whole and record['accepted']==accepted,'ActualNewPIDRoundClosure')
        boundary=whole
    tracked=old_gate.metadata_closure(d,[dispatch_pin,controller_pin,hp,ps_pin,pe_pin,runner_pin,old_gate_pin,
        plan_pin,d.pin(Path(__file__).resolve()),*prior.values(),*resource_pins])
    count=old_gate.process_count(v['legacy'],v['eight_handoff']['completed_rounds'],h['completed_rounds'])
    for p in tracked.values():d.bound(p)
    return dict(status='PASS_EXACT_FOUR_EPOCH_EXECUTION_UNION_NOT_NUMERICAL_ACCEPTANCE',
        actual_rounds=204,legacy_rounds=66,successor_rounds=138,
        predecessor_closed_rounds=v['eight_handoff']['completed_rounds'],eight_worker_closed_rounds=h['completed_rounds'],
        incremental_rounds=v['eight_handoff']['completed_rounds']-66,
        eight_worker_rounds=h['completed_rounds']-v['eight_handoff']['completed_rounds'],
        eleven_worker_rounds=204-h['completed_rounds'],legacy_processes=198,successor_processes=count-198,
        legacy_entry=v['context']['legacy_R'],incremental_entry=h['new_R'],fit_entry=h['fit_R'],merge_entry=h['merge_R'],
        eight_fit_entry=v['eight_handoff']['fit_R'],eight_merge_entry=v['eight_handoff']['merge_R'],finish_entry=h['finish_R'],
        handoff=hp,pins=list(tracked.values()),scientific_acceptance=False,old_controller_success_asserted=False,
        scientific_actor_epochs=4,eleven_controller_attempts=2,eleven_first_attempt_closed_through=FIRST_UNUSED-1,
        eleven_stopped_original_actual_parent_exit_code=1,eleven_continuation_actual_parent_exit_code=0,
        eleven_original_controller=d.read(prior['controller_started']['path']),eleven_continuation_controller=controller,
        resource_scope='Authenticated recorded admission, not current resource remeasurement',
        predecessor_actual_native_exit_code=None,predecessor_actual_parent_exit_code=-15)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('driver','dispatch','parent_started','parent_exit'):p.add_argument(key,type=Path)
    args=p.parse_args()
    try:
        import json
        print(json.dumps(verify(load_fixed(args.driver,DRIVER_SHA,'saved_original_eleven_driver'),
            args.dispatch,args.parent_started,args.parent_exit),sort_keys=True,allow_nan=False))
    except Exception as exc:
        raise SystemExit('SAVED_SAME_EPOCH_GATE_STOP_CLASS_ONLY:'+type(exc).__name__)
