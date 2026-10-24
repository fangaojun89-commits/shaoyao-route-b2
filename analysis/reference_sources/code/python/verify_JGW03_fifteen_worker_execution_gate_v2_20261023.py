"""Saved-only five-actor execution gate after an honest childless 11-to-15 handoff.

Original eleven attempt88..92 retains genuine parent1; resumed93..N retains
genuine parent-15. New fifteen owns N+1..204 and must have genuine parent0.
No fits, resource measurement, scientific RDS decoding or inferred zero exits.
"""
import argparse
import hashlib
from pathlib import Path
import re
import types

DRIVER_SHA='880b85ee7905d82b560c655b6dd0be6e2242ebffc2713e8b1996b69591988c0c'
COORDINATOR_SHA='b37b3a51bd82ca5854b6273013fbf9f4f04caeaf2b5912709d5f7724e2a09baa'
OLD_DRIVER_SHA='718d394a97b6052c427855169fd24dab68402ed96ef59cf5338eb2259151e22c'
RUNNER_SHA='cea59028f4e94d6b366d32c1006c6cdaf0a33b8fa29683e71fa4c35b91d4c998'
SAME_GATE_SHA='ae29f10783aaa54392e7f0aab532a11a20e491e267103d192d5d9152c6cfb9da'
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



def check_retired_parent(d,controller,ps,pe,start_pin,runner,original_driver,prior):
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
         pe['dispatch_exit'] is None and type(pe['actual_native_exit_code']) is int and
         pe['actual_native_exit_code']==-15 and pe['observer_failure_class_only'] is None and
         pe['raw_streams_persisted'] is False and pe['old_controller_success_asserted'] is False and
         pe['scientific_acceptance'] is False,'ActualRetiredContinuationParentMinus15')
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
        '--parent-exit':prior['parent_exit']['path'],'--attempt':Path(start_pin['path']).parent.name,
        '--next-round':str(FIRST_UNUSED)})
    for stream in ('stdout','stderr'):
        need(d.integer(pe[stream+'_bytes']) and isinstance(pe[stream+'_sha256'],str) and
             re.fullmatch('[a-f0-9]{64}',pe[stream+'_sha256']),'DigestOnlyActualParentStreams')


def exact_remaining_records(records,remaining):
    need(len(records)==len(remaining) and all(
        set(record)=={'number','whole','accepted'} and type(record['number']) is int and record['number']==n
        for record,(n,_) in zip(records,remaining)),'AllAndOnlyUnfinishedRoundsInOrder')


def coverage(legacy,predecessor_closed,eight_closed,eleven_closed):
    need(type(predecessor_closed) is int and 67<=predecessor_closed<eight_closed and
         type(eight_closed) is int and eight_closed<eleven_closed and
         type(eleven_closed) is int and 93<=eleven_closed<=200,'ExactFiveActorBoundaries')
    count=3*predecessor_closed+11*(eight_closed-predecessor_closed)
    for number,jobs in enumerate(legacy.plan(),1):
        if number<=eight_closed:continue
        if number<=201:
            n=jobs[0]['last']-jobs[0]['first']+1
            count+=3+sum(min(k,n) for k in ((4,3,4) if number<=eleven_closed else (6,4,5)))
        else:count+=len(jobs)
    need(count==3615-8*predecessor_closed-3*eight_closed-4*eleven_closed,
         'ExactFiniteFiveActorProcessCount')
    return count


def check_retirement(d,h,controller,record):
    n=h['completed_rounds']
    need(record['schema']=='jgw03-fifteen-worker-predecessor-terminal-v1' and
         record['status']=='ACTUAL_PROCESS_EXIT_OBSERVED_FOR_PLANNED_HANDOFF' and
         record['pid']==controller['pid'] and record['start_ticks']==controller['start_ticks'] and
         record['actual_native_exit_code'] is None and record['pidfd_exit_observed'] is True and
         record['exit_observed_after_signal'] is True and record['termination_signal_sent']==15 and
         type(record['termination_signal_sent']) is int and
         record['pending_TERM_delivered_via_SIGCONT'] is True and record['no_active_fit_children'] is True and
         record['driver']==h['old_driver'] and record['runner']==h['old_resume_runner'] and
         record['started']==h['old_resume_controller_started'] and
         record['completed_rounds']==n and type(record['completed_rounds']) is int and
         record['next_round']==n+1 and type(record['next_round']) is int and
         record['previous_boundary']==h['previous_boundary'] and record['initial_previous']==h['initial_previous'] and
         record['coordinator_source']==h['handoff_coordinator'] and
         d.bound(record['coordinator_source'])['sha256']==COORDINATOR_SHA,
         'ActualChildlessRetirementNotAZeroNativeExit')
    resource=d.read(d.bound(record['handoff_resource_admission'])['path'])
    need(resource['resources_admitted'] is True and resource['maximum_workers']==15 and
         type(resource['maximum_workers']) is int and resource['direct_memory_gate_bytes']==53*d.GIB,
         'ActualFifteenResourceAdmissionBeforeRetirement')

def check_fifteen_parent(d, started, terminal, controller, driver, handoff, dispatch, started_pin):
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
         started['coordinator_identity'] == terminal['coordinator_identity'], 'ActualFifteenWorkerDirectParentZero')
    owner = started['coordinator_identity']
    need(all(d.integer(owner[k], 1) for k in ('pid', 'parent', 'birth', 'start_ticks')) and
         owner['birth'] == owner['start_ticks'] and owner['pid'] == controller['parent'],
         'SameActualCoordinatorPIDAndBirth')
    for stream in ('stdout', 'stderr'):
        need(d.integer(terminal[stream+'_bytes']) and
             isinstance(terminal[stream+'_sha256'], str) and
             re.fullmatch('[0-9a-f]{64}', terminal[stream+'_sha256']), 'ParentDigestOnlyStreams')



MIRROR_ROOT=Path('/root/autodl-tmp/JGW03_full_reference_linux_20261017r03/windows_metadata_20261023_r01')
MIRROR_SHA='11b9a80c3892d04693369410758d8527379e0884eb08f367314f2b44bf6d3755'

def mirror_bindings(d,root,expected_sha):
    """Authenticate exact logical identities and one-to-one physical metadata."""
    from pathlib import PureWindowsPath
    mp=d.pin(root/'MANIFEST.json')
    need(mp['sha256']==expected_sha,'ExactMirrorManifestSHA')
    manifest=d.read(mp['path'])
    need(set(manifest)=={'schema','role','rows'} and
         manifest['schema']=='jgw03-exact-windows-metadata-mirror-v1' and
         manifest['role']=='historical_metadata_not_new_scientific_evidence',
         'ExactMirrorManifestSchema')
    need(isinstance(manifest['rows'],list) and len(manifest['rows'])==65,'ExactMirrorInventory')
    mapping={};physical_paths=set();physical_pins=[mp]
    for row in manifest['rows']:
        need(set(row)=={'logical','physical'},'ExactMirrorRow')
        logical=row['logical'];physical=row['physical']
        need(isinstance(logical,dict) and set(logical)=={'path','bytes','sha256'} and
             isinstance(logical['path'],str) and
             re.fullmatch('[CR]:[\\\\/].+',logical['path']) is not None and
             PureWindowsPath(logical['path']).is_absolute() and
             PureWindowsPath(logical['path']).suffix.lower() in ('.json','.py','.r','.md') and
             type(logical['bytes']) is int and 0<=logical['bytes']<=16*1024**2 and
             isinstance(logical['sha256'],str) and re.fullmatch('[a-f0-9]{64}',logical['sha256']),
             'ExactLogicalMetadataIdentity')
        need(isinstance(physical,dict) and set(physical)=={'path','bytes','sha256'},'ExactPhysicalRecord')
        p=Path(physical['path'])
        need(p.is_absolute() and p.parent==root/'files' and
             p.name==logical['sha256']+PureWindowsPath(logical['path']).suffix.lower() and
             physical['bytes']==logical['bytes'] and physical['sha256']==logical['sha256'],
             'ContainedContentAddressedMirror')
        need(logical['path'] not in mapping and str(p) not in physical_paths,'OneToOneMirror')
        d.bound(physical)
        mapping[logical['path']]=row;physical_paths.add(str(p));physical_pins.append(physical)
    need({str(p) for p in root.rglob('*') if p.is_file()}==set(physical_paths)|{mp['path']},
         'ExactMirrorMembersNoExtras')
    return mapping,physical_pins

def cross_platform_metadata_closure(d,records):
    """Same recursive metadata closure; foreign records require exact bindings."""
    mapping,mirror_pins=mirror_bindings(d,MIRROR_ROOT,MIRROR_SHA)
    logical_tracked={};physical_tracked={p['path']:p for p in mirror_pins}
    pending=list(records)
    def discover(value):
        if isinstance(value,dict):
            if set(value)=={'path','bytes','sha256'}:
                # PureWindowsPath also handles original backslash-form suffixes.
                from pathlib import PureWindowsPath
                if PureWindowsPath(value['path']).suffix.lower() in ('.json','.py','.r','.md'):
                    pending.append(value)
            else:
                for item in value.values():discover(item)
        elif isinstance(value,list):
            for item in value:discover(item)
    while pending:
        item=pending.pop();path=item['path']
        if path in logical_tracked:
            need(logical_tracked[path]==item,'MetadataReferenceIdentity');continue
        if Path(path).is_absolute():
            physical=d.bound(item)
        else:
            need(path in mapping,'UnknownForeignMetadata')
            row=mapping[path]
            need(item==row['logical'],'ExactForeignMetadataIdentity')
            physical=d.bound(row['physical'])
        logical_tracked[path]=item
        if physical['path'] in physical_tracked:
            need(physical_tracked[physical['path']]==physical,'PhysicalMetadataIdentity')
        physical_tracked[physical['path']]=physical
        if Path(physical['path']).suffix.lower()=='.json':discover(d.read(physical['path']))
    for pin in physical_tracked.values():d.bound(pin)
    return physical_tracked


def verify(d,dispatch_path,parent_start_path,parent_exit_path):
    dispatch_pin=d.pin(dispatch_path);dispatch=d.read(dispatch_path);attempt=Path(dispatch_path).parent
    need(attempt.is_relative_to(d.ROOT/'provenance') and dispatch_path.name=='DISPATCH_EXIT.json',
         'RunBoundActualFifteenDispatch')
    controller_pin=d.pin(attempt/'CONTROLLER_STARTED.json');controller=d.read(controller_pin['path'])
    hp=d.bound(controller['handoff']);h=d.read(hp['path'])
    need(h['schema']=='jgw03-fifteen-worker-shard-handoff-v1' and
         h['status']=='ACCEPTED_CLOSED_ELEVEN_BOUNDARY_STAR_COMPLETE' and
         h['actual_driver']==d.pin(d.__file__) and h['actual_driver']['sha256']==DRIVER_SHA and
         controller['schema']=='jgw03-fifteen-worker-controller-v1' and
         controller['actual_driver']==h['actual_driver'] and
         controller['assigned_CPUs']==h['worker_cpu_ids'] and controller['resource_epoch']==h['resource_epoch'] and
         controller['old_controller_success_asserted'] is False,'HonestNewFifthActorController')
    n=h['completed_rounds']
    need(type(n) is int and 93<=n<=200 and type(h['next_round']) is int and h['next_round']==n+1,
         'ExactClosedElevenAndFirstUnusedFifteenBoundary')
    # Authenticate unchanged original88..92 before the separately retired resume.
    d.bound(h['old_driver'],OLD_DRIVER_SHA);d.bound(h['old_resume_runner'],RUNNER_SHA)
    old=load_fixed(h['old_driver']['path'],OLD_DRIVER_SHA,'saved_original_eleven_driver')
    runner=load_fixed(h['old_resume_runner']['path'],RUNNER_SHA,'saved_original_resume_runner')
    raw=h['old_resume_args']
    need(set(raw)=={'driver','dispatch','parent_started','parent_exit','attempt','next_round'} and
         raw['driver']==h['old_driver']['path'] and type(raw['next_round']) is int and
         raw['next_round']==FIRST_UNUSED,'ExactOriginalResumeAuthenticationArguments')
    args=argparse.Namespace(**{k:(Path(v) if k in ('driver','dispatch','parent_started','parent_exit') else v)
                              for k,v in raw.items()})
    oldh,oldhp,oldv,previous,boundary,oldremaining,prior=runner.authenticate(old,args)
    oldattempt=old.ROOT/'provenance'/args.attempt
    need(d.bound(h['old_resume_controller_started'])==old.pin(oldattempt/'CONTROLLER_STARTED.json') and
         d.bound(h['old_resume_parent_started'])==old.pin(oldattempt/'ACTUAL_PARENT_STARTED.json') and
         d.bound(h['old_resume_parent_exit'])==old.pin(oldattempt/'ACTUAL_PARENT_NATIVE_EXIT.json') and
         not (oldattempt/'DISPATCH_EXIT.json').exists(),'ActualRetiredResumeOwnedPathsWithoutPretendedDispatch')
    oldctl=old.read(h['old_resume_controller_started']['path'])
    oldps,oldpe=old.read(h['old_resume_parent_started']['path']),old.read(h['old_resume_parent_exit']['path'])
    check_retired_parent(old,oldctl,oldps,oldpe,h['old_resume_parent_started'],
                         h['old_resume_runner'],h['old_driver'],prior)
    need(oldctl['handoff']==oldhp and oldctl['resource_epoch']==oldh['resource_epoch'] and
         oldctl['assigned_CPUs']==oldh['worker_cpu_ids'] and oldv['execution_authorized'] is True and
         oldv['old_handoff']['completed_rounds']==66,'UnchangedFourHistoricalActors')
    original_records=old.read(args.dispatch)['completed_remaining_rounds']
    oldplan_pin=old.pin(oldattempt/'PLAN.json');oldplan=old.read(oldplan_pin['path'])
    need(oldplan==dict(remaining=[[i,j] for i,j in oldremaining],first_unused_round=FIRST_UNUSED,
        same_scientific_epoch=True,maximum_workers=11,original_driver=oldh['actual_driver'],
        inherited_closed_rounds=original_records),'ExactOriginalResumePlanNotTruncated')
    oldresources=[]
    for number,_ in oldremaining:
        if number>n:break
        rp=old.pin(oldattempt/f'RESOURCES_{number:04d}.json');oldresources.append(rp);resource=old.read(rp['path'])
        need(resource['resource_guard_source']==h['old_resume_runner'] and
             resource['resource_guard_policy']==GUARD and resource['direct_memory_gate_bytes']==41*old.GIB and
             resource['resources_admitted'] is True,'HistoricalUnchangedElevenAdmissionAndQCGuard')
        previous,whole,accepted=old.validate_saved_round(oldattempt/f'round_{number:04d}',number,
            oldh,previous,boundary,oldv,oldctl['pid'])
        boundary=whole
    need(h['previous_boundary']==boundary and set(h['initial_previous'])==set(old.CELLS) and
         all(h['initial_previous'][c]['prefix']==previous[c] and h['initial_previous'][c]['through']==n*50
             for c in old.CELLS),'ExactNewBoundaryJoinsActualResumedPrefix')
    check_retirement(d,h,oldctl,d.read(d.bound(h['predecessor_retirement'])['path']))
    # The new driver separately validates the source-bound complete-STAR/resource
    # contract and its own handoff. It never produces a scientific draw here.
    v=d.validate_handoff(h,pre_retirement=False)
    need(dispatch['status']=='REMAINING_SCOPE_ACTUAL_NATIVE_ZERO_NOT_FULL_SCIENCE_ACCEPTANCE' and
         dispatch['actual_driver']==h['actual_driver'] and dispatch['handoff']==hp and
         dispatch['completed_before_this_epoch']==n and dispatch['old_controller_success_asserted'] is False and
         dispatch['scientific_acceptance'] is False,'ActualEntireFifteenScopeNotPreparedOrPartial')
    ps_pin,pe_pin=d.pin(parent_start_path),d.pin(parent_exit_path)
    check_fifteen_parent(d,d.read(ps_pin['path']),d.read(pe_pin['path']),controller,h['actual_driver'],
                        hp,dispatch_pin,ps_pin)
    need(d.read(ps_pin['path'])['coordinator_source']==h['handoff_coordinator'],
         'SameReviewedRetirementAndFifteenActualParent')
    remaining=oldv['predecessor'].remaining_plan(oldv['legacy'],n)
    plan_pin=d.pin(attempt/'PLAN.json');plan=d.read(plan_pin['path'])
    need(plan==dict(remaining=[[number,jobs] for number,jobs in remaining],maximum_workers=15,
         per_worker_memory_bytes=d.MEMORY_PER_WORKER,ordered_merge_workers=3,
         exact_and_finish_workers=3,completed_before_this_epoch=n,scientific_acceptance=False),
         'ExactUnusedFifteenPlan')
    records=dispatch['completed_remaining_rounds'];exact_remaining_records(records,remaining);resources=[]
    for (number,_),record in zip(remaining,records):
        rp=d.pin(attempt/f'RESOURCES_{number:04d}.json');resources.append(rp);resource=d.read(rp['path'])
        need(resource['resources_admitted'] is True and resource['maximum_workers']==15 and
             type(resource['maximum_workers']) is int and resource['direct_memory_gate_bytes']==53*d.GIB and
             resource['global_CPU_floor']==5,'RecordedFifteenAdmissionNotLowered')
        previous,whole,accepted=d.validate_saved_round(attempt/f'round_{number:04d}',number,h,
            previous,boundary,v,controller['pid'])
        need(record['whole']==whole and record['accepted']==accepted,'ActualFifteenPIDRoundClosure')
        boundary=whole
    old_gate_path=Path(oldh['actual_driver']['path']).parent/'verify_JGW03_eleven_worker_execution_gate_20261022.py'
    old_gate_pin=d.pin(old_gate_path)
    old_gate=load_fixed(old_gate_path,OLD_GATE_SHA,'saved_historical_metadata_closure')
    same_gate_path=Path(__file__).with_name('verify_JGW03_same_eleven_epoch_execution_gate_20261022.py')
    same_gate_pin=d.pin(same_gate_path)
    need(same_gate_pin['sha256']==SAME_GATE_SHA,'PreservedOriginalV7Source')
    tracked=cross_platform_metadata_closure(d,[dispatch_pin,controller_pin,hp,ps_pin,pe_pin,plan_pin,
        d.pin(Path(__file__).resolve()),old_gate_pin,same_gate_pin,oldplan_pin,
        *prior.values(),*oldresources,*resources,h['old_driver'],h['old_resume_runner'],
        h['old_resume_controller_started'],h['old_resume_parent_started'],h['old_resume_parent_exit']])
    p=oldv['eight_handoff']['completed_rounds'];e=oldh['completed_rounds']
    count=coverage(oldv['legacy'],p,e,n)
    for item in tracked.values():d.bound(item)
    return dict(status='PASS_EXACT_FIVE_ACTOR_EXECUTION_UNION_NOT_NUMERICAL_ACCEPTANCE',
        actual_rounds=204,legacy_rounds=66,successor_rounds=138,
        predecessor_closed_rounds=p,eight_worker_closed_rounds=e,eleven_worker_closed_rounds=n,
        incremental_rounds=p-66,eight_worker_rounds=e-p,eleven_worker_rounds=n-e,fifteen_worker_rounds=204-n,
        legacy_processes=198,successor_processes=count-198,legacy_entry=oldv['context']['legacy_R'],
        incremental_entry=oldh['new_R'],fit_entry=oldh['fit_R'],merge_entry=oldh['merge_R'],
        eight_fit_entry=oldv['eight_handoff']['fit_R'],eight_merge_entry=oldv['eight_handoff']['merge_R'],
        fifteen_fit_entry=h['fit_R'],fifteen_merge_entry=h['merge_R'],finish_entry=h['finish_R'],
        handoff=hp,pins=list(tracked.values()),scientific_acceptance=False,old_controller_success_asserted=False,
        scientific_actor_epochs=5,eleven_controller_attempts=2,eleven_first_attempt_closed_through=92,
        eleven_stopped_original_actual_parent_exit_code=1,eleven_continuation_actual_parent_exit_code=-15,
        fifteen_actual_parent_exit_code=0,eleven_original_controller=old.read(prior['controller_started']['path']),
        eleven_continuation_controller=oldctl,fifteen_controller=controller,
        historical_metadata_mirror_manifest=d.pin(MIRROR_ROOT/'MANIFEST.json'),
        resource_scope='Authenticated recorded admission, not current resource remeasurement',
        predecessor_actual_native_exit_code=None,predecessor_actual_parent_exit_code=-15)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('driver','dispatch','parent_started','parent_exit'):p.add_argument(key,type=Path)
    args=p.parse_args()
    try:
        import json
        print(json.dumps(verify(load_fixed(args.driver,DRIVER_SHA,'saved_fifteen_driver'),
            args.dispatch,args.parent_started,args.parent_exit),sort_keys=True,allow_nan=False))
    except Exception as exc:
        raise SystemExit('SAVED_FIVE_ACTOR_GATE_STOP_CLASS_ONLY:'+type(exc).__name__)





