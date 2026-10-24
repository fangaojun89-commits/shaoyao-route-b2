"""Fifteen-worker unused-label successor, plan-only by default.
Five disjoint shards per cell; unchanged scientific callbacks and saved ordered merge.
The previous controller may be retired only by the separately reviewed coordinator.
This executor contains no signals, retries, cache advice or completed scientific fit.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime as dt
import hashlib
import os
import json
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import threading
import time
import types

ROOT=Path('/root/autodl-tmp/JGW03_full_reference_linux_20261017r03')
STAR=Path('/root/autodl-tmp/GSE205906_matched_STAR_Ensembl101_20261015r01')
CELLS=('GSE289929','GSE227494','GSE299063')
GIB=1024**3
OLD_DRIVER_SHA='718d394a97b6052c427855169fd24dab68402ed96ef59cf5338eb2259151e22c'
OLD_RESUME_SHA='cea59028f4e94d6b366d32c1006c6cdaf0a33b8fa29683e71fa4c35b91d4c998'
FIT_SHA='7eabb266f29fb73fe172f08e090de8e1d9c1984b7ec9dec34104d0fa84d51450'
MERGE_SHA='92945fe0dba456f6d2594e138a37115e7cfa8925738e1d640a5dde4683a76347'
CONTEXT_SHA='452acfa10dbbaa29631ad41877feab120ba421dec6852bd42ab5180ef9b1b56c'
FINISH_SHA='69ba3dd31a7854462a22e03c0121d825ea8b228a0278d35d5862f1ac8a7d5613'
HELPER_SHA='74bb8dc0d3cdb16bb9373fdf5a9e6153af7298d5c104276e58c3247fc41b972e'
STAR_SHA='7fe64fcfcee554e1500f60da2cd7b4a3ce8a84fcd4e46b84c3c7e14ec1fd907d'
QUEUE_SHA='8c4e33bd9e7ecb3489b121540a77ab66b4f27392e2b46f5ada6bd1ea97d3c5cf'
OLD_PID,OLD_BIRTH=57056,2616685761
MEMORY_PER_WORKER=3*GIB
AFFINITY_EXEC='import os,sys;os.sched_setaffinity(0,{int(sys.argv[1])});os.execv(sys.argv[2],sys.argv[2:])'
def need(ok, label):
    if not ok: raise RuntimeError(label)

def integer(value, minimum=0):
    return type(value) is int and value >= minimum

def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def pin(path):
    p = Path(path)
    need(p.is_absolute() and not any(q.is_symlink() for q in (p, *p.parents)), 'AbsoluteNonsymlinkPin')
    s = p.stat(); need(stat.S_ISREG(s.st_mode), 'RegularPin')
    digest = hashlib.sha256()
    with p.open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''): digest.update(block)
    t = p.stat(); fields = ('st_dev', 'st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns')
    need(all(getattr(s, x) == getattr(t, x) for x in fields), 'PinStable')
    return dict(path=str(p), bytes=t.st_size, sha256=digest.hexdigest())

def bound(record, sha=None):
    need(isinstance(record, dict) and set(record) == {'path', 'bytes', 'sha256'} and
         integer(record['bytes']) and isinstance(record['sha256'], str) and
         re.fullmatch('[a-f0-9]{64}', record['sha256']), 'FullPin')
    need(sha is None or record['sha256'] == sha, 'FixedSourceSHA')
    need(pin(record['path']) == record, 'ExactCurrentPin'); return record

def read(path):
    need(Path(path).stat().st_size <= 16*1024**2, 'SmallJSONOnly')
    def unique(pairs):
        result = {}
        for k, v in pairs: need(k not in result, 'DuplicateJSONKey'); result[k] = v
        return result
    return json.loads(Path(path).read_text(encoding='utf-8'), object_pairs_hook=unique,
        parse_constant=lambda _: (_ for _ in ()).throw(ValueError('NonfiniteJSON')))

def write_new(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n'); stream.flush(); os.fsync(stream.fileno())
    return pin(path)

def load_module(record, name):
    bound(record); raw = Path(record['path']).read_bytes()
    need(len(raw) == record['bytes'] and hashlib.sha256(raw).hexdigest() == record['sha256'], 'ExecutedBufferPin')
    m = types.ModuleType(name); m.__file__ = record['path']
    exec(compile(raw, record['path'], 'exec'), m.__dict__); return m

def proc(pid):
    directory = Path('/proc') / str(pid)
    f = (directory/'stat').read_text().rsplit(') ', 1)[1].split()
    return dict(pid=pid, parent=int(f[1]), start_ticks=int(f[19]), state=f[0],
        argv=[os.fsdecode(x) for x in (directory/'cmdline').read_bytes().rstrip(b'\0').split(b'\0')])

def partitions(cell, first, last):
    need(cell in CELLS and integer(first, 3) and integer(last, first) and last <= 10001 and
         last-first+1 <= 50, 'OneUnusedBoundedRound')
    n = last-first+1; count = min(5, n)
    base, extra = divmod(n, count); result = []; start = first
    for i in range(count):
        end = start+base+(i < extra)-1; result.append((start,end)); start = end+1
    return result

def fit_jobs(jobs):
    result = []
    for original in jobs:
        need(original['mode'] == 'run' and original['space'] == 'four_group_MC', 'MCOnlyFitShards')
        for i, (first, last) in enumerate(partitions(original['accession'], original['first'], original['last']), 1):
            result.append(dict(mode='fit_only', accession=original['accession'], space='four_group_MC',
                first=first, last=last, round_first=original['first'], round_last=original['last'], shard_index=i))
    need(len(result) <= 15, 'MaximumFifteenFitProcesses'); return result


def validate_STAR_closure(h):
    c=read(bound(h['STAR_closure'])['path'])
    need(c['schema']=='jgw03-completed-STAR-queue-closure-v1','ExplicitSTARClosure')
    original=load_module(bound(c['source'],STAR_SHA),'jgw15_original_STAR')
    queue=load_module(bound(c['queue_source'],QUEUE_SHA),'jgw15_original_queue')
    need(c['source']['path']==str(STAR/'source_staging_20261015_01/run_gse205906_matched_STAR_20261015.py') and
         c['queue_source']['path']==str(STAR/'remaining_serial_source_20261022_r01/code/python/continue_GSE205906_STAR_remaining_serial_20261022.py'),
         'OriginalSTARAndQueueSourcePaths')
    expected=STAR/'provenance/remaining_serial_queue_20261022_r01'
    parent=STAR/'provenance/remaining_serial_queue_parent_20261022_r01'
    for field,path in [('queue_started',expected/'QUEUE_STARTED.json'),('all_nine',expected/'ALL_NINE_LIBRARIES_COMPLETE.json'),
                       ('queue_parent_started',parent/'STARTED.json'),('queue_parent_exit',parent/'NATIVE_EXIT.json')]:
        need(bound(c[field])['path']==str(path),'FixedActualSTARReceiptPath')
    qs,ps,pe=read(c['queue_started']['path']),read(c['queue_parent_started']['path']),read(c['queue_parent_exit']['path'])
    queue.native_zero(ps,pe)
    need(ps['actual_pid']==31592 and ps['identity']['start_ticks']==2614845177 and
         all(qs['identity'][k]==ps['identity'][k] for k in ('pid','start_ticks','argv')) and
         qs['source']==c['queue_source'] and ps['argv']==['/usr/bin/python3','-I','-B',c['queue_source']['path'],
             '--execute','--attempt','remaining_serial_queue_20261022_r01'],'FixedActualSTARQueueIdentity')
    allnine=read(c['all_nine']['path']); pins=[]
    for lib in original.metadata()['libraries']:
        folder=STAR/'libraries'/lib['SRA_run']; d=original.completed(folder,STAR_SHA)
        need(d['technical_reservoir_only'] is False and d['input_pairs']==lib['registered_pairs'] and
             len(d['events'])==8 and all(type(e['actual_exit']) is int and e['actual_exit']==0 for e in d['events']),
             'ActualNineFullNativeLibraries')
        pins.append(pin(folder/'COMPLETE.json'))
    need(c['library_closures']==allnine['checkpoints']==pins and
         allnine['scope']=='ORIGINAL_LIBRARY_CHECKPOINTS_ONLY','ExactNineSTARClosures')
    try: live=proc(31592)
    except (FileNotFoundError,ProcessLookupError): live=None
    need(live is None or live['start_ticks']!=2614845177 or live['state'] in ('Z','X'),'STARQueueActuallyGone')
    return c

def old_history(h,allow_active_future=False):
    old=load_module(bound(h['old_driver'],OLD_DRIVER_SHA),'jgw15_saved_eleven')
    resume=load_module(bound(h['old_resume_runner'],OLD_RESUME_SHA),'jgw15_saved_resume93')
    raw=dict(h['old_resume_args'])
    need(set(raw)=={'driver','dispatch','parent_started','parent_exit','attempt','next_round'} and
         raw['next_round']==93 and raw['attempt']=='eleven_checkpoint_resume93_01','ExactOriginalResumeArguments')
    for k in ('driver','dispatch','parent_started','parent_exit'): raw[k]=Path(raw[k])
    args=argparse.Namespace(**raw)
    need(pin(args.driver)==h['old_driver'],'OriginalDriverArgument')
    oh,ohp,ov,previous,boundary,remaining,pins=resume.authenticate(old,args)
    attempt=ROOT/'provenance'/raw['attempt']
    need(h['old_resumed_attempt']==str(attempt) and
         bound(h['old_resume_controller_started'])['path']==str(attempt/'CONTROLLER_STARTED.json') and
         bound(h['old_resume_parent_started'])['path']==str(attempt/'ACTUAL_PARENT_STARTED.json'),'OriginalResumeReceiptPaths')
    controller=read(h['old_resume_controller_started']['path'])
    need(controller['pid']==OLD_PID and controller['start_ticks']==OLD_BIRTH and
         controller['continuation_runner']==h['old_resume_runner'] and controller['actual_driver']==h['old_driver'] and
         controller['handoff']==ohp and controller['first_unused_round']==93 and controller['prior_receipts']==pins,
         'ActualSameElevenEpochController')
    history=h['completed_predecessor_rounds']
    need(len(history)==h['completed_rounds']-92,'AllResumedClosedRoundsPresent')
    for n,item in enumerate(history,93):
        need(set(item)=={'number','whole','accepted'} and type(item['number']) is int and item['number']==n,
             'OrderedOriginalResumedRoundHistory')
        previous,whole,accepted=old.validate_saved_round(attempt/f'round_{n:04d}',n,oh,previous,boundary,ov,OLD_PID)
        need(item['whole']==whole and item['accepted']==accepted,'ExactOriginalResumedRoundPins')
        boundary=whole
    if not allow_active_future:
        for n,_ in remaining:
            if n>h['completed_rounds']: need(not (attempt/f'round_{n:04d}').exists(),'NoOldPartialFutureRound')
    need(h['previous_boundary']==boundary,'ActualLatestOldRoundBoundary')
    for cell in CELLS:
        p=h['initial_previous'][cell]
        need(p['prefix']==previous[cell] and p['through']==h['completed_rounds']*50 and
             p['scientific_identity_sha256']==oh['initial_previous'][cell]['scientific_identity_sha256'],
             'ExactCurrentPrefixAndPreservedScientificIdentity')
        bound(p['prefix'])
    return old,resume,oh,ov,previous,boundary,controller

def validate_old_parent(h,controller,pre_retirement):
    ps=read(bound(h['old_resume_parent_started'])['path'])
    need(ps['schema']=='jgw03-same-eleven-epoch-parent-started-v1' and
         ps['actual_pid']==controller['pid']==OLD_PID and ps['start_ticks']==controller['start_ticks']==OLD_BIRTH and
         ps['runner']==h['old_resume_runner'] and ps['original_driver']==h['old_driver'] and
         ps['actual_native_exit_not_yet_observed'] is True and ps['argv']==controller['argv'] and
         ps['observer']['pid']==controller['parent'],'ExactOriginalResumeActualParent')
    if pre_retirement:
        need('old_resume_parent_exit' not in h and 'predecessor_retirement' not in h,'NoInventedFutureParentExit')
        live=proc(OLD_PID)
        need(live['start_ticks']==OLD_BIRTH and live['argv']==controller['argv'] and live['parent']==controller['parent'],
             'SameLiveControllerBeforeRetirement')
        return ps
    need(bound(h['old_resume_parent_exit'])['path']==str(Path(h['old_resumed_attempt'])/'ACTUAL_PARENT_NATIVE_EXIT.json'),
         'FixedRealResumeNativeExitPath')
    pe=read(h['old_resume_parent_exit']['path'])
    need(pe['schema']=='jgw03-same-eleven-epoch-parent-native-exit-v1' and
         type(pe['actual_native_exit_code']) is int and pe['actual_native_exit_code']==-15 and
         pe['actual_pid']==ps['actual_pid'] and pe['start_ticks']==ps['start_ticks'] and pe['argv']==ps['argv'] and
         pe['started']==h['old_resume_parent_started'] and pe['runner']==ps['runner'] and
         pe['original_driver']==ps['original_driver'] and pe['observer']==ps['observer'] and
         pe['observer_failure_class_only'] is None and pe['raw_streams_persisted'] is False and
         pe['dispatch_exit'] is None and pe['scientific_acceptance'] is False and
         pe['old_controller_success_asserted'] is False,'ActualParentMinus15NotInventedSuccess')
    for stream in ('stdout','stderr'):
        need(integer(pe[stream+'_bytes']) and re.fullmatch('[0-9a-f]{64}',pe[stream+'_sha256']),'DigestOnlyActualParentStreams')
    t=read(bound(h['predecessor_retirement'])['path'])
    need(t['schema']=='jgw03-fifteen-worker-predecessor-terminal-v1' and
         t['status']=='ACTUAL_PROCESS_EXIT_OBSERVED_FOR_PLANNED_HANDOFF' and t['pid']==OLD_PID and
         t['start_ticks']==OLD_BIRTH and t['actual_native_exit_code'] is None and
         t['pidfd_exit_observed'] is True and t['exit_observed_after_signal'] is True and
         t['termination_signal_sent']==15 and t['pending_TERM_delivered_via_SIGCONT'] is True and
         t['no_active_fit_children'] is True and t['driver']==h['old_driver'] and
         t['runner']==h['old_resume_runner'] and t['started']==h['old_resume_controller_started'] and
         t['completed_rounds']==h['completed_rounds'] and t['next_round']==h['next_round'] and
         t['previous_boundary']==h['previous_boundary'] and t['initial_previous']==h['initial_previous'] and
         t['coordinator_source']==h['handoff_coordinator'],'GuardedClosedChildlessRetirementEvidence')
    r=read(bound(t['handoff_resource_admission'])['path'])
    need(r['resources_admitted'] is True and r['policy']=='POST_STAR_FIFTEEN_WORKER_EXACT_ACTOR_EPOCH_V1' and
         r['actual_driver']==h['actual_driver'] and r['direct_memory_gate_bytes']==53*GIB,'ActualRetirementResourceGate')
    try: live=proc(OLD_PID)
    except (FileNotFoundError,ProcessLookupError): live=None
    need(live is None or live['start_ticks']!=OLD_BIRTH or live['state'] in ('Z','X'),'OldControllerActuallyGone')
    return ps

def validate_handoff(h,pre_retirement=False,allow_active_future=False):
    need(type(allow_active_future) is bool and (not allow_active_future or pre_retirement is True),
         'ActiveFutureAllowanceOnlyForNonexecutingPreflight')
    need(type(pre_retirement) is bool and h['schema']=='jgw03-fifteen-worker-shard-handoff-v1' and
         h['status']==('MATERIALS_PREPARED_NOT_HANDOFF_NOT_EXECUTION_ACCEPTANCE' if pre_retirement else
                       'ACCEPTED_CLOSED_ELEVEN_BOUNDARY_STAR_COMPLETE') and h['run_root']==str(ROOT) and
         integer(h['completed_rounds'],93) and h['completed_rounds']<=200 and
         type(h['next_round']) is int and h['next_round']==h['completed_rounds']+1,'ExplicitFifteenBoundary')
    need(bound(h['actual_driver'])==pin(Path(__file__).resolve()),'ActualFifteenSource')
    for k,sha in [('fit_R',FIT_SHA),('merge_R',MERGE_SHA),('new_R',CONTEXT_SHA),('finish_R',FINISH_SHA),('increment_helper',HELPER_SHA)]:
        bound(h[k],sha); need(Path(h[k]['path']).is_relative_to(ROOT),'RunBoundFifteenActors')
    old,resume,oh,ov,previous,boundary,controller=old_history(h,allow_active_future=allow_active_future)
    for k in ('new_R','finish_R','increment_helper','config','runtime','legacy_acceptance','keys_sidecar','legacy_R'):
        need(h[k]==ov['context'][k],'UnchangedOriginalScientificContext')
    for k in ('review','resource_epoch','handoff_coordinator'): bound(h[k])
    review=read(h['review']['path'])
    need(review['schema']=='jgw03-fifteen-worker-source-review-v1' and review['independent_review_accepted'] is True and
         all(h[k] in review['sources'] for k in ('actual_driver','fit_R','merge_R','handoff_coordinator')),
         'ExplicitIndependentSourceReview')
    epoch=read(h['resource_epoch']['path']);cpus=h['worker_cpu_ids']
    need(epoch['schema']=='jgw03-post-STAR-fifteen-worker-resource-v1' and
         epoch['instance']=='62c8pzy536-50fb405a' and epoch['review']==h['review'] and
         epoch['worker_cpu_ids']==cpus and len(cpus)==len(set(cpus))==15 and all(integer(x) for x in cpus) and
         epoch['direct_memory_gate_bytes']==53*GIB and epoch['global_CPU_floor']==5,'ExplicitFifteenResourceEpoch')
    validate_STAR_closure(h); validate_old_parent(h,controller,pre_retirement)
    context=dict(ov['context'],actual_driver=h['actual_driver'],fit_R=h['fit_R'],merge_R=h['merge_R'],
                 review=h['review'],completed_rounds=h['completed_rounds'],initial_previous=h['initial_previous'])
    return dict(scope='FIVE_ACTOR_EPOCH_BOUNDARY_NOT_FULL_SCIENCE_ACCEPTANCE',execution_authorized=not pre_retirement,
                predecessor=ov['predecessor'],legacy=ov['legacy'],context=context,keys=ov['keys'],epoch=epoch,
                epoch_validator=ov['epoch_validator'],old_handoff=oh,old_validated=ov,old_driver=old,old_resume=resume)

def build_boundary_material(template,next_round):
    need(integer(next_round,94) and next_round<=201,'BoundedNextUnusedFifteenRound')
    h=json.loads(json.dumps(template));h['completed_rounds']=next_round-1;h['next_round']=next_round
    h['status']='MATERIALS_PREPARED_NOT_HANDOFF_NOT_EXECUTION_ACCEPTANCE'
    h.pop('old_resume_parent_exit',None);h.pop('predecessor_retirement',None)
    attempt=Path(h['old_resumed_attempt']);h['completed_predecessor_rounds']=[]
    for n in range(93,next_round):
        folder=attempt/f'round_{n:04d}'
        h['completed_predecessor_rounds'].append(dict(number=n,whole=pin(folder/'merge/ROUND_EXIT.json'),
                                                      accepted=pin(folder/'ROUND_ACCEPTED.json')))
    h['previous_boundary']=h['completed_predecessor_rounds'][-1]['whole']
    for cell in CELLS:
        h['initial_previous'][cell]['through']=50*(next_round-1)
        h['initial_previous'][cell]['prefix']=pin(ROOT/'derived'/cell/f'four_group_MC_TAIL_PREFIX_{50*(next_round-1):05d}.rds')
    return h

def actor_for(h,validated,cell,key):
    item=validated['keys']['cells'][cell];mc=item['four_group_MC']
    need(key in mc or key in item['high_companion_exact'],'KnownFrozenKeyOnly')
    if key not in mc: return validated['context']['legacy_R']
    i=mc.index(key)+1
    if i<=h['completed_rounds']*50:
        return validated['old_driver'].actor_for(validated['old_handoff'],validated['old_validated'],cell,key)
    return h['fit_R']
def current_credit(h, validated):
    """Same logical-byte credit rules, with exact four-epoch actor ownership.

    No scientific payload is opened or hashed. Pending/exclusive new leaves get
    no credit. Old overhead/25% variation remain in the frozen calculate().
    """
    total = 0; inodes = set()
    for cell in CELLS:
        root = ROOT/'derived'/cell; identity = h['initial_previous'][cell]['scientific_identity_sha256']
        prepared = (ROOT/'intermediate'/cell/'PREPARED.json').stat().st_mtime
        for leaf in root.iterdir():
            if not leaf.is_dir(): continue
            key = leaf.name; expected = actor_for(h, validated, cell, key)
            need(not leaf.is_symlink(), 'NoSymlinkArchiveLeaf')
            if not (leaf/'COMPLETE.json').exists(): continue  # no credit, never remove/retry
            start, done = read(leaf/'STARTED.json'), read(leaf/'COMPLETE.json')
            need(start['status'] == 'ALLOCATION_OWNED_BEFORE_NATIVE_CALL' and
                 done['status'] == 'IMMUTABLE_COMPLETE_COMPACT_CHECKPOINT' and
                 start['allocation_key'] == done['allocation_key'] == key and
                 start['scientific_identity_sha256'] == done['scientific_identity_sha256'] == identity and
                 start['execution_entry'] == expected and type(done['new_native_calls']) is int and
                 done['new_native_calls'] in (0, 1) and type(start['new_fit_intended']) is bool and
                 int(start['new_fit_intended']) == done['new_native_calls'], 'ExactClosedArchiveCreditIdentity')
            path = leaf/'COMPACT.rds'; s = path.lstat()
            need(stat.S_ISREG(s.st_mode) and s.st_nlink == 1 and s.st_mtime >= prepared and
                 s.st_size == done['archive']['bytes'] and type(done['archive']['bytes']) is int and
                 re.fullmatch('[0-9a-f]{64}', done['archive']['sha256']) and
                 (s.st_dev, s.st_ino) not in inodes, 'StableSingleLinkLogicalArchiveCredit')
            inodes.add((s.st_dev, s.st_ino)); total += s.st_size
    return total


def resource_arithmetic(limit,current,free,cpu_quota,cpus,credit,compact_max):
    need(all(integer(x) for x in (limit,current,free,cpu_quota,credit)) and limit==96*GIB and current<=limit and
         cpu_quota==20 and len(cpus)==len(set(cpus))==15 and all(integer(c) for c in cpus),
         'MeasuredPostSTARFifteenDomain')
    base=sum(compact_max)*10020
    compact=(base*5+3)//4
    remaining=max(0,compact-credit)+40*GIB
    required=remaining+8*GIB+50*GIB
    categories=[]
    if limit-current<53*GIB: categories.append('WAIT_FIFTEEN_WORKER_DIRECT_MEMORY_HEADROOM')
    if free<required: categories.append('WAIT_POST_STAR_DISK_SCENARIO')
    return dict(policy='POST_STAR_FIFTEEN_WORKER_EXACT_ACTOR_EPOCH_V1',resources_admitted=not categories,
        wait_categories=categories,maximum_workers=15,per_worker_memory_bytes=3*GIB,worker_memory_bytes=45*GIB,
        memory_floor_bytes=8*GIB,direct_memory_gate_bytes=53*GIB,cgroup_limit_bytes=limit,
        cgroup_current_bytes=current,direct_headroom_bytes=limit-current,cgroup_vCPU_quota=cpu_quota,
        assigned_CPU_ids=cpus,global_CPU_floor=5,remaining_JGW_disk_scenario_bytes=remaining,
        credited_materialized_compact_logical_bytes=credit,original_finite_disk_allowance_bytes=8*GIB,
        one_global_disk_floor_bytes=50*GIB,required_available_disk_bytes=required,available_cloud_bytes=free,
        STAR_active_workers=0,STAR_future_alignment_or_count_growth_bytes=0,G128_active_workers=0,
        disk_scenario_not_hard_stop_or_size_upper_bound=True,planning_is_not_aggregate_RSS_hard_limit=True,
        no_other_project_processes_terminated_or_modified=True)

def inspect_resources(h,validated):
    bound(h['resource_epoch']);validated['epoch_validator'].no_g128_executors()
    # Rebind closed small receipts. Full original count/QC validation occurs at
    # handoff and final saved acceptance, not repeatedly during each memory wait.
    closure=read(bound(h['STAR_closure'])['path'])
    for k in ('source','queue_source','queue_started','queue_parent_started','queue_parent_exit','all_nine'):
        bound(closure[k])
    for p in closure['library_closures']: bound(p)
    try: live=proc(31592)
    except (FileNotFoundError,ProcessLookupError): live=None
    need(live is None or live['start_ticks']!=2614845177 or live['state'] in ('Z','X'),'STARCannotBeRunningAgain')
    c=Path('/sys/fs/cgroup');quota,period=(c/'cpu.max').read_text().split()
    need(quota.isdigit() and period.isdigit() and int(quota)==20*int(period),'ActualTwentyCPUQuota')
    need(set(h['worker_cpu_ids'])<=set(os.sched_getaffinity(0)),'ActualAllowedFifteenCPUs')
    admission=validated['predecessor'].load_module(validated['context']['admission'],'jgw15_original_disk')
    value=resource_arithmetic(int((c/'memory.max').read_text()),int((c/'memory.current').read_text()),
                              shutil.disk_usage(ROOT).free,20,h['worker_cpu_ids'],current_credit(h,validated),admission.COMPACT_MAX)
    return dict(value,utc=utc(),resource_epoch=h['resource_epoch'],actual_driver=h['actual_driver'],
                STAR_closure=h['STAR_closure'],STAR_complete_validation_at_handoff=True)

def predecessor_gone():
    try: value=proc(OLD_PID)
    except (FileNotFoundError,ProcessLookupError): return True
    need(value['start_ticks']!=OLD_BIRTH or value['state'] in ('Z','X'),'PredecessorStillLive')
    return True
def permission(job, h, previous, boundary, resource, output, shards=None):
    common = dict(resources_admitted=True, independent_wiring_review_accepted=True,
        previous_boundary_authenticated_by_caller=True, per_worker_memory_bytes=MEMORY_PER_WORKER,
        legacy_scientific_entry=h['legacy_R'], successor_context_source=h['new_R'],
        config=h['config'], legacy_acceptance=h['legacy_acceptance'], output_leaf=str(output),
        evidence=dict(review=h['review'], resources=resource, previous_boundary=boundary))
    if job['mode'] == 'finish':
        return dict(common, schema='jgw03-authenticated-saved-prefix-finish-acceptance-v1',
            status='ACCEPTED_FOR_ONE_SAVED_PREFIX_FINISH', maximum_parallel_single_core_workers=3,
            actual_finish=h['finish_R'], increment_helper=h['increment_helper'],
            job=dict(accession=job['accession'], space=job['space'], previous_prefix=previous))
    common.update(maximum_parallel_single_core_workers=15,
        single_owner_disjoint_partitions_authenticated_by_caller=True)
    if job['mode'] == 'merge':
        return dict(common, schema='jgw03-fifteen-worker-saved-fit-shard-merge-acceptance-v1',
            status='ACCEPTED_FOR_ONE_ORDERED_SAVED_SHARD_MERGE', actual_merge=h['merge_R'],
            fit_executor=h['fit_R'], increment_helper=h['increment_helper'], shards=shards,
            job=dict(accession=job['accession'], space=job['space'], first=job['first'], last=job['last'], previous_prefix=previous))
    need(job['mode'] == 'fit_only', 'OneExplicitFitShard')
    return dict(common, schema='jgw03-fifteen-worker-fit-only-shard-acceptance-v1',
        status='ACCEPTED_FOR_ONE_UNUSED_FIT_ONLY_SHARD', actual_executor=h['fit_R'],
        job={**{k: v for k, v in job.items() if k != 'mode'}, 'previous_prefix': previous})

def actual_entry(job, h):
    return h[{'fit_only': 'fit_R', 'merge': 'merge_R', 'finish': 'finish_R'}.get(job['mode'], 'legacy_R')]

def command(job, h, permit, previous, m):
    if job['mode'] in ('fit_only', 'merge'):
        args = ['--execute', h['config']['path'], job['accession'], job['space'], str(job['first']), str(job['last']),
            previous['path'], str(previous['bytes']), previous['sha256'], h['legacy_acceptance']['path'],
            h['legacy_acceptance']['sha256'], h['new_R']['path'], h['new_R']['sha256'], permit['path'], permit['sha256']]
        need(len(args) == 15, 'ExactFitOrMergeRArguments')
        return [h['runtime']['path'], '--vanilla', actual_entry(job, h)['path'], *args]
    return m.command(job, h, permit, previous)

def run_child(job, folder, cpu, h, permit, previous, m):
    """Same single-core child monitor, bounded stream hashes, always real wait.

    No automatic retry or fit signal. 3 GiB is monitored RSS, not RLIMIT_AS.
    """
    r_call = command(job, h, permit, previous, m); actor = actual_entry(job, h)
    call = [sys.executable, '-I', '-B', '-c', AFFINITY_EXEC, str(cpu), *r_call]
    env = dict(os.environ, R_LIBS_USER='/root/autodl-tmp/GSE128623_R_overlay_20260920/library:/root/autodl-tmp/route_b2_public_analysis/R-library',
        OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', LANG='C', LC_ALL='C', TMPDIR=str(ROOT/'tmp'))
    cmd = write_new(folder/'COMMAND.json', dict(argv=call, R_argv=r_call, job=job, assigned_CPU=cpu, utc=utc(),
        driver=h['actual_driver'], Rscript=h['runtime'], config_sha256=h['config']['sha256'], entry_sha256=actor['sha256'],
        legacy_scientific_entry=h['legacy_R'], executor_acceptance=permit, legacy_acceptance=h['legacy_acceptance'],
        fixed_thread_environment={k: env[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','LANG','LC_ALL','R_LIBS_USER','TMPDIR')}))
    child = subprocess.Popen(call, cwd=ROOT, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    started = time.monotonic(); peak = 0; observer = None; identity = None; start_pin = None; streams = {}; errors = []
    def drain(name, stream):
        digest = hashlib.sha256(); n = 0
        try:
            for block in iter(lambda: stream.read(65536), b''): digest.update(block); n += len(block)
        except BaseException as exc: errors.append(type(exc).__name__)
        finally:
            streams[name] = dict(bytes=n, sha256=digest.hexdigest())
            try: stream.close()
            except BaseException as exc: errors.append(type(exc).__name__)
    threads = []
    for name in ('stdout', 'stderr'):
        try:
            thread = threading.Thread(target=drain, args=(name, getattr(child, name)), daemon=True)
            thread.start(); threads.append(thread)
        except BaseException as exc:
            errors.append(type(exc).__name__)
            try: getattr(child, name).close()
            except BaseException as close_exc: errors.append(type(close_exc).__name__)
    try:
        identity = proc(child.pid)
        start_pin = write_new(folder/'STARTED.json', {k: identity[k] for k in ('pid','parent','start_ticks')} |
            dict(utc=utc(), job=job, config_sha256=h['config']['sha256'], entry_sha256=actor['sha256']))
        while child.poll() is None:
            try:
                text = (Path('/proc')/str(child.pid)/'status').read_text()
                peak = max([peak]+[int(x[1])*1024 for x in re.finditer(r'^Vm(?:HWM|RSS):\s+(\d+) kB$', text, re.M)])
            except (FileNotFoundError, ProcessLookupError): pass
            time.sleep(1)
    except BaseException as exc: observer = type(exc).__name__
    finally:
        actual = None
        while actual is None:
            try: actual = child.wait()
            except BaseException as exc: errors.append(type(exc).__name__)
        for thread in threads:
            joined = False
            while not joined:
                try: thread.join(); joined = True
                except BaseException as exc: errors.append(type(exc).__name__)
    if errors: observer = observer or errors[0]
    import resource
    child_max = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss*1024
    receipt = None
    if child.returncode == 0 and job['mode'] != 'run':
        try: receipt = pin(folder/'producer/R_COMPLETE.json')
        except BaseException as exc: observer = observer or type(exc).__name__
    result = dict(schema='jgw03-fifteen-worker-child-native-exit-v1',
        status='ACTUAL_BOUNDED_EXIT_NOT_INDEPENDENT_NUMERICAL_ACCEPTANCE', job=job,
        actual_pid=child.pid, actual_R_pid=child.pid, actual_native_exit_code=child.returncode,
        execution_kind='AFFINITY_BOOTSTRAP_EXECV_SAME_PID',
        actual_R_exit_code=0 if child.returncode == 0 else None, actual_executor=actor,
        command=cmd, started=start_pin, R_argv=r_call, r_complete=receipt,
        start_ticks=identity['start_ticks'] if identity else None, finished_utc=utc(),
        elapsed_seconds=time.monotonic()-started, sampled_child_peak_RSS_bytes=peak,
        controller_children_max_RSS_bytes=child_max, memory_budget_met=max(peak,child_max) <= MEMORY_PER_WORKER,
        observer_failure_class_only=observer, raw_streams_persisted=False, automatic_native_retry=False,
        stdout_bytes=streams.get('stdout',{}).get('bytes'), stderr_bytes=streams.get('stderr',{}).get('bytes'),
        stdout_sha256=streams.get('stdout',{}).get('sha256'), stderr_sha256=streams.get('stderr',{}).get('sha256'))
    write_new(folder/'PROCESS_EXIT.json', result); return result

def check_fit_receipt(result, job, h, previous, permit, keys):
    value = read(bound(result['r_complete'])['path'])
    need(value['schema'] == 'jgw03-fifteen-worker-fit-only-shard-result-v1' and
         value['status'] == 'FIT_ONLY_SHARD_COMPLETE_NOT_PREFIX_OR_INDEPENDENT_ACCEPTANCE' and
         value['actual_executor'] == h['fit_R'] and value['legacy_scientific_entry'] == h['legacy_R'] and
         value['accession'] == job['accession'] and value['space'] == job['space'] and
         all(type(value[k]) is int and value[k] == job[k] for k in ('first','last','round_first','round_last')) and
         value['scientific_identity_sha256'] == h['initial_previous'][job['accession']]['scientific_identity_sha256'] and
         value['previous_prefix'] == previous and value['executor_acceptance'] == permit and
         value['prefix_written'] is False and value['actual_native_exit_observed_here'] is False,
         'ExactFitOnlyProducerReceipt')
    records = value['allocation_records']
    need(len(records) == job['last']-job['first']+1, 'CompleteFitShardIDs')
    for i, row in zip(range(job['first'],job['last']+1), records):
        need(type(row['allocation_id']) is int and row['allocation_id'] == i and
             row['allocation_key'] == keys['cells'][job['accession']]['four_group_MC'][i-1] and
             row['execution_entry'] == h['fit_R'] and row['prior_leaf_existed_at_preflight'] is False and
             type(row['new_native_calls_recorded_in_archive']) is int and row['new_native_calls_recorded_in_archive'] == 1,
             'OrderedUnusedFitShardRecords')
        bound(row['completion_metadata'])
    return result['r_complete']

def check_merge_receipt(result, job, h, previous, permit):
    value = read(bound(result['r_complete'])['path'])
    need(value['schema'] == 'jgw03-fifteen-worker-saved-fit-shard-merge-result-v1' and
         value['status'] == 'ORDERED_SAVED_SHARD_PREFIX_COMPLETE_NOT_INDEPENDENT_ACCEPTANCE' and
         value['actual_merge'] == h['merge_R'] and value['fit_executor'] == h['fit_R'] and
         value['legacy_scientific_entry'] == h['legacy_R'] and value['previous_prefix'] == previous and
         value['executor_acceptance'] == permit and value['accession'] == job['accession'] and
         value['space'] == job['space'] and value['first'] == job['first'] and value['last'] == job['last'] and
         value['scientific_identity_sha256'] == h['initial_previous'][job['accession']]['scientific_identity_sha256'] and
         type(value['new_fits']) is int and value['new_fits'] == 0 and value['scientific_acceptance'] is False,
         'OrderedSavedMergeNotNewFitsOrScienceAcceptance')
    p = bound(value['new_prefix'])
    need(p['path'] == str(ROOT/'derived'/job['accession']/f"four_group_MC_TAIL_PREFIX_{job['last']:05d}.rds"), 'CanonicalMergedPrefix')
    return p

def run_batch(jobs, folder, cpus, h, previous, boundary, resource, validated, shard_pins=None, runner=run_child):
    need(1 <= len(jobs) <= 15 and len(cpus) >= len(jobs) and len(set(cpus)) == len(cpus),
         'EveryJobHasDistinctExplicitCPU')
    folder.mkdir(); permits = {}; leaves = []
    for job in jobs:
        key = job['accession']+(f"__shard_{job['shard_index']}" if job['mode']=='fit_only' else '')
        leaf = folder/key; leaf.mkdir(); leaves.append(leaf)
        if job['mode'] != 'run':
            permits[key] = write_new(leaf/'EXECUTOR_ACCEPTANCE.json', permission(job, h,
                previous[validated['predecessor'].prefix_key(job)], boundary, resource, leaf/'producer',
                shard_pins.get(job['accession']) if shard_pins else None))
    results = []
    with ThreadPoolExecutor(max_workers=len(jobs)) as pool:
        futures = [pool.submit(runner, job, leaf, cpu, h, permits.get(leaf.name),
            previous.get(validated['predecessor'].prefix_key(job)), validated['predecessor']) for job,leaf,cpu in zip(jobs,leaves,cpus)]
        for future in futures:
            try: results.append(future.result())
            except BaseException as exc:
                results.append(dict(actual_native_exit_code=None, actual_R_exit_code=None, memory_budget_met=False,
                    observer_failure_class_only=type(exc).__name__))
    accepted = all(validated['predecessor'].strict_child(x, validated['legacy']) for x in results)
    whole = write_new(folder/'ROUND_EXIT.json', dict(utc=utc(), results=results,
        all_actual_exits_zero_and_observed=accepted, actual_driver=h['actual_driver']))
    need(accepted, 'ActualNativeBatchFailedNoRetry')
    return results, leaves, permits, whole

def validate_saved_round(folder, number, h, previous, boundary, validated, controller_pid):
    """Metadata-only counterpart of execute_round, shared with the final gate.

    It never reads an allocation RDS or recomputes a tail. The independent R
    checker still owns the one full mathematical archive verification.
    """
    context = validated['context']; m = validated['predecessor']; legacy = validated['legacy']
    originals = legacy.plan()[number-1]; resource = pin(Path(folder).parent/f'RESOURCES_{number:04d}.json')
    need(read(resource['path'])['resources_admitted'] is True, 'ActualAdmittedRoundResourceReceipt')
    def batch(jobs, subfolder, previous_boundary, shard_pins=None):
        values = []; pins = {}; receipts = {}; actual_cpus = []
        for job in jobs:
            key = job['accession']+(f"__shard_{job['shard_index']}" if job['mode']=='fit_only' else '')
            leaf = subfolder/key
            cp, sp, np = (pin(leaf/name) for name in ('COMMAND.json','STARTED.json','PROCESS_EXIT.json'))
            c, s, n = (read(p['path']) for p in (cp,sp,np)); permit = None
            if job['mode'] != 'run':
                permit = pin(leaf/'EXECUTOR_ACCEPTANCE.json')
                need(read(permit['path']) == permission(job, context, previous[m.prefix_key(job)], previous_boundary,
                    resource, leaf/'producer', shard_pins.get(job['accession']) if shard_pins else None), 'ExactSavedJobPermit')
            expected = command(job,context,permit,previous.get(m.prefix_key(job)),m); actor = actual_entry(job,context)
            need(c['job'] == s['job'] == n['job'] == job and
                 all(type(c['job'].get(k)) is type(v) for k,v in job.items()) and
                 c['R_argv'] == n['R_argv'] == expected and c['driver'] == h['actual_driver'] and
                 c['argv'] == [validated['context']['legacy_python'], '-I','-B','-c',
                              AFFINITY_EXEC,str(c['assigned_CPU']),*expected] and
                 c['entry_sha256'] == s['entry_sha256'] == actor['sha256'] and
                 c['Rscript'] == context['runtime'] and integer(c['assigned_CPU']) and
                 c['assigned_CPU'] in h['worker_cpu_ids'] and
                 n['schema'] == 'jgw03-fifteen-worker-child-native-exit-v1' and
                 n['execution_kind'] == 'AFFINITY_BOOTSTRAP_EXECV_SAME_PID' and
                 n['actual_executor'] == actor and n['command'] == cp and n['started'] == sp and
                 integer(s['pid'],1) and integer(s['start_ticks'],1) and type(s['parent']) is int and
                 s['parent'] == controller_pid and s['pid'] == n['actual_pid'] == n['actual_R_pid'] and
                 s['start_ticks'] == n['start_ticks'] and m.strict_child(n,legacy), 'ActualSavedChildSourceNativeAndParent')
            actual_cpus.append(c['assigned_CPU']); values.append(n); pins[key] = np
            if job['mode'] == 'fit_only':
                receipts[key] = check_fit_receipt(n,job,context,previous[job['accession']],permit,validated['keys'])
            elif job['mode'] == 'merge':
                receipts[key] = (check_merge_receipt(n,job,context,previous[job['accession']],permit),n['r_complete'])
            elif job['mode'] == 'finish':
                receipts[key] = m.validate_finish_receipt(leaf/'producer/R_COMPLETE.json',job,context,
                    previous[m.prefix_key(job)],permit)
        whole = pin(subfolder/'ROUND_EXIT.json'); w = read(whole['path'])
        need(len(set(actual_cpus)) == len(jobs) and w['results'] == values and
             w['all_actual_exits_zero_and_observed'] is True and w['actual_driver'] == h['actual_driver'], 'SavedWholeEqualsActualChildren')
        return whole,pins,receipts
    updated = dict(previous); bindings = {}
    if number <= 201:
        jobs = fit_jobs(originals); fit_whole,natives,receipts = batch(jobs,Path(folder)/'fit',boundary)
        shards = {cell:[] for cell in CELLS}
        for job in jobs:
            key = job['accession']+f"__shard_{job['shard_index']}"
            shards[job['accession']].append(dict(r_complete=receipts[key],native_exit=natives[key]))
        whole,_,receipts = batch([dict(j,mode='merge') for j in originals],Path(folder)/'merge',fit_whole,shards)
        for cell in CELLS: updated[cell],bindings[cell] = receipts[cell]
    else:
        whole,_,receipts = batch(originals,Path(folder)/'native',boundary)
        for job in originals:
            if job['mode'] == 'finish': bindings[job['accession']] = receipts[job['accession']]
            else:
                need(number == 202 and job['space'] == 'high_companion_exact', 'SavedOriginalExact20Only')
                updated[m.prefix_key(job)] = pin(ROOT/'derived'/job['accession']/'high_companion_exact_TAIL_PREFIX_00020.rds')
    acceptance = pin(Path(folder)/'ROUND_ACCEPTED.json'); a = read(acceptance['path'])
    need(a == dict(scope='ACTUAL_NATIVE_ZERO_AND_ORDERED_PREFIX_BINDING_ONLY',whole=whole,
         producer_receipts=bindings,next_previous=updated,scientific_acceptance=False,old_controller_success_asserted=False),
         'SavedRoundAcceptanceExactScope')
    return updated,whole,acceptance

def execute_round(number, original_jobs, folder, h, previous, boundary, resource, validated):
    folder.mkdir(); context = validated['context']; m = validated['predecessor']
    updated = dict(previous); receipts = {}; cpus = h['worker_cpu_ids']
    if number <= 201:
        jobs = fit_jobs(original_jobs)
        result, leaves, permits, fit_whole = run_batch(jobs, folder/'fit', cpus, context, previous, boundary, resource, validated)
        shards = {c: [] for c in CELLS}
        for job, leaf, native in zip(jobs, leaves, result):
            r = check_fit_receipt(native, job, context, previous[job['accession']], permits[leaf.name], validated['keys'])
            shards[job['accession']].append(dict(r_complete=r, native_exit=pin(leaf/'PROCESS_EXIT.json')))
        merge_jobs = [dict(j, mode='merge') for j in original_jobs]
        result, leaves, permits, whole = run_batch(merge_jobs, folder/'merge', cpus[:3], context, previous,
            fit_whole, resource, validated, shards)
        for job, leaf, native in zip(merge_jobs, leaves, result):
            updated[job['accession']] = check_merge_receipt(native, job, context, previous[job['accession']], permits[leaf.name])
            receipts[job['accession']] = native['r_complete']
    else:
        result, leaves, permits, whole = run_batch(original_jobs, folder/'native', cpus[:3], context,
            previous, boundary, resource, validated)
        for job, leaf, native in zip(original_jobs, leaves, result):
            if job['mode'] == 'finish':
                receipts[job['accession']] = m.validate_finish_receipt(leaf/'producer/R_COMPLETE.json', job, context,
                    previous[m.prefix_key(job)], permits[leaf.name])
            else:
                need(number == 202 and job['space'] == 'high_companion_exact', 'OriginalFullExactOnly')
                updated[m.prefix_key(job)] = pin(ROOT/'derived'/job['accession']/'high_companion_exact_TAIL_PREFIX_00020.rds')
    accepted = write_new(folder/'ROUND_ACCEPTED.json', dict(scope='ACTUAL_NATIVE_ZERO_AND_ORDERED_PREFIX_BINDING_ONLY',
        whole=whole, producer_receipts=receipts, next_previous=updated, scientific_acceptance=False,
        old_controller_success_asserted=False))
    return updated, whole, accepted

def run(handoff_path, handoff_sha, attempt):
    need(sys.platform.startswith('linux') and re.fullmatch('[A-Za-z0-9_-]{1,80}', attempt), 'LinuxFreshAttemptOnly')
    hp = pin(handoff_path); need(hp['sha256'] == handoff_sha and Path(hp['path']).is_relative_to(ROOT/'provenance'), 'ExternalHandoff')
    h = read(handoff_path); v = validate_handoff(h); folder = ROOT/'provenance'/attempt
    need(not folder.exists(), 'FreshAttemptNoAutomaticReplay')
    import fcntl
    with (ROOT/'provenance/FULL_REFERENCE.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB); predecessor_gone(); folder.mkdir()
        identity = proc(os.getpid())
        write_new(folder/'CONTROLLER_STARTED.json', dict(schema='jgw03-fifteen-worker-controller-v1', **{k: identity[k] for k in ('pid','parent','start_ticks')}, **dict(actual_driver=h['actual_driver'], handoff=hp, argv=sys.argv, executable=sys.executable,
                assigned_CPUs=h['worker_cpu_ids'], resource_epoch=h['resource_epoch'], old_controller_success_asserted=False, utc=utc())))
        write_new(folder/'PLAN.json', dict(remaining=v['predecessor'].remaining_plan(v['legacy'],h['completed_rounds']),
            maximum_workers=15, per_worker_memory_bytes=MEMORY_PER_WORKER, ordered_merge_workers=3,
            exact_and_finish_workers=3, completed_before_this_epoch=h['completed_rounds'], scientific_acceptance=False))
        previous = {c: h['initial_previous'][c]['prefix'] for c in CELLS}
        boundary = h['previous_boundary']; completed = []
        try:
            for number, jobs in v['predecessor'].remaining_plan(v['legacy'], h['completed_rounds']):
                while True:
                    bound(hp); predecessor_gone()
                    for k in ('actual_driver','fit_R','merge_R','new_R','finish_R','increment_helper','review','resource_epoch'):
                        bound(h[k])
                    for k in ('legacy_driver','admission','legacy_R','config','runtime','legacy_acceptance','keys_sidecar'):
                        bound(v['context'][k])
                    state = inspect_resources(h, v)
                    if state['resources_admitted'] is True: break
                    v['predecessor'].progress(folder, dict(status='WAIT_RESOURCE_WINDOW_NO_NATIVE_RESTART', next_round=number, resource_state=state))
                    time.sleep(10)
                resource = write_new(folder/f'RESOURCES_{number:04d}.json',state)
                v['predecessor'].progress(folder, dict(status='RUNNING_BOUND_FIFTEEN_WORKER_ROUND', round=number, completed_remaining_rounds=completed))
                previous,boundary,accepted = execute_round(number,jobs,folder/f'round_{number:04d}',h,previous,boundary,resource,v)
                completed.append(dict(number=number,whole=boundary,accepted=accepted))
            result = dict(status='REMAINING_SCOPE_ACTUAL_NATIVE_ZERO_NOT_FULL_SCIENCE_ACCEPTANCE',
                completed_remaining_rounds=completed, completed_before_this_epoch=h['completed_rounds'],
                actual_driver=h['actual_driver'], handoff=hp, old_controller_success_asserted=False, scientific_acceptance=False, utc=utc())
            write_new(folder/'DISPATCH_EXIT.json', result); v['predecessor'].progress(folder,result); return 0
        except BaseException as exc:
            write_new(folder/'DISPATCH_EXIT.json',dict(status='STOP_RETAINED_NO_AUTOMATIC_NATIVE_RETRY',
                exception_class_only=type(exc).__name__,completed_remaining_rounds=completed,
                old_controller_success_asserted=False,scientific_acceptance=False,raw_error_persisted=False,utc=utc()))
            return 1

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute', action='store_true'); parser.add_argument('--handoff',type=Path)
    parser.add_argument('--handoff-sha256'); parser.add_argument('--attempt'); args = parser.parse_args()
    if not args.execute:
        print(json.dumps(dict(status='PLAN_ONLY_NO_SIGNALS_OR_FITS',maximum_workers=15,
            per_worker_memory_bytes=MEMORY_PER_WORKER, source=pin(Path(__file__).resolve())))); return 0
    need(all((args.handoff,args.handoff_sha256,args.attempt)), 'ExplicitFrozenHandoffAndFreshAttempt')
    return run(args.handoff,args.handoff_sha256,args.attempt)

if __name__=='__main__':
    try: raise SystemExit(main())
    except Exception as exc:
        print(json.dumps(dict(status='NO_NEW_FIFTEEN_DISPATCH',exception_class_only=type(exc).__name__,raw_error_persisted=False)))
        raise SystemExit(2)
