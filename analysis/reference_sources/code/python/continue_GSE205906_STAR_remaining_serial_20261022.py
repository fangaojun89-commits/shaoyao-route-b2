"""Observe the existing SRR19631206, then six unused libraries in fixed order.

Plan-only by default. No retry, cache advice, download, model or matrix stage.
Each new child is the unchanged reviewed one-library frontend, actually waited
by its direct parent. Existing scientific processes are never signalled.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import socket
import sys
import time
import types

ROOT = Path('/root/autodl-tmp/GSE205906_matched_STAR_Ensembl101_20261015r01')
HOST = 'autodl-container-62c8pzy536-50fb405a'
FRONTEND = ROOT/'one_library_source_20261022_r01/code/python/run_gse205906_STAR_one_remaining_library_20261022.py'
FRONTEND_SHA = '05fb686b5401ea11a3ed1c248b13cf0fe42fd0464c8ecf5946cc37031ef78685'
OBSERVER = ROOT/'post_G128_source_20261022_r02/code/python/start_GSE205906_STAR_post_G128_prefix_v2_20261022.py'
OBSERVER_SHA = '835552fbe28150c38bb0f022d8c22b8d99270f5beed4ba02a611d36f92b0509e'
ORIGINAL = ROOT/'source_staging_20261015_01/run_gse205906_matched_STAR_20261015.py'
ORIGINAL_SHA = '7fe64fcfcee554e1500f60da2cd7b4a3ce8a84fcd4e46b84c3c7e14ec1fd907d'
PREVIOUS = ROOT/'provenance/remaining_SRR19631206_20261022_r01'
PREVIOUS_START_SHA = '15056741e631357cfba86a60d979884f893201b24ac209bc903503a0353d723f'
LIBRARIES = ('SRR19631205','SRR19631204','SRR19631203','SRR19631201','SRR19631200','SRR19631199')
RESERVE = 270287983845
BUDGET = 56*1024**3
MEMORY = 40*1024**3

def need(ok, category):
    if not ok: raise ValueError(category)

def module(path, sha, name):
    need(path.is_file() and not path.is_symlink(), 'RegularReviewedSource')
    raw=path.read_bytes();need(hashlib.sha256(raw).hexdigest()==sha,'ReviewedSourceSHA')
    m=types.ModuleType(name);m.__file__=str(path);exec(compile(raw,str(path),'exec'),m.__dict__)
    return m

def native_zero(start, end):
    need(type(end['actual_native_exit_code']) is int and end['actual_native_exit_code']==0,
         'ActualPreviousNativeZero')
    need(end['observer_failure_class_only'] is None and end['raw_streams_persisted'] is False,
         'ActualPreviousObserverComplete')
    need(start['actual_pid']==end['actual_pid']==start['identity']['pid']==end['identity']['pid']
         and start['identity']['start_ticks']==end['identity']['start_ticks']
         and start['argv']==end['argv']==start['identity']['argv']==end['identity']['argv'],
         'ActualPreviousStartExitIdentity')
    need(set(end['stream_summaries'])=={'stdout','stderr'} and
         all(x['complete'] is True for x in end['stream_summaries'].values()),'CompleteDigestStreams')

def serial(ops):
    """Small scheduling core; callbacks are invented in pure tests."""
    ops.wait_previous()
    for run in LIBRARIES:
        ops.unused(run)
        resources=ops.wait_resources(run)
        ops.unused(run)
        ops.once(run, resources)
    ops.finish()

class Actual:
    def __init__(self, out, observer, original):
        self.out,self.o,self.a=out,observer,original

    def complete(self, run):
        d=self.a.completed(ROOT/'libraries'/run,ORIGINAL_SHA)
        registered=next(x for x in self.a.metadata()['libraries'] if x['SRA_run']==run)
        need(d['SRA_run']==run and d['input_pairs']==registered['registered_pairs'] and
             d['technical_reservoir_only'] is False and len(d['events'])==8 and
             all(type(x['actual_exit']) is int and x['actual_exit']==0 for x in d['events']),
             'EightOriginalNativeCommandsZero')
        return d

    def wait_previous(self):
        need(self.o.pin(PREVIOUS/'STARTED.json')['sha256']==PREVIOUS_START_SHA,'FixedExistingLibraryStart')
        start=self.a.load(PREVIOUS/'STARTED.json');gone_since=None
        while not (PREVIOUS/'NATIVE_EXIT.json').exists():
            try: current=self.o.proc(start['actual_pid'])
            except (FileNotFoundError,ProcessLookupError): current=None
            if current and current['state'] not in ('Z','X'):
                need(current['start_ticks']==start['identity']['start_ticks'] and current['argv']==start['argv'],
                     'ExistingLibraryProcessIdentity');gone_since=None
            else:
                if gone_since is None: gone_since=time.monotonic()
                need(time.monotonic()-gone_since<=30,'MissingActualPreviousTerminalNoRestart')
            time.sleep(2)
        native_zero(start,self.a.load(PREVIOUS/'NATIVE_EXIT.json'));self.complete('SRR19631206')
        self.o.put(self.out/'PREDECESSOR_ACCEPTED.json',{'native_exit':self.o.pin(PREVIOUS/'NATIVE_EXIT.json'),
            'checkpoint':self.o.pin(ROOT/'libraries/SRR19631206/COMPLETE.json'),
            'existing_library_restarted':False,'utc':self.o.utc()})

    def unused(self, run):
        need(not (ROOT/'libraries'/run).exists(),'ExistingPartialOrCompletedTargetNoReplay')
        need(not (ROOT/'libraries'/run).is_symlink(),'NoTargetSymlink')

    def wait_resources(self, run):
        args=argparse.Namespace(quiet_window_confirmed=True,reserved_growth_bytes=RESERVE,disk_budget_bytes=BUDGET)
        noted=False
        while True:
            self.a.unchanged(ORIGINAL_SHA)
            need(shutil.disk_usage(ROOT).free>=RESERVE+BUDGET,'FullFixedDiskReservationUnavailable')
            maximum=Path('/sys/fs/cgroup/memory.max').read_text().strip()
            need(maximum!='max','FiniteMemoryLimit')
            current=int(Path('/sys/fs/cgroup/memory.current').read_text())
            if int(maximum)-current>=MEMORY: return self.a.resources(args)
            if not noted:
                self.o.put(self.out/(run+'_WAIT_RESOURCE.json'),{'utc':self.o.utc(),
                    'library':run,'waiting_only_no_scientific_call':True,'required_direct_headroom_bytes':MEMORY})
                noted=True
            time.sleep(15)

    def once(self, run, resources):
        for path,sha in ((FRONTEND,FRONTEND_SHA),(OBSERVER,OBSERVER_SHA),(ORIGINAL,ORIGINAL_SHA)):
            need(self.o.pin(path)['sha256']==sha,'UnchangedBeforeNextLibrary')
        leaf=self.out/run;leaf.mkdir()
        cmd=['/usr/bin/python3','-I','-B',str(FRONTEND),run,'--execute','--quiet-window-confirmed',
             '--reserved-growth-bytes',str(RESERVE),'--disk-budget-bytes',str(BUDGET)]
        self.o.put(leaf/'INTENT.json',{'argv':cmd,'assigned_CPUs':[203,204],'resource_before':resources,
            'frontend':self.o.pin(FRONTEND),'original':self.o.pin(ORIGINAL),'observer':self.o.pin(OBSERVER),
            'automatic_retry':False,'original_CLI_execute_invoked':False,'utc':self.o.utc()})
        end=self.o.observe(cmd,[203,204],leaf);self.o.put(leaf/'NATIVE_EXIT.json',end)
        native_zero(self.a.load(leaf/'STARTED.json'),end);self.complete(run)
        self.o.put(leaf/'ACCEPTED.json',{'status':'ORIGINAL_LIBRARY_COMPLETE_AND_ACTUAL_PARENT_ZERO',
            'checkpoint':self.o.pin(ROOT/'libraries'/run/'COMPLETE.json'),
            'native_exit':self.o.pin(leaf/'NATIVE_EXIT.json'),'utc':self.o.utc()})

    def finish(self):
        pins=[]
        for lib in self.a.metadata()['libraries']:
            self.complete(lib['SRA_run']);pins.append(self.o.pin(ROOT/'libraries'/lib['SRA_run']/'COMPLETE.json'))
        self.o.put(self.out/'ALL_NINE_LIBRARIES_COMPLETE.json',{'scope':'ORIGINAL_LIBRARY_CHECKPOINTS_ONLY',
            'checkpoints':pins,'new_matrix_executed':False,'new_model_executed':False,
            'matched_comparison_complete':False,'utc':self.o.utc()})

def execute(attempt):
    need(sys.platform.startswith('linux') and socket.gethostname()==HOST,'FixedShaoyaoInstance')
    need(re.fullmatch('[A-Za-z0-9_-]{1,80}',attempt) is not None,'SimpleFreshAttempt')
    o=module(OBSERVER,OBSERVER_SHA,'unchanged_parent_observer')
    a=module(ORIGINAL,ORIGINAL_SHA,'unchanged_original_library')
    module(FRONTEND,FRONTEND_SHA,'unchanged_one_library_frontend')
    a.unchanged(ORIGINAL_SHA)
    import fcntl
    with (ROOT/'provenance/REMAINING_SERIAL_QUEUE.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        out=ROOT/'provenance'/attempt;out.mkdir()
        ops=Actual(out,o,a)
        try:
            for run in ('SRR19631202','SRR19631207'): ops.complete(run)
            for run in LIBRARIES: ops.unused(run)
            o.put(out/'QUEUE_STARTED.json',{'identity':o.proc(os.getpid()),'source':o.pin(Path(__file__)),
                'existing_predecessor_start':o.pin(PREVIOUS/'STARTED.json'),'fixed_remaining':list(LIBRARIES),
                'automatic_retry':False,'new_models':0,'cache_advice':False,'utc':o.utc()})
            serial(ops)
            return 0
        except Exception as exc:
            o.put(out/'QUEUE_STOP.json',{'exception_class_only':type(exc).__name__,
                'raw_error_text_persisted':False,'automatic_retry':False,'utc':o.utc()})
            return 1

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--execute',action='store_true')
    p.add_argument('--attempt');args=p.parse_args()
    if not args.execute:
        print(json.dumps({'status':'PLAN_ONLY_NO_IO_OR_SCIENTIFIC_CALL','existing_wait_only':'SRR19631206',
                          'fixed_remaining':LIBRARIES,'automatic_retry':False}));return 0
    need(args.attempt is not None,'ExplicitFreshAttempt')
    return execute(args.attempt)

if __name__=='__main__':
    try: raise SystemExit(main())
    except Exception as exc:
        print(json.dumps({'status':'STOP_NO_AUTOMATIC_RETRY','exception_class_only':type(exc).__name__}));raise SystemExit(2)
