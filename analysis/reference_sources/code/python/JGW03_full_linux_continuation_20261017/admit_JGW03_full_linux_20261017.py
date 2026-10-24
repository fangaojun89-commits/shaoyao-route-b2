"""Saved-evidence/full-resource admission only; never launch R, transfer or delete.

Default CLI is inspect. Acceptance requires the explicit assemble subcommand
and --create-acceptance. Source/config identities remain the frozen v3 identities;
this helper and the bounded driver are separately pinned execution tools.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import re
import stat

GIB = 1024**3
ROOT = Path('/root/autodl-tmp/JGW03_full_reference_linux_20261017r03')
G128 = Path('/root/autodl-tmp/gse128623_label_reference_20261015r01')
RELEASE = Path('/root/autodl-tmp/GSE205906_SAM_archive_receipts_20261017r01/release_01')
SOURCE_PARENT = '/root/autodl-tmp/GSE205906_full_read_counts_fasp_v18e_20261005/formal/intermediate'
ARCHIVE = 'R:/Shaoyao_storage/archives/GSE205906_cloud_SAM_preservation_20261017_r01'
CELLS = ('GSE289929', 'GSE227494', 'GSE299063')
FAMILIES = {'GSE289929': (31, 12), 'GSE227494': (29, 11), 'GSE299063': (29, 11)}
ENTRY = 'code/R/run_JGW03_full_reference_linux_v3_20261017.R'
CHECKER = 'code/R/check_JGW03_linux_native_compatibility_v3_20261017.R'
CONFIG_SHA = 'eb1a72bbecd61878e61a073b025480a0683c6be76c135257ec2e7943f30bf01c'
ENTRY_SHA = 'c0c610da98b6511689b46bb5e3a98e87b64b02fa1ed978718db67b412591dd2d'
CHECKER_SHA = '63fe74ca1529e3d4a9fa020b47e1a8b9550d2d3a448b1acdabda4e2a81ad848d'
SIX_SHA = '43f610bf42cee00f009f3e436c6a39d2fd3618956a86f5d9fad7a605be27874a'
RESOURCE_SHA = 'c27f7db86584f0c8ea2c344b675bc6849a77d531e5678eca69059366308d2ec7'
PLAN_SHA = '35c5cb4d790529854756a778db46704cc6013bac5ed7eab5997ea73bb43ae647'
OLD_RESOURCE_SHA = 'f388ecf2574d935865fb46d6cd0fa926284ed0fb8690ec2a043d5c0318d85e16'
SELECTION_SHA = 'b285ad6d35679a03df2fb190c0eb6e33316d145b36b7adcb62b1f5beddb6c386'
APPROVAL_SHA = '907365c737505c9ce9396c10716a0bbc42ff681403a7bcb97eedff4e446790d0'
COPY_SHA = '47a9e0471b25975076e76868549409ea5b88e1cf160886c076b8f617e6eab374'
RELEASE_SHA = '7bcd64669267cfe082572fa22707c240814deae9cc5e410320beadbd29d5931c'
SAM_BYTES = (25745781978, 22500689884, 27146492496, 28859470133, 26317228454,
             26782082390, 25995660079, 26720732929, 26093300590)
COMPACT_MAX = (5437884, 6025095, 6052566)
RUNTIME = Path('/root/autodl-tmp/route_b2_public_analysis/runtime/bin/Rscript')


def need(ok, category):
    if not ok:
        raise ValueError(category)


def integer(x, minimum=0):
    return type(x) is int and x >= minimum


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def clean(path):
    p = Path(path)
    need(p.is_absolute() and not any(q.is_symlink() for q in (p, *p.parents)), 'PathIdentity')
    return p


def pin(path):
    p = clean(path)
    need(p.name.lower() != 'checked_payload.rds', 'ForbiddenPayload')
    s = p.stat()
    need(stat.S_ISREG(s.st_mode), 'RegularFileRequired')
    before = (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
    h = hashlib.sha256()
    with p.open('rb') as stream:
        # Native admission is Linux. Windows invented fixtures have differing
        # path-stat/fd-stat ctime semantics; do not weaken the Linux comparison.
        fields = ('st_dev', 'st_ino', 'st_size', 'st_mtime_ns')
        if os.name != 'nt':
            fields += ('st_ctime_ns',)
        need(tuple(getattr(os.fstat(stream.fileno()), k) for k in fields) ==
             tuple(getattr(s, k) for k in fields), 'PinChangedBeforeRead')
        for block in iter(lambda: stream.read(1024**2), b''):
            h.update(block)
    s = p.stat()
    need((s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns) == before, 'PinChangedDuringRead')
    return {'path': str(p), 'bytes': s.st_size, 'sha256': h.hexdigest()}


def bound(record):
    need(isinstance(record, dict) and integer(record.get('bytes')) and
         isinstance(record.get('sha256'), str) and re.fullmatch('[0-9a-f]{64}', record['sha256']), 'FileRecordSchema')
    actual = pin(record['path'])
    need(actual['bytes'] == record['bytes'] and actual['sha256'] == record['sha256'], 'FilePinMismatch')
    return actual


def read(path):
    p = clean(path)
    need(p.suffix == '.json' and p.stat().st_size <= 8 * 1024**2, 'SmallJSONRequired')
    def unique(pairs):
        result = {}
        for key, value in pairs:
            need(key not in result, 'DuplicateJSONKey')
            result[key] = value
        return result
    return json.loads(p.read_text(encoding='utf-8'), object_pairs_hook=unique,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError('NonfiniteJSON')))


def configuration():
    cfg_pin = pin(ROOT / 'CONFIG.json')
    need(cfg_pin['sha256'] == CONFIG_SHA, 'FrozenConfigIdentity')
    cfg = read(cfg_pin['path'])
    need(cfg['schema'] == 'jgw03-separate-linux-reference-config-v2' and
         cfg['run_root'] == str(ROOT) and cfg['project_root'] == str(ROOT / 'project') and
         set(cfg['cells']) == set(CELLS) and len(cfg['source_pins']) == 22, 'FrozenConfigNamespace')
    for rel, expected in cfg['source_pins'].items():
        need(not Path(rel).is_absolute() and '..' not in Path(rel).parts, 'SourceRelativePath')
        actual = pin(ROOT / 'project' / rel)
        need({k: actual[k] for k in ('bytes', 'sha256')} == expected, 'FrozenSourceIdentity')
    need(cfg['source_pins'][ENTRY]['sha256'] == ENTRY_SHA and
         cfg['source_pins'][CHECKER]['sha256'] == CHECKER_SHA and
         cfg['source_pins']['code/python/run_JGW03_linux_bridge_20261017.py']['sha256'] == OLD_RESOURCE_SHA,
         'FrozenEntryCheckerResourceSource')
    return cfg, cfg_pin


def calculate(workers=3, g128_remaining=0, materialized_compact_bytes=0):
    need(integer(workers, 1) and workers <= 3 and integer(g128_remaining) and
         integer(materialized_compact_bytes), 'PlanningIntegerDomain')
    base = sum(COMPACT_MAX) * 10020
    compact = (base * 5 + 3) // 4
    overhead = (8 + 8 + 8 + 16) * GIB
    remaining = max(0, compact - materialized_compact_bytes) + overhead
    return {'planning_scenario_not_size_upper_bound': True, 'workers': workers,
        'single_CPU_per_worker': True, 'maximum_workers': 3, 'compact_base_bytes': base,
        'compact_with_25_percent_variation_bytes': compact, 'metadata_prefix_bytes': 8*GIB,
        'QA_bytes': 8*GIB, 'capsule_scratch_bytes': 8*GIB, 'safety_bytes': 16*GIB,
        'initial_JGW_disk_scenario_bytes': compact + overhead,
        'credited_materialized_compact_logical_bytes': materialized_compact_bytes,
        'remaining_JGW_disk_scenario_bytes': remaining, 'original_finite_disk_allowance_bytes': 8*GIB,
        'remaining_original_G128_growth_and_floor_bytes': g128_remaining,
        'required_available_disk_bytes': remaining + 8*GIB + g128_remaining,
        'worker_memory_bytes': workers*3*GIB, 'memory_floor_bytes': 8*GIB}


def release_evidence(expected=None):
    exit_path = RELEASE / 'RELEASE_EXIT.json'
    rec = bound(expected) if expected is not None else pin(exit_path)
    need(rec['path'] == str(exit_path), 'ExactReleaseReceipt')
    result = read(exit_path)
    intent = read(RELEASE / 'ALL_NINE_RELEASE_INTENT.json')
    need(result['status'] == 'EXACT_NINE_CLOUD_SOURCES_RELEASED' and result['exit_code'] == 0 and
         result['released_files'] == 9 and result['released_logical_bytes'] == sum(SAM_BYTES) and
         result['original_manifests_modified'] is False and result['Atlas_modified'] is False,
         'AllNineActualReleaseRequired')
    need(intent['approval_sha256'] == APPROVAL_SHA and intent['selection_sha256'] == SELECTION_SHA and
         intent['copy_source_sha256'] == COPY_SHA and intent['release_source_sha256'] == RELEASE_SHA and
         intent['archive_root'] == ARCHIVE and intent['cloud_receipt'] == str(RELEASE) and
         re.fullmatch('[0-9a-f]{64}', intent['all_nine_acceptance_sha256']) and
         len(intent['files']) == len(intent['rows']) == len(intent['R_fresh_full_readbacks']) == len(result['released']) == 9,
         'ReleaseAcceptanceChain')
    pins = [rec, pin(RELEASE / 'ALL_NINE_RELEASE_INTENT.json')]
    for i, (accepted, row, back, released) in enumerate(zip(intent['files'], intent['rows'],
               intent['R_fresh_full_readbacks'], result['released'])):
        cid = f'B2a-{i+1:02d}'; source = f'{SOURCE_PARENT}/SRR{19631199+i}/full_rat.sam'
        r_path = ARCHIVE + '/payload' + source
        need(accepted['candidate_id'] == row['candidate_id'] == back['candidate_id'] == released['candidate_id'] == cid
             and accepted['source_path'] == row['source_absolute_path'] == released['source_path'] == source
             and accepted['R_path'] == released['R_path'] == r_path
             and accepted['bytes'] == back['bytes'] == released['bytes'] == int(row['bytes']) == SAM_BYTES[i]
             and accepted['sha256'] == back['sha256'] == released['sha256']
             and re.fullmatch('[0-9a-f]{64}', accepted['sha256'])
             and accepted['native_transfer_exit'] == 0 and accepted['complete_closed_file_readback'] is True,
             'ExactNinePerFileSHAChain')
        before, after = accepted['source_before'], accepted['source_after']
        need(before['ok'] is True and after['ok'] is True and before['sha256'] == accepted['sha256'] and
             before['read_bytes'] == SAM_BYTES[i] and before['stat'] == after['stat'] and
             before['ctime_ns'] == after['ctime_ns'], 'OriginalSourceFullReadbackChain')
        one_intent = read(RELEASE / (cid + '_INTENT.json'))
        one_done = read(RELEASE / (cid + '_RELEASED.json'))
        need(one_done == released and one_intent['source_path'] == source and
             one_intent['R_path'] == r_path and one_intent['bytes'] == SAM_BYTES[i] and
             one_intent['sha256'] == accepted['sha256'] and
             one_intent['archive_acceptance_sha256'] == intent['all_nine_acceptance_sha256'] and
             not os.path.lexists(source), 'ActualPerFileReleaseAndSourceAbsence')
        pins += [pin(RELEASE / (cid + '_INTENT.json')), pin(RELEASE / (cid + '_RELEASED.json'))]
    return {'all_nine_actual_sources_absent': True, 'released_files': 9,
        'released_logical_bytes': sum(SAM_BYTES), 'receipt_pins': pins,
        'R_acceptance_checked_by_frozen_release_source_not_live_Linux_R_mount': True,
        'all_nine_acceptance_sha256': intent['all_nine_acceptance_sha256']}


def proc(pid):
    p = Path('/proc') / str(pid)
    try:
        text = (p / 'stat').read_text(); fields = text[text.rfind(')')+2:].split()
        return {'parent': int(fields[1]), 'start_ticks': int(fields[19]),
                'cmd': (p / 'cmdline').read_bytes().decode('utf-8', 'strict')}
    except (FileNotFoundError, ProcessLookupError):
        return None


def g128_controller(plan):
    fields = {'key', 'first', 'last', 'worker', 'job', 'job_sha256'}
    need(len(plan['tasks']) == 918 and len({t['worker'] for t in plan['tasks']}) == 918 and
         all(set(t) == fields and Path(t['worker']).name == t['worker'] and
             integer(t['first'], 1) and integer(t['last'], t['first']) for t in plan['tasks']),
         'OriginalG128TaskSchemaAndUniqueWorkers')
    controller = proc(522800)
    if controller is not None:
        need(controller['start_ticks'] == 2556696464 and str(G128) in controller['cmd'], 'G128ControllerIdentity')
        return {'state': 'BOUND_ORIGINAL_CONTROLLER_LIVE', 'pid': 522800, 'birth_ticks': 2556696464}
    folder = G128 / 'dispatch/main_02'
    result = read(folder / 'DISPATCH_EXIT.json')
    started = read(folder / 'STARTED.json')
    need(started['pid'] == 522800 and result['status'] == 'ALL_REQUESTED_SHARDS_EXITED_ZERO_NEEDS_INDEPENDENT_CHECK'
         and result['started_shards'] == result['finished_shards'] == result['planned_shards'] == len(plan['tasks']) == 918,
         'MissingG128ControllerWithoutCompleteTerminal')
    for task in plan['tasks']:
        leaf = folder / task['worker']
        exit_code = read(leaf / 'PROCESS_EXIT.json')['exit_code']
        need(read(leaf / 'STARTED.json')['task'] == task and type(exit_code) is int and
             exit_code == 0, 'G128FullDispatchTaskTerminal')
    return {'state': 'ALL_918_ACTUAL_SHARD_TERMINALS_ZERO_NOT_SCIENCE_ACCEPTANCE',
            'terminal': pin(folder / 'DISPATCH_EXIT.json'), 'original_budget_still_reserved': True}


def allocated_file(path, inodes):
    p = clean(path); s = p.lstat()
    need(stat.S_ISREG(s.st_mode) and s.st_nlink == 1 and s.st_blocks*512 >= s.st_size and
         (s.st_dev, s.st_ino) not in inodes, 'ArchiveDiskAccountingIdentity')
    inodes.add((s.st_dev, s.st_ino))
    return s


def g128_publication_pending(path):
    """Zero credit only for the pinned producer's stable two-link publication."""
    p = clean(path)
    need(p.is_relative_to(G128 / 'allocations') and
         re.fullmatch(r'allocation_[0-9]{5}[.]rds', p.name), 'OriginalPublicationPath')
    first = p.lstat()
    if not stat.S_ISREG(first.st_mode) or first.st_nlink != 2:
        return False
    partial = clean(Path(str(p) + '.transport.partial'))
    fields = ('st_dev', 'st_ino', 'st_mode', 'st_nlink', 'st_size',
              'st_mtime_ns', 'st_ctime_ns')
    identity = lambda value: tuple(getattr(value, key) for key in fields)
    try:
        linked = partial.lstat()
        final_after, partial_after = p.lstat(), partial.lstat()
    except FileNotFoundError:
        return False  # Publication may just have completed; strict fresh stat follows.
    return (stat.S_ISREG(linked.st_mode) and linked.st_nlink == 2 and
            (first.st_dev, first.st_ino) == (linked.st_dev, linked.st_ino) and
            identity(first) == identity(linked) == identity(final_after) == identity(partial_after))


def disk_credits(decision, plan):
    inodes = set(); original = current = 0
    pending_publications = 0
    entry = pin(ROOT / 'project' / ENTRY)
    cutoff = dt.datetime.fromisoformat(decision['recorded_utc']).timestamp()
    need(len(plan['counts']) == 24 and sum(plan['counts'].values()) == 44792 and len(plan['tasks']) == 918,
         'OriginalG128CompletePlan')
    for key, count in plan['counts'].items():
        folder = clean(G128 / 'allocations' / key)
        need(folder.is_dir(), 'OriginalArchiveParent')
        for p in folder.iterdir():
            m = re.fullmatch(r'allocation_([0-9]{5})\.rds', p.name)
            if not m:
                continue
            need(1 <= int(m[1]) <= count and p.lstat().st_mtime >= cutoff,
                 'OriginalArchiveCountOrDate')
            if g128_publication_pending(p):
                pending_publications += 1
                continue  # No credit: retain the entire original remaining-growth allowance.
            s = allocated_file(p, inodes)
            need(1 <= int(m[1]) <= count and s.st_mtime >= cutoff, 'OriginalArchiveCountOrDate')
            original += s.st_size
    for accession in CELLS:
        folder = ROOT / 'derived' / accession
        if not folder.exists():
            continue
        clean(folder); need(folder.is_dir(), 'CurrentArchiveParent')
        identity = None
        complete_n = 0
        prepared_mtime = (ROOT / 'intermediate' / accession / 'PREPARED.json').stat().st_mtime
        for leaf in folder.iterdir():
            if not leaf.is_dir():
                continue  # Prefix/summary overhead stays in its separate 8 GiB reserve.
            key = leaf.name
            need(re.fullmatch('[0-3]{12}', key) and all(key.count(c) == 3 for c in '0123'), 'CurrentArchiveLabelKey')
            start = read(leaf / 'STARTED.json'); complete = read(leaf / 'COMPLETE.json')
            source = start['execution_entry']
            need(start['status'] == 'ALLOCATION_OWNED_BEFORE_NATIVE_CALL' and
                 complete['status'] == 'IMMUTABLE_COMPLETE_COMPACT_CHECKPOINT' and
                 start['allocation_key'] == complete['allocation_key'] == key and
                 re.fullmatch('[0-9a-f]{64}', start['scientific_identity_sha256']) and
                 start['scientific_identity_sha256'] == complete['scientific_identity_sha256'] and
                 source == entry,
                 'CurrentArchiveMetadataCoupling')
            if identity is None:
                identity = start['scientific_identity_sha256']
            need(identity == start['scientific_identity_sha256'], 'CurrentCellScientificIdentityChanged')
            p = leaf / 'COMPACT.rds'; s = allocated_file(p, inodes); rec = complete['archive']
            need(integer(rec['bytes']) and s.st_size == rec['bytes'] and re.fullmatch('[0-9a-f]{64}', rec['sha256']) and
                 s.st_mtime >= prepared_mtime and integer(complete['new_native_calls']) and
                 complete['new_native_calls'] in (0, 1) and type(start['new_fit_intended']) is bool and
                 int(start['new_fit_intended']) == complete['new_native_calls'], 'CurrentArchiveStatMetadata')
            current += s.st_size
            complete_n += 1
        need(complete_n <= 10020, 'CurrentCompleteSpaceArchiveCount')
    floor = math.ceil(decision['minimum_available_disk_bytes'])
    remaining = max(floor, math.ceil(decision['complete_additional_disk_scenario_and_floor_bytes']) - original)
    return {'original_publication_pending_files_zero_credit': pending_publications,
            'original_credited_logical_bytes': original, 'original_remaining_including_floor_bytes': remaining,
            'original_floor_bytes': floor, 'current_credited_logical_bytes': current,
            'credit_is_stat_and_saved_metadata_not_RDS_numeric_validation': True}


def inspect_resources(workers=3):
    """Read-only barrier inspection; resource/window waits are not identity failures."""
    calculate(workers)  # strict caller/type check before any filesystem observation
    configuration()
    release = release_evidence()
    rp, pp = pin(G128 / 'RESOURCE_DECISION.json'), pin(G128 / 'dispatch/main_02/PLAN.json')
    need(rp['sha256'] == RESOURCE_SHA and pp['sha256'] == PLAN_SHA, 'OriginalG128BudgetOrPlanIdentity')
    decision, plan = read(rp['path']), read(pp['path'])
    controller = g128_controller(plan)
    credits = disk_credits(decision, plan)
    calculation = calculate(workers, credits['original_remaining_including_floor_bytes'],
                            credits['current_credited_logical_bytes'])
    waits = set()
    own = {os.getpid(), os.getppid()}
    for p in Path('/proc').iterdir():
        if not p.name.isdigit() or int(p.name) in own:
            continue
        item = proc(int(p.name))
        if item is None:
            continue
        exe = item['cmd'].split('\0')[0]
        if 'STAR' in Path(exe).name or 'atlas' in item['cmd'].lower():
            waits.add('WAIT_ATLAS_OR_STAR_WINDOW')
        if ('/R' in exe or 'Rscript' in exe) and str(G128) not in item['cmd']:
            waits.add('WAIT_OTHER_SCIENTIFIC_R_WINDOW')
    limit = int(Path('/sys/fs/cgroup/memory.max').read_text())
    used = int(Path('/sys/fs/cgroup/memory.current').read_text())
    quota, period = map(int, Path('/sys/fs/cgroup/cpu.max').read_text().split())
    need(decision['maximum_workers_accepted'] == 12, 'OriginalG128CPUReservation')
    cpu = cpu_budget(quota, period, len(os.sched_getaffinity(0)), workers)
    waits.update(cpu['wait_categories'])
    v = os.statvfs(ROOT); free = v.f_bavail*v.f_frsize
    waits.update(resource_waits(calculation, decision['complete_memory_scenario_bytes'], limit, used, free))
    return {'schema': 'jgw03-full-resource-inspection-v1', 'utc': utc(), 'resources_admitted': not waits,
        'wait_categories': sorted(waits), 'calculation': calculation, 'disk_credits': credits,
        'free_disk_bytes': free, 'memory_limit_bytes': limit, 'memory_current_bytes': used,
        'original_G128_complete_memory_scenario_bytes': decision['complete_memory_scenario_bytes'],
        'CPU': cpu, 'original_G128_controller_evidence': controller,
        'original_budget': rp, 'original_plan': pp, 'SAM_release_evidence': release,
        'window_observation_not_future_idle_guarantee': True, 'scientific_models_launched': False}


def resource_waits(calculation, original_memory, limit, used, free):
    """Pure threshold arithmetic; does not certify an actual machine snapshot."""
    need(all(integer(x) for x in (original_memory, limit, used, free)), 'ResourceIntegerDomain')
    waits = set()
    extra = calculation['worker_memory_bytes'] + calculation['memory_floor_bytes']
    if limit-used < extra:
        waits.add('WAIT_DIRECT_MEMORY_HEADROOM')
    if limit < original_memory + extra:
        waits.add('WAIT_FULL_G128_MEMORY_SCENARIO')
    if free < calculation['required_available_disk_bytes']:
        waits.add('WAIT_FULL_DISK_SCENARIO')
    return waits


def cpu_budget(quota, period, affinity_count, workers=3):
    need(all(integer(x, 1) for x in (quota, period, affinity_count, workers)) and workers <= 3,
         'CPUQuotaAffinitySchema')
    effective = min(quota/period, affinity_count)
    required = 12 + workers + 5
    return {'cgroup_quota_CPU': quota/period, 'affinity_CPU_count': affinity_count,
        'effective_CPU_capacity': effective, 'original_G128_reserved_CPU': 12,
        'new_single_CPU_workers': workers, 'other_CPU_floor': 5, 'required_total_CPU': required,
        'wait_categories': [] if effective >= required else ['WAIT_FULL_CPU_CAPACITY']}


def actual_process(group, entry, cfg_pin, mode, accession, allocation, output=None):
    records = {name: bound(group[name]) for name in ('command', 'started', 'terminal')}
    command, started, terminal = (read(records[name]['path']) for name in ('command', 'started', 'terminal'))
    need(command['entry'] == {k: entry[k] for k in ('bytes', 'sha256')} and
         command['config'] == {k: cfg_pin[k] for k in ('bytes', 'sha256')} and
         command['mode'] == mode and command['accession'] == accession and command['allocation'] == allocation,
         'ActualCommandSourceConfigAssociation')
    argv = command['argv']
    expected = [str(RUNTIME), '--vanilla', entry['path'], cfg_pin['path']]
    expected += ['bridge', accession, str(allocation)] if mode == 'bridge' else [accession, str(allocation), output]
    need(argv == expected and integer(started['actual_R_pid'], 1) and integer(started['start_ticks'], 1) and
         started['actual_R_pid'] == terminal['actual_R_pid'] == group['actual_R_pid'] and
         terminal['actual_R_exit_code'] == group['actual_R_exit_code'] == 0 and
         started['config_sha256'] == terminal['config_sha256'] == cfg_pin['sha256'] and
         terminal['entry_sha256'] == entry['sha256'] and terminal['bridge_memory_budget_met'] is True and
         terminal['observer_failure_class_only'] is None and terminal['peak_child_RSS_bytes'] <= 3*GIB,
         'ActualProcessTerminalOrIdentity')
    return records


def six_evidence(record, cfg, cfg_pin):
    six_pin = bound(record)
    need(six_pin['path'] == str(ROOT / 'provenance/SIX_FIXED_NATIVE_COMPATIBILITY_01.json') and
         six_pin['sha256'] == SIX_SHA, 'FrozenSixSummaryIdentity')
    summary = read(six_pin['path'])
    need(summary['status'] == 'SIX_FIXED_VECTOR_CHECKS_PASS_NOT_FULL_RESOURCE_ACCEPTANCE' and
         summary['six_fixed_checks_passed'] is True and summary['new_native_calls'] == 5 and
         summary['config'] == cfg_pin and summary['original_actual_exit_code'] == 2 and
         summary['full_resources_admitted'] is False and set(summary['cells']) == set(CELLS), 'SixSummaryScope')
    original = cfg['reused_first_PC12']
    original_records = {k: bound(original[k]) for k in ('process_exit', 'process_command', 'process_started', 'native')}
    old_exit, old_command, old_start = (read(original_records[k]['path']) for k in
                                       ('process_exit', 'process_command', 'process_started'))
    need(original_records['process_exit'] == bound(summary['original_actual_exit']) and
         old_exit['actual_R_exit_code'] == 2 and old_start['actual_R_pid'] == old_exit['actual_R_pid'] and
         integer(old_start['start_ticks'], 1) and old_start['config_sha256'] == original['config']['sha256'] and
         old_command['config']['sha256'] == original['config']['sha256'] and
         old_command['entry']['sha256'] == '57e3643a06b5c17f66c46eb658d150665c9bc55a1b4cc025ceb3e21d6d2cdc8f' and
         old_command['mode'] == 'bridge' and old_command['accession'] == 'GSE289929' and old_command['allocation'] == 1,
         'OriginalPC12ActualExitTwoChain')
    entry, checker = pin(ROOT / 'project' / ENTRY), pin(ROOT / 'project' / CHECKER)
    cells, processes = {}, {}
    compact_maxima = []
    for accession in CELLS:
        items = summary['cells'][accession]
        need(set(items) == {'1', '2'}, 'SixFixedItems')
        prepared_meta = read(ROOT / 'intermediate' / accession / 'PREPARED.json')
        prep = bound(prepared_meta['prepared'])
        need(prep['path'] == str(ROOT / 'intermediate' / accession / 'PREPARED.rds') and
             prepared_meta['config_sha256'] == cfg_pin['sha256'] and
             prepared_meta['entry_sha256'] == entry['sha256'], 'PreparedAssociation')
        cell = {'prepared': prep, 'bridges': {}, 'checks': {}, 'producer_processes': {}}
        process_cell = {}; sizes = []
        for label in ('1', '2'):
            i = int(label); item = items[label]
            bridge_pin = bound(item['producer_receipt']); check_pin = bound(item['independent_receipt'])
            bridge, check = read(bridge_pin['path']), read(check_pin['path'])
            need(bridge_pin['path'] == str(ROOT / 'intermediate' / accession / ('bridge_' + label) / 'BRIDGE_COMPLETE.json') and
                 bridge['status'] == 'LINUX_NATIVE_COMPATIBILITY_PASS_V2' and bridge['compatibility_all_fields_pass'] is True and
                 bridge['config_sha256'] == check['config_sha256'] == cfg_pin['sha256'] and
                 bridge['entry_sha256'] == check['producer_entry_sha256'] == entry['sha256'] and
                 bridge['prepared_sha256'] == check['prepared_sha256'] == prep['sha256'] and
                 bridge['accession'] == check['accession'] == accession and bridge['allocation_id'] == check['allocation_id'] == i and
                 check['status'] == 'PASS_SAVED_NATIVE_COMPATIBILITY_V2_NOT_EQUIVALENCE' and
                 check['all_hard_invariants_pass'] is True and check['sources_inputs_native_unchanged'] is True and
                 check['checker_sha256'] == checker['sha256'] and check['native_sha256'] == bridge['native']['sha256'],
                 'SixNativeIndependentAssociation')
            for name in ('native', 'compact', 'checks', 'semantic_checks', 'aligned_strict_checks', 'checks_saved', 'session_info'):
                actual = bound(bridge[name])
                if name in item['products']:
                    need(actual == item['products'][name], 'SummaryProductIdentity')
            if i == 1:
                bound(bridge['domain'])
            for product in check['outputs']:
                bound(product)
            reused = accession == 'GSE289929' and i == 1
            need(bridge['new_native_calls'] == item['new_native_calls'] == (0 if reused else 1) and
                 bridge['native_generation_kind'] == item['native_generation_kind'] ==
                 ('REUSED_ORIGINAL_LINUX_NATIVE' if reused else 'NEW_LINUX_NATIVE_V2'), 'SixNativeCallAccounting')
            if reused:
                need(bound(bridge['native']) == original_records['native'] and
                     bridge['origin']['actual_original_R_exit_code'] == 2, 'OriginalNativeReuseNotRefit')
            registry = check['registry']
            need(check['this_object_full_MC_families_checked'] == item['MC_families'] == FAMILIES[accession][0] and
                 check['this_cell_source_exact_families'] == item['exact_family_definition_n'] == FAMILIES[accession][1] and
                 registry['MC_rows'] == registry['MC_unique'] == 10001 and
                 registry['companion_rows'] == registry['companion_unique'] == 20 and registry == item['registry'],
                 'FrozenCompleteFamilyAndRegistryScope')
            producer_process = actual_process(item['producer_process'], entry, cfg_pin, 'bridge', accession, i)
            checker_process = actual_process(item['checker_process'], checker, cfg_pin, 'check', accession, i,
                                             str(Path(check_pin['path']).parent))
            cell['bridges'][label] = bridge_pin
            cell['checks'][label] = {'receipt': check_pin, 'process': checker_process['terminal']}
            cell['producer_processes'][label] = producer_process['terminal']
            process_cell[label] = {'producer': producer_process, 'checker': checker_process}
            sizes.append(bridge['compact']['bytes'])
        compact_maxima.append(max(sizes)); cells[accession] = cell; processes[accession] = process_cell
    need(tuple(compact_maxima) == COMPACT_MAX, 'CompactPlanningSourceSizes')
    return cells, {'six_summary': six_pin, 'actual_processes': processes,
                   'original_PC12_exit_two': original_records, 'model_fits_executed_by_assembler': False}


def write_new(path, value):
    with Path(path).open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n'); stream.flush(); os.fsync(stream.fileno())
    return pin(path)


def assemble(request_json, attempt):
    """Explicit one-shot construction only; no dispatcher or scientific process."""
    need(re.fullmatch('[A-Za-z0-9_-]{1,80}', attempt), 'FreshAttemptName')
    request_pin = pin(request_json); request = read(request_json)
    need(request['schema'] == 'jgw03-full-admission-request-v1' and request['run_root'] == str(ROOT), 'RequestSchemaRoot')
    cfg, cfg_pin = configuration()
    need(bound(request['config']) == cfg_pin, 'RequestConfigIdentity')
    tools = request['execution_tools']
    need(set(tools) == {'driver', 'admission'}, 'ExplicitTwoExecutionTools')
    for name, record in tools.items():
        actual = bound(record)
        need(Path(actual['path']).parent == ROOT / 'execution_tools_02' and Path(actual['path']).suffix == '.py', 'ExecutionToolLocation')
    need(bound(tools['admission']) == pin(Path(__file__).resolve()) and
         tools['driver']['path'] != tools['admission']['path'], 'ActualAdmissionToolBinding')
    cells, evidence = six_evidence(request['six_checks'], cfg, cfg_pin)
    release = release_evidence(request['sam_release'])
    resource = inspect_resources(workers=3)
    need(resource['resources_admitted'] is True, 'ResourcesNotAdmittedNoAcceptanceWritten')
    need(resource['SAM_release_evidence'] == release, 'ReleaseChangedDuringAssembly')
    folder = clean(ROOT / 'provenance') / attempt
    need(not folder.exists() and not folder.is_symlink(), 'ImmutableAttemptAlreadyExists')
    folder.mkdir()
    saved_request = write_new(folder / 'REQUEST.json', request)
    resource_pin = write_new(folder / 'RESOURCE.json', dict(resource, execution_tools=tools))
    evidence_pin = write_new(folder / 'PROCESS_EVIDENCE.json', dict(evidence,
        source_request=request_pin, saved_request=saved_request, execution_tools=tools))
    acceptance = {'schema': 'jgw03-linux-full-execution-acceptance-v2', 'status': 'ACCEPTED_FOR_FULL_LINUX_REFERENCE_V2',
        'config_sha256': cfg_pin['sha256'], 'entry_sha256': ENTRY_SHA, 'original_PC12_native_actual_exit_code': 2,
        'remaining_five_producer_actual_exits_zero': True, 'six_independent_checks_actual_exits_zero': True,
        'independent_bridge_acceptance': True, 'resources_admitted': True,
        'full_MC_family_n': 89, 'full_companion_family_n': 34, 'MC_nonidentity_per_cell': 10000, 'companion_per_cell': 20,
        'cells': cells, 'independent_check': bound(request['six_checks']), 'resource_receipt': resource_pin,
        'process_evidence': evidence_pin, 'execution_tools': tools, 'request': saved_request,
        'scientific_configuration_modified': False, 'scientific_processes_launched': False,
        'planning_scenario_not_worst_case_guarantee': True, 'utc': utc()}
    # Closure: pins are re-read before the only acceptance file is published.
    configuration(); bound(request['six_checks']); release_evidence(request['sam_release'])
    for record in tools.values():
        bound(record)
    write_new(folder / 'ACCEPTANCE.json', acceptance)
    return pin(folder / 'ACCEPTANCE.json')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', nargs='?', choices=('calculate', 'inspect', 'assemble'), default='inspect')
    parser.add_argument('--workers', type=int, default=3)
    parser.add_argument('--request', type=Path)
    parser.add_argument('--attempt')
    parser.add_argument('--create-acceptance', action='store_true')
    args = parser.parse_args()
    try:
        if args.mode == 'calculate':
            value = calculate(args.workers)
        elif args.mode == 'inspect':
            value = inspect_resources(args.workers)
        else:
            need(args.create_acceptance and args.request is not None and args.attempt is not None,
                 'ExplicitAcceptanceCreationRequired')
            value = assemble(args.request, args.attempt)
        print(json.dumps(value, sort_keys=True, allow_nan=False))
        return 0
    except Exception as exc:
        print(json.dumps({'status': 'ADMISSION_STOP_NO_MODEL_NO_RETRY',
                          'exception_class_only': type(exc).__name__, 'raw_error_persisted': False}))
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
