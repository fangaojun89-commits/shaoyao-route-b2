"""Separate 21 October round54 bounded after G128 completion cache operation, never delete/read payload or launch R.

This new fixed pool is disjoint from every prior processed batch.
The extra 24-GiB advice cap applies to this operation once, not per retry.

Default is metadata/residency inspection. --execute changes ONLY cache state
for fixed closed raw_all30 IDs3..10000 and raw_predefined29 IDs3..2500, excluding fixed sentinels and prior pool7501..7900. All earlier processed batch paths are excluded. It is not a read-only operation.
Original17/40-GiB admission rules, processes and scientific files stay unchanged.
No global cache operation, fsync, configuration write or payload hash is used.
Advice may initiate ordinary asynchronous writeback; it does not certify that
every target page was clean, and previously issued advice is not rolled back.
The live finalizer can read these files later; only present consumers are checked.
Future reload I/O is possible. No claim of permanent nonuse is made.
"""

import argparse
import csv
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import socket
import stat
import subprocess


ROOT = Path('/root/autodl-tmp/gse128623_label_reference_20261015r01')
JGW = Path('/root/autodl-tmp/JGW03_full_reference_linux_20261017r03/provenance/full_reference_02')
CGROUP = Path('/sys/fs/cgroup')
CELLS = ('raw_all30', 'raw_predefined29')
PAGE = os.sysconf('SC_PAGE_SIZE')
MAX_ADVISED = 24 * 1024**3
TARGET_HEADROOM = 18 * 1024**3


def need(ok, label):
    if not ok:
        raise RuntimeError(label)


def read(path):
    return json.loads(path.read_text())


def signature(value):
    return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)


def anchored_open(path):
    """Reject symlinks at every ancestor, retaining a directory-FD anchor."""
    directory = os.open('/', os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        for component in path.parts[1:-1]:
            child = os.open(component, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                            dir_fd=directory)
            os.close(directory)
            directory = child
        descriptor = os.open(path.name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=directory)
    finally:
        os.close(directory)
    try:
        need(os.readlink('/proc/self/fd/' + str(descriptor)) == str(path), 'FDPathIdentity')
    except BaseException:
        os.close(descriptor)
        raise
    return descriptor


def completed_g128():
    data = (ROOT / 'dispatch/main_02/DISPATCH_EXIT.json').read_bytes()
    need(hashlib.sha256(data).hexdigest() ==
         '4ec5049eb50086ff1cfe14ffda2e5a2f8cf6d93e16afa5ed286786e0d3aedf11',
         'ExactG128ActualTerminal')
    terminal = json.loads(data)
    need(terminal['status'] == 'ALL_REQUESTED_SHARDS_EXITED_ZERO_NEEDS_INDEPENDENT_CHECK'
         and terminal['whole_run_independently_accepted'] is False
         and all(type(terminal[k]) is int and terminal[k] == 918
                 for k in ('planned_shards', 'started_shards', 'finished_shards')),
         'G128Actual918CompletedNotScientificAcceptance')
    return '4ec5049eb50086ff1cfe14ffda2e5a2f8cf6d93e16afa5ed286786e0d3aedf11'


def snapshot():
    counters = dict(line.split() for line in (CGROUP / 'memory.stat').read_text().splitlines())
    current = int((CGROUP / 'memory.current').read_text())
    maximum = int((CGROUP / 'memory.max').read_text())
    processes = {}
    for pid, birth in ((659622, 2576253226), (896255, 2602926060)):
        fields = (Path('/proc') / str(pid) / 'stat').read_text().rsplit(') ', 1)[1].split()
        need(int(fields[19]) == birth and fields[0] not in ('Z', 'X'), 'BoundJGWAndFinalizerContinuity')
        processes[str(pid)] = {'birth': birth, 'state': fields[0]}
    status = read(JGW / 'STATUS.json')
    need(status.get('owner_pid') == 659622, 'ExactOriginalControllerOwner')
    if status.get('status') == 'WAIT_RESOURCE_WINDOW_NO_NATIVE_PROCESS_RESTART':
        need(status.get('next_round') == 54, 'ExactRound54OriginalControllerWait')
    return {'utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'current': current, 'maximum': maximum, 'headroom': maximum - current,
            'file': int(counters['file']), 'anon': int(counters['anon']),
            'dirty': int(counters['file_dirty']), 'writeback': int(counters['file_writeback']),
            'pressure': (CGROUP / 'memory.pressure').read_text().splitlines(),
            'bound_live_processes': processes, 'completed_G128_terminal_sha256': completed_g128(), 'JGW_status': status['status']}


def no_consumers(keys):
    mapped = {(os.major(dev), os.minor(dev), inode) for dev, inode in keys}
    checked = 0
    for proc in Path('/proc').iterdir():
        if not proc.name.isdigit() or int(proc.name) == os.getpid():
            continue
        try:
            for descriptor in (proc / 'fd').iterdir():
                try:
                    value = descriptor.stat()
                except FileNotFoundError:
                    continue
                need((value.st_dev, value.st_ino) not in keys, 'CurrentFDConsumer')
            lines = (proc / 'maps').read_text().splitlines()
        except FileNotFoundError:
            continue  # A process/descriptor that genuinely disappeared is not an unreadable live one.
        for line in lines:
            fields = line.split(None, 5)
            if len(fields) >= 5:
                major, minor = fields[3].split(':')
                need((int(major, 16), int(minor, 16), int(fields[4])) not in mapped,
                     'CurrentMappedConsumer')
        checked += 1
    # Permission errors are deliberately not suppressed: do not advise blind.
    return checked


def selected(aid, cell='raw_all30'):
    if cell == 'raw_predefined29':
        return 3 <= aid <= 2500 and aid not in (101, 1001)
    return 3 <= aid <= 10000 and aid not in (101, 1001, 5001) and not 7501 <= aid <= 7900


def candidates():
    need(socket.gethostname() == 'autodl-container-62c8pzy536-50fb405a', 'ExactCurrentInstance')
    filesystem = subprocess.run(['findmnt', '-T', str(ROOT), '-n', '-o', 'FSTYPE'],
                                capture_output=True, text=True, check=True).stdout.strip()
    need(filesystem == 'xfs', 'ReviewedDiskBackedFilesystem')
    plan = read(ROOT / 'dispatch/main_02/PLAN.json')
    config = ROOT / 'DISPATCH_CONFIG_v2.json'
    need(hashlib.sha256(config.read_bytes()).hexdigest() == plan['config_sha256'] ==
         '847a049b410828b0b920a9f355ef93e34d9d51d30b68cb4559bf1059f20d44a8', 'FrozenPlanConfig')
    controller = ROOT / 'source/dispatch_GSE128623_reference_20261015_v2.py'
    need(hashlib.sha256(controller.read_bytes()).hexdigest() == plan['controller_sha256'] ==
         'ba55d2ec1d7aa742ccf0c310802271d8a4e72d224db235c90a75de2db5434a40', 'FrozenController')
    completed_g128()
    need(len(plan['tasks']) == 918 and len({t['worker'] for t in plan['tasks']}) == 918,
         'ActualCompletePlan')
    for task in plan['tasks']:
        folder = ROOT / 'dispatch/main_02' / task['worker']
        ended, started = read(folder / 'PROCESS_EXIT.json'), read(folder / 'STARTED.json')
        need(type(ended['exit_code']) is int and ended['exit_code'] == 0
             and started['task'] == task, 'All918ActualZeroExitTaskBindings')
    finalizer = ROOT / 'finalization_preparation_20261015_r01/finalize_GSE128623_label_reference_20261015.R'
    contract = ROOT / 'full_finalization_execution_20261020_r01/CONTRACT.json'
    need(hashlib.sha256(finalizer.read_bytes()).hexdigest() ==
         '2d3298699771c82a983edf4f599c0cfaeb0e210f52f824f95934d9e0220c2240'
         and hashlib.sha256(contract.read_bytes()).hexdigest() ==
         '35297c8beb13e50e571888f488f7e6a1419606e682cdbf3af00f86c5f325bbb5',
         'ExactAlreadyRunningFinalizerAndContract')
    argv = Path('/proc/896255/cmdline').read_bytes().decode().rstrip(chr(0)).split(chr(0))
    need('--file=' + str(finalizer) in argv and '--args' in argv
         and argv[argv.index('--args') + 1:] == [str(ROOT), str(contract),
             '35297c8beb13e50e571888f488f7e6a1419606e682cdbf3af00f86c5f325bbb5',
             str(ROOT / 'full_reference_finalization_20261020_r01')],
         'ActualFinalizerExactArguments')
    wanted = {}
    now = datetime.datetime.now(datetime.timezone.utc).timestamp()
    for task in plan['tasks']:
        if task['key'].split('/')[0] not in CELLS or not task['key'].endswith('/six_group') or not any(selected(aid, task['key'].split('/')[0]) for aid in range(task['first'], task['last'] + 1)):
            continue
        folder = ROOT / 'dispatch/main_02' / task['worker']
        ended, started = read(folder / 'PROCESS_EXIT.json'), read(folder / 'STARTED.json')
        need(type(ended['exit_code']) is int and ended['exit_code'] == 0 and started['task'] == task,
             'ExactClosedTaskBinding')
        finished = datetime.datetime.fromisoformat(ended['finished_utc']).timestamp()
        need(now - finished >= 3600, 'ClosedForAtLeastOneHour')
        with (ROOT / 'workers' / task['worker'] / 'ARCHIVE_MANIFEST.tsv').open(newline='') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        need([int(row['allocation_id']) for row in rows] == list(range(task['first'], task['last'] + 1)),
             'CompleteFixedTaskManifest')
        cell = task['key'].split('/')[0]
        need(cell in CELLS, 'FixedCell')
        for row in rows:
            aid = int(row['allocation_id'])
            if selected(aid, cell):
                path = ROOT / 'allocations' / cell / 'six_group' / ('allocation_%05d.rds' % aid)
                need(str(path) == row['file'] and path.resolve(strict=True) == path, 'ExactClosedOutputPath')
                value = path.lstat()
                listed_size = float(row['bytes'])
                need(stat.S_ISREG(value.st_mode) and value.st_nlink == 1 and value.st_uid == os.geteuid()
                     and math.isfinite(listed_size) and listed_size.is_integer()
                     and value.st_size == listed_size, 'RegularOwnedSingleLinkManifestSize')
                wanted[(cell, aid)] = {'path': path, 'signature': signature(value),
                                       'manifest_sha256': row['sha256']}
    need(len(wanted) == 12091, 'Fixed12091PoolBeforePriorProcessedPathExclusion')
    previous = ROOT / 'closed_output_cache_advice_20261020_r01' / 'EXECUTION_STDOUT.json'
    previous_bytes = previous.read_bytes()
    need(hashlib.sha256(previous_bytes).hexdigest() ==
         'a89b2db4180b03336260f9d20faf732a92dbf005609d5d400e5e52a4a7222d98', 'PreviousOperationExactReceipt')
    previous_record = json.loads(previous_bytes)
    need(previous_record['schema'] == 'GSE128623_CLOSED_OUTPUT_CACHE_ADVICE_20261020/v1'
         and previous_record['source_sha256'] ==
         '5c8d2ac82ce4e50fe377e43e1fadf47804a8698989a0270fb7eaf8b08e443f89'
         and previous_record['execute'] is True
         and previous_record['limit_bytes'] == 20 * 1024**3, 'PreviousOperationExactIdentity')
    need(previous_record['status'] == 'BOUNDED_CACHE_STATE_OPERATION_COMPLETE', 'PreviousOperationClosed')
    previous_paths = {item['path'] for batch in previous_record['batches'] for item in batch['files']}
    prior34 = ROOT / 'closed_output_cache_advice_round34_20261020_r01' / 'EXECUTION_STDOUT.json'
    prior34_bytes = prior34.read_bytes()
    need(hashlib.sha256(prior34_bytes).hexdigest() ==
         '82f55a5d039befcdaf3f46ef9c538fa91b5f61b37c66c2dacd8199bcb2181a0d', 'Round34ExactReceipt')
    prior34_record = json.loads(prior34_bytes)
    need(prior34_record['schema'] == 'GSE128623_CLOSED_OUTPUT_CACHE_ADVICE_ROUND34_20261020/v1'
         and prior34_record['source_sha256'] ==
         'bc13b41999429ea791ccc5c7197395825b143f0be7cef64279f92f958efd236f'
         and prior34_record['execute'] is True
         and prior34_record['limit_bytes'] == 12 * 1024**3
         and prior34_record['status'] == 'BOUNDED_CACHE_STATE_OPERATION_COMPLETE', 'Round34ExactIdentity')
    previous_paths.update(item['path'] for batch in prior34_record['batches'] for item in batch['files'])
    prior35 = ROOT / 'closed_output_cache_advice_round35_20261020_r01' / 'EXECUTION_STDOUT.json'
    prior35_bytes = prior35.read_bytes()
    need(hashlib.sha256(prior35_bytes).hexdigest() ==
         '881f8c017a3135a7b80ddce814c563f6eddd98250338dac0fbbe9e23a8446446', 'Round35ExactReceipt')
    prior35_record = json.loads(prior35_bytes)
    need(prior35_record['schema'] == 'GSE128623_CLOSED_OUTPUT_CACHE_ADVICE_ROUND35_20261020/v1'
         and prior35_record['source_sha256'] ==
         'b7dda09643f2892b17e461b058cd9150f6e3af82c95fbd01e39ef784be707060'
         and prior35_record['execute'] is True
         and prior35_record['limit_bytes'] == 4 * 1024**3
         and prior35_record['status'] == 'BOUNDED_CACHE_STATE_OPERATION_COMPLETE', 'Round35ExactIdentity')
    previous_paths.update(item['path'] for batch in prior35_record['batches'] for item in batch['files'])
    prior37 = ROOT / 'closed_output_cache_advice_round37_20261020_r01' / 'EXECUTION_STDOUT.json'
    prior37_bytes = prior37.read_bytes()
    need(hashlib.sha256(prior37_bytes).hexdigest() ==
         '35ff7f9031852642d8aa1e3747b9e5c261fabc836daf910ccea5b80c61ecb3b9', 'Round37ExactReceipt')
    prior37_record = json.loads(prior37_bytes)
    need(prior37_record['schema'] == 'GSE128623_CLOSED_OUTPUT_CACHE_ADVICE_ROUND37_20261020/v1'
         and prior37_record['source_sha256'] ==
         '5b89956c7bd3bac5d97c3b0ab447618a20406c2e4f32eff1df83f54a8daa10f7'
         and prior37_record['execute'] is True
         and prior37_record['limit_bytes'] == 8 * 1024**3
         and prior37_record['status'] == 'BOUNDED_CACHE_STATE_OPERATION_COMPLETE', 'Round37ExactIdentity')
    previous_paths.update(item['path'] for batch in prior37_record['batches'] for item in batch['files'])
    prior39 = ROOT / 'closed_output_cache_advice_round39_20261020_r01' / 'EXECUTION_STDOUT.json'
    prior39_bytes = prior39.read_bytes()
    need(hashlib.sha256(prior39_bytes).hexdigest() ==
         '6b1e6e2138e6e62dbdd0ad49c2db3401a135724733468ae64ec4a895e81f44fa', 'Round39ExactReceipt')
    prior39_record = json.loads(prior39_bytes)
    need(prior39_record['schema'] == 'GSE128623_CLOSED_OUTPUT_CACHE_ADVICE_ROUND39_20261020/v1'
         and prior39_record['source_sha256'] ==
         '64556af01b27cf9dbc75ae6bda006396676f8e2001c6b67d4645274d63834c81'
         and prior39_record['execute'] is True
         and prior39_record['limit_bytes'] == 6 * 1024**3
         and prior39_record['status'] == 'BOUNDED_CACHE_STATE_OPERATION_COMPLETE', 'Round39ExactIdentity')
    previous_paths.update(item['path'] for batch in prior39_record['batches'] for item in batch['files'])
    prior41 = ROOT / 'closed_output_cache_advice_round41_20261020_r01' / 'EXECUTION_STDOUT.json'
    prior41_bytes = prior41.read_bytes()
    need(hashlib.sha256(prior41_bytes).hexdigest() ==
         'a33e2897d358f9dbf9ce98f2ce64a17be4e1e315c407c4361f53f57e6f8004f6', 'Round41ExactReceipt')
    prior41_record = json.loads(prior41_bytes)
    need(prior41_record['schema'] == 'GSE128623_CLOSED_OUTPUT_CACHE_ADVICE_ROUND41_20261020/v1'
         and prior41_record['source_sha256'] ==
         'e11b3ff82033adcbee1644568be020e6a5fea8dba76f435dbe32420b49070bfa'
         and prior41_record['execute'] is True
         and prior41_record['limit_bytes'] == 6 * 1024**3
         and prior41_record['status'] == 'BOUNDED_CACHE_STATE_OPERATION_COMPLETE', 'Round41ExactIdentity')
    previous_paths.update(item['path'] for batch in prior41_record['batches'] for item in batch['files'])
    prior46 = ROOT / 'closed_output_cache_advice_round46_20261020_r01' / 'EXECUTION_STDOUT.json'
    prior46_bytes = prior46.read_bytes()
    need(hashlib.sha256(prior46_bytes).hexdigest() ==
         '12977378f787101ea8e1f3ad589699721fee7630851e222da6569dd4316d2bfd', 'Round46ExactReceipt')
    prior46_record = json.loads(prior46_bytes)
    need(prior46_record['schema'] == 'GSE128623_CLOSED_OUTPUT_CACHE_ADVICE_ROUND46_20261020/v1'
         and prior46_record['source_sha256'] ==
         '013038ca7dc490a677c7239c892ea784bcb7f177c0d74cb267100eeeeb4c003b'
         and prior46_record['execute'] is True
         and prior46_record['limit_bytes'] == 8 * 1024**3
         and prior46_record['status'] == 'BOUNDED_CACHE_STATE_OPERATION_COMPLETE', 'Round46ExactIdentity')
    previous_paths.update(item['path'] for batch in prior46_record['batches'] for item in batch['files'])
    prior19 = ROOT / 'closed_output_cache_advice_20261019_r01/advise_GSE128623_closed_output_cache_20261019.py'
    need(hashlib.sha256(prior19.read_bytes()).hexdigest() ==
         'a9c5a7382a9b6126a9f97458f4fdaccc694266406bd373ef4a33939162815737',
         'Prior19FixedSourceOnly7501Through7900Disjoint')
    prior47 = ROOT / 'closed_output_cache_advice_round47_20261020_r01' / 'EXECUTION_STDOUT.json'
    prior47_bytes = prior47.read_bytes()
    need(hashlib.sha256(prior47_bytes).hexdigest() ==
         '688c21a47698998fcfc5779cf45faa2c0a9457a9ce5b511b095381995475d929', 'Round47ExactReceipt')
    prior47_record = json.loads(prior47_bytes)
    need(prior47_record['schema'] == 'GSE128623_CLOSED_OUTPUT_CACHE_ADVICE_ROUND47_20261020/v1'
         and prior47_record['source_sha256'] ==
         '20531c01383726691ebe62f88ba153ed61f4c468412c8fd9bfc5560b13dbd1e4'
         and prior47_record['execute'] is True
         and prior47_record['limit_bytes'] == 16 * 1024**3
         and prior47_record['status'] == 'BOUNDED_CACHE_STATE_OPERATION_COMPLETE', 'Round47ExactIdentity')
    previous_paths.update(item['path'] for batch in prior47_record['batches'] for item in batch['files'])
    prior48 = ROOT / 'closed_output_cache_advice_round48_20261021_r01' / 'EXECUTION_STDOUT.json'
    prior48_bytes = prior48.read_bytes()
    need(hashlib.sha256(prior48_bytes).hexdigest() ==
         'ab4846b8aacdb26c26968ecddb3c66e18f726d1632af21b2f71bc2b0408565fa', 'Round48ExactReceipt')
    prior48_record = json.loads(prior48_bytes)
    need(prior48_record['schema'] == 'GSE128623_CLOSED_OUTPUT_CACHE_ADVICE_ROUND48_20261021/v1'
         and prior48_record['source_sha256'] ==
         'd19cb39f5823248bd35716126ddf7b81353135ac9b0613da0eb1b2e312900b64'
         and prior48_record['execute'] is True
         and prior48_record['limit_bytes'] == 16 * 1024**3
         and prior48_record['status'] == 'BOUNDED_CACHE_STATE_OPERATION_COMPLETE', 'Round48ExactIdentity')
    previous_paths.update(item['path'] for batch in prior48_record['batches'] for item in batch['files'])
    prior49 = ROOT / 'closed_output_cache_advice_round49_20261021_r01' / 'EXECUTION_STDOUT.json'
    prior49_bytes = prior49.read_bytes()
    need(hashlib.sha256(prior49_bytes).hexdigest() ==
         '1f075e488b5251b2d6e4c8d9d7ca43c6955174199c104297e318fcfadf7f4671', 'Round49ExactReceipt')
    prior49_record = json.loads(prior49_bytes)
    need(prior49_record['schema'] == 'GSE128623_CLOSED_OUTPUT_CACHE_ADVICE_ROUND49_20261021/v1'
         and prior49_record['source_sha256'] ==
         '5d6609d0fa5f70f91ff149fe3cc82349b35e013da8a7b7085a9340ed5ad94514'
         and prior49_record['execute'] is True
         and prior49_record['limit_bytes'] == 16 * 1024**3
         and prior49_record['status'] == 'BOUNDED_CACHE_STATE_OPERATION_COMPLETE', 'Round49ExactIdentity')
    previous_paths.update(item['path'] for batch in prior49_record['batches'] for item in batch['files'])
    prior53 = ROOT / 'closed_output_cache_advice_round53_20261021_r01' / 'EXECUTION_STDOUT.json'
    prior53_bytes = prior53.read_bytes()
    need(hashlib.sha256(prior53_bytes).hexdigest() ==
         'ba17bf6f4c8ba0ace620a2ff84f94211e36aa1eb72384ddfb53371cc49021288', 'Round53ExactReceipt')
    prior53_record = json.loads(prior53_bytes)
    need(prior53_record['schema'] == 'GSE128623_CLOSED_OUTPUT_CACHE_ADVICE_ROUND53_20261021/v1'
         and prior53_record['source_sha256'] ==
         'ab46f955d9785802092785c4952ceceafd967fa2500d4f513afae83d13889973'
         and prior53_record['execute'] is True
         and prior53_record['limit_bytes'] == 16 * 1024**3
         and prior53_record['status'] == 'BOUNDED_CACHE_STATE_OPERATION_COMPLETE', 'Round53ExactIdentity')
    previous_paths.update(item['path'] for batch in prior53_record['batches'] for item in batch['files'])
    # Raw-all-30's producer/independent phase is already closed. Other targets
    # are completed model outputs; present consumers are checked before advice.
    # Raw-predefined-29 may be reread later by the original finalizer.
    output = ROOT / 'full_reference_finalization_20261020_r01/raw_all30__six_group'
    need(all((output / name).is_file() for name in (
        'PRODUCER_TAIL.rds', 'INDEPENDENT_ARITHMETIC.rds',
        'TRANSPORT_LINEAGE.tsv', 'COMPACT_DIAGNOSTICS_SCOPE.tsv',
        'SELECTED_OUTER_MC_SIDECAR.rds')), 'RawAll30AlreadyWrittenByOriginalFinalizer')
    wanted = {key: item for key, item in wanted.items() if str(item['path']) not in previous_paths}
    need(len(wanted) > 0, 'NonemptyNeverPreviouslyProcessedPool')
    need(not ({str(item['path']) for item in wanted.values()} & previous_paths), 'NoPriorProcessedBatchOverlap')
    return [wanted[key] for key in sorted(wanted, key=lambda key: (CELLS.index(key[0]), key[1]))]


def validate_residency(rows, batch, fd_paths):
    need(isinstance(rows, list) and len(rows) == len(batch), 'CompleteResidencyQuery')
    for row, item, fd_path in zip(rows, batch, fd_paths):
        need(isinstance(row, dict) and row.get('file') == fd_path, 'ResidencyHeldFDPath')
        need(all(type(row.get(key)) is int for key in ('size', 'res', 'pages')),
             'ResidencyIntegerFields')
        size = item['signature'][2]
        need(row['size'] == size and 0 <= row['pages'] <= (size + PAGE - 1) // PAGE
             and row['res'] == row['pages'] * PAGE, 'ResidencyHeldFDSizeAndResidentPages')


def run(execute, limit, result):
    need(type(limit) is int and 0 < limit <= MAX_ADVISED, 'FixedMaximumTwentyFourGiB')
    pool = candidates()
    result.update(kernel=os.uname().release, mount='xfs', candidate_files=len(pool),
                  before=snapshot(), batches=[], advised_bytes=0, advised_files=0)
    for offset in range(0, len(pool), 32):
        state = snapshot()
        if state['JGW_status'] != 'WAIT_RESOURCE_WINDOW_NO_NATIVE_PROCESS_RESTART':
            result['stop_reason'] = 'OriginalJGWControllerLeftResourceWait'
            break
        if state['headroom'] >= TARGET_HEADROOM:
            result['stop_reason'] = 'MeasuredEighteenGiBHeadroomReached'
            break
        if result['advised_bytes'] >= limit:
            result['stop_reason'] = 'PresetAdviceVolumeReached'
            break
        need(state['dirty'] + state['writeback'] <= 16 * 1024**2, 'NoLargeDirtyWritebackWindow')
        batch = pool[offset:offset + 32]
        descriptors = []
        try:
            for item in batch:
                descriptor = anchored_open(item['path'])
                descriptors.append(descriptor)
                need(signature(os.fstat(descriptor)) == item['signature'], 'HeldFDIdentityBefore')
            keys = {(os.fstat(fd).st_dev, os.fstat(fd).st_ino) for fd in descriptors}
            checked = no_consumers(keys)
            fd_paths = ['/proc/%d/fd/%d' % (os.getpid(), fd) for fd in descriptors]
            arguments = ['/usr/bin/fincore', '--json', '--bytes', '--output', 'RES,PAGES,SIZE,FILE']
            residency = json.loads(subprocess.run(arguments + fd_paths, capture_output=True,
                                                  text=True, check=True).stdout)['fincore']
            validate_residency(residency, batch, fd_paths)
            advised = 0
            for item, fd, resident, fd_path in zip(batch, descriptors, residency, fd_paths):
                length = min(item['signature'][2] // PAGE * PAGE,
                             (limit - result['advised_bytes']) // PAGE * PAGE)
                if execute and resident['res'] > 0 and length > 0:
                    os.posix_fadvise(fd, 0, length, os.POSIX_FADV_DONTNEED)
                    result['advised_bytes'] += length
                    result['advised_files'] += 1
                    advised += length
                need(signature(os.fstat(fd)) == item['signature'] == signature(item['path'].lstat()),
                     'StableFileIdentityAfter')
            after_residency = json.loads(subprocess.run(arguments + fd_paths, capture_output=True,
                                                        text=True, check=True).stdout)['fincore']
            validate_residency(after_residency, batch, fd_paths)
            result['batches'].append({'files': [{'path': str(item['path']),
                'signature': item['signature'], 'manifest_sha256': item['manifest_sha256']}
                for item in batch], 'visible_processes_checked': checked,
                'resident_before': sum(row['res'] for row in residency),
                'resident_after': sum(row['res'] for row in after_residency),
                'advised_bytes': advised, 'cgroup_after': snapshot()})
        finally:
            for descriptor in descriptors:
                os.close(descriptor)
        if not execute:
            result['stop_reason'] = 'ReadOnlyFirst32ResidencyInspection'
            break
    result['after'] = snapshot()
    result['status'] = 'BOUNDED_CACHE_STATE_OPERATION_COMPLETE' if execute else 'READ_ONLY_INSPECTION_COMPLETE'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--limit-bytes', type=int, default=256 * 1024**2)
    args = parser.parse_args()
    result = {'schema': 'GSE128623_CLOSED_OUTPUT_CACHE_ADVICE_ROUND54_20261021/v1', 'execute': args.execute,
              'limit_bytes': args.limit_bytes, 'helper_pid': os.getpid(),
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'scientific_payload_read': False, 'scientific_content_hash_rechecked': False,
              'scientific_files_opened_readonly_only': True, 'model_processes_launched': 0,
              'global_cache_operations': 0, 'configuration_or_admission_rules_changed': False,
              'cache_mutation_is_not_read_only': args.execute,
              'system_wide_future_nonuse_or_memory_safety_certified': False}
    try:
        run(args.execute, args.limit_bytes, result)
        code = 0
    except Exception as error:
        result.update(status='STOP_RETAINED_NO_RETRY', exception_class_only=type(error).__name__,
                      diagnostic_sha256=hashlib.sha256(str(error).encode()).hexdigest(),
                      prior_cache_advice_not_rolled_back=result.get('advised_files', 0) > 0,
                      scientific_files_retained=True)
        code = 1
    print(json.dumps(result))
    return code


if __name__ == '__main__':
    raise SystemExit(main())
