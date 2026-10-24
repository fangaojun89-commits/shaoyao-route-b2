"""One fixed remaining library through the unchanged original library function.

Import/default CLI performs no file or scientific IO. Execution never invokes
the old remaining CLI, clears cache, retries, changes parameters or schedules
other jobs. An external real parent owns process exits and boundary admission.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import socket
import sys
import types

ROOT = Path('/root/autodl-tmp/GSE205906_matched_STAR_Ensembl101_20261015r01')
SOURCE = ROOT/'source_staging_20261015_01/run_gse205906_matched_STAR_20261015.py'
SOURCE_SHA = '7fe64fcfcee554e1500f60da2cd7b4a3ce8a84fcd4e46b84c3c7e14ec1fd907d'
SOURCE_BYTES = 34422
HOST = 'autodl-container-62c8pzy536-50fb405a'
LARGEST = 'SRR19631202'
REMAINING = tuple('SRR'+str(i) for i in range(19631207, 19631198, -1) if i != 19631202)


class GuardError(ValueError):
    def __init__(self, category):
        self.category = category
        super().__init__(category)


def need(ok, category):
    if not ok:
        raise GuardError(category)


def load_original():
    need(SOURCE.is_file() and not SOURCE.is_symlink(), 'FixedRegularOriginalSource')
    raw = SOURCE.read_bytes()
    need(len(raw) == SOURCE_BYTES and hashlib.sha256(raw).hexdigest() == SOURCE_SHA,
         'FixedOriginalSourceBytesAndSHA')
    module = types.ModuleType('unchanged_gse205906_STAR_library')
    module.__file__ = str(SOURCE)
    exec(compile(raw, str(SOURCE), 'exec'), module.__dict__)
    need(module.RUN == ROOT and module.META == ROOT/'metadata_preparation_01/METADATA_AND_RESERVOIR_BINDING.json'
         and module.INDEX == ROOT/'index_01/index', 'FixedOriginalRootAndInputs')
    return module


def validate_completed(source, library, folder):
    done = source.completed(folder, SOURCE_SHA)
    need(done['SRA_run'] == library['SRA_run'] and
         done['input_pairs'] == library['registered_pairs'], 'CompletedLibraryIdentity')
    return done


def result(library, reused, projection=None, resources=None):
    return dict(schema='GSE205906_STAR_ONE_REMAINING_LIBRARY_FRONTEND/v1',
                frontend_source_path=str(Path(__file__)),
                status='REUSED_ORIGINAL_LIBRARY_COMPLETE_NO_NATIVE_CALL' if reused else
                'ORIGINAL_LIBRARY_FUNCTION_RETURNED_COMPLETE', library=library,
                original_source=dict(path=str(SOURCE), bytes=SOURCE_BYTES, sha256=SOURCE_SHA),
                original_CLI_execute_invoked=False, original_library_function_calls=0 if reused else 1,
                full_pending_projection_bytes=projection, resource_before=resources,
                new_models=0, automatic_retry=False, cache_advice=False,
                actual_parent_exit_not_created_here=True)


def execute(args):
    need(sys.platform.startswith('linux') and socket.gethostname() == HOST, 'FixedLinuxInstance')
    need(args.library in REMAINING, 'FixedEightRemainingLibraries')
    need(ROOT.is_dir() and not ROOT.is_symlink() and not (ROOT/'libraries').is_symlink(),
         'FixedIndependentRoot')
    source = load_original()
    source.unchanged(SOURCE_SHA)
    meta = source.metadata()
    remaining = [x for x in meta['libraries'] if x['SRA_run'] != LARGEST]
    need(tuple(x['SRA_run'] for x in remaining) == REMAINING, 'FrozenRemainingMetadataOrder')
    selected = next(x for x in remaining if x['SRA_run'] == args.library)
    target = ROOT/'libraries'/args.library
    need(not target.is_symlink(), 'OrdinaryLibraryFolder')
    if (target/'COMPLETE.json').exists():
        validate_completed(source, selected, target)
        source.unchanged(SOURCE_SHA)
        return result(args.library, True)
    need(args.quiet_window_confirmed and type(args.reserved_growth_bytes) is int and
         args.reserved_growth_bytes >= 0 and type(args.disk_budget_bytes) is int and
         args.disk_budget_bytes > 0, 'ExplicitOriginalResourceContract')
    largest = next(x for x in meta['libraries'] if x['SRA_run'] == LARGEST)
    admission = source.completed(ROOT/'reference_admission_01', SOURCE_SHA)
    for tool in admission['tools'].values():
        for pin in [tool['binary'], *tool['linked_libraries']]:
            need(source.pin(Path(pin['path'])) == pin, 'AdmittedToolOrLinkedLibraryDrift')
    reference = admission['reference']
    for pin in reference['reference_stat']:
        need(source.stat_file(Path(pin['path'])) == pin, 'AdmittedFullReferenceStatDrift')
    index = source.completed(ROOT/'index_01', SOURCE_SHA)
    for pin in index['index_stat']:
        need(source.stat_file(Path(pin['path'])) == pin, 'AdmittedIndexStatDrift')
    source.completed(ROOT/'reservoir_01', SOURCE_SHA)
    large_done = validate_completed(source, largest, ROOT/'libraries'/LARGEST)
    pending = []
    for library in remaining:
        folder = ROOT/'libraries'/library['SRA_run']
        need(not folder.is_symlink(), 'OrdinaryLibraryFolder')
        if (folder/'COMPLETE.json').exists():
            validate_completed(source, library, folder)
        else:
            need(not folder.exists(), 'PartialOrFailedLibraryPreservedNoRetry')
            pending.append(library)
    projection = math.ceil(large_done['retained_stage_bytes'] *
        sum(x['registered_pairs'] for x in pending) / largest['registered_pairs'] +
        max(0, large_done['peak_stage_bytes_observed'] - large_done['retained_stage_bytes']))
    need(args.disk_budget_bytes >= projection, 'FullPendingRemainingMeasuredProjection')
    resources = source.resources(args)
    source.unchanged(SOURCE_SHA)
    need(pending and args.library == pending[0]['SRA_run'], 'FirstPendingLibraryOnly')
    selected = pending[0]
    source.library(meta, selected, target, args, reference, SOURCE_SHA)
    validate_completed(source, selected, target)
    source.unchanged(SOURCE_SHA)
    return result(args.library, False, projection, resources)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('library', choices=REMAINING)
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--quiet-window-confirmed', action='store_true')
    parser.add_argument('--reserved-growth-bytes', type=int)
    parser.add_argument('--disk-budget-bytes', type=int)
    args = parser.parse_args(argv)
    if not args.execute:
        print(json.dumps(dict(status='PLAN_ONLY_NO_IO_OR_NATIVE_CALL', library=args.library,
            original_source_path=str(SOURCE), original_source_sha256=SOURCE_SHA,
            invocation='original library() only, not original remaining CLI',
            pending_order_not_read_or_admitted=True)))
        return 0
    try:
        print(json.dumps(execute(args), allow_nan=False))
        return 0
    except Exception as exc:
        print(json.dumps(dict(status='STOP_RETAINED_NO_AUTOMATIC_RETRY',
            exception_class_only=type(exc).__name__,
            category=exc.category if isinstance(exc, GuardError) else 'OriginalFunctionOrIOFailure',
            raw_error_text_persisted=False)))
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
