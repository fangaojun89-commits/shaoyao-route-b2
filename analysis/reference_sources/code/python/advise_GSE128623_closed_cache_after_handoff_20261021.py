"""Finite cache maintenance after the authenticated no-fit serial handoff.

Retains the existing finite utility and all consumer/path/cap checks. A real
negative native exit with accepted saved outputs replaces only serial liveness.
"""
import argparse
import ast
import datetime
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path('/root/autodl-tmp/gse128623_label_reference_20261015r01')
PREVIOUS = ROOT / 'cache_window_after_compact_source_20261021_r02/advise_GSE128623_closed_cache_after_compact_20261021.py'
PREVIOUS_SHA = '26d278011122663e75982e5700010c2827f39f2f2c02ce47d621cf61782d2f03'
PINS = {
    'serial_saved_arithmetic_handoff_20261021_r01/HANDOFF_INTENT.json':
        '107514e99a84851d2eb599abb5955c24cb073d234ad00fb93e855400c1759d68',
    'full_finalization_execution_20261020_r01/FINALIZER_PROCESS_EXIT.json':
        '061045f02458972570207060caff33e95412740b1887094010240e018ac0039a',
    'saved_sampled_output_QA_20261021_r01/raw_all30/SAVED_OUTPUT_QA.json':
        'd176cdee49b8690c0988f4d787468586fb8b49604298eb8c5d892cf631e12333',
    'saved_sampled_output_QA_20261021_r01/raw_predefined29/SAVED_OUTPUT_QA.json':
        '4b4f7e103d11b85aa122277414cf578f32383acfe0037ae4c788e4ad5cfb0667',
}


def accepted_handoff(q, root=ROOT):
    records = {}
    for relative, expected in PINS.items():
        p = root / relative
        q.need(hashlib.sha256(p.read_bytes()).hexdigest() == expected,
               'ExactSavedOutputHandoffPin')
        records[relative] = q.read(p)
    e = records['full_finalization_execution_20261020_r01/FINALIZER_PROCESS_EXIT.json']
    q.need(type(e['actual_R_exit_code']) is int and e['actual_R_exit_code'] == -15
           and e['pid'] == 896255 and e['start_ticks'] == 2602926060
           and e['aggregate'] is None, 'TruthfulSerialHandoffNotCLISuccess')
    intent = records['serial_saved_arithmetic_handoff_20261021_r01/HANDOFF_INTENT.json']
    q.need(intent['pid'] == 896255 and intent['birth'] == 2602926060
           and intent['original_CLI_success_asserted'] is False
           and type(intent['models_restarted']) is int and intent['models_restarted'] == 0,
           'NoModelRestartHandoff')
    for cell in ('raw_all30', 'raw_predefined29'):
        r = records[f'saved_sampled_output_QA_20261021_r01/{cell}/SAVED_OUTPUT_QA.json']
        q.need(r['schema'] == 'gse128623-saved-sampled-output-qa-v1'
               and r['status'] == 'PASS_SELECTED_SAVED_SAMPLE_JOB_OUTPUTS_AND_PINNED_METADATA_ONLY'
               and r['cell_id'] == cell
               and all(type(r[k]) is int and r[k] == value for k, value in
                       (('checked_allocation_n', 10001), ('checked_family_n', 67), ('model_calls', 0)))
               and r['original_finalizer_CLI_success_asserted'] is False
               and r['final_scientific_acceptance'] is False
               and r['producer_independent_saved_tails_agree'] is True,
               'BothSavedSampledJobsAccepted')
    return dict(pid=896255, birth=2602926060,
                state='REAL_NATIVE_MINUS15_AFTER_ACCEPTED_SAVED_OUTPUT_HANDOFF',
                original_CLI_success=False, receipt_pins=PINS)


def handoff_state(q, utility, previous):
    fields = Path('/proc/659622/stat').read_text().rsplit(') ', 1)[1].split()
    q.need(int(fields[19]) == 2576253226 and fields[0] not in ('Z', 'X'),
           'ExactOriginalJGWControllerContinuity')
    c = Path('/sys/fs/cgroup')
    counters = dict(line.split() for line in (c / 'memory.stat').read_text().splitlines())
    used, limit = (int((c / name).read_text()) for name in ('memory.current', 'memory.max'))
    status = q.read(q.JGW / 'provenance/full_reference_02/STATUS.json')
    q.need(status['owner_pid'] == 659622, 'OriginalControllerOwner')
    return dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                headroom=limit-used, current=used, maximum=limit,
                file=int(counters['file']), anon=int(counters['anon']),
                dirty=int(counters['file_dirty']), writeback=int(counters['file_writeback']),
                JGW_status=status['status'], next_round=status.get('next_round'),
                bound_processes={'659622': dict(birth=2576253226, state=fields[0]),
                                 '896255': accepted_handoff(q),
                                 '947489': previous.accepted_compact_exit(q)},
                completed_G128_terminal_sha256=utility.completed_g128(),
                continuity_adapter_source=dict(path=str(Path(__file__)),
                    sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-round', type=int, required=True)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    assert hashlib.sha256(PREVIOUS.read_bytes()).hexdigest() == PREVIOUS_SHA
    spec = importlib.util.spec_from_file_location('accepted_compact_continuity', PREVIOUS)
    previous = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(previous)
    assert hashlib.sha256(previous.UTILITY.read_bytes()).hexdigest() == previous.UTILITY_SHA
    spec = importlib.util.spec_from_file_location('unchanged_finite_cache_window', previous.UTILITY)
    utility = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(utility)
    original_dependencies = utility.dependencies

    def dependencies():
        q, module = original_dependencies()
        replace_legacy_argv_check(q, module)
        q.state = lambda m: handoff_state(q, m, previous)
        return q, module

    utility.dependencies = dependencies
    result, code = utility.run(args.execute, args.expected_round)
    print(json.dumps({k: result.get(k) for k in
                     ('status', 'advised_bytes', 'advised_files', 'stop_reason', 'exception_class_only')}))
    raise SystemExit(code)


def replace_legacy_argv_check(q, module):
    """Replace only the obsolete live argv check in the pinned legacy function."""
    source = Path(module.__file__)
    q.need(hashlib.sha256(source.read_bytes()).hexdigest() ==
           'df1f143ea8a73479fd12fd7e409d3149a46e807c3b056cd0f1ddc18640976c89',
           'PinnedLegacyCandidateImplementation')
    tree = ast.parse(source.read_text())
    functions = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'candidates']
    q.need(len(functions) == 1, 'ExactLegacyCandidateFunction')
    fn = functions[0]
    positions = [i for i, n in enumerate(fn.body) if isinstance(n, ast.Assign)
                 and len(n.targets) == 1 and isinstance(n.targets[0], ast.Name)
                 and n.targets[0].id == 'argv']
    q.need(len(positions) == 1, 'OneObsoleteArgvAssignment')
    i = positions[0]
    check = fn.body[i + 1]
    q.need(isinstance(check, ast.Expr) and isinstance(check.value, ast.Call)
           and isinstance(check.value.func, ast.Name) and check.value.func.id == 'need'
           and len(check.value.args) == 2 and isinstance(check.value.args[1], ast.Constant)
           and check.value.args[1].value == 'ActualFinalizerExactArguments',
           'OnlyObsoleteLiveArgvGuardReplaced')
    replacement = ast.parse('_accepted_serial_handoff()').body
    fn.body[i:i + 2] = replacement
    isolated = ast.fix_missing_locations(ast.Module(body=[fn], type_ignores=[]))
    namespace = dict(module.__dict__)
    namespace['_accepted_serial_handoff'] = lambda: accepted_handoff(q)
    exec(compile(isolated, str(source) + ':terminal_handoff_guard', 'exec'), namespace)
    module.candidates = namespace['candidates']


if __name__ == '__main__':
    main()
