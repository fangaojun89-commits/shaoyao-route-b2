"""Validate retained anchor outputs with original R checks and a separate gate readback."""
import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]
PROVENANCE = 'data/anchor_saved_outputs/PROVENANCE.json'


def identity(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            digest.update(block)
    return {'bytes': path.stat().st_size, 'sha256': digest.hexdigest()}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def rows(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def number(value):
    return float('nan') if value in ('NA', '') else float(value)


def same_number(left, right):
    return (math.isnan(left) and math.isnan(right)) or math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-12)


def boolean(value):
    if value == 'NA':
        return None
    require(value in ('TRUE', 'FALSE'), 'Unexpected logical value')
    return value == 'TRUE'


def sign(value):
    return None if not math.isfinite(value) else 1 if value > 1e-12 else -1 if value < -1e-12 else 0


def equal_sign(left, right):
    return None if left is None or right is None else left == right


def logical_and(values):
    return False if False in values else None if None in values else True


def indexed(records, fields):
    result = {}
    for record in records:
        key = tuple(record[field] for field in fields)
        require(key not in result, 'Duplicate saved-result key: ' + str(key))
        result[key] = record
    return result


def check_gate(project):
    """Re-evaluate recorded gate arithmetic; never fit or resample."""
    program = project / 'results/GSE16710/program_level'
    diagnostic = project / 'results/GSE16710/diagnostics'
    machine = indexed(rows(program / 'S1_GATE.tsv'), ['collection', 'feature_id'])
    permutations = indexed(rows(diagnostic / 'program_Freedman_Lane_10000.tsv'), ['family', 'feature_id_original'])
    expected = {'Reactome': 1656, 'PROGENy': 14, 'DoRothEA': 101}
    require(dict(collections.Counter(key[0] for key in machine)) == expected, 'Gate family membership differs')
    checked = set()
    numeric_fields = 0
    for family, count in expected.items():
        tables = rows(program / ('S1_gate_' + family + '.tsv'))
        features = indexed(tables, ['feature_id'])
        require(len(features) == count, 'Family feature count differs')
        omitted = collections.defaultdict(list)
        loo = rows(diagnostic / (family.lower() + '_leave_one_animal_out_long.tsv.gz'))
        require(len(loo) == count * 26, 'Incomplete saved omission table')
        indexed(loo, ['feature_id', 'omitted_sample_id'])
        for record in loo:
            omitted[record['feature_id']].append(number(record['loo_estimate']))
        for (feature,), record in features.items():
            key = (family, feature)
            require(key in machine and key in permutations and feature in omitted, 'Saved result join is incomplete')
            summary = machine[key]
            effect = number(record['activity_estimate'])
            low, high = number(record['activity_ci_low']), number(record['activity_ci_high'])
            direction = sign(effect)
            ci = None if not all(map(math.isfinite, (low, high))) else (low > 0 and high > 0) or (low < 0 and high < 0)
            original = equal_sign(sign(number(record['original_annotation_estimate'])), direction)
            all27 = equal_sign(sign(number(record['all27_estimate'])), direction)
            estimates = omitted[feature]
            require(len(estimates) == 26 and all(map(math.isfinite, estimates)), 'Omission estimates incomplete or nonfinite')
            full = number(record['full_estimate'])
            flips = sum(sign(value) != sign(full) for value in estimates)
            fraction = sum(sign(value) == sign(full) for value in estimates) / 26
            no_flip = flips == 0
            calculated = {'loo_min': min(estimates), 'loo_max': max(estimates),
                          'loo_direction_fraction': fraction, 'loo_sign_flips': flips,
                          'max_relative_change': max(abs(value-full)/(abs(full)+1e-12) for value in estimates)}
            for field, value in calculated.items():
                require(same_number(number(record[field]), value), 'Omission arithmetic differs: ' + field)
                numeric_fields += 1
            require(boolean(record['no_single_animal_sign_flip']) == no_flip, 'Omission sign gate differs')
            require(boolean(record['passes_direction_fraction_0_80']) == (fraction >= .8), 'Omission fraction gate differs')
            require(boolean(record['activity_ci_excludes_zero']) == ci, 'CI gate differs')
            require(boolean(record['original_annotation_direction_stable']) == original, 'Original-annotation sign differs')
            require(boolean(record['all27_direction_stable']) == all27, 'All-27 sign differs')
            require(same_number(number(record['activity_sign']), float('nan') if direction is None else direction), 'Activity sign differs')
            q_name = 'camera_q_value_BH' if family == 'Reactome' else 'activity_q_value_BH'
            q = number(record[q_name])
            primary = None if not math.isfinite(q) else q <= .05
            if family == 'Reactome':
                camera = {'Up': 1, 'Down': -1}.get(record['camera_direction'])
                agreement = equal_sign(camera, direction)
                require(boolean(record['primary_camera_pass']) == primary, 'CAMERA threshold differs')
                require(boolean(record['camera_activity_direction_agreement']) == agreement, 'CAMERA/activity sign differs')
                passing = logical_and([primary, ci, agreement, original, all27, no_flip])
                feature_name = record['term_name']
            else:
                require(boolean(record['primary_activity_pass']) == primary, 'Activity threshold differs')
                passing = logical_and([primary, ci, original, all27, no_flip])
                feature_name = feature
            require(boolean(record['S1_pass']) == passing == boolean(summary['S1_pass']), 'Final admission differs')
            require(summary['feature_name'] == feature_name, 'Feature label differs')
            require(boolean(summary['original_annotation_direction_stable']) == original, 'Machine annotation flag differs')
            require(boolean(summary['all27_direction_stable']) == all27, 'Machine all-27 flag differs')
            comparisons = {'primary_effect': effect, 'primary_ci_low': low, 'primary_ci_high': high,
                           'primary_q_value': q, 'loo_direction_fraction': fraction,
                           'direction': float('nan') if direction is None else direction,
                           'permutation_q_value_BH': number(permutations[key]['empirical_q_value_BH'])}
            for field, value in comparisons.items():
                require(same_number(number(summary[field]), value), 'Machine value differs: ' + field)
                numeric_fields += 1
            require(same_number(number(record['permutation_p_value']), number(permutations[key]['empirical_p_value'])), 'Saved permutation P join differs')
            require(same_number(number(record['permutation_q_value_BH']), comparisons['permutation_q_value_BH']), 'Saved permutation q join differs')
            checked.add(key)
    require(checked == set(machine) == set(permutations), 'Saved gate/permutation key spaces differ')
    return dict(status='PASS_SAVED_GATE_AND_OMISSION_ARITHMETIC', family_features=expected,
                programs_checked=len(checked), admitted=sum(boolean(r['S1_pass']) is True for r in machine.values()),
                omission_estimates_checked=1771*26, numeric_field_comparisons=numeric_fields,
                model_fits=0, activity_scores_recomputed=0, new_permutations=0,
                scope='Consistency of retained gate, omission estimates and stored permutation P/q joins; no model or phenotype validation')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rscript', required=True)
    parser.add_argument('--r-library', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output must be a new directory')
    library = args.r_library.resolve(strict=True)
    require(library.is_dir(), 'R library must be a directory')
    manifest = json.loads((REPO / 'MANIFEST.json').read_text(encoding='utf-8'))
    contract = json.loads((REPO / PROVENANCE).read_text(encoding='utf-8'))
    dependencies = ['pipeline/anchor_saved.py', 'analysis/anchor_saved/run.R', PROVENANCE,
                    contract['validator']] + [record['public_path'] for record in contract['files']]
    pins = {}
    for name in dict.fromkeys(dependencies):
        actual = identity(REPO / name)
        require(actual == manifest['files'][name], 'Release source/input mismatch: ' + name)
        pins[name] = actual
    for record in contract['files']:
        require(pins[record['public_path']] == record['public_identity'], 'Saved input provenance differs')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    receipts = output / 'receipts'
    receipts.mkdir()
    project = output / 'project'
    for record in contract['files']:
        target = project / record['runtime_path']
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / record['public_path'], target)
    code = project / 'code/R'
    code.mkdir(parents=True)
    shutil.copyfile(REPO / contract['validator'], code / 'validate_GSE16710_outputs.R')
    shutil.copyfile(REPO / 'analysis/anchor_saved/run.R', code / 'run.R')
    write_json(receipts / 'INPUT_PINS.json', pins)
    command = [args.rscript, '--vanilla', (code / 'run.R').as_posix()]
    env = os.environ.copy()
    locale = 'English_United States.utf8' if os.name == 'nt' else 'C.UTF-8'
    env.update(R_LIBS_USER=library.as_posix(), RENV_CONFIG_AUTOLOADER_ENABLED='FALSE',
               GSE16710_ADJUDICATION_SOURCE='CANONICAL_POST_QC_RECORD',
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               LANG=locale, LC_ALL=locale, LC_CTYPE=locale)
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    with (receipts / 'stdout.txt').open('wb') as stdout, (receipts / 'stderr.txt').open('wb') as stderr:
        result = subprocess.run(command, cwd=output, env=env, stdout=stdout, stderr=stderr, check=False)
    independent = check_gate(project) if result.returncode == 0 else {'status': 'NOT_RUN_AFTER_R_FAILURE'}
    write_json(receipts / 'INDEPENDENT_GATE_QA.json', independent)
    same = all(identity(REPO / name) == pin for name, pin in pins.items())
    staged_same = all(identity(project / record['runtime_path']) == record['public_identity'] for record in contract['files'])
    summary_path = project / 'results/GSE16710/diagnostics/GSE16710_output_validation_summary.json'
    summary = json.loads(summary_path.read_text(encoding='utf-8')) if summary_path.exists() else {}
    success = (result.returncode == 0 and same and staged_same and summary.get('overall_status') == 'PASS'
               and summary.get('checks_total') == 51 and independent.get('status') == 'PASS_SAVED_GATE_AND_OMISSION_ARITHMETIC')
    products = {path.relative_to(output).as_posix(): identity(path) for path in output.rglob('*') if path.is_file()}
    write_json(receipts / 'PARENT.json', dict(status='PASS' if success else 'FAIL', command=command,
        invocation=sys.argv, python=sys.version, platform=platform.platform(),
        release_version=manifest['version'], release_manifest_identity=identity(REPO / 'MANIFEST.json'),
        started_utc=started, finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        returncode=result.returncode, source_and_public_inputs_unchanged=same, staged_inputs_unchanged=staged_same,
        original_checks=summary, independent_gate=independent, outputs=products,
        model_fits=0, scores_recomputed=0, new_allocations=0,
        scope='Retained-output validation; not raw-to-model reconstruction or independent biological validation'))
    print(json.dumps({'status': 'PASS' if success else 'FAIL', 'original_checks': summary.get('checks_pass'),
                      'gate_programs_checked': independent.get('programs_checked'), 'output': str(output)}))
    if not success:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
