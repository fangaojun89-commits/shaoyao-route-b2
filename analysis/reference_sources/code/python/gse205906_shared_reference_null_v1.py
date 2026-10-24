"""Analytical and synthetic shared-reference diagnostic; no study expression inputs."""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import random
import sys

ROOT = Path(__file__).resolve().parents[2]
REGISTRY = 'config/GSE205906_shared_reference_null_design_v1_20260909.json'
CONFIG = 'config/GSE205906_bulk_model_v2_r4.json'
CONFIG_SHA = 'd8304d8713593f98d80119dc2f0d5d3d39f6b7eb096a3ccee0bc7b13cab9f16c'
OUTPUT = 'results/diagnostics/GSE205906_shared_reference_null_v1_20260909.json'

def need(ok, label):
    if not ok: raise ValueError(label)

def contrast_moments(group_sizes, first, second):
    need(len(group_sizes) == len(first) == len(second) and len(group_sizes) > 0, 'DIMENSION')
    need(all(type(n) is int and n > 0 for n in group_sizes), 'GROUP_SIZE')
    need(all(type(c) is int for c in (*first, *second)), 'COEFFICIENT')
    covariance = sum((Fraction(a*b, n) for a,b,n in zip(first,second,group_sizes)), Fraction(0))
    var1 = sum((Fraction(a*a,n) for a,n in zip(first,group_sizes)), Fraction(0))
    var2 = sum((Fraction(b*b,n) for b,n in zip(second,group_sizes)), Fraction(0))
    need(var1 > 0 and var2 > 0, 'ZERO_VARIANCE')
    return {'covariance': str(covariance), 'variance_first': str(var1), 'variance_second': str(var2),
            'correlation': float(covariance)/math.sqrt(float(var1*var2))}

def simulate(sizes, seed, draws):
    need(type(draws) is int and draws > 1, 'DRAWS')
    rng = random.Random(seed)
    scales = [1/math.sqrt(n) for n in sizes]
    records = {'shared': [0.0]*5+[0], 'independent_reference': [0.0]*5+[0]}
    for _ in range(draws):
        a,b,c = [rng.gauss(0, scale) for scale in scales]
        independent_b = rng.gauss(0, scales[1])
        x = b-a
        for label,y in (('shared',c-b), ('independent_reference',c-independent_b)):
            sums = records[label]
            for i,value in enumerate((x,y,x*x,y*y,x*y)): sums[i] += value
            sums[5] += int(x*y < 0)
    output = {}
    for label,(sx,sy,sxx,syy,sxy,opposite) in records.items():
        correlation = (sxy-sx*sy/draws)/math.sqrt((sxx-sx*sx/draws)*(syy-sy*sy/draws))
        output[label] = {'pearson_correlation': correlation, 'opposite_sign_n': opposite,
                         'opposite_sign_fraction': opposite/draws}
    return output

def audit():
    raw = (ROOT/CONFIG).read_bytes()
    need(hashlib.sha256(raw).hexdigest() == CONFIG_SHA, 'FROZEN_CONFIG_CHANGED')
    cfg = json.loads(raw)
    registry_raw = (ROOT/REGISTRY).read_bytes(); registry = json.loads(registry_raw)
    need(registry['source_config_sha256'] == CONFIG_SHA and registry['seed'] == 20260909 and
         registry['draws'] == 50000 and registry['production_go'] is False, 'REGISTRY_CHANGED')
    design = cfg['design']; groups = design['group_levels']; sizes = [design['expected_group_n'][g] for g in groups]
    need(groups == ['sham','SNI','PLM'] and sizes == [3,3,3], 'DESIGN_CHANGED')
    need(len(design['sample_ids_in_order']) == len(set(design['sample_ids_in_order'])) == 9 and
         design['groups_in_order'] == ['sham']*3+['SNI']*3+['PLM']*3, 'UNIT_MAPPING')
    first = design['contrasts']['local_disease']['vector']; second = design['contrasts']['intervention']['vector']
    need(first == [-1,1,0] and second == [0,-1,1], 'CONTRAST_CHANGED')
    analytical = contrast_moments(sizes,first,second)
    loao = {group: contrast_moments([n-int(j==i) for j,n in enumerate(sizes)],first,second)
            for i,group in enumerate(groups)}
    return {'schema': 'gse205906-shared-reference-null-result-v1', 'source_config_sha256': CONFIG_SHA,
            'registry_sha256': hashlib.sha256(registry_raw).hexdigest(), 'registry': registry,
            'analytical_full': analytical, 'analytical_one_animal_omitted_by_group': loao,
            'independent_reference_expected_correlation': 0.0,
            'gaussian_expected_opposite_sign_fraction_shared': 0.5-math.asin(analytical['correlation'])/math.pi,
            'simulation': simulate(sizes, registry['seed'], registry['draws']),
            'python': sys.version, 'argv': sys.argv, 'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'real_expression_data_read': False, 'real_count_model_fitted': False,
            'biological_negative_control_completed': False, 'production_go': False,
            'overall_objective': 'FULL_OBJECTIVE_NOT_COMPLETE'}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--write-report', action='store_true')
    args = parser.parse_args(); result = audit()
    if args.write_report:
        with (ROOT/OUTPUT).open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(result,stream,indent=2); stream.write('\n')
    print(json.dumps({k: result[k] for k in ('analytical_full','analytical_one_animal_omitted_by_group','simulation','production_go')},indent=2))
