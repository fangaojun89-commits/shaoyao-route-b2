"""Descriptive paired comparison of retained PC12 fits; no model estimation."""
from pathlib import Path
import argparse, hashlib, json, math
import numpy as np
import pandas as pd

def rank_correlations(a, b):
    """Descriptive Spearman and tie-adjusted Kendall, without significance tests."""
    rho = float(a.rank(method='average').corr(b.rank(method='average')))
    av, bv = a.to_numpy(), b.to_numpy()
    numerator = 0
    for i in range(len(av)-1):
        numerator += int(np.sum(np.sign(av[i+1:]-av[i]) * np.sign(bv[i+1:]-bv[i])))
    n0 = len(av)*(len(av)-1)//2
    ties_a = sum(int(n)*(int(n)-1)//2 for n in a.value_counts())
    ties_b = sum(int(n)*(int(n)-1)//2 for n in b.value_counts())
    tau = numerator / math.sqrt((n0-ties_a)*(n0-ties_b))
    return rho, tau

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data/pc12_dose_saved'
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir', required=True, type=Path)
args = parser.parse_args()
OUT = args.output_dir.resolve()
if OUT.exists():
    raise SystemExit('Output directory already exists; select a new directory.')
manifest = json.loads((BASE / 'PROVENANCE.json').read_text(encoding='utf-8'))
for item in manifest['tables']:
    p = BASE / item['path']
    if hashlib.sha256(p.read_bytes()).hexdigest() != item['sha256']:
        raise SystemExit('Input hash mismatch: ' + item['path'])
OUT.mkdir(parents=True)
BULK = OUT / 'paired'
BULK.mkdir()
STEM = 'Route_B2_PC12_dose_encoding_comparison_20261025'
specs = [
 ('edgeR_genes', 'genes/edgeR_QL_ordered_trend_all_genes.tsv.gz', 'genes/edgeR_QL_logdose_trend_all_genes.tsv.gz', 'raw_gene_id', 'logFC', None, None),
 ('Reactome_score', 'programs/Reactome97_score_trend.tsv', 'programs/Reactome97_score_logtrend.tsv', 'feature', 'effect', 'CI_low', 'CI_high'),
 ('QuickGO_SynGO_score', 'programs/QuickGO_SynGO_sensitivity_score_trend.tsv', 'programs/QuickGO_SynGO_sensitivity_score_logtrend.tsv', 'feature', 'effect', 'CI_low', 'CI_high'),
 ('Reactome_CAMERA', 'programs/Reactome97_camera_ordered_trend.tsv', 'programs/Reactome97_camera_logdose_trend.tsv', 'term_id', None, None, None),
 ('QuickGO_SynGO_CAMERA', 'programs/QuickGO_SynGO_sensitivity_camera_ordered_trend.tsv', 'programs/QuickGO_SynGO_sensitivity_camera_logdose_trend.tsv', 'term_id', None, None, None),
 ('PROGENy', 'activities/PROGENy_trend.tsv', 'activities/PROGENy_logtrend.tsv', 'feature', 'effect', 'CI_low', 'CI_high'),
 ('DoRothEA', 'activities/DoRothEA_trend.tsv', 'activities/DoRothEA_logtrend.tsv', 'feature', 'effect', 'CI_low', 'CI_high'),
]
inputs = {}
def read(rel):
    p = BASE / rel
    inputs[str(p.relative_to(ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
    return pd.read_csv(p, sep='\t')
key = read('dose_response/dose_key.tsv')
assert key.ordered_dose.tolist() == [0, 1, 2, 3]
assert key.dose_ug_ml.tolist() == [0, 4, 20, 100]
scale_a = float(key.ordered_dose.max() - key.ordered_dose.min())
scale_b = math.log10(101)
rows, pairs = [], {}
for family, afile, bfile, idcol, eff, lo, hi in specs:
    a, b = read(afile), read(bfile)
    assert a[idcol].notna().all() and b[idcol].notna().all()
    assert a[idcol].is_unique and b[idcol].is_unique
    assert set(a[idcol]) == set(b[idcol]), family
    a = a.set_index(idcol).sort_index(); b = b.set_index(idcol).reindex(a.index)
    x = pd.DataFrame(index=a.index)
    x['FDR_ordinal'] = a.FDR; x['FDR_logdose'] = b.FDR
    x['P_ordinal'] = a.PValue; x['P_logdose'] = b.PValue
    assert np.isfinite(x.to_numpy()).all(), family
    assert ((x >= 0) & (x <= 1)).all().all()
    sa, sb = x.FDR_ordinal.lt(.05), x.FDR_logdose.lt(.05)
    union, both = sa | sb, sa & sb
    if eff:
        for tag, frame, scale in [('ordinal', a, scale_a), ('logdose', b, scale_b)]:
            x['coefficient_'+tag] = frame[eff]
            x['endpoint_effect_'+tag] = frame[eff] * scale
            if lo is not None:
                x['endpoint_CI_low_'+tag] = frame[lo] * scale
                x['endpoint_CI_high_'+tag] = frame[hi] * scale
        assert np.isfinite(x.to_numpy()).all(), family
        da, db = np.sign(a[eff]), np.sign(b[eff])
        rho, tau = rank_correlations(a[eff], b[eff])
        k = min(20, len(a))
        # IDs are sorted before stable sorting to make boundary ties deterministic.
        topa = set(a[eff].abs().sort_values(ascending=False, kind='stable').index[:k])
        topb = set(b[eff].abs().sort_values(ascending=False, kind='stable').index[:k])
        overlap = len(topa & topb)
        absdiff = (x.endpoint_effect_logdose-x.endpoint_effect_ordinal).abs()
        median_diff, p95_diff = float(absdiff.median()), float(absdiff.quantile(.95))
        ci_switch = None
        if lo is not None:
            zero_a = (x.endpoint_CI_low_ordinal <= 0) & (x.endpoint_CI_high_ordinal >= 0)
            zero_b = (x.endpoint_CI_low_logdose <= 0) & (x.endpoint_CI_high_logdose >= 0)
            ci_switch = int((zero_a != zero_b).sum())
    else:
        assert set(a.Direction) <= {'Up','Down'} and set(b.Direction) <= {'Up','Down'}
        da, db = a.Direction.map({'Up':1,'Down':-1}), b.Direction.map({'Up':1,'Down':-1})
        rho = tau = k = overlap = median_diff = p95_diff = ci_switch = None
    x['direction_ordinal'] = da; x['direction_logdose'] = db
    x['significant_ordinal'] = sa; x['significant_logdose'] = sb
    same = da.eq(db); opposite = (da * db).lt(0)
    x.to_csv(BULK / (family+'_paired.tsv.gz'), sep='\t', compression='gzip')
    pairs[family] = x
    rows.append(dict(family=family, n_pairs=len(x), missing_numeric_pairs=0,
        direction_same=int(same.sum()), direction_opposite=int(opposite.sum()),
        direction_same_fraction=float(same.mean()), effect_spearman=rho, effect_kendall_tau_b=tau,
        ordinal_FDR05=int(sa.sum()), logdose_FDR05=int(sb.sum()), both_FDR05=int(both.sum()),
        union_FDR05=int(union.sum()), discovery_jaccard=float(both.sum()/union.sum()) if union.any() else None,
        union_direction_opposite=int((union & opposite).sum()),
        top_absolute_effect_k=k, top_absolute_effect_overlap=overlap,
        endpoint_abs_difference_median=median_diff, endpoint_abs_difference_p95=p95_diff,
        CI_zero_inclusion_changes=ci_switch))
summary = pd.DataFrame(rows)
summary.to_csv(OUT/(STEM+'.tsv'), sep='\t', index=False, na_rep='NA')
high = read('genes/edgeR_QL_high_vs_vehicle_all_genes.tsv.gz').set_index('raw_gene_id')
g = pairs['edgeR_genes']
assert high.index.is_unique and set(high.index)==set(g.index)
high = high.reindex(g.index)
coherent = {}
for tag in ['ordinal', 'logdose']:
    yes = high.FDR.lt(.05) & g['FDR_'+tag].lt(.05) & np.sign(high.logFC).eq(g['direction_'+tag])
    coherent[tag] = set(g.index[yes])
coh = dict(ordinal=len(coherent['ordinal']), logdose=len(coherent['logdose']),
           shared=len(coherent['ordinal'] & coherent['logdose']),
           ordinal_only=len(coherent['ordinal']-coherent['logdose']),
           logdose_only=len(coherent['logdose']-coherent['ordinal']))
assert coh['ordinal'] == 40
receipt = dict(scientific_fits_run=0, scope='canonical PC12 saved-fit descriptive comparison',
    primary='ordinal 0,1,2,3', alternative='log10(dose_ug_ml+1)',
    endpoint_scales=dict(ordinal=scale_a, logdose=scale_b),
    high_and_trend_coherent=coh, inputs_sha256=inputs,
    full_pairs_directory=str(BULK), original_candidate_gates_changed=False,
    unsupported_edgeR_approximate_intervals_excluded=True)
for rel, digest in inputs.items():
    assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest() == digest
receipt['source_files_unchanged'] = True
(OUT/(STEM+'_sources.json')).write_text(json.dumps(receipt, indent=2), encoding='utf-8')

print(summary.to_string(index=False))
print(json.dumps(coh))
