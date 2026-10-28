"""Descriptive all-program comparison; never refits or pools study-level effects."""
from pathlib import Path
import argparse, datetime, hashlib, json
import numpy as np
import pandas as pd

PROGRAMS = ['Androgen', 'EGFR', 'Estrogen', 'Hypoxia', 'JAK-STAT', 'MAPK',
            'NFkB', 'p53', 'PI3K', 'TGFb', 'TNFa', 'Trail', 'VEGF', 'WNT']
NUMERIC_SOURCE = ['estimate', 'ci_low', 'ci_high', 'PValue', 'FDR']


def numeric(series):
    result = pd.to_numeric(series.replace({'NA': None, '': None}), errors='raise')
    if np.isinf(result.to_numpy(dtype=float)).any():
        raise ValueError('Infinite input is not an interpretable program effect')
    return result


def compare(basis, contrasts, coverage):
    """Keep every frozen cell; summarize only finite matched effects, without P values."""
    if basis[['comparator_id', 'program_id']].duplicated().any():
        raise ValueError('Duplicate retained comparator/program')
    ids = list(dict.fromkeys(basis.comparator_id))
    if len(ids) != 19 or len(basis) != 266:
        raise ValueError('Expected the complete frozen 19 by 14 comparator panel')
    for cid in ids:
        if set(basis.loc[basis.comparator_id == cid, 'program_id']) != set(PROGRAMS):
            raise ValueError('Incomplete frozen program panel')
    cov = coverage.loc[coverage.collection == 'PROGENy'].copy()
    new = contrasts.loc[contrasts.collection == 'PROGENy'].copy()
    if cov.term_id.duplicated().any() or new.term_id.duplicated().any():
        raise ValueError('Duplicate new-source program')
    if set(cov.term_id) != set(PROGRAMS) or not set(new.term_id).issubset(PROGRAMS):
        raise ValueError('New-source coverage must enumerate all 14 programs')
    cov = cov.set_index('term_id').loc[PROGRAMS].reset_index()
    cov['eligible'] = cov.eligible.map(str).str.lower().map({'true': True, 'false': False})
    if cov.eligible.isna().any():
        raise ValueError('Unrecognized eligibility value')
    cov['nonzero_expressed_targets'] = numeric(cov.nonzero_expressed_targets)
    if not ((cov.nonzero_expressed_targets >= 0) &
            (cov.nonzero_expressed_targets % 1 == 0)).all():
        raise ValueError('Invalid expressed-target count')
    if not (cov.eligible == cov.nonzero_expressed_targets.between(10, 500)).all():
        raise ValueError('Eligibility differs from the frozen 10–500 expressed-target rule')
    if set(new.term_id) != set(cov.loc[cov.eligible, 'term_id']):
        raise ValueError('Every eligible program requires a result row; ineligible programs remain absent')
    for col in NUMERIC_SOURCE:
        new[col] = numeric(new[col])
    for col in ['PValue', 'FDR']:
        if not new[col].dropna().between(0, 1).all():
            raise ValueError('Invalid source probability')
    if not new.estimate.notna().all():
        raise ValueError('Eligible program has no finite effect')
    if not new.status.isin(['TESTED', 'NUMERICALLY_DEGENERATE_NO_INFERENCE']).all():
        raise ValueError('Unexpected signed-activity status')
    tested = new.status == 'TESTED'
    if new.loc[tested, NUMERIC_SOURCE].isna().any().any():
        raise ValueError('Tested program lacks its source uncertainty or P value')
    if not ((new.loc[tested, 'ci_low'] <= new.loc[tested, 'estimate']) &
            (new.loc[tested, 'estimate'] <= new.loc[tested, 'ci_high'])).all():
        raise ValueError('Invalid source confidence interval')
    sci = cov[['term_id', 'eligible', 'nonzero_expressed_targets']].merge(
        new[['term_id', 'status'] + NUMERIC_SOURCE], on='term_id', how='left', validate='one_to_one')
    sci.loc[~sci.eligible, 'status'] = 'INELIGIBLE_EXPRESSED_TARGET_COVERAGE'
    sci = sci.rename(columns={k: 'SCI_' + k for k in sci if k != 'term_id'})
    sci = sci.rename(columns={'term_id': 'program_id', 'SCI_estimate': 'SCI_effect'})
    long = basis.drop(columns=['SCI_effect']).merge(sci, on='program_id', how='left', validate='many_to_one', sort=False)
    long['direction_relation'] = 'MISSING_EFFECT'
    summaries = []
    for cid in ids:
        take = long.comparator_id == cid
        old = numeric(long.loc[take, 'effect'])
        drug = long.loc[take, 'SCI_effect']
        paired = old.notna() & drug.notna()
        zero = paired & ((old == 0) | (drug == 0))
        same = paired & ~zero & (np.sign(old) == np.sign(drug))
        opposite = paired & ~zero & ~same
        for mask, label in [(zero, 'AT_LEAST_ONE_EXACT_ZERO'), (same, 'SAME_SIGN'), (opposite, 'OPPOSITE_SIGN')]:
            long.loc[mask.index[mask], 'direction_relation'] = label
        rho = None
        reason = 'INSUFFICIENT_FINITE_PAIRS'
        if paired.sum() >= 2:
            x, y = old[paired].rank(method='average'), drug[paired].rank(method='average')
            if x.nunique() > 1 and y.nunique() > 1:
                rho = float(x.corr(y)); reason = 'DESCRIPTIVE_ONLY_NO_CORRELATION_TEST'
            else:
                reason = 'CONSTANT_EFFECT_VECTOR'
        assert same.sum()+opposite.sum()+zero.sum() == paired.sum()
        first = long.loc[take].iloc[0]
        summaries.append(dict(comparator_id=cid, panel_role=first.panel_role,
            estimand=first.estimand, score_definition=first.score_definition,
            direction_interpretation=first.direction_interpretation,
            programs_total=14, finite_pairs=int(paired.sum()), same_sign=int(same.sum()),
            opposite_sign=int(opposite.sum()), exact_zero_pairs=int(zero.sum()),
            missing_pairs=int((~paired).sum()), descriptive_spearman=rho,
            spearman_status=reason, inference='No pooled effect, agreement P value or efficacy inference'))
    return long, pd.DataFrame(summaries)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--basis-dir', type=Path, required=True)
    parser.add_argument('--contrasts', type=Path, required=True)
    parser.add_argument('--coverage', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args=parser.parse_args()
    basis_path=args.basis_dir/'retained_PROGENy_comparator_basis.tsv'
    frozen=json.loads((args.basis_dir/'receipt.json').read_text())
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    if digest(basis_path) != frozen['output_sha256']:
        raise ValueError('Retained comparator basis differs from frozen staging record')
    read=lambda p:pd.read_csv(p,sep='\t',dtype=str,keep_default_na=False)
    long,summary=compare(read(basis_path),read(args.contrasts),read(args.coverage))
    args.out.mkdir(parents=True,exist_ok=False)
    long.to_csv(args.out/'complete_program_comparisons.tsv',sep='\t',index=False,na_rep='NA')
    summary.to_csv(args.out/'descriptive_comparison_summary.tsv',sep='\t',index=False,na_rep='NA')
    receipt=dict(status='COMPLETE_DESCRIPTIVE_COMPARISON',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        inputs={str(p):digest(p) for p in [basis_path,args.contrasts,args.coverage]},
        program_rows=len(long),comparators=len(summary),biological_models_run=0,
        source_q_values_recomputed=False,significance_filtering=False,
        source_wide_outputs_replaced=False,protected_holdout=False)
    (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':
    main()
