"""Readback verification only: independent BH and score-contrast arithmetic."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

BASE = Path('R:/Shaoyao_storage/data/GSE298414_preparation/analysis_20261025_03')
OUT = BASE / 'program_and_stratified_extension_01'


def bh(p):
    x = np.asarray(p, dtype=float)
    assert np.isfinite(x).all() and ((x >= 0) & (x <= 1)).all()
    order = np.argsort(x, kind='stable')
    q = np.minimum.accumulate((x[order] * len(x) / np.arange(1, len(x) + 1))[::-1])[::-1]
    ans = np.empty(len(x)); ans[order] = np.minimum(q, 1)
    return ans


def check_bh(df, grouping, p, q):
    for _, g in df.groupby(grouping, sort=False):
        assert np.allclose(bh(g[p]), g[q], atol=1e-12, rtol=1e-10)


def main():
    receipt = json.loads((OUT / 'receipt.json').read_text())
    for x in receipt['source_manifest']:
        assert hashlib.sha256(Path(x['path']).read_bytes()).hexdigest() == x['sha256']
    inventory = pd.read_csv(OUT / 'contrast_inventory.tsv', sep='\t')
    samples = pd.read_csv(BASE / 'full_72/sample_qc.tsv', sep='\t')
    assert len(inventory) == 36 and not inventory.contrast_id.duplicated().any()
    assert inventory.groupby('contrast_class').size().tolist() == [12, 12, 12]
    gene_rows = {}
    for fn, p, q in [('new_24_edgeR_genes.tsv.gz', 'PValue', 'FDR'),
                     ('new_24_voom_genes.tsv.gz', 'P.Value', 'adj.P.Val')]:
        df = pd.read_csv(OUT / fn, sep='\t')
        assert len(df) == 13761 * 24 and not df.duplicated(['contrast_id', 'gene_id']).any()
        check_bh(df, ['contrast_id'], p, q)
        check_bh(df, ['contrast_class'], p, 'FDR_class')
        gene_rows[fn] = len(df)
    camera = pd.read_csv(OUT / 'camera_all_programs.tsv.gz', sep='\t')
    check_bh(camera, ['collection', 'contrast_id'], 'PValue', 'FDR')
    check_bh(camera, ['collection', 'contrast_class'], 'PValue', 'FDR_class')
    assert not camera.duplicated(['collection', 'contrast_id', 'term_id']).any()
    assert camera.NGenes.between(10, 500).all()
    full = pd.read_csv(OUT / 'signed_activity_all_programs.tsv', sep='\t')
    influence = pd.read_csv(OUT / 'signed_activity_seven_omissions.tsv.gz', sep='\t')
    all_signed = pd.concat([full, influence], ignore_index=True)
    check_bh(all_signed, ['analysis', 'collection', 'contrast_id'], 'PValue', 'FDR')
    check_bh(all_signed, ['analysis', 'collection', 'contrast_class'], 'PValue', 'FDR_class')
    assert not all_signed.duplicated(['analysis', 'collection', 'contrast_id', 'term_id']).any()
    assert (all_signed.CI_low <= all_signed.estimate).all() and (all_signed.estimate <= all_signed.CI_high).all()
    checked = 0
    for collection in ['PROGENy', 'DoRothEA']:
        scores = pd.read_csv(OUT / (collection + '_scores.tsv'), sep='\t').set_index('term_id')
        assert list(scores.columns) == list(samples.sample_id)
        for analysis, tab in all_signed[all_signed.collection == collection].groupby('analysis'):
            ss = samples if analysis == 'full_72' else samples[samples.sample_id != analysis.removeprefix('omit_')]
            assert len(ss) == (72 if analysis == 'full_72' else 71)
            for r in inventory.itertuples():
                pos = ss.loc[ss.group == r.positive_group, 'sample_id']
                neg = ss.loc[ss.group == r.negative_group, 'sample_id']
                expected = scores[pos].mean(axis=1) - scores[neg].mean(axis=1)
                t = tab[tab.contrast_id == r.contrast_id].set_index('term_id')
                assert set(t.index) == set(scores.index)
                assert np.allclose(t.estimate, expected.loc[t.index], atol=1e-10, rtol=1e-9)
                checked += len(t)
    families = pd.read_csv(OUT / 'program_family_summary.tsv', sep='\t')
    assert len(families) == 144
    for r in families.itertuples():
        data = camera if r.method == 'CAMERA_membership' else full
        t = data[(data.collection == r.collection) & (data.contrast_id == r.contrast_id)]
        assert len(t) == r.tested and (t.FDR <= .05).sum() == r.within_FDR05
        assert (t.FDR_class <= .05).sum() == r.class_FDR05
    cov = pd.read_csv(OUT / 'unsigned_program_coverage.tsv', sep='\t')
    assert not cov.duplicated(['collection', 'term_id']).any()
    for name, g in cov.groupby('collection'):
        assert set(g.loc[g.eligible, 'term_id']) == set(camera.loc[camera.collection == name, 'term_id'])
    result = dict(status='INDEPENDENT_BH_AND_SCORE_CONTRAST_READBACK_PASSED',
                  new_gene_rows=gene_rows, camera_rows=len(camera), signed_primary_rows=len(full),
                  signed_effects_checked_directly_from_library_means=checked,
                  program_families=144, omission_analyses=7, new_models_run=0,
                  previous_12_voom_results_bridged_by_extension=True)
    (OUT / 'review_receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
