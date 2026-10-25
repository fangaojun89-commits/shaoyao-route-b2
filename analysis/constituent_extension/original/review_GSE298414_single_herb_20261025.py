"""Read retained outputs only; independent BH, coverage, summaries and influence."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd

BASE = Path('R:/Shaoyao_storage/data/GSE298414_preparation/single_herb_20261025_01')
OUT = BASE / 'programs_01'

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

def rho(a, b):
    return np.corrcoef(pd.Series(a).rank(), pd.Series(b).rank())[0, 1]

def genes():
    receipt = json.loads((BASE / 'gene_receipt.json').read_text())
    inv = pd.read_csv(BASE / 'contrast_inventory.tsv', sep='\t')
    samples = pd.read_csv(BASE / 'samples.tsv', sep='\t')
    assert len(samples) == 126 and samples.groupby('group').size().eq(3).all()
    assert samples.group.nunique() == 42 and len(inv) == 60
    assert inv.groupby('contrast_class').size().eq(30).all()
    old = pd.read_csv(BASE.parent / 'analysis_20261025_03/program_and_stratified_extension_01/contrast_inventory.tsv', sep='\t')
    assert not set(inv.contrast_id) & set(old.contrast_id)
    totals = pd.read_csv(BASE / 'import_totals.tsv', sep='\t')
    assert np.allclose(totals.raw_expected_total, totals.adjusted_total, rtol=1e-12, atol=1e-7)
    stored = {}
    manifests = []
    for phase, n in [('full_126', 126), ('omit_GSM9013891_125', 125)]:
        folder = BASE / phase
        complete = json.loads((folder / 'COMPLETE.json').read_text())
        assert complete['libraries'] == n and complete['groups'] == 42
        assert complete['residual_df'] == n - 42 and complete['contrasts'] == 60
        filt = pd.read_csv(folder / 'gene_filter.tsv', sep='\t')
        expected = set(filt.loc[filt.retained, 'gene_id'])
        assert len(expected) == complete['genes'] == receipt['primary_genes']
        assert filt.loc[filt.retained, 'length_eligible'].all()
        tabs = []
        for fn, p, q in [('edgeR_all_genes.tsv.gz', 'PValue', 'FDR'), ('voom_all_genes.tsv.gz', 'P.Value', 'adj.P.Val')]:
            path = folder / fn
            df = pd.read_csv(path, sep='\t')
            assert len(df) == len(expected) * 60 and not df.duplicated(['contrast_id', 'gene_id']).any()
            assert set(df.contrast_id) == set(inv.contrast_id)
            for id, t in df.groupby('contrast_id', sort=False):
                assert set(t.gene_id) == expected
                assert (t.contrast_class == inv.set_index('contrast_id').loc[id, 'contrast_class']).all()
            assert np.isfinite(df.logFC).all()
            check_bh(df, ['contrast_id'], p, q)
            check_bh(df, ['contrast_class'], p, 'FDR_class')
            if 'CI.L' in df:
                assert (df['CI.L'] <= df.logFC).all() and (df.logFC <= df['CI.R']).all()
            manifests.append(dict(path=str(path), bytes=path.stat().st_size, sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
            tabs.append(df[['contrast_id', 'contrast_class', 'gene_id', 'logFC', q, 'FDR_class']].copy())
            del df
        a, b = tabs
        assert a[['contrast_id','gene_id']].equals(b[['contrast_id','gene_id']])
        summary = pd.read_csv(folder / 'contrast_summary.tsv', sep='\t').set_index('contrast_id')
        assert len(summary) == 60
        for id, x in a.groupby('contrast_id', sort=False):
            y = b.loc[x.index]; r = summary.loc[id]
            assert len(x) == r.tested_genes
            for field, count in dict(edgeR_within05=(x.FDR <= .05).sum(), edgeR_class05=(x.FDR_class <= .05).sum(),
                                     voom_within05=(y['adj.P.Val'] <= .05).sum(), voom_class05=(y.FDR_class <= .05).sum(),
                                     shared_class05=((x.FDR_class <= .05) & (y.FDR_class <= .05)).sum()).items():
                assert count == r[field], (phase, id, field)
            assert np.isclose(rho(x.logFC, y.logFC), r.effect_spearman, atol=1e-12)
            assert np.isclose(np.median(abs(x.logFC)), r.median_abs_edgeR_logFC, atol=1e-12)
        stored[phase] = a[['contrast_id','gene_id','logFC','FDR_class']].copy()
        del tabs, a, b
    a, b = stored.values()
    assert a[['contrast_id','gene_id']].equals(b[['contrast_id','gene_id']])
    impact = []
    for id, x in a.groupby('contrast_id', sort=False):
        y = b.loc[x.index]; hit = x.FDR_class <= .05; retained = hit & (y.FDR_class <= .05)
        impact.append(dict(contrast_id=id, tested_genes=len(x), full_edgeR_class05=int(hit.sum()),
                           omission_edgeR_class05=int((y.FDR_class <= .05).sum()),
                           full_discoveries_retained=int(retained.sum()),
                           retained_same_direction=int((retained & (np.sign(x.logFC)==np.sign(y.logFC))).sum()),
                           effect_spearman=rho(x.logFC,y.logFC),
                           median_abs_full_logFC=float(np.median(abs(x.logFC))),
                           median_abs_omission_logFC=float(np.median(abs(y.logFC)))))
    pd.DataFrame(impact).to_csv(BASE/'gene_omission_comparison.tsv',sep='\t',index=False)
    result=dict(status='GENE_READBACK_PASSED',libraries=126,contrasts=60,genes=receipt['primary_genes'],
                phases=2,within_and_class_BH_checked=True,summary_arithmetic_checked=True,
                prior_36_contrast_overlap=0,all_gene_effect_omission_comparisons=60,models_run=0,files=manifests)
    (BASE/'gene_review_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='files'},indent=2))

def programs():
    receipt=json.loads((OUT/'receipt.json').read_text())
    for x in receipt['source_manifest']:
        assert hashlib.sha256(Path(x['path']).read_bytes()).hexdigest()==x['sha256']
    inv=pd.read_csv(BASE/'contrast_inventory.tsv',sep='\t')
    samples=pd.read_csv(BASE/'samples.tsv',sep='\t')
    camera=pd.read_csv(OUT/'camera_all_programs.tsv.gz',sep='\t')
    check_bh(camera,['collection','contrast_id'],'PValue','FDR')
    check_bh(camera,['collection','contrast_class'],'PValue','FDR_class')
    assert not camera.duplicated(['collection','contrast_id','term_id']).any()
    assert camera.NGenes.between(10,500).all()
    full=pd.read_csv(OUT/'signed_activity_all_programs.tsv',sep='\t')
    omit=pd.read_csv(OUT/'signed_activity_three_omissions.tsv.gz',sep='\t')
    signed=pd.concat([full,omit],ignore_index=True)
    assert set(omit.analysis)=={'omit_GSM9013891','omit_GSM9013914','omit_GSM9013917'}
    check_bh(signed,['analysis','collection','contrast_id'],'PValue','FDR')
    check_bh(signed,['analysis','collection','contrast_class'],'PValue','FDR_class')
    assert not signed.duplicated(['analysis','collection','contrast_id','term_id']).any()
    assert (signed.CI_low<=signed.estimate).all() and (signed.estimate<=signed.CI_high).all()
    checked=0
    for collection in ['PROGENy','DoRothEA']:
        scores=pd.read_csv(OUT/(collection+'_scores.tsv'),sep='\t').set_index('term_id')
        assert list(scores.columns)==list(samples.sample_id)
        for analysis,tab in signed[signed.collection==collection].groupby('analysis'):
            ss=samples if analysis=='full_126' else samples[samples.sample_id!=analysis.removeprefix('omit_')]
            assert len(ss)==(126 if analysis=='full_126' else 125)
            assert tab.residual_df.eq(len(ss)-42).all()
            for r in inv.itertuples():
                pos=ss.loc[ss.group==r.positive_group,'sample_id'];neg=ss.loc[ss.group==r.negative_group,'sample_id']
                expected=scores[pos].mean(axis=1)-scores[neg].mean(axis=1)
                t=tab[tab.contrast_id==r.contrast_id].set_index('term_id')
                assert set(t.index)==set(scores.index)
                assert np.allclose(t.estimate,expected.loc[t.index],atol=1e-10,rtol=1e-9)
                checked+=len(t)
    families=pd.read_csv(OUT/'program_family_summary.tsv',sep='\t')
    assert len(families)==240
    for r in families.itertuples():
        data=camera if r.method=='CAMERA_membership' else full
        t=data[(data.collection==r.collection)&(data.contrast_id==r.contrast_id)]
        assert len(t)==r.tested and (t.FDR<=.05).sum()==r.within_FDR05 and (t.FDR_class<=.05).sum()==r.class_FDR05
    for cov_file,data in [('unsigned_program_coverage.tsv',camera),('signed_program_coverage.tsv',full)]:
        cov=pd.read_csv(OUT/cov_file,sep='\t')
        assert not cov.duplicated(['collection','term_id']).any()
        for name,g in cov.groupby('collection'):
            terms=set(g.loc[g.eligible,'term_id'])
            for id in inv.contrast_id:
                assert terms==set(data.loc[(data.collection==name)&(data.contrast_id==id),'term_id'])
    result=dict(status='PROGRAM_READBACK_PASSED',camera_rows=len(camera),signed_primary_rows=len(full),
                signed_effects_checked_from_library_means=checked,program_families=240,omission_analyses=3,models_run=0)
    (OUT/'review_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    {'genes':genes,'programs':programs}[sys.argv[1]]()
