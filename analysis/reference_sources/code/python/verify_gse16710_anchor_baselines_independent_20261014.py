"""Independent cached-score anchor baseline checks, not a science worker.

Standard library only. Reconstructs four-cell sufficient-statistic OLS and
restricted means directly, with an independent incomplete-beta Student tail.
Prepare mode reads existing metadata only, never nonexistent/new results.
Verify mode is used only after the main executor confirms its final output.
All diagnostic copies/receipts go to a separate R-drive directory.
"""
import argparse
from collections import Counter, defaultdict
import csv
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys

FAMILIES={'Reactome':1656,'PROGENy':14,'DoRothEA':101}
GROUPS=['Sham21','Sham60','Spi21','Spi60']
C=[1,-1,-1,1]
CONTRAST='spasticity_progression_interaction'
CHECKS=Counter()
MAX_ABS=defaultdict(float)

def pin(path):
    p=Path(path);h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return {'path':p.as_posix(),'bytes':p.stat().st_size,'sha256':h.hexdigest()}

def read(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def number(s):
    try:
        x=float(s);return x if math.isfinite(x) else None
    except (ValueError,TypeError):return None

def avg(x):return math.fsum(x)/len(x)

def qtile(x,p):
    a=sorted(x)
    if not a:return None
    z=(len(a)-1)*p;i=int(z);f=z-i
    return a[i]+f*(a[min(i+1,len(a)-1)]-a[i])

def ranks(x):
    positions=defaultdict(list)
    for i,v in enumerate(x):positions[v].append(i)
    result=[None]*len(x);previous=0
    for v in sorted(positions):
        rr=previous+(len(positions[v])+1)/2
        for i in positions[v]:result[i]=rr
        previous+=len(positions[v])
    return result

def spearman(x,y):
    pair=[(a,b) for a,b in zip(x,y) if a is not None and b is not None]
    if len(pair)<2:return None
    a,b=zip(*pair);a=ranks(a);b=ranks(b);middle=(len(a)+1)/2
    aa=[v-middle for v in a];bb=[v-middle for v in b]
    denominator=math.sqrt(math.fsum(v*v for v in aa)*math.fsum(v*v for v in bb))
    return math.fsum(v*w for v,w in zip(aa,bb))/denominator if denominator else None

def beta_fraction(a,b,x):
    # Lentz continued fraction, independent of R pt/qt or worker functions.
    qab=a+b;qap=a+1;qam=a-1;floor=1e-300
    c=1.0;d=1-qab*x/qap
    if abs(d)<floor:d=floor
    d=1/d;h=d
    for m in range(1,501):
        m2=2*m;aa=m*(b-m)*x/((qam+m2)*(a+m2))
        d=1+aa*d;c=1+aa/c
        if abs(d)<floor:d=floor
        if abs(c)<floor:c=floor
        d=1/d;h*=d*c
        aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
        d=1+aa*d;c=1+aa/c
        if abs(d)<floor:d=floor
        if abs(c)<floor:c=floor
        d=1/d;change=d*c;h*=change
        if abs(change-1)<=4e-15:return h
    raise ArithmeticError('Incomplete beta did not converge')

def regularized_beta(x,a,b):
    if x<=0:return 0.0
    if x>=1:return 1.0
    coefficient=math.exp(math.lgamma(a+b)-math.lgamma(a)-math.lgamma(b)+a*math.log(x)+b*math.log1p(-x))
    if x<(a+1)/(a+b+2):return coefficient*beta_fraction(a,b,x)/a
    return 1-coefficient*beta_fraction(b,a,1-x)/b

def t_two_sided(t,df=22):
    if t==0:return 1.0
    return regularized_beta(df/(df+t*t),df/2,.5)

def t_critical_95(df=22):
    lo=0.;hi=64.
    for _ in range(100):
        mid=(lo+hi)/2
        if t_two_sided(mid,df)>.05:lo=mid
        else:hi=mid
    return (lo+hi)/2

def bh(p,n):
    ordered=sorted((v,i) for i,v in enumerate(p) if v is not None)
    result=[None]*len(p);minimum=1.
    for j in range(len(ordered)-1,-1,-1):
        value,i=ordered[j];minimum=min(minimum,n*value/(j+1));result[i]=minimum
    return result

def exact_rank(a):
    if not a:return 0
    m=[[Fraction(v) for v in r] for r in a];pivot=0
    for c in range(len(m[0])):
        k=next((i for i in range(pivot,len(m)) if m[i][c]),None)
        if k is None:continue
        m[k],m[pivot]=m[pivot],m[k]
        for i in range(pivot+1,len(m)):
            if m[i][c]:
                f=m[i][c]/m[pivot][c];m[i]=[u-f*v for u,v in zip(m[i],m[pivot])]
        pivot+=1
        if pivot==len(m):break
    return pivot

def metadata(root):
    ledger=sorted(read(root/'data/processed/GSE16710/sample_ledger_verified.tsv'),key=lambda r:int(r['sample_order']))
    qc=read(root/'data/processed/GSE16710/qc_inclusion_decision.tsv');qmap={r['sample_id']:r for r in qc}
    assert len(ledger)==len(qc)==len(qmap)==len({r['sample_id'] for r in ledger})==27
    assert set(qmap)=={r['sample_id'] for r in ledger}
    retained=[r for r in ledger if qmap[r['sample_id']]['primary_analysis_decision']=='retain']
    excluded=[r['sample_id'] for r in ledger if qmap[r['sample_id']]['primary_analysis_decision']!='retain']
    assert excluded==['GSM418682'] and len(retained)==26
    assert all(r['experimental_arm'] in ['Sham','Spi'] and r['time_day'] in ['21','60'] for r in retained)
    samples=[r['sample_id'] for r in retained];groups=[r['experimental_arm']+r['time_day'] for r in retained]
    assert Counter(groups)==Counter({'Sham21':6,'Sham60':4,'Spi21':8,'Spi60':8})
    design=read(root/'results/GSE16710/diagnostics/design_matrix_primary.tsv')
    columns=['(Intercept)','injurySpi','time_factor60','injurySpi:time_factor60']
    assert list(design[0])==['sample_id']+columns and [r['sample_id'] for r in design]==samples
    expected=[[1,int(r['experimental_arm']=='Spi'),int(r['time_day']=='60'),int(r['experimental_arm']=='Spi' and r['time_day']=='60')] for r in retained]
    assert [[int(r[c]) for c in columns] for r in design]==expected
    assert exact_rank(expected)==4 and exact_rank([r[:3] for r in expected])==3
    missing={key:sum(r[key] in ['','NA','N/A','unknown'] for r in retained) for key in ['batch','biosample','sra_experiment','biological_unit_id','statistical_unit','experimental_arm','time_day','sample_id']}
    warning={'batch_missing_n':missing['batch'],'biosample_missing_n':missing['biosample'],'sra_experiment_missing_n':missing['sra_experiment'],
        'biological_unit_id_equals_GSM_n':sum(r['biological_unit_id']==r['sample_id'] for r in retained),
        'statistical_unit_labels':dict(Counter(r['statistical_unit'] for r in retained)),
        'interpretation':'GSM-derived animal labels are not independently authenticated subject linkage; batch/independence unresolved; not a train/test task'}
    return samples,groups,dict(metadata_rows_before=27,metadata_rows_after=27,arrays_in_cached_analysis=26,
        historical_QC_excluded=excluded,records_deleted_by_validator=0,group_counts=dict(Counter(groups)),missing_metadata=missing,warnings=warning,
        metadata_key_alignment='PASS',broader_independence_validation='NEEDS REVIEW')

def reconstruct(y,samples,groups):
    fields=['estimate','se','t','p','ci_low','ci_high','mse','sse_full','sse_reduced','nested_f','group_effect_check','se_formula_check',
        'rank_effect','rank_unique','rank_tied_observations','max_abs_studentized','max_cook','max_cook_sample','group_variance_ratio']
    out={k:None for k in fields};out['status']='nonfinite';out['negative_sse_roundoff']=False
    for g in GROUPS:
        for prefix in ['mean_','reduced_mean_','variance_','rank_mean_']:out[prefix+g]=None
    if any(v is None for v in y):return out,[]
    ix={g:[i for i,k in enumerate(groups) if k==g] for g in GROUPS}
    n=[len(ix[g]) for g in GROUPS];means=[avg([y[i] for i in ix[g]]) for g in GROUPS]
    delta=math.fsum(c*m for c,m in zip(C,means));kappa=math.fsum(1/z for z in n)
    restricted=[m-c*delta/(z*kappa) for m,c,z in zip(means,C,n)]
    residual=[y[i]-means[GROUPS.index(groups[i])] for i in range(26)]
    sse=math.fsum(r*r for r in residual);sse0=math.fsum((y[i]-restricted[GROUPS.index(groups[i])])**2 for i in range(26))
    total=math.fsum((z-avg(y))**2 for z in y);degenerate=len(set(y))==1 or sse<=64*sys.float_info.epsilon*total
    rr=ranks(y);rmeans=[avg([rr[i] for i in ix[g]]) for g in GROUPS]
    variances=[math.fsum((y[i]-m)**2 for i in ix[g])/(len(ix[g])-1) for g,m in zip(GROUPS,means)]
    out.update(status='numerically_degenerate' if degenerate else 'ok_nominal_working_model',estimate=delta,mse=sse/22,sse_full=sse,sse_reduced=sse0,
        rank_effect=math.fsum(c*m for c,m in zip(C,rmeans)),rank_unique=len(set(rr)),rank_tied_observations=sum(Counter(y)[z]>1 for z in y),
        group_effect_check=0.,group_variance_ratio=max(variances)/min(variances) if min(variances)>0 else None)
    for i,g in enumerate(GROUPS):
        for prefix,values in [('mean_',means),('reduced_mean_',restricted),('variance_',variances),('rank_mean_',rmeans)]:out[prefix+g]=values[i]
    cook=[]
    if not degenerate:
        mse=sse/22;se=math.sqrt(mse*kappa);t=delta/se;critical=t_critical_95()
        student=[r/math.sqrt(mse*(1-1/n[GROUPS.index(g)])) for r,g in zip(residual,groups)]
        cook=[u*u*(1/n[GROUPS.index(g)])/(4*(1-1/n[GROUPS.index(g)])) for u,g in zip(student,groups)]
        maximum=max(cook)
        out.update(se=se,t=t,p=t_two_sided(t),ci_low=delta-critical*se,ci_high=delta+critical*se,
            nested_f=max(0,sse0-sse)/mse,se_formula_check=0.,max_abs_studentized=max(abs(u) for u in student),
            max_cook=maximum,max_cook_sample=samples[cook.index(maximum)])
    return out,cook

def compare(actual,expected,label):
    CHECKS[label]+=1
    if expected is None:assert actual in [None,'NA','','NaN'],(label,actual,'unavailable expected')
    elif isinstance(expected,bool):assert str(actual).upper()==str(expected).upper(),(label,actual,expected)
    elif isinstance(expected,str):assert actual==expected,(label,actual,expected)
    elif isinstance(expected,int):assert int(actual)==expected,(label,actual,expected)
    else:
        value=number(actual);assert value is not None,(label,actual,expected)
        error=abs(value-expected);MAX_ABS[label]=max(MAX_ABS[label],error)
        atol,rtol=(5e-13,1e-10) if label.endswith('.p') or label.endswith('q_value_BH_nominal') else (1e-10,1e-9)
        assert error<=atol+rtol*abs(expected),(label,actual,expected,error,atol+rtol*abs(expected))

def overlap(p1,p2):
    ix=[i for i,(a,b) in enumerate(zip(p1,p2)) if a is not None and b is not None];n=len(ix)
    if not n:return dict(paired_n=0,rho=None,top_k=0,top_left=0,top_right=0,top_overlap=0)
    k=min(20,n);ba=sorted(p1[i] for i in ix)[k-1];bb=sorted(p2[i] for i in ix)[k-1]
    a={i for i in ix if p1[i]<=ba};b={i for i in ix if p2[i]<=bb}
    return dict(paired_n=n,rho=spearman(p1,p2),top_k=k,top_left=len(a),top_right=len(b),top_overlap=len(a&b))

def compare_structure(actual,expected,label):
    if isinstance(expected,dict):
        assert set(actual)==set(expected),(label,'keys',set(actual)^set(expected))
        for k,v in expected.items():compare_structure(actual[k],v,label+'.'+k)
    else:compare(actual,expected,label)

def math_test():
    for t in [.01,.1,1,5,20]:assert abs(t_two_sided(t,2)-(1-t/math.sqrt(2+t*t)))<3e-14
    critical=t_critical_95();assert abs(critical-2.0738730679040147)<1e-12
    assert ranks([1,1,3])==[1.5,1.5,3.0]
    assert bh([.01,None,.04],3)==[.03,None,.06]
    return {'status':'PASS_SMALL_INDEPENDENT_MATH_ONLY','student22_critical95':critical,'no_existing_scores_used':True}

def verify(root,out,diagnostics):
    assert out.as_posix().startswith('R:/Shaoyao_storage/cache/') and diagnostics.as_posix().startswith('R:/Shaoyao_storage/cache/')
    destination=diagnostics/'INDEPENDENT_NUMERIC_CHECK.json';assert not destination.exists()
    summary=json.loads((out/'RUN_SUMMARY.json').read_text());reg=json.loads((out/'REGISTRY.json').read_text())
    before=read(out/'INPUT_MANIFEST.tsv');after=read(out/'INPUT_MANIFEST_AFTER.tsv');assert before==after
    assert len(before)==21 and len({r['role'] for r in before})==21
    source_files=[]
    for r in before:
        path=(root/r['relative_path']).resolve();assert path.is_relative_to(root.resolve())
        item=pin(path);assert item['sha256']==r['sha256'] and item['bytes']==int(r['bytes']);source_files.append(item)
    assert pin(out/'INPUT_MANIFEST.tsv')['sha256']==reg['source_manifest_sha256']
    assert pin(out/'FEATURE_REGISTRY.tsv')['sha256']==reg['feature_registry_sha256']
    assert pin(out/'SAMPLE_REGISTRY.tsv')['sha256']==reg['sample_registry_sha256']
    expected_product_pins=[]
    for p in summary['outputs']:
        current=pin(out/p['file']);assert current['bytes']==p['bytes'] and current['sha256']==p['sha256'];expected_product_pins.append(current)
    samples,groups,validation=metadata(root);sample_rows=read(out/'SAMPLE_DESIGN.tsv')
    assert [r['sample_id'] for r in sample_rows]==samples and [r['group'] for r in sample_rows]==groups
    assert read(out/'SAMPLE_REGISTRY.tsv')==[{'sample_id':s,'group':g} for s,g in zip(samples,groups)]
    for r in sample_rows:compare(r['leverage'],1/Counter(groups)[r['group']],'leverage')
    scores=read(out/'SCORE_EXPORT_FOR_VERIFICATION.tsv');results=read(out/'BASELINE_RESULTS.tsv');features=read(out/'FEATURE_REGISTRY.tsv')
    assert len(scores)==len(results)==len(features)==1771
    assert list(scores[0])==['collection','feature_id']+samples
    keys=[(r['collection'],r['feature_id']) for r in scores]
    assert len(set(keys))==1771 and keys==[(r['collection'],r['feature_id']) for r in results]==[(r['collection'],r['feature_id']) for r in features]
    assert Counter(r['collection'] for r in scores)==Counter(FAMILIES)
    reconstructed=[];cook_alternate_top=0
    for s,r in zip(scores,results):
        y=[number(s[k]) for k in samples];expected,cooks=reconstruct(y,samples,groups)
        for field,value in expected.items():
            if field=='max_cook_sample' and cooks and r[field]!=value:
                assert r[field] in samples
                actual_cook=cooks[samples.index(r[field])]
                assert abs(actual_cook-max(cooks))<=1e-10*(1+max(cooks)),('Cook maximum sample mismatch',r['collection'],r['feature_id'])
                cook_alternate_top+=1
            elif field=='negative_sse_roundoff':
                actualnegative=str(r[field]).upper()=='TRUE'
                if actualnegative:assert number(r['sse_full'])-number(r['sse_reduced'])<=1e-10*(1+number(r['sse_reduced']))
                CHECKS[field]+=1
            else:compare(r[field],value,'baseline.'+field)
        expected['collection']=r['collection'];expected['feature_id']=r['feature_id'];reconstructed.append(expected)
        assert r['no_interaction_fixed_effect']=='0' and r['rank_inference_status']=='descriptive_only_no_P_q_CI'
        assert r['nominal_inference_status']=='conditional_homoskedastic_Gaussian_array_independence_not_verified'
    roles={r['role']:root/r['relative_path'] for r in before}
    native_fields=['estimate','standard_error','ci_low','ci_high','moderated_t','df_total','p_value','q_value_BH']
    fl=read(roles['observed_t']);camera=[r for r in read(roles['camera']) if r['contrast']==CONTRAST]
    assert len({r['feature_id'] for r in camera})==len(camera)
    cam={r['feature_id']:r for r in camera};family_summaries={}
    for family,n in FAMILIES.items():
        ix=[i for i,r in enumerate(results) if r['collection']==family]
        native=[r for r in read(roles['native_'+family]) if r['contrast']==CONTRAST]
        nmap={r['feature_id']:r for r in native};assert len(nmap)==n and set(nmap)=={results[i]['feature_id'] for i in ix}
        fmap={r['feature_id_original']:r for r in fl if r['family']==family};assert len(fmap)==n and set(fmap)==set(nmap)
        rows=[results[i] for i in ix];derived=[reconstructed[i] for i in ix];p=[r['p'] for r in derived];q=bh(p,n)
        nativep=[];nativeq=[];delta=[];rankdelta=[];tdiff=[];fdiff=[];ndiff=[]
        for j,(r,d) in enumerate(zip(rows,derived)):
            key=r['feature_id'];nr=nmap[key]
            compare(r['q_value_BH_nominal'],q[j],'baseline.q_value_BH_nominal')
            for field in native_fields:compare(r['retained_moderated_'+field],number(nr[field]),'retained.'+field)
            compare(r['retained_FL_observed_t'],number(fmap[key]['ordinary_t_observed']),'retained_FL_t')
            difference=number(r['estimate'])-number(nr['estimate']);compare(r['native_effect_difference'],difference,'native_effect_difference')
            ndiff.append(abs(difference));nativep.append(number(nr['p_value']));nativeq.append(number(nr['q_value_BH']))
            delta.append(d['estimate']);rankdelta.append(d['rank_effect'])
            if d['t'] is not None:
                tdiff.append(abs(number(r['t'])-number(fmap[key]['ordinary_t_observed'])))
                fdiff.append(abs(number(r['nested_f'])-number(r['t'])**2))
                assert tdiff[-1]<=1e-7*(1+abs(number(fmap[key]['ordinary_t_observed'])))
                assert fdiff[-1]<=1e-7*(1+number(r['t'])**2)
            if family=='Reactome':
                cm=cam.get(key)
                for outfield,sourcefield in [('camera_direction','direction'),('camera_p_value','p_value'),('camera_q_value_BH','q_value_BH')]:
                    value=None if cm is None else cm[sourcefield] if sourcefield=='direction' else number(cm[sourcefield])
                    compare(r[outfield],value,'camera.'+outfield)
            else:
                for field in ['camera_direction','camera_p_value','camera_q_value_BH']:compare(r[field],None,'camera_not_applicable')
        directions=[i for i,(a,b) in enumerate(zip(delta,rankdelta)) if a is not None and b is not None and a!=0 and b!=0]
        actual_summary=summary['collections'][family]
        cooked=[d['max_cook'] for d in derived if d['max_cook'] is not None];variances=[d['group_variance_ratio'] for d in derived if d['group_variance_ratio'] is not None]
        expected=dict(n=n,status_counts=dict(Counter(d['status'] for d in derived)),
            native_effect_compared_n=sum(d['estimate'] is not None and number(nmap[r['feature_id']]['estimate']) is not None for r,d in zip(rows,derived)),
            ordinary_t_compared_n=len(tdiff),native_effect_max_abs_difference=max(ndiff,default=None),ordinary_t_max_abs_difference=max(tdiff,default=None),
            f_t_squared_max_abs_difference=max(fdiff,default=None),moderated_vs_OLS=overlap(p,nativep),
            OLS_q_eligible_n=sum(v is not None for v in q),retained_q_eligible_n=sum(v is not None for v in nativeq),
            q_common_eligible_n=sum(a is not None and b is not None for a,b in zip(q,nativeq)),OLS_nominal_q05_n=sum(v is not None and v<=.05 for v in q),
            retained_moderated_q05_n=sum(v is not None and v<=.05 for v in nativeq),q05_overlap_n=sum(a is not None and b is not None and a<=.05 and b<=.05 for a,b in zip(q,nativeq)),
            rank_effect_paired_n=len(directions),rank_effect_direction_agree_n=sum((delta[i]>0)==(rankdelta[i]>0) for i in directions),
            rank_effect_zero_n=sum(a is not None and b is not None and (a==0 or b==0) for a,b in zip(delta,rankdelta)),
            rank_effect_missing_n=sum(a is None or b is None for a,b in zip(delta,rankdelta)),rank_effect_zero_or_missing_n=n-len(directions),
            absolute_native_rank_effect_rho=spearman([None if d is None else abs(d) for d in delta],[None if d is None else abs(d) for d in rankdelta]),
            variance_ratio_quantiles={key:qtile(variances,pct) for key,pct in [('0%',0),('50%',.5),('90%',.9),('100%',1)]},
            max_cook_quantiles={key:qtile(cooked,pct) for key,pct in [('0%',0),('50%',.5),('90%',.9),('100%',1)]},
            tied_score_feature_n=sum(d['rank_tied_observations'] is not None and d['rank_tied_observations']>0 for d in derived),
            negative_sse_roundoff_n=sum(str(r['negative_sse_roundoff']).upper()=='TRUE' for r in rows),camera_vs_OLS=None)
        if family=='Reactome':
            cp=[number(cam[r['feature_id']]['p_value']) if r['feature_id'] in cam else None for r in rows]
            ok=[i for i,r in enumerate(rows) if r['feature_id'] in cam and cam[r['feature_id']]['direction'] in ['Up','Down'] and delta[i] is not None and delta[i]!=0]
            expected['camera_vs_OLS']={**overlap(p,cp), 'direction_paired_n':len(ok),
                'direction_agree_n':sum((delta[i]>0)==(cam[rows[i]['feature_id']]['direction']=='Up') for i in ok),
                'missing_keys':sum(r['feature_id'] not in cam for r in rows),'extra_keys':sum(k not in nmap for k in cam),'null':'competitive, not same score-level null'}
        compare_structure(actual_summary,expected,'summary.'+family);family_summaries[family]=expected
    for old in source_files:
        now=pin(old['path']);assert now==old
    for old in expected_product_pins:assert pin(old['path'])==old
    diagnostics.mkdir(parents=True,exist_ok=True)
    receipt=dict(status='PASS_INDEPENDENT_CACHED_SCORE_NUMERIC_CHECK',verifier=pin(__file__),science_output=out.as_posix(),registry=pin(out/'REGISTRY.json'),run_summary=pin(out/'RUN_SUMMARY.json'),
        source_pins=source_files,products=expected_product_pins,source_and_product_pins_unchanged=True,metadata_validation=validation,features=1771,scores_per_feature=26,
        scalar_checks=dict(CHECKS),maximum_absolute_differences=dict(MAX_ABS),cook_max_sample_roundoff_ties=cook_alternate_top,
        independent_math=math_test(),collections=family_summaries,original_science_model_run=False,R_worker_sourced=False,
        RDS_value_bridge='RDS bytes checked; calculations use complete score export; no independent RDS decoder in standard library',
        working_assumptions='nominal OLS only; rank descriptive; no new primary/candidate/gate or independent validation')
    with destination.open('x',encoding='utf-8') as f:json.dump(receipt,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({'status':receipt['status'],'receipt':pin(destination),'features':1771,'metadata':validation}))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['prepare','verify']);parser.add_argument('--root',type=Path,default=Path(r'<PROJECT_ROOT>'))
    parser.add_argument('--output',type=Path);parser.add_argument('--diagnostics',type=Path,default=Path(r'R:/Shaoyao_storage/cache/GSE16710_anchor_baseline_20261014_review_r01'))
    args=parser.parse_args();assert sys.flags.isolated and sys.flags.dont_write_bytecode
    if args.mode=='prepare':
        _,_,validation=metadata(args.root);value={'status':'PREPARED_NO_REAL_BASELINE_RESULTS_READ','verifier':pin(__file__),'math':math_test(),'metadata_validation':validation}
        assert args.diagnostics.as_posix().startswith('R:/Shaoyao_storage/cache/');args.diagnostics.mkdir(parents=True,exist_ok=True)
        with (args.diagnostics/'PREPARATION.json').open('x',encoding='utf-8') as f:json.dump(value,f,indent=2);f.write('\n')
        print(json.dumps(value))
    else:verify(args.root,args.output,args.diagnostics)
