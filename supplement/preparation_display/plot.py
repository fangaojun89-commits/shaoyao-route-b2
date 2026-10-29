"""Display all five retained descriptive intersections from the full420-pair table."""
from pathlib import Path
import hashlib,json
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
source=Path(__file__).resolve().parent/'all_420_primary_pairs.tsv'
out=args.output.resolve()
if out.exists():raise FileExistsError('Use a new output directory; previous figures are preserved')
out.mkdir(parents=True)
allpairs=pd.read_csv(source,sep='\t')
assert len(allpairs)==420 and not allpairs.duplicated(['positive_group','term_id']).any()
flag=(allpairs.FDR_class_PR<=.05)&(allpairs.FDR_class_GR<=.05)&(np.sign(allpairs.estimate_PR)==np.sign(allpairs.estimate_GR))
d=allpairs.loc[flag].sort_values(['positive_group','term_id']).copy()
assert list(zip(d.positive_group,d.term_id))==[('JGE02_L3','Androgen'),('JGE02_L3','EGFR'),('JGE07_L3','MAPK'),('JGE07_L3','PI3K'),('JGE08_L3','MAPK')]
assert all(d.solvent=='E') and all(d.ordinal_level==3)
for tag in ['PR','GR']:
 assert np.isfinite(d[[f'estimate_{tag}',f'CI_low_{tag}',f'CI_high_{tag}']]).all().all()
 assert ((d[f'CI_low_{tag}']<=d[f'estimate_{tag}'])&(d[f'CI_high_{tag}']>=d[f'estimate_{tag}'])).all()
d.to_csv(out/'plotted_data_all_five_pairs.tsv',sep='\t',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'ps.fonttype':42})
fig=plt.figure(figsize=(8.8,6.2));ax=fig.add_axes([.305,.34,.43,.47]);qax=fig.add_axes([.765,.34,.22,.47])
y=np.arange(5)[::-1]
labels={'JGE02_L3':'2:1 co-extract†','JGE07_L3':'1:1 separate mixture','JGE08_L3':'1:2 separate mixture'}
for tag,offset,color,marker,name in [('PR',.12,'#0072B2','o','Mixture − peony'),('GR',-.12,'#D55E00','D','Mixture − licorice')]:
 est=d[f'estimate_{tag}'].to_numpy();low=d[f'CI_low_{tag}'].to_numpy();high=d[f'CI_high_{tag}'].to_numpy()
 ax.errorbar(est,y+offset,xerr=[est-low,high-est],fmt=marker,color=color,markersize=5,elinewidth=1.3,capsize=2,label=name)
ax.set(yticks=y,yticklabels=[labels[g]+'\n'+t for g,t in zip(d.positive_group,d.term_id)],xlim=(-.95,1.2),ylim=(-.6,4.6),xlabel='Within-model score difference (nominal 95% CI)')
ax.axvline(0,color='#777',linewidth=.8,zorder=0);ax.grid(axis='x',color='#e5e5e5',linewidth=.6)
ax.tick_params(axis='y',length=0,pad=9)
for s in ['top','right','left']:ax.spines[s].set_visible(False)
fig.legend(*ax.get_legend_handles_labels(),loc='lower left',bbox_to_anchor=(.04,.19),ncol=2,frameon=False,fontsize=9.5)
qax.set(xlim=(0,1),ylim=ax.get_ylim());qax.axis('off')
qax.text(.5,5.10,'Original class-adjusted q',ha='center',fontsize=9.5)
qax.text(.20,4.70,'vs peony',ha='center',fontsize=9.5);qax.text(.78,4.70,'vs licorice',ha='center',fontsize=9.5)
for yy,(_,r) in zip(y,d.iterrows()):
 for xx,tag in [(.2,'PR'),(.78,'GR')]:qax.text(xx,yy,format(r[f'FDR_class_{tag}'],'.3g'),ha='center',va='center',fontsize=9.5)
fig.text(.04,.95,'Preparation responses relative to each single herb',fontsize=13,weight='bold')
fig.text(.04,.89,'Complete five-pair intersection from all 420 preparation–program pairs',fontsize=9.5)
fig.text(.04,.85,'All displayed conditions: 70% ethanol, high ordinal exposure level',fontsize=9.5,color='#444')
fig.text(.04,.13,'Shared mixture libraries; comparator-specific original adjustment families.',fontsize=9.5)
fig.text(.04,.05,'† Both co-extract pairs lose the intersection after GSM9013891 omission.\nThe other three pairs remain; that omission is outside their comparison groups.',fontsize=9.5,linespacing=1.4)
fig.savefig(out/'Figure_preparation_single_herb_effects.pdf');fig.savefig(out/'Figure_preparation_single_herb_effects.png',dpi=260)
plt.close(fig)
rec={'status':'GENERATED_FROM_COMPLETE_FIVE_OF_420_PAIR_INTERSECTION_VISUAL_REVIEW_PENDING','source_bytes':source.stat().st_size,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'all_pairs':420,'displayed_pairs':5,'displayed_estimates':10,'models_refitted':0,'joint_inference_performed':False,'retrospective_illustration':True,'original_intervals_and_adjustment_families_preserved':True}
(out/'SOURCE_AND_DISPLAY.json').write_text(json.dumps(rec,indent=2),encoding='utf-8');print(json.dumps(rec))
