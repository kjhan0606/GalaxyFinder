from pathlib import Path
import json, math
import numpy as np
import pandas as pd
from scipy.stats import binomtest
ROOT=Path(__file__).resolve().parents[1];PREV=ROOT.parent/'PGalF_opus_revision_20261013';OUT=ROOT/'results'
def stats(a):
 x=np.array(a,float);x=x[np.isfinite(x)];x.sort();n=len(x)
 return dict(n=n,median=float(np.median(x)),p16_p84=np.quantile(x,[.16,.84]).tolist(),median_ci96=[float(x[2]),float(x[-3])] if n==12 else None,positive=int((x>0).sum()),negative=int((x<0).sum()),sign_p=float(binomtest(int((x>0).sum()),int((x!=0).sum()),.5).pvalue) if (x!=0).any() else 1.)
n=pd.read_csv(OUT/'outer_cut.tsv',sep='\t');r=pd.read_csv(OUT/'population_controls.tsv',sep='\t');comp=pd.read_csv(OUT/'components.tsv',sep='\t')
s={'outer_cut':{k:stats(n[k]) for k in ['removed_fraction','actual_change','null_change','residual_pp']},'population':{},'components':{}}
for (name,width),g in r.groupby(['control','width_dex']):
 s['population'][name+'_'+str(width)]={k:stats(g[k]) for k in ['age_fractional','age_delta','Z_fractional','logZ_delta','coverage','radius_fractional']}
for name,g in comp.groupby('component'):s['components'][name]={k:stats(g[k]) for k in ['age_Gyr','Z_Zsun','rmean_kpc','sigma_kms','central_mass_fraction','small_galaxy_mass_fraction','background_mass_fraction']}
c=pd.read_csv(OUT/'recentred_curves.tsv',sep='\t',dtype={'step':str});p=pd.read_csv(OUT/'baseline_curves.tsv',sep='\t',dtype={'step':str});p=p[(p.variant=='primary')&(p.step=='01001')]
joined=c.merge(p,on=['step','nstar','band','los','b'],suffixes=('_recentered','_original'),validate='one_to_one');joined['delta_iou']=joined.iou_recentered-joined.iou_original
s['recentering']={}
for band,g in joined.groupby('band'):
 rec={}
 for bv in [.1,.11,.2]:
  a=g[g.b==bv];h=a.groupby('nstar')[['iou_original','iou_recentered','delta_iou']].median()
  rec[str(bv)]=dict(original_median=float(h.iou_original.median()),recentered_median=float(h.iou_recentered.median()),median_paired_gain=float(h.delta_iou.median()),original_failed_views=int((a.iou_original<.9).sum()),recentered_failed_views=int((a.iou_recentered<.9).sum()),original_failed_entries=int((h.iou_original<.9).sum()),recentered_failed_entries=int((h.iou_recentered<.9).sum()))
 for suffix in ['original','recentered']:
  view=g.loc[g.groupby(['nstar','los'])['iou_'+suffix].idxmax()].copy();fit=view.groupby('nstar').filter(lambda a:((a.b>.02)&(a.b<.2)).all())
  rec[suffix+'_n_fitted']=fit.nstar.nunique();rec[suffix+'_median_selected_b']=float(fit.groupby('nstar').b.median().median())
 s['recentering'][band]=rec
joined.to_csv(OUT/'recentred_paired_curves.tsv',sep='\t',index=False)
(OUT/'control_summary.json').write_text(json.dumps(s,indent=2)+'\n')
for k,v in s['outer_cut'].items():print(k,v['median'],v['p16_p84'])
for k,v in s['population'].items():print(k,'age',v['age_fractional']['median'],'Z',v['Z_fractional']['median'],'logZ',v['logZ_delta']['median'],'coverage',v['coverage']['median'],'signs',v['age_fractional']['positive'],v['Z_fractional']['negative'])
print('RECENTER',json.dumps(s['recentering']))
