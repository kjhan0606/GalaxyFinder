from pathlib import Path
import json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];PREV=ROOT.parent/'PGalF_opus_revision_20261013';PHYS=ROOT.parent/'PGalF_BCG_physics_20261013'
c=pd.read_csv(PREV/'results/photometric_curves.tsv',sep='\t',dtype={'step':str})
rows=[]
for (variant,band),g in c[(c.b==.1)&c.variant.isin(['primary','light_peak'])].groupby(['variant','band']):
 h=g.groupby(['step','nstar']).iou.median()
 rows.append(dict(variant=variant,band=band,n_entries=len(h),median_iou=h.median(),failed_entries=int((h<.9).sum()),n_views=len(g),failed_views=int((g.iou<.9).sum())))
pd.DataFrame(rows).to_csv(ROOT/'results/light_peak_all_entries.tsv',sep='\t',index=False)
# Explain different maximization/selection statistics at z=.227, on identical samples.
h=pd.read_csv(PREV/'results/photometric_halo_metrics.tsv',sep='\t',dtype={'step':str})
h=h[(h.step=='00967')&(h.variant=='primary')&(h.weighting=='position')]
r=[]
for band,g in h.groupby('band'):
 for subset,q in [('all',g),('interior_only',g[g.interior==1])]:
  r.append(dict(band=band,subset=subset,n=len(q),view_maxima_median=q.fitted_view_iou.median(),best_common_median=q.best_common_b_iou.median(),selected_common_median=q.selected_halo_iou.median(),median_paired_loss=(q.fitted_view_iou-q.selected_halo_iou).median()))
pd.DataFrame(r).to_csv(ROOT/'results/z0227_score_definitions.tsv',sep='\t',index=False)
props=[]
for d in sorted((PHYS/'results').glob('01001_*')):
 if not d.is_dir() or not (d/'properties.tsv').exists():continue
 t=pd.read_csv(d/'properties.tsv',sep='\t',dtype={'step':str});t=t[(t.aperture_kpc==0)&t.b.isin([.1,.2])].copy()
 for band in ['V','g']:
  t['Rhalf_'+band+'_median_kpc']=t[[f'Rhalf_{band}_{los}_kpc' for los in 'xyz']].median(axis=1)
  t['flux_'+band+'_uJy']=t['L_'+band]*3631e6
 cols=['step','nstar','b','mass_msun','rhalf_mass_kpc','Rhalf_V_median_kpc','Rhalf_g_median_kpc','flux_V_uJy','flux_g_uJy','age_mass_Gyr','age_V_Gyr','age_g_Gyr','metal_mass_Zsun','sigma_mass_kms']
 props.append(t[cols])
pd.concat(props).to_csv(ROOT/'results/absolute_properties.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).to_string(index=False));print(pd.DataFrame(r).to_string(index=False))
