from pathlib import Path
import pandas as pd,numpy as np
R=Path(__file__).resolve().parents[1];P=R.parent/'PGalF_BCG_physics_20261013';O=R.parent/'PGalF_opus_revision_20261013'
a=pd.read_csv(R/'results/absolute_properties.tsv',sep='\t',dtype={'step':str})
f=pd.read_csv(O/'results/failure_view_catalogue.tsv',sep='\t',dtype={'step':str});f=f[(f.step=='01001')&(f.weighting=='position')]
t=pd.read_csv(R/'results/recentred_targets.tsv',sep='\t');scale=(t.req_proper_kpc/t.req_arcsec).median()
a['mass_removed_fraction']=1-a.mass_msun/a.nstar.map(a[a.b==.2].set_index('nstar').mass_msun)
for band in ['V','g']:
 q=f[f.band==band].groupby('nstar').req_arcsec.median()*scale
 a['isophote_req_'+band+'_kpc']=a.nstar.map(q)
 a['ABmag_'+band]=-2.5*np.log10(a['flux_'+band+'_uJy']/(3631e6))
rows=[]
for n in a.nstar.unique():
 c=pd.read_csv(P/'refinement'/f'01001_{n}'/'centres.tsv',sep='\t').set_index('centre')
 off=c.loc['catalogue','offset_from_peak_kpc'];row=dict(nstar=int(n),offset_3D_kpc=off)
 for b in [.1,.2]:row['offset_over_rhalf_b'+str(b)]=off/a[(a.nstar==n)&(a.b==b)].rhalf_mass_kpc.iloc[0]
 for band in ['V','g']:
  q=f[(f.nstar==n)&(f.band==band)]
  row['median_projected_offset_over_isophote_'+band]=float((q.seed_offset_arcsec/q.req_arcsec).median())
 rows.append(row)
a.to_csv(R/'results/absolute_properties.tsv',sep='\t',index=False)
pd.DataFrame(rows).to_csv(R/'results/normalized_centre_offsets.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).median(numeric_only=True).to_string())
