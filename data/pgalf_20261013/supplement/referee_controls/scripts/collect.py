from pathlib import Path
import json
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
dirs=sorted(d for d in (ROOT/'results').glob('01001_*') if (d/'DONE.json').exists())
assert len(dirs)==12, f'Only {len(dirs)}/12 complete'
for name in ['outer_cut','population_controls','components','recentred_curves','recentred_targets']:
 data=pd.concat([pd.read_csv(d/(name+'.tsv'),sep='\t',dtype={'step':str}) for d in dirs],ignore_index=True)
 if name=='population_controls':
  # Fixed-support radius limits must exclude particles outside the frozen support,
  # even when the coarse bin containing a frozen fine bin also spans other stars.
  for key,g in data[data.control=='fixed_radial_support'].groupby(['step','nstar']):
   fine=g[g.width_dex==.025].iloc[0]
   for col in ['rmin_kpc','rmax_kpc']:data.loc[g.index,col]=fine[col]
 data.to_csv(ROOT/'results'/(name+'.tsv'),sep='\t',index=False)
print('COLLECTED 12')
