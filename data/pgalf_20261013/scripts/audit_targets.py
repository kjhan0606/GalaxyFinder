"""Read cached seed indices and catalogue counts, without loading particles."""
import csv,json
from collections import defaultdict
from pathlib import Path
import numpy as np
PAPER=Path('/tmp/overleaf-pgalf')
PREVIOUS=Path('/gpfs/kjhan/Others/Jaehyun/2026.Apr/PGalF_major_revision_20261013')
SUPP=PAPER/'supplement'
def read(name):
 return list(csv.DictReader([s for s in (SUPP/name).read_text().splitlines() if s and not s.startswith('#')],delimiter='\t'))
gal=defaultdict(list)
for r in read('catalogue_galaxies.tsv'):gal[(r['step'],int(r['nstar']))].append(r)
rows=[]
for key,rs in gal.items():
 rs.sort(key=lambda r:int(r['galaxy_index']));central=next(r for r in rs if int(r['selected']))
 with np.load(PREVIOUS/'results'/f'{key[0]}_{key[1]}'/'membership.npz') as d:seed=int(d['seed']);n=int(d['n_stars'])
 offset=0;seedgal=-1
 for r in rs:
  if int(r['selected']) or int(r['n_stars'])<100:
   count=int(r['n_stars'])
   if offset<=seed<offset+count:seedgal=int(r['galaxy_index'])
   offset+=count
 rows.append(dict(step=key[0],nstar=key[1],central_index=central['galaxy_index'],seed_index=seed,n_input=n,seed_catalogue_galaxy=seedgal,seed_in_selected_central=int(seedgal==int(central['galaxy_index'])),seed_in_background=int(seedgal<0)))
with (SUPP/'target_membership.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
print(json.dumps(dict(n=len(rows),seed_in_central=sum(r['seed_in_selected_central'] for r in rows),seed_in_background=sum(r['seed_in_background'] for r in rows))))
