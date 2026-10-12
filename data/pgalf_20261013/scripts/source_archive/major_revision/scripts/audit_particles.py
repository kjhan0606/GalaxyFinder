"""Finite read-only stellar-population audit on Lageunha; no new FoF runs."""
import os
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1')
import csv,json,sys,socket,subprocess,concurrent.futures
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OLD=ROOT.parent/'PGalF_major_revision_20261013'
sys.path.insert(0,str(OLD/'scripts/source'))
import sb_vs_b as sb
import sb_revise as rv
SUPP=OLD/'inputs/supplement'
def read(name):return list(csv.DictReader([s for s in (SUPP/name).read_text().splitlines() if s and not s.startswith('#')],delimiter='\t'))
def one(step,nstar):
 out=ROOT/'results'/f'clip_{step}_{nstar}.json'
 if out.exists():return
 meta=next(r for r in read('halo_metadata.tsv') if r['step']==step and int(r['nstar'])==nstar)
 a,z,gpath,bpath,tag=rv.catalogue_paths(step)
 with np.load(sb.OUT/f'offsets_{tag}_bg.npz') as d:offset=int(d['off'][int(meta['catalogue_row'])])
 pieces=sb.load_halo(gpath,int(meta['byte_offset']));bg=sb.load_halo(bpath,offset,background=True)[0]
 cat,centre,nsat,nobj=sb.parent_stars(pieces,bg);del pieces,bg
 with np.load(SUPP/f'fsps_z{z:.3f}.npz') as d:grids=(*sb.build_time(a),d['age'],d['logz'],d['flux'])
 ages=sb.proper_age(np.clip(np.interp(cat['tp'],grids[1],grids[0]),.02,a),a,grids[0],grids[2])
 masks={'young':ages<grids[3][0],'old':ages>grids[3][-1],'metal_low':cat['zp']<1e-4,'metal_high':cat['zp']>.05}
 _,_,_,_,flux=rv.unwrap_flux(cat,centre,a,grids)
 result=dict(step=step,nstar=nstar,n_input=len(ages),z=z,scope='Non-satellite input stars',counts={},flux_fractions={})
 for key,m in masks.items():result['counts'][key]=int(m.sum());result['flux_fractions'][key]=(flux[m].sum(axis=0)/flux.sum(axis=0)).tolist()
 out.write_text(json.dumps(result,indent=2)+'\n')
 print('DONE',step,nstar,flush=True)
def main():
 if socket.gethostname().lower()!='lageunha':raise SystemExit('Particle audit is restricted to Lageunha.')
 (ROOT/'results').mkdir(exist_ok=True)
 if len(sys.argv)>1:one(sys.argv[1],int(sys.argv[2]));return
 rows=read('halo_metadata.tsv')
 def run(r):
  dest=ROOT/'logs'/f"clip_{r['step']}_{r['nstar']}.log"
  with dest.open('w') as f:
   subprocess.run(['bash','-c','ulimit -v 67108864; exec "$@"','pgalf',sys.executable,str(Path(__file__).resolve()),r['step'],r['nstar']],stdout=f,stderr=subprocess.STDOUT,check=True)
  print('AUDITED',r['step'],r['nstar'],flush=True)
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(run,rows))
 print('COMPLETE',len(rows),flush=True)
if __name__=='__main__':main()
