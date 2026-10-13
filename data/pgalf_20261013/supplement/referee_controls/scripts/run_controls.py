from pathlib import Path
import os,json,sys,socket,subprocess,concurrent.futures
ROOT=Path(__file__).resolve().parents[1]
if socket.gethostname().lower()!='lageunha':raise SystemExit('Lageunha only')
pairs=sorted([d.name for d in (ROOT.parent/'PGalF_BCG_physics_20261013/refinement').glob('01001_*') if (d/'DONE.json').exists()])
assert len(pairs)==12
state={'total':12,'workers':2,'completed':[],'status':'running'}
def save():
 p=ROOT/'results/status.json';tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(state,indent=2));tmp.replace(p)
def one(key):
 step,nstar=key.split('_')
 with (ROOT/'logs'/f'{key}.log').open('a') as f:subprocess.run(['bash','-c','ulimit -v 67108864; exec "$@"','controls',sys.executable,'-u',str(ROOT/'scripts/measure_controls.py'),step,nstar],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=14400)
 return key
save()
try:
 pilot=pairs[0];one(pilot);state['completed'].append(pilot);save();print('PILOT_COMPLETE',pilot,flush=True)
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
  for fut in concurrent.futures.as_completed([pool.submit(one,k) for k in pairs[1:]]):state['completed'].append(fut.result());save();print('COMPLETED',len(state['completed']),flush=True)
 state['status']='complete';save()
except BaseException as e:state.update(status='failed',error=str(e));save();raise
