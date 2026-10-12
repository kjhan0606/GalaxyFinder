#!/usr/bin/env python3
"""Primary calibration and validation from consistent-distance saved curves.

Training uses every entry in the other outputs after a positional exclusion.
Evaluation includes all 76 entries, irrespective of individual fitted endpoints.
Archived finer-grid and thinning results are controls, never combined with the
corrected-distance grid. No particle calculations are performed by this script.
"""
import csv,json
from collections import defaultdict
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]; SUPP=ROOT/'supplement'
B=np.array([.02,.03,.05,.08,.10,.12,.15,.18,.20]); BOX=44.171267

def read(name):
 return list(csv.DictReader([s for s in (SUPP/name).read_text().splitlines() if s.strip() and not s.startswith('#')],delimiter='\t'))
def write(name,rows):
 with (SUPP/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def q(a):
 a=np.asarray(a,float);a=a[np.isfinite(a)]
 return np.quantile(a,[.16,.5,.84]).tolist() if len(a) else []
def main():
 curves=defaultdict(list)
 for r in read('recomputed_curves.tsv'):
  if r['variant']=='consistent_distance':curves[(r['step'],int(r['nstar']),r['los'],r['band'])].append(r)
 views=[]
 for key,rs in curves.items():
  rs.sort(key=lambda r:float(r['b'])); scores=np.array([float(r['iou']) for r in rs]);j=int(np.argmax(scores))
  assert np.allclose([float(r['b']) for r in rs],B)
  views.append(dict(key=key,rows=rs,j=j,b=float(B[j]),edge=int(j in (0,len(B)-1))))
 keys=sorted(set(v['key'][:2] for v in views));lookup={v['key']:v for v in views}
 assert len(keys)==76 and len(views)==456
 good={(*k,band):all(not lookup[(*k,los,band)]['edge'] for los in 'xyz') for k in keys for band in ['V','g']}
 byband={}; outputs=[];viewtable=[]
 for v in views:
  r=v['rows'][v['j']]
  viewtable.append(dict(step=v['key'][0],nstar=v['key'][1],los=v['key'][2],band=v['key'][3],b_best=v['b'],edge=v['edge'],area_arcsec2=r['area_arcsec2'],empty=int(float(r['area_arcsec2'])==0),one_pixel=int(np.isclose(float(r['area_arcsec2']),.04)),iou=r['iou'],completeness=r['completeness'],precision=r['precision']))
 for band in ['V','g']:
  kk=[k for k in keys if good[(*k,band)]]
  vv=[[lookup[(*k,los,band)] for los in 'xyz'] for k in kk]
  def stat(field,node=None):
   return q([np.median([float(v['rows'][v['j'] if node is None else int(np.argmin(abs(B-node)))][field]) for v in arr]) for arr in vv])
  matched=[np.median([v['b'] for v in arr]) for arr in vv]
  widths=[np.median([float(np.ptp(B[np.array([float(r['iou']) for r in v['rows']])>=float(v['rows'][v['j']]['iou'])-.02])) for v in arr]) for arr in vv]
  ranges=[np.ptp([v['b'] for v in arr]) for arr in vv]
  one=[]; losses=[]; n_view_different=0
  for arr,bb in zip(vv,matched):
   j=int(np.argmin(abs(B-bb)));score=np.median([float(v['rows'][j]['iou']) for v in arr]); fit=np.median([float(v['rows'][v['j']]['iou']) for v in arr]);one.append(score);losses.append(fit-score);n_view_different+=int(np.ptp([v['b'] for v in arr])>1e-8)
  fixed={}
  for b in [.08,.10,.20]:
   j=int(np.argmin(abs(B-b)))
   ious=[np.median([float(v['rows'][j]['iou']) for v in arr]) for arr in vv]
   absdiff=[np.median([abs(float(v['rows'][j]['L_group'])/float(v['rows'][j]['L_inside'])-1) for v in arr]) for arr in vv]
   fixviews=[v['rows'][j] for arr in vv for v in arr]
   fixed[str(b)]=dict(iou=q(ious),absolute_luminosity_error=q(absdiff),completeness_view=q([float(r['completeness']) for r in fixviews]),precision_view=q([float(r['precision']) for r in fixviews]),n_below_090=sum(a<.9 for a in ious))
  gain=[];improved=0;errimproved=0
  for arr in vv:
   a=np.median([float(v['rows'][4]['iou']) for v in arr]);c=np.median([float(v['rows'][-1]['iou']) for v in arr]);gain.append(a-c)
   improved+=sum(float(v['rows'][4]['iou'])>float(v['rows'][-1]['iou']) for v in arr)
   e=lambda j:np.median([abs(float(v['rows'][j]['L_group'])/float(v['rows'][j]['L_inside'])-1) for v in arr])
   errimproved+=int(e(4)<e(-1))
  boundary={str(b):{m:stat(m,b) for m in ['ring_iou','ring_completeness','ring_precision']} for b in [.08,.10,.20]}
  boundary['selected']={m:stat(m) for m in ['ring_iou','ring_completeness','ring_precision']}
  byband[band]=dict(boundary_defined_counts={str(b):{m:sum(np.isfinite(np.median([float(v['rows'][int(np.argmin(abs(B-b)))][m]) for v in arr])) for arr in vv) for m in ['ring_iou','ring_completeness','ring_precision']} for b in B},n=len(kk),matched_b=q(matched),fitted_iou=stat('iou'),view_range=q(ranges),width=q(widths),single_halo_iou=q(one),single_halo_loss=q(losses),max_single_halo_loss=max(losses),n_different_views=n_view_different,fixed=fixed,gain_b010_over_b020=q(gain),n_improved_views=improved,n_luminosity_improved=errimproved,boundary=boundary,n_fitted_below_090=sum(np.median([float(v['rows'][v['j']]['iou']) for v in arr])<.9 for arr in vv))
  for step in sorted(set(k[0] for k in kk)):
   inds=[i for i,k in enumerate(kk) if k[0]==step]
   outputs.append(dict(step=step,z=lookup[(*kk[inds[0]],'x',band)]['rows'][0]['z'],band=band,n=len(inds),b_median=float(np.median([matched[i] for i in inds])),iou_median=float(np.median([np.median([float(v['rows'][v['j']]['iou']) for v in vv[i]]) for i in inds])),view_range=float(np.median([ranges[i] for i in inds])),width=float(np.median([widths[i] for i in inds]))))
 positions={(r['step'],int(r['nstar'])):np.array([float(r[c]) for c in ['x_hinv_cMpc','y_hinv_cMpc','z_hinv_cMpc']]) for r in read('halo_positions.tsv')}
 halo={(*k,band):np.median([[float(r['iou']) for r in lookup[(*k,los,band)]['rows']] for los in 'xyz'],axis=0) for k in keys for band in ['V','g']}
 folds=[];predictions=[];validation={};perf=[]
 for step in sorted(set(k[0] for k in keys)):
  held=[k for k in keys if k[0]==step];train=[];near=0
  for k in keys:
   if k[0]==step:continue
   sep=np.array([positions[k]-positions[h] for h in held]);sep-=BOX*np.rint(sep/BOX)
   if np.min(np.linalg.norm(sep,axis=1))<.5:near+=1
   else:train.append(k)
  for band in ['V','g']:
   objective=np.median([halo[(*k,band)] for k in train],axis=0);j=int(np.argmax(objective))
   folds.append(dict(held_out_step=step,band=band,n_train=len(train),n_neighbours_removed=near,n_held=len(held),selected_b=float(B[j])))
   for k in held:
    a=halo[(*k,band)];predictions.append(dict(step=k[0],nstar=k[1],band=band,interior=int(good[(*k,band)]),selected_b=float(B[j]),held_out_iou=float(a[j]),iou_b010=float(a[4]),iou_b020=float(a[-1]),best_common_b_iou=float(a.max())))
 for band in ['V','g']:
  rs=[r for r in predictions if r['band']==band]
  validation[band]=dict(n=len(rs),held_out_iou=q([r['held_out_iou'] for r in rs]),gain_over_b020=q([r['held_out_iou']-r['iou_b020'] for r in rs]),n_below_090=sum(r['held_out_iou']<.9 for r in rs),fold_b=sorted(set(r['selected_b'] for r in rs)),conditional_75_iou=q([r['held_out_iou'] for r in rs if r['interior']]))
  for step in sorted(set(r['step'] for r in rs)):
   ss=[r for r in rs if r['step']==step];gain=[r['held_out_iou']-r['iou_b020'] for r in ss]
   perf.append(dict(step=step,z=next(r['z'] for r in outputs if r['step']==step),band=band,n=len(ss),selected_b=ss[0]['selected_b'],iou_median=float(np.median([r['held_out_iou'] for r in ss])),gain_median=float(np.median(gain)),n_improved=sum(g>0 for g in gain),n_below_090=sum(r['held_out_iou']<.9 for r in ss)))
 thinning=[]
 for tag in ['thin2','thin4','thin8']:
  for band in ['V','g']:
   groups=defaultdict(list)
   for line in (SUPP/'sightline_refined.tsv').read_text().splitlines():
    p=line.split()
    if not p or p[0]!=tag or p[5]!=band:continue
    k=(p[1],int(p[3]));a=np.array(p[6:],float).reshape(9,3); j=int(np.argmax(a[:,0]));groups[k].append((j,a[j,0],a[4,0]))
   base=[k for k in keys if good[(*k,band)]]
   nlower=sum(any(j==0 for j,_,_ in groups[k]) for k in base);nupper=sum(any(j==8 for j,_,_ in groups[k]) for k in base);nret=sum(all(j not in [0,8] for j,_,_ in groups[k]) for k in base)
   thinning.append(dict(kind=tag,band=band,n_original=len(base),n_retained=nret,n_lower_halos=nlower,n_upper_halos=nupper,n_below090_at_b010=sum(np.median([x[2] for x in groups[k]])<.9 for k in base),iou_b010_median=float(np.median([np.median([x[2] for x in groups[k]]) for k in base]))))
 for name,rows in [('primary_view_selection.tsv',viewtable),('primary_output_summary.tsv',outputs),('validation_corrected_folds.tsv',folds),('validation_corrected_predictions.tsv',predictions),('validation_output_summary.tsv',perf),('thinning_all_entries.tsv',thinning)]:write(name,rows)
 report=dict(photometry='consistent_distance',grid=B.tolist(),selection=byband,validation=validation,thinning=thinning,training='All other-output entries, excluding positional neighbours. No fitted-interior filter.')
 (SUPP/'primary_summary.json').write_text(json.dumps(report,indent=2,default=lambda value:value.item())+'\n');print(json.dumps(report,indent=2,default=lambda value:value.item()))
if __name__=='__main__':main()
