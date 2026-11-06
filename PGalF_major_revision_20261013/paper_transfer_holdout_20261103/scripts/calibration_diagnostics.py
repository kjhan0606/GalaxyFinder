#!/usr/bin/env python3
"""Compute held-out output scores from saved tables.

No particle calculations. Training excludes all entries from the held-out
output and any remaining halo within 0.5 h^-1 cMpc of a held-out halo.
Positional separation cannot replace a merger tree.
"""
import csv
import json
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SUPP=ROOT/'supplement'
B=np.array([.02,.03,.05,.06,.07,.08,.09,.10,.11,.12,.13,.14,.15,.18,.20])
BOX=44.171267


def read_table(name):
    lines=(SUPP/name).read_text().splitlines()
    return list(csv.DictReader([s for s in lines if s.strip() and not s.startswith('#')],delimiter='\t'))


def quantiles(a):
    return np.quantile(a,[.16,.5,.84]).tolist()


def main():
    meta=read_table('halo_metadata.tsv')
    keys=[(r['step'],int(r['nstar'])) for r in meta]
    good={(r['step'],int(r['nstar'])):bool(int(r['interior_V']) and int(r['interior_g'])) for r in meta}
    positions={(r['step'],int(r['nstar'])):np.array([float(r[k]) for k in
                ['x_hinv_cMpc','y_hinv_cMpc','z_hinv_cMpc']]) for r in read_table('halo_positions.tsv')}
    curves={}
    for line in (SUPP/'sightline_refined.tsv').read_text().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        p=line.split()
        if p[0]!='POS':
            continue
        key=(p[1],int(p[3]),p[4],p[5])
        curves[key]=np.asarray(p[6:],float).reshape(len(B),3)
    halo={}
    for key in keys:
        for band in ['V','g']:
            halo[(*key,band)]=np.median([curves[(*key,los,band)][:,0] for los in 'xyz'],axis=0)
    folds=[]
    predictions=[]
    for step in sorted({k[0] for k in keys}):
        held=[k for k in keys if k[0]==step]
        train=[]
        n_near=0
        for k in keys:
            if k[0]==step or not good[k]:
                continue
            sep=np.array([positions[k]-positions[h] for h in held])
            sep-=BOX*np.rint(sep/BOX)
            if np.min(np.linalg.norm(sep,axis=1))<.5:
                n_near+=1
            else:
                train.append(k)
        for band in ['V','g']:
            objective=np.median([halo[(*k,band)] for k in train],axis=0)
            selected=int(np.argmax(objective))
            folds.append(dict(held_out_step=step,band=band,n_train=len(train),
                              n_neighbours_removed=n_near,n_held=len(held),
                              selected_b=float(B[selected])))
            for k in held:
                a=halo[(*k,band)]
                predictions.append(dict(step=k[0],nstar=k[1],band=band,
                                        baseline_interior=int(good[k]),selected_b=float(B[selected]),
                                        held_out_iou=float(a[selected]),
                                        iou_b010=float(a[np.argmin(abs(B-.10))]),
                                        iou_b020=float(a[-1]),
                                        fitted_iou=float(a.max())))
    results={}
    for band in ['V','g']:
        rs=[r for r in predictions if r['band']==band and r['baseline_interior']]
        results[band]=dict(n=len(rs),held_out_iou=quantiles([r['held_out_iou'] for r in rs]),
                           gain_over_b020=quantiles([r['held_out_iou']-r['iou_b020'] for r in rs]),
                           shortfall_from_fitted=quantiles([r['fitted_iou']-r['held_out_iou'] for r in rs]),
                           fold_b=sorted({r['selected_b'] for r in rs}))
        all_rows=[r for r in predictions if r['band']==band]
        results[band]['all_entries']=dict(
            n=len(all_rows),
            held_out_iou=quantiles([r['held_out_iou'] for r in all_rows]),
            gain_over_b020=quantiles([r['held_out_iou']-r['iou_b020'] for r in all_rows]),
            n_below_090=sum(r['held_out_iou']<.9 for r in all_rows))
    for name,rs in [('validation_folds.tsv',folds),('validation_predictions.tsv',predictions)]:
        with (SUPP/name).open('w') as fh:
            writer=csv.DictWriter(fh,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n')
            writer.writeheader()
            writer.writerows(rs)
    (SUPP/'validation_summary.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))


if __name__=='__main__':
    main()
