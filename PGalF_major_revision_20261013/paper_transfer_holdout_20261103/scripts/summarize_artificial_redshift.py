#!/usr/bin/env python3
"""Summarize raw static-population redshift measurements per observed z.

Usage: python summarize_artificial_redshift.py INPUT OUTPUT.csv
INPUT may be the results directory containing 75 versioned run folders, or
the archived combined per-view file artificial_redshift_metrics_v6.tsv.gz.
All views receive equal weight in the detection fraction and mean overlap
across all views. The conditional mean overlap includes detected views only.
Completeness is undefined for non-detections and excluded from its mean.
"""
from pathlib import Path
import sys
import pandas as pd

if len(sys.argv)!=3:
    raise SystemExit(__doc__)
source=Path(sys.argv[1]); out=Path(sys.argv[2])
raw_by_entry={}
if source.is_file():
    archive=pd.read_csv(source,sep='\t',compression='infer')
    for (step,nstar),raw in archive.groupby(['step','nstar'],sort=True):
        key=f'{int(step):05d}_{int(nstar)}'
        raw_by_entry[key]=raw
    if len(raw_by_entry)!=75 or len(archive)!=27000 or any(len(q)!=360 for q in raw_by_entry.values()):
        raise RuntimeError('Archived raw file must contain 75 entries with 360 rows each (27,000 total)')
    entry_names={key:key for key in raw_by_entry}
else:
    runs={}
    for d in source.iterdir():
        if not d.is_dir() or '_static_population_redshift' not in d.name: continue
        key=d.name.split('_static_population_redshift')[0]
        rank=lambda x: int(x.name.split('_v',1)[1].split('_',1)[0]) if '_v' in x.name else 1
        old=runs.get(key)
        if old is None or rank(d)>rank(old): runs[key]=d
    if len(runs)!=75:
        raise RuntimeError(f'Expected 75 locked entries, found {len(runs)}')
    if any('v6_properfixed_simdistance' not in d.name for d in runs.values()):
        raise RuntimeError('Refusing to mix older geometry or photometry runs into the final summary')
    for key,directory in runs.items():
        f=directory/'metrics.tsv'
        if not f.is_file(): raise RuntimeError(f'Missing per-view metrics: {f}')
        raw_by_entry[key]=pd.read_csv(f,sep='\t')
    entry_names={key:directory.name for key,directory in runs.items()}
if len(raw_by_entry)!=75:
    raise RuntimeError(f'Expected 75 locked entries, found {len(raw_by_entry)}')
records=[]
for key,raw in sorted(raw_by_entry.items()):
    entry_name=entry_names[key]
    for z in [.190,.227,.500,1.002,1.990]:
        for mode in ['observer_full','observer_no_tolman_factor','rest_sed_with_distance','rest_sed_no_dimming']:
          for band in ['V','g']:
            q=raw[(raw.flux_mode==mode) & (raw.band==band) &
                  (pd.to_numeric(raw.b,errors='coerce').round(2)==.10) &
                  (raw.grid.round(3)==z)]
            if len(q)!=3:
                raise RuntimeError(f'{entry_name}: expected 3 views for z={z}, {mode}, {band}, found {len(q)}')
            detected=q.area_arcsec2>0
            records.append(dict(entry=key,z_obs=z,flux_mode=mode,band=band,
                n_views=len(q),n_detected=int(detected.sum()),
                detected_iou_sum=float(q.loc[detected,'iou'].sum()),
                detected_completeness_sum=float(q.loc[detected,'completeness'].sum()),
                all_view_iou_sum=float(q.iou.fillna(0).sum())))
df=pd.DataFrame(records)
rows=[]
for (z,mode,band),q in df.groupby(['z_obs','flux_mode','band'],sort=True):
    nd=int(q.n_views.sum()); ndet=int(q.n_detected.sum())
    rows.append(dict(z_obs=z,flux_mode=mode,band=band,n_entries=len(q),n_views=nd,n_detected_views=ndet,
       detection_fraction=ndet/nd,
       mean_iou_detected_views=q.detected_iou_sum.sum()/ndet if ndet else float('nan'),
       mean_completeness_detected_views=q.detected_completeness_sum.sum()/ndet if ndet else float('nan'),
       mean_iou_all_views=q.all_view_iou_sum.sum()/nd,
       arithmetic_error=(q.all_view_iou_sum.sum()/nd -
           (ndet/nd)*(q.detected_iou_sum.sum()/ndet) if ndet else q.all_view_iou_sum.sum()/nd)))
if len(rows)!=40 or any(r['n_entries']!=75 or r['n_views']!=225 for r in rows):
    raise RuntimeError('Expected 40 summaries, each from 75 entries and 225 views')
if max(abs(r['arithmetic_error']) for r in rows)>1e-12:
    raise RuntimeError('All-view mean does not match the detection fraction times the mean among detections')
out.parent.mkdir(parents=True,exist_ok=True);pd.DataFrame(rows).to_csv(out,index=False)
print(f'Wrote {len(rows)} redshift-band summaries from {len(raw_by_entry)} entries to {out}')
