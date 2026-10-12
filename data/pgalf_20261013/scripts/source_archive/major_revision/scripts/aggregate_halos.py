"""Aggregate only after every planned halo has completed baseline reproduction.

Outputs remain in the revision directory for inspection before updating the paper.
"""
import csv
import json
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SUPP=ROOT/'inputs/supplement'


def read(path):
    return list(csv.DictReader([s for s in path.read_text().splitlines()
                               if s.strip() and not s.startswith('#')],delimiter='\t'))


def write(path,rows):
    with path.open('w') as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def q(values):
    a=np.asarray(values,float)
    finite=a[np.isfinite(a)]
    return dict(n=len(a),n_finite=len(finite),
                p16_median_p84=np.quantile(finite,[.16,.5,.84]).tolist() if len(finite) else [])


def main():
    meta=read(SUPP/'halo_metadata.tsv')
    legacy=read(ROOT/'inputs/boundary_ring_legacy.tsv')
    original={(r['step'],int(r['nstar']),r['los'],r['band']):r for r in legacy}
    baseline_good={(r['step'],int(r['nstar'])):bool(int(r['interior_V']) and int(r['interior_g'])) for r in meta}
    allcurves=[]
    allcatalogues=[]
    for r in meta:
        dest=ROOT/'results'/f"{r['step']}_{r['nstar']}"
        if not (dest/'DONE').exists():
            raise RuntimeError(f'Incomplete halo {dest.name}. No partial sample summary permitted.')
        allcurves.extend(read(dest/'curves.tsv'))
        allcatalogues.append(json.loads((dest/'catalogue_audit.json').read_text()))
    curves=defaultdict(list)
    for r in allcurves:
        key=(r['step'],int(r['nstar']),r['los'],r['band'],r['variant'])
        curves[key].append(r)
    rings=[]
    chosen=[]
    for key,rs in curves.items():
        rs.sort(key=lambda r:float(r['b']))
        scores=np.array([float(r['iou']) for r in rs])
        finite=np.isfinite(scores)
        best=int(np.argmax(np.where(finite,scores,-np.inf)))
        selected=rs[best]
        chosen.append(dict(step=key[0],nstar=key[1],los=key[2],band=key[3],variant=key[4],
                           b_best=float(selected['b']),iou_best=float(selected['iou']),
                           edge=int(not finite.any() or best in [0,len(rs)-1])))
        if key[4]!='baseline':
            continue
        ref=original[key[:4]]
        out={k:ref[k] for k in ['step','z','nstar','los','band','area_arcsec2','b_best','iou_best','edge']}
        for name,b in [('best',float(ref['b_best'])),('0.10',.1),('0.20',.2)]:
            row=next(r for r in rs if np.isclose(float(r['b']),b))
            for metric in ['ring_iou','ring_completeness','ring_precision']:
                out[f'{metric}_{name}']=float(row[metric])
        rings.append(out)
    write(ROOT/'results/corrected_boundary_ring.tsv',rings)
    write(ROOT/'results/recomputed_curves.tsv',allcurves)
    write(ROOT/'results/selected_parameters.tsv',chosen)
    write(ROOT/'results/catalogue_audit.tsv',allcatalogues)
    summary={}
    for band in ['V','g']:
        good=[r for r in rings if r['band']==band and baseline_good[(r['step'],int(r['nstar']))]]
        group=defaultdict(list)
        for r in good:
            group[(r['step'],int(r['nstar']))].append(r)
        stats={}
        for metric in ['ring_iou_best','ring_completeness_best','ring_precision_best']:
            # An undefined view remains NaN. Do not silently drop a view.
            stats[metric]=q([np.median([r[metric] for r in rs]) for rs in group.values()])
        normalized=defaultdict(list)
        for r in chosen:
            if r['variant']=='consistent_distance' and r['band']==band:
                normalized[(r['step'],int(r['nstar']))].append(r)
        db=[]
        matches=[]
        newvalid=[]
        selected_b=[]
        delta_iou=[]
        for key,rs in normalized.items():
            if len(rs)!=3 or any(r['edge'] for r in rs):
                continue
            newvalid.append(key)
            selected_b.append(np.median([r['b_best'] for r in rs]))
            if baseline_good[key]:
                old=np.median([float(original[(*key,los,band)]['b_best']) for los in 'xyz'])
                db.append(np.median([r['b_best'] for r in rs])-old)
                current=np.median([r['iou_best'] for r in rs])
                matches.append(current)
                before=[]
                for los in 'xyz':
                    chosen_b=float(original[(*key,los,band)]['b_best'])
                    row=next(r for r in curves[(*key,los,band,'baseline')]
                             if np.isclose(float(r['b']),chosen_b))
                    before.append(float(row['iou']))
                delta_iou.append(current-np.median(before))
        stats['consistent_distance']=dict(n_interior=len(newvalid),halo_b=q(selected_b),
                                          n_originally_retained_lost=75-len(db),
                                          n_newly_retained=len(newvalid)-len(db),
                                          n_paired_b_changed=int(np.count_nonzero(np.abs(db)>1e-8)),
                                          paired_delta_b=q(db),paired_delta_iou=q(delta_iou),paired_iou=q(matches))
        summary[band]=stats
    summary['catalogue']=dict(n_halos=len(allcatalogues),
                              count_reconstruction_matches=sum(r['nstar']==r['reconstructed_count'] for r in allcatalogues),
                              central_is_brightest_V=sum(r['selected_index']==r['brightest_V'] for r in allcatalogues),
                              central_is_brightest_g=sum(r['selected_index']==r['brightest_g'] for r in allcatalogues),
                              central_is_most_massive=sum(r['selected_index']==r['most_massive'] for r in allcatalogues))
    (ROOT/'results/revision_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
