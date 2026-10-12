#!/usr/bin/env python3
"""Summarize luminosity differences from archived radial flux sums."""
import csv
import json
from collections import defaultdict
from pathlib import Path
import numpy as np

SUPP=Path(__file__).resolve().parents[1]/'supplement'


def read(name):
    lines=(SUPP/name).read_text().splitlines()
    return list(csv.DictReader([s for s in lines if s.strip() and not s.startswith('#')],delimiter='\t'))


def q(a):
    return np.quantile(a,[.16,.5,.84]).tolist()


def main():
    selected={(r['step'],int(r['nstar'])) for r in read('halo_metadata.tsv')
              if int(r['interior_V']) and int(r['interior_g'])}
    views=defaultdict(list)
    for r in read('radial_luminosity.tsv'):
        key=(r['step'],int(r['nstar']),r['los'],r['band'],float(r['b']))
        if key[:2] in selected:
            views[key].append(r)
    halos=defaultdict(list)
    for (step,ns,los,band,b),rs in views.items():
        assert len(rs)==4
        li=sum(float(r['L_inside']) for r in rs)
        lg=sum(float(r['L_group']) for r in rs)
        # The outer annuli begin at the equivalent radius, not the contour.
        out=[r for r in rs if int(r['annulus'])>=2]
        outer_excess=sum(float(r['L_group'])-float(r['L_both']) for r in out)
        total_excess=sum(float(r['L_group'])-float(r['L_both']) for r in rs)
        missed=sum(float(r['L_inside'])-float(r['L_both']) for r in rs)
        halos[(step,ns,band,b)].append([lg/li-1.,outer_excess/li,total_excess/li,missed/li,abs(lg/li-1.)])
    summaries={}
    for band in ['V','g']:
        summaries[band]={}
        paired={}
        for b in [.1,.2]:
            vals={key[:2]:np.median(rs,axis=0) for key,rs in halos.items() if key[2:]==(band,b)}
            assert len(vals)==75
            a=np.array(list(vals.values()))
            summaries[band][str(b)]=dict(n=len(a),fractional_luminosity_difference=q(a[:,0]),
                                        absolute_luminosity_difference=q(a[:,4]),
                                        exterior_beyond_equivalent_radius=q(a[:,1]),
                                        all_exterior_light=q(a[:,2]),missed_interior_light=q(a[:,3]))
            paired[b]=vals
        improvements=[paired[.2][k][4]-paired[.1][k][4] for k in paired[.1]]
        summaries[band]['paired_absolute_error_reduction']=q(improvements)
        summaries[band]['n_absolute_error_improved']=int(np.count_nonzero(np.array(improvements)>0))
    (SUPP/'luminosity_summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
    print(json.dumps(summaries,indent=2))


if __name__=='__main__':
    main()
