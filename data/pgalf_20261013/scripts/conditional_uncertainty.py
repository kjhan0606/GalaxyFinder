"""Conditional bootstrap over entries grouped by the ten known positional pairs."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
from pathlib import Path
import csv,json,argparse
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'results'
def read(name):
    return list(csv.DictReader((DATA/name).open(),delimiter='\t'))
def main():
    global DATA
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,default=ROOT/'results')
    parser.add_argument('--out',type=Path,default=ROOT/'results')
    args=parser.parse_args();DATA=args.data;args.out.mkdir(parents=True,exist_ok=True)
    rows=read('photometric_halo_metrics.tsv')
    keys=sorted({(r['step'],r['nstar']) for r in rows})
    parent={k:k for k in keys}
    for r in read('known_positional_pairs.tsv'):
        a=(r['step_a'],r['nstar_a']);b=(r['step_b'],r['nstar_b']);parent[b]=a
    groups={}
    for k in keys:groups.setdefault(parent[k],[]).append(k)
    groups=list(groups.values());rng=np.random.default_rng(20261013)
    draws=[rng.integers(0,len(groups),len(groups)) for _ in range(2000)]
    report=dict(n_entries=len(keys),n_groups=len(groups),n_known_pairs=10,n_draws=len(draws),seed=20261013,
        scope='Conditional resampling of one-volume output entries. Positional pairs are grouped; unrecognized progenitors remain possible. Intervals are not independent-population confidence intervals.',bands={})
    for band in ['V','g']:
        report['bands'][band]={}
        for weighting in ['position','PSF']:
            table={(r['step'],r['nstar']):r for r in rows if r['variant']=='primary' and r['band']==band and r['weighting']==weighting}
            values={'median_selected_b':[],'median_iou_b010':[],'median_gain_b010_over_b020':[],'median_absolute_luminosity_error_b010':[],'median_absolute_luminosity_error_reduction_b010_vs_b020':[]}
            for draw in draws:
                rr=[table[k] for i in draw for k in groups[i]]
                values['median_selected_b'].append(np.median([float(r['b']) for r in rr if int(r['interior'])]))
                values['median_iou_b010'].append(np.median([float(r['iou_b0.10']) for r in rr]))
                values['median_gain_b010_over_b020'].append(np.median([float(r['iou_b0.10'])-float(r['iou_b0.20']) for r in rr]))
                values['median_absolute_luminosity_error_b010'].append(np.median([float(r['error_b0.10']) for r in rr]))
                values['median_absolute_luminosity_error_reduction_b010_vs_b020'].append(np.median([float(r['error_b0.20'])-float(r['error_b0.10']) for r in rr]))
            report['bands'][band][weighting]={k:np.quantile(v,[.025,.5,.975]).tolist() for k,v in values.items()}
    (args.out/'conditional_uncertainty.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Conditional uncertainty summary written to',args.out)
if __name__=='__main__':main()
