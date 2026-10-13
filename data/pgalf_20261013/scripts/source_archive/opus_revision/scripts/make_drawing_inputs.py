"""Export bounded corrected map cutouts and replace the cluster contour."""
from pathlib import Path
import csv,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'results'
CASES=[('00475',5557744,'g','z',5),('00967',17494810,'V','x',60),('00967',17494810,'g','z',60),
       ('00178',7062464,'V','x',12),('00178',7062464,'V','y',12),('00178',7062464,'g','z',12)]
def main():
    curve_file=DATA/'photometric_curves.tsv'
    if curve_file.is_file():rows=list(csv.DictReader(curve_file.open(),delimiter='\t'))
    else:
        rows=[]
        for step,nstar in sorted({(c[0],c[1]) for c in CASES}):
            rows+=list(csv.DictReader((DATA/f'{step}_{nstar}'/'fine_curves.tsv').open(),delimiter='\t'))
    arrays={};metadata=[]
    for i,(step,nstar,band,los,half) in enumerate(CASES):
        d=DATA/f'{step}_{nstar}';geom=json.loads((d/'drawing_geometry.json').read_text());n=geom['n_pixels']
        with np.load(d/'images.npz') as image, np.load(d/'drawing_membership.npz') as member:
            p=int(round(half/.2));lo=n//2-p;hi=n//2+p
            arr=image[f'image_{los}_{band}'][lo:hi,lo:hi]
            mu=np.full_like(arr,40.);pos=arr>0;mu[pos]=-2.5*np.log10(arr[pos]/.2**2)
            mask=np.unpackbits(image[f'mask_{los}_{band}'],count=n*n).reshape(n,n)[lo:hi,lo:hi]
            keep=member['b0.10'].astype(bool)
            with np.load(d/'drawing_positions.npz') as positions:
                x=positions['x_'+los][keep];y=positions['y_'+los][keep]
            on=(abs(x)<=half)&(abs(y)<=half)
            arrays[f'mu_{i}']=mu;arrays[f'mask_{i}']=mask;arrays[f'x_{i}']=x[on];arrays[f'y_{i}']=y[on]
        row=next(r for r in rows if r['variant']=='primary' and r['step']==step and int(r['nstar'])==nstar and r['band']==band and r['los']==los and float(r['b'])==.1)
        metadata.append(dict(step=step,nstar=nstar,band=band,los=los,half_arcsec=half,z=float(row['z']),iou=float(row['iou']),
                             seed_offsets_arcsec=geom['seed_offsets_arcsec'][los]))
    np.savez_compressed(DATA/'failure_drawing.npz',**arrays)
    (DATA/'failure_drawing.json').write_text(json.dumps(metadata,indent=2)+'\n')
    old=ROOT/'inputs/supplement/cluster_illustration.npz'
    with np.load(old) as d:cluster={k:d[k] for k in d.files}
    with np.load(DATA/'01001_38857678/images.npz') as d:
        cluster['region']=np.unpackbits(d['mask_y_g'],count=1800**2).reshape(1800,1800).astype(bool)
    np.savez_compressed(DATA/'cluster_illustration.npz',**cluster)
    print('Exported corrected compact figure inputs')
if __name__=='__main__':main()
