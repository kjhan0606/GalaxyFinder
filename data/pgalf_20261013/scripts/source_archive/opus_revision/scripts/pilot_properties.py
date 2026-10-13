"""Supplement the initial pilot using saved images/memberships, without new FoF."""
import socket
import sys
import json
import numpy as np
from measure_halo import ROOT,OLD,load,sb,rv,FINE,fraction,stellar_properties,write

def main():
    if socket.gethostname().lower()!='lageunha':raise SystemExit('Lageunha only')
    step=sys.argv[1] if len(sys.argv)>1 else '00081'
    nstar=int(sys.argv[2]) if len(sys.argv)>2 else 4902218
    dest=ROOT/'results'/f'{step}_{nstar}'
    h=load(step,nstar);scale,n=rv.angular_scale(h['a'],h['z']);coords=[h[k] for k in ['dx','dy','dz']]
    orig=np.load(OLD/'results'/f'{step}_{nstar}'/'membership.npz');fine=np.load(dest/'fine_membership.npz');images=np.load(dest/'images.npz')
    selections={}
    geometry=dict(step=step,nstar=nstar,z=h['z'],scale_arcsec_per_hinv_cMpc=scale,n_pixels=n,seed_pixels={},seed_offsets_arcsec={})
    for ia,los in enumerate(rv.LOS):
        j,k=sb.LOS_AXES[ia];pixels=rv.pixels(coords[j],coords[k],scale,n)
        seed=int(orig['seed']);geometry['seed_pixels'][los]=[int(pixels[0][seed]),int(pixels[1][seed])]
        geometry['seed_offsets_arcsec'][los]=[float(coords[j][seed]*scale),float(coords[k][seed]*scale)]
        for ib,band in enumerate(rv.BAND):
            region=np.unpackbits(images[f'mask_{los}_{band}'],count=n*n).astype(bool).reshape(n,n)
            inside=rv.sample_inside(region,*pixels)
            kernel=sb.psf_kernel(sb.PIX_AS,sb.FWHM_AS,sb.MOFFAT_BETA,sb.WING_A[ib],sb.KERNEL_AS)[0]
            q=fraction(region,kernel,pixels,len(h['mass']));selections[(los,band)]=(inside,q,ib)
    rows=[]
    for b in FINE:
        tag=f'b{b:.2f}';cache=orig if tag in orig else fine
        keep=np.unpackbits(cache[tag],count=len(h['mass'])).astype(bool)
        for (los,band),(inside,q,ib) in selections.items():
            rows.append(dict(step=step,nstar=nstar,z=h['z'],los=los,band=band,b=float(b),**stellar_properties(h,keep,inside,q,ib)))
    write(dest/'stellar_properties.tsv',rows)
    (dest/'drawing_geometry.json').write_text(json.dumps(geometry,indent=2)+'\n')
    print('PILOT_PROPERTIES_COMPLETE',step,nstar,flush=True)

if __name__=='__main__':main()
