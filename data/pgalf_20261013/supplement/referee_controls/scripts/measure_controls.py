"""Finite existing-data controls. Read raw catalogues only on Lageunha."""
import os
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1')
from pathlib import Path
import sys,json,csv,socket,time
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
PHYS=ROOT.parent/'PGalF_BCG_physics_20261013'
sys.path.insert(0,str(PHYS/'scripts'))
import measure_bcg as mb
pr,sb,rv=mb.prior,mb.sb,mb.rv

def radial_compare(r,m,age,z,ret,extra,bins,support=None):
    nbin=int(bins.max())+1
    E=(~ret)&extra; R=ret&extra
    ec=np.bincount(bins[E],minlength=nbin);rc=np.bincount(bins[R],minlength=nbin)
    good=(ec>=100)&(rc>=100) if support is None else support
    ew=np.bincount(bins[E],weights=m[E],minlength=nbin)
    rw=np.bincount(bins[R],weights=m[R],minlength=nbin)
    good=good&(ew>0)&(rw>0)
    total=ew[good].sum()
    row={'n_bins':int(good.sum()),'coverage':float(total/m[~ret].sum()),'eligible_coverage':float(total/m[E].sum()) if E.any() else np.nan}
    incommon=good[bins]&extra
    row.update(rmin_kpc=float(r[incommon].min()) if incommon.any() else np.nan,rmax_kpc=float(r[incommon].max()) if incommon.any() else np.nan)
    for key,x in [('age',age),('Z',z),('logZ',np.log10(np.maximum(z,1e-12))),('radius',r)]:
        ea=np.divide(np.bincount(bins[E],weights=m[E]*x[E],minlength=nbin),ew,out=np.zeros(nbin),where=ew>0)
        ra=np.divide(np.bincount(bins[R],weights=m[R]*x[R],minlength=nbin),rw,out=np.zeros(nbin),where=rw>0)
        v1=np.sum(ea[good]*ew[good])/total if total else np.nan
        v0=np.sum(ra[good]*ew[good])/total if total else np.nan
        row.update({key+'_excluded':float(v1),key+'_retained':float(v0),key+'_delta':float(v1-v0),key+'_fractional':float(v1/v0-1) if v0 else np.nan})
    return row,good

def run(step,nstar):
    if socket.gethostname().lower()!='lageunha':raise SystemExit('Lageunha only')
    dest=ROOT/'results'/f'{step}_{nstar}';dest.mkdir(exist_ok=True)
    if (dest/'DONE.json').exists():return
    start=time.monotonic();sb.read_stars=mb.read_stars
    meta=next(r for r in pr.read('halo_metadata.tsv') if r['step']==step and int(r['nstar'])==nstar)
    a,z,gpath,bpath,tag=rv.catalogue_paths(step)
    with np.load(sb.OUT/f'offsets_{tag}_bg.npz') as d:off=int(d['off'][int(meta['catalogue_row'])])
    pieces=sb.load_halo(gpath,int(meta['byte_offset']));bg=sb.load_halo(bpath,off,background=True)[0]
    central=sb.brightest(pieces);centre=central[1];chunks=[];orig=[]
    for p in pieces:
        if p[2] is None or (p is not central and p[0]>=100):continue
        chunks.append(p[2]);orig.append(np.full(p[0],0 if p is central else 1,dtype=np.uint8))
    if bg[2] is not None:chunks.append(bg[2]);orig.append(np.full(bg[0],2,dtype=np.uint8))
    cat={k:np.concatenate([c[k] for c in chunks]) for k in chunks[0]};origin=np.concatenate(orig)
    del pieces,bg,chunks,orig,central
    with np.load(mb.OLD/'inputs/supplement'/f'fsps_z{z:.3f}.npz') as d:grids=(*sb.build_time(a),d['age'],d['logz'],d['flux'])
    dx,dy,dz,mass,flux=rv.unwrap_flux(cat,centre,a,grids)
    audit=json.loads((mb.PREV/'inputs/supplement/photometry_audit.json').read_text());flux*=next(x['flux_factor'] for x in audit['records'] if abs(x['z']-z)<1e-6)
    birth=np.clip(np.interp(cat['tp'],grids[1],grids[0]),.02,a)
    age=np.interp(a,grids[0],grids[2])-np.interp(birth,grids[0],grids[2]);Z=cat['zp']/sb.ZSUN_SIM
    velocity=np.column_stack([cat[k] for k in ['vx','vy','vz']]); del cat
    old=np.load(mb.OLD/'results'/f'{step}_{nstar}'/'membership.npz');fine=np.load(mb.PREV/'results'/f'{step}_{nstar}'/'fine_membership.npz')
    seed=int(old['seed']);masks={}
    for b in pr.FINE:
        key=f'b{b:.2f}';cache=old if key in old else fine;masks[float(b)]=np.unpackbits(cache[key],count=len(mass)).astype(bool)
    old.close();fine.close()
    base=masks[.2];ret=masks[.1][base]
    assert not (masks[.1]&~base).any()
    scale=a*1000/(sb.H0/100);peak=np.array([dx[seed],dy[seed],dz[seed]])
    pos=(np.column_stack([dx[base],dy[base],dz[base]])-peak)*scale
    m=mass[base]/(sb.H0/100); ages=age[base];met=Z[base];origins=origin[base];v=velocity[base];r=np.linalg.norm(pos,axis=1)
    del age,Z,velocity,origin
    ref=pd.read_csv(PHYS/'results'/f'{step}_{nstar}'/'properties.tsv',sep='\t')
    assert np.isclose(m.sum(),ref[(ref.b==.2)&(ref.aperture_kpc==0)].mass_msun.iloc[0],rtol=1e-10)
    order=np.argsort(r);cw=np.cumsum(m[order]);quant=lambda f:float(np.interp(f*cw[-1],cw,r[order]))
    retained_fraction=m[ret].sum()/m.sum();rbase=quant(.5);rnull=quant(.5*retained_fraction);ractual=mb.weighted_quantile(r[ret],m[ret],[.5])[0]
    null=dict(step=step,nstar=nstar,mass_b020=m.sum(),mass_b010=m[ret].sum(),removed_fraction=1-retained_fraction,rhalf_b020=rbase,rhalf_b010=ractual,rhalf_outer_cut=rnull,actual_change=ractual/rbase-1,null_change=rnull/rbase-1,residual_pp=100*(ractual-rnull)/rbase)
    pd.DataFrame([null]).to_csv(dest/'outer_cut.tsv',sep='\t',index=False)
    # Fine bins define one shared support. Every coarser bin uses exactly the same particles.
    finebins=np.clip(np.floor((np.log10(np.maximum(r,1e-4))+4)/.025).astype(int),0,399)
    allmask=np.ones(len(r),bool); rad=[]
    fineout,support=radial_compare(r,m,ages,met,ret,allmask,finebins)
    common=support[finebins]
    for factor in [1,2,4,8]:
        out,_=radial_compare(r,m,ages,met,ret,common,finebins//factor)
        rad.append(dict(control='fixed_radial_support',width_dex=.025*factor,**out))
    for name,cut in [('all_stars',allmask),('age_ge_0p1Gyr',ages>=.1),('age_ge_1Gyr',ages>=1),('exclude_small_catalogue_galaxies',origins!=1)]:
        out,_=radial_compare(r,m,ages,met,ret,cut,finebins)
        rad.append(dict(control=name,width_dex=.025,**out))
    # Joint age/radius conditioning: each cell must contain >=100 stars per population.
    agebin=np.clip(np.floor(ages/.5).astype(int),0,39)
    out,_=radial_compare(r,m,ages,met,ret,allmask,(finebins//2)*40+agebin)
    rad.append(dict(control='radius_0p05dex_age_0p5Gyr',width_dex=.05,**out))
    # Ellipsoidal shells from fixed b=.20 stellar second moment; volume-preserving axes.
    cov=(pos*m[:,None]).T@pos/m.sum();vals,vec=np.linalg.eigh(cov);axes=np.sqrt(np.maximum(vals,1e-10));axes/=np.prod(axes)**(1/3)
    rell=np.linalg.norm((pos@vec)/axes,axis=1);ellbins=np.clip(np.floor((np.log10(np.maximum(rell,1e-4))+4)/.025).astype(int),0,399)
    out,_=radial_compare(rell,m,ages,met,ret,allmask,ellbins);rad.append(dict(control='ellipsoidal_radius',width_dex=.025,**out))
    pd.DataFrame(rad).assign(step=step,nstar=nstar).to_csv(dest/'population_controls.tsv',sep='\t',index=False)
    comp=[]
    for name,q in [('retained',ret),('excluded',~ret),('excluded_matched',(~ret)&common),('excluded_unmatched',(~ret)&~common)]:
        w=m[q];vv=v[q];mu=np.average(vv,weights=w,axis=0)
        comp.append(dict(component=name,mass_msun=w.sum(),age_Gyr=np.average(ages[q],weights=w),Z_Zsun=np.average(met[q],weights=w),rmean_kpc=np.average(r[q],weights=w),sigma_kms=np.sqrt(np.average(np.sum((vv-mu)**2,axis=1),weights=w)/3),central_mass_fraction=w[origins[q]==0].sum()/w.sum(),small_galaxy_mass_fraction=w[origins[q]==1].sum()/w.sum(),background_mass_fraction=w[origins[q]==2].sum()/w.sum()))
    pd.DataFrame(comp).assign(step=step,nstar=nstar).to_csv(dest/'components.tsv',sep='\t',index=False)
    del pos,m,ages,met,origins,v,r,order,cw,finebins,allmask,common,agebin,rell,ellbins,ret
    print('POPULATIONS_DONE',step,nstar,flush=True)
    centres=pd.read_csv(PHYS/'refinement'/f'{step}_{nstar}'/'centres.tsv',sep='\t');c=centres[centres.centre=='shrinking_from_peak'].iloc[0]
    shift=peak+np.array([c.x_kpc,c.y_kpc,c.z_kpc])/scale
    coords=[dx-shift[0],dy-shift[1],dz-shift[2]]
    ang,n=rv.angular_scale(a,z);raw,pixels=rv.project_raw(*coords,flux,ang,n)
    rows=[];targets=[]
    for ia,los in enumerate(rv.LOS):
        for ib,band in enumerate(rv.BAND):
            kernel=sb.psf_kernel(sb.PIX_AS,sb.FWHM_AS,sb.MOFFAT_BETA,sb.WING_A[ib],sb.KERNEL_AS)[0]
            image=rv.convolve_band(raw[ia][ib],kernel)
            region,inside,wide,flags=pr.image_selection(image,sb.THRESHOLDS[ib],pixels[ia],seed)
            q=pr.fraction(region,kernel,pixels[ia],len(mass))
            j,k=sb.LOS_AXES[ia];R=np.hypot(coords[j],coords[k])*scale
            isoR=mb.weighted_quantile(R[inside],flux[inside,ib],[.5])[0]
            targets.append(dict(step=step,nstar=nstar,los=los,band=band,**flags,req_proper_kpc=flags['req_arcsec']/(ang/scale),Rhalf_isophote_kpc=isoR))
            for b,keep in masks.items():
                score=pr.sums(flux[:,ib],keep,inside);psf=pr.sums(flux[:,ib],keep,q)
                rows.append(dict(step=step,nstar=nstar,z=z,band=band,los=los,b=b,**score,psf_iou=psf['iou']))
            print('IMAGE_DONE',step,nstar,los,band,flush=True)
    pd.DataFrame(rows).to_csv(dest/'recentred_curves.tsv',sep='\t',index=False)
    pd.DataFrame(targets).to_csv(dest/'recentred_targets.tsv',sep='\t',index=False)
    (dest/'DONE.json').write_text(json.dumps(dict(step=step,nstar=nstar,seconds=time.monotonic()-start),indent=2))
    print('CONTROL_DONE',step,nstar,flush=True)
if __name__=='__main__':run(sys.argv[1],int(sys.argv[2]))
