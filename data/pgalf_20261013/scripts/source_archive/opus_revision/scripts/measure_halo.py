"""Consistent-distance photometry and inverse-probability particle sampling.

Original catalogues and archived membership caches are read-only. All new
outputs are confined to this task. This is not a hydrodynamic convergence run.
"""
import os
os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                  MKL_NUM_THREADS='1', NUMBA_NUM_THREADS='1')
import csv
import gc
import json
import socket
import sys
from pathlib import Path
import numpy as np
from scipy import ndimage
from scipy.signal import fftconvolve

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT.parent / 'PGalF_major_revision_20261013'
SUPP = ROOT / 'inputs/supplement'
sys.path.insert(0, str(ROOT / 'scripts/source'))
import sb_vs_b as sb
import sb_revise as rv
FINE = np.array([.02,.03,.05,.06,.07,.08,.09,.10,.11,.12,.13,.14,.15,.18,.20])
ORIGINAL = np.array([.02,.03,.05,.08,.10,.12,.15,.18,.20])
EXTENDED = np.array([.25,.30,.40,.50])


def read(name):
    return list(csv.DictReader([s for s in (SUPP/name).read_text().splitlines()
                               if s.strip() and not s.startswith('#')], delimiter='\t'))


def write(path, rows):
    tmp = path.with_suffix(path.suffix+'.tmp')
    with tmp.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)
    tmp.replace(path)


def load(step, nstar):
    meta = next(r for r in read('halo_metadata.tsv')
                if r['step']==step and int(r['nstar'])==nstar)
    a,z,gpath,bpath,tag = rv.catalogue_paths(step)
    offset_cache = sb.OUT / f'offsets_{tag}_bg.npz'
    if not offset_cache.is_file():
        raise RuntimeError('Missing offsets; implicit catalogue scan is forbidden.')
    with np.load(offset_cache) as d:
        offset = int(d['off'][int(meta['catalogue_row'])])
    pieces = sb.load_halo(gpath, int(meta['byte_offset']))
    bg = sb.load_halo(bpath, offset, background=True)[0]
    if sum(p[0] for p in pieces)+bg[0] != nstar:
        raise RuntimeError('Catalogue count mismatch')
    cat,centre,nsat,nobj = sb.parent_stars(pieces,bg)
    del pieces,bg
    with np.load(OLD/'inputs/supplement'/f'fsps_z{z:.3f}.npz') as d:
        grids = (*sb.build_time(a),d['age'],d['logz'],d['flux'])
    dx,dy,dz,mass,flux = rv.unwrap_flux(cat,centre,a,grids)
    birth=np.clip(np.interp(cat['tp'],grids[1],grids[0]),.02,a)
    ages=(np.interp(a,grids[0],grids[2])-np.interp(birth,grids[0],grids[2])).astype(np.float32)
    del birth
    del cat
    audit=json.loads((SUPP/'photometry_audit.json').read_text())
    factor=next(r['flux_factor'] for r in audit['records'] if abs(r['z']-z)<1e-6)
    flux *= factor
    return dict(dx=dx,dy=dy,dz=dz,mass=mass,flux=flux,ages=ages,a=a,z=z,
                centre=centre,nsat=nsat,nobj=nobj,factor=factor)


def fraction(region, kernel, pixels, nstars):
    ix,iy,on = pixels
    out=np.zeros(nstars, np.float32)
    if region.any():
        image=fftconvolve(region.astype(float),kernel,mode='same')
        out[on]=image[iy[on],ix[on]]
        if out.max()>1+1e-5:
            raise RuntimeError('Invalid PSF fraction')
        np.clip(out,0,1,out=out)
    return out


def sums(w, keep, selection):
    lg=float(w[keep].sum())
    if selection.dtype==np.bool_:
        li=float(w[selection].sum()); lb=float(w[keep&selection].sum())
    else:
        li=float(np.sum(w*selection)); lb=float(np.sum(w[keep]*selection[keep]))
    den=lg+li-lb
    return dict(iou=lb/den if den>0 else np.nan,
                completeness=lb/li if li>0 else np.nan,
                precision=lb/lg if lg>0 else np.nan,
                L_group=lg,L_target=li,L_shared=lb,
                luminosity_error=lg/li-1 if li>0 else np.nan)


def stellar_properties(h,keep,inside,q,ib):
    """Current mass and luminosity-weighted true stellar age, not clipped SSP age."""
    w=h['flux'][:,ib];m=h['mass'];ages=h['ages']
    def mean_age(weight):
        den=float(weight.sum())
        return float(np.sum(weight*ages))/den if den>0 else np.nan
    group_age=mean_age(w*keep)
    position_age=mean_age(w*inside)
    aperture_age=mean_age(w*q)
    target_mass=float(m[inside].sum());group_mass=float(m[keep].sum())
    return dict(group_mass_hinv_msun=group_mass,target_mass_hinv_msun=target_mass,
                mass_error=group_mass/target_mass-1 if target_mass>0 else np.nan,
                age_group_Gyr=group_age,age_position_Gyr=position_age,age_aperture_Gyr=aperture_age,
                age_difference_position_Gyr=group_age-position_age,
                age_difference_aperture_Gyr=group_age-aperture_age)


def image_selection(image, threshold, pixels, seed, peak=False):
    mu,mask,labels,nlab=rv.mu_labels(image,threshold)
    if peak:
        py,px=np.unravel_index(int(np.argmax(image)),image.shape)
    else:
        py=px=image.shape[0]//2
    cid,fallback,off=rv.component_at(labels,mask,py,px)
    region,area=rv.region_of(labels,cid,sb.PIX_AS)
    ix,iy,on=pixels
    seed_inside=bool(on[seed] and region[iy[seed],ix[seed]])
    if region.any():
        distances=ndimage.distance_transform_edt(region) - ndimage.distance_transform_edt(~region)
        wide={str(width):rv.sample_inside(np.abs(distances)*sb.PIX_AS<=width,ix,iy,on)
              for width in [1.,2.]}
    else:
        wide={str(width):np.zeros(len(ix),bool) for width in [1.,2.]}
    inside=rv.sample_inside(region,ix,iy,on)
    offset=float(np.hypot(ix[seed]-image.shape[1]//2,iy[seed]-image.shape[0]//2)*sb.PIX_AS)
    flags=dict(area_arcsec2=area,req_arcsec=np.sqrt(area/np.pi),fallback=int(fallback),
               empty=int(area==0),one_pixel=int(np.isclose(area,.04)),
               seed_inside=int(seed_inside),seed_offset_arcsec=offset,n_components=nlab)
    return region,inside,wide,flags


def run(step,nstar):
    dest=ROOT/'results'/f'{step}_{nstar}';dest.mkdir(exist_ok=True)
    if (dest/'DONE').is_file():
        print('ALREADY_DONE',step,nstar,flush=True);return
    print('LOAD',step,nstar,flush=True)
    h=load(step,nstar); coords=[h[k] for k in ['dx','dy','dz']]
    ninput=len(h['mass']);scale,n=rv.angular_scale(h['a'],h['z'])
    cache=OLD/'results'/f'{step}_{nstar}'/'membership.npz'
    if not cache.is_file():raise RuntimeError('Original membership cache is required')
    original=np.load(cache);seed=int(original['seed'])
    kernels=[sb.psf_kernel(sb.PIX_AS,sb.FWHM_AS,sb.MOFFAT_BETA,a,sb.KERNEL_AS)[0] for a in sb.WING_A]
    wingoff=sb.psf_kernel(sb.PIX_AS,sb.FWHM_AS,sb.MOFFAT_BETA,0.,sb.KERNEL_AS)[0]
    raw,pixels=rv.project_raw(*coords,h['flux'],scale,n)
    geometry=dict(step=step,nstar=nstar,z=h['z'],scale_arcsec_per_hinv_cMpc=scale,n_pixels=n,
                  seed_pixels={los:[int(pixels[ia][0][seed]),int(pixels[ia][1][seed])] for ia,los in enumerate(rv.LOS)},
                  seed_offsets_arcsec={los:[float(coords[j][seed]*scale),float(coords[k][seed]*scale)]
                                       for los,(j,k) in zip(rv.LOS,sb.LOS_AXES)})
    (dest/'drawing_geometry.json').write_text(json.dumps(geometry,indent=2)+'\n')
    masks={};targetrows=[]; illustration={}
    ref={(r['los'],r['band'],float(r['b'])):r for r in read('recomputed_curves.tsv')
         if r['variant']=='consistent_distance' and r['step']==step and int(r['nstar'])==nstar}
    for ia,los in enumerate(rv.LOS):
        for ib,band in enumerate(rv.BAND):
            image=rv.convolve_band(raw[ia][ib],kernels[ib])
            for variant,threshold,img,peak in [
                ('primary',sb.THRESHOLDS[ib],image,False),
                ('threshold_bright',sb.THRESHOLDS[ib]-.5,image,False),
                ('threshold_faint',sb.THRESHOLDS[ib]+.5,image,False),
                ('light_peak',sb.THRESHOLDS[ib],image,True),
                ('wing_off',sb.THRESHOLDS[ib],rv.convolve_band(raw[ia][ib],wingoff),False)]:
                region,inside,wide,flags=image_selection(img,threshold,pixels[ia],seed,peak)
                item=dict(inside=inside,flags=flags,ib=ib)
                if variant=='primary':
                    expected=float(ref[(los,band,.1)]['area_arcsec2'])
                    if abs(flags['area_arcsec2']-expected)>5e-4:raise RuntimeError('Primary area not reproduced')
                    q=fraction(region,kernels[ib],pixels[ia],ninput)
                    photo=float(image[region].sum());pred=float(np.sum(h['flux'][:,ib]*q))
                    if abs(photo-pred)>1e-5*max(photo,pred,1e-30):raise RuntimeError('PSF aperture flux mismatch')
                    item.update(q=q,wide=wide)
                    flags['aperture_flux']=photo
                    illustration[f'image_{los}_{band}']=image.astype(np.float32)
                    illustration[f'mask_{los}_{band}']=np.packbits(region)
                masks[(variant,los,band)]=item
                targetrows.append(dict(step=step,nstar=nstar,z=h['z'],variant=variant,los=los,band=band,**flags))
            print('IMAGE',los,band,flush=True)
    del raw;gc.collect()
    write(dest/'targets.tsv',targetrows)
    # Drawing arrays are bounded and contain no full particle catalogue.
    sample=np.linspace(0,ninput-1,min(ninput,12000),dtype=int)
    for ia,los in enumerate(rv.LOS):
        illustration[f'ix_{los}']=pixels[ia][0][sample];illustration[f'iy_{los}']=pixels[ia][1][sample]
    illustration['sample_indices']=sample
    illustration['seed']=np.array(seed)
    np.savez_compressed(dest/'images.npz',**illustration);del illustration
    massq=np.quantile(h['mass'],[.16,.5,.84]);linkq=sb.link02(massq,.1)*1e3
    (dest/'mass_scale.json').write_text(json.dumps(dict(step=step,nstar=nstar,z=h['z'],n_input=ninput,
        mass_hinv_msun=massq.tolist(),link_b010_hinv_ckpc=linkq.tolist(),
        link_b010_proper_kpc=(linkq*h['a']/(sb.H0/100)).tolist(),
        link_b010_cells=(linkq/(44.171267*1e3/2**17)).tolist()),indent=2)+'\n')
    rows=[];properties=[];newpacked={};dense_samples={}
    extra_cache=dest/'fine_membership.npz'
    cached=np.load(extra_cache) if extra_cache.is_file() else None
    for b in FINE:
        tag=f'b{b:.2f}'
        if tag in original:
            keep=np.unpackbits(original[tag],count=ninput).astype(bool)
        elif cached is not None and tag in cached:
            keep=np.unpackbits(cached[tag],count=ninput).astype(bool)
        else:
            roots=sb.fof_largest(*coords,sb.link02(h['mass'],b));keep=roots==roots[seed];del roots
            newpacked[tag]=np.packbits(keep)
        if b in [.08,.09,.1,.2]:dense_samples[tag]=keep[sample]
        print('FINE',float(b),int(keep.sum()),flush=True)
        for (variant,los,band),item in masks.items():
            ib=item['ib'];w=h['flux'][:,ib];inside=item['inside']
            score=sums(w,keep,inside)
            if variant=='primary' and b in ORIGINAL:
                expected=ref[(los,band,float(b))]
                if abs(score['iou']-float(expected['iou']))>1e-8:raise RuntimeError('Original primary score not reproduced')
            row=dict(step=step,nstar=nstar,z=h['z'],variant=variant,los=los,band=band,b=float(b),**score)
            if variant=='primary':
                row.update({f'psf_{k}':v for k,v in sums(w,keep,item['q']).items()})
                properties.append(dict(step=step,nstar=nstar,z=h['z'],los=los,band=band,b=float(b),
                                       **stellar_properties(h,keep,inside,item['q'],ib)))
                for width,ring in item['wide'].items():
                    row.update({f'band{width}_{k}':v for k,v in sums(w[ring],keep[ring],inside[ring]).items()})
            else:
                row.update({f'psf_{k}':np.nan for k in score})
                for width in ['1.0','2.0']:row.update({f'band{width}_{k}':np.nan for k in score})
            rows.append(row)
        del keep
    if cached is not None:cached.close()
    if newpacked:np.savez(extra_cache,**newpacked)
    np.savez(dest/'drawing_membership.npz',**dense_samples)
    write(dest/'fine_curves.tsv',rows)
    write(dest/'stellar_properties.tsv',properties)
    original.close();del rows;gc.collect()
    rng=np.random.default_rng(np.random.SeedSequence([20261013,int(step),nstar]))
    draw=rng.random(ninput);thinrows=[];normalization=[]
    for factor in [2,4,8]:
        selected=draw<1/factor
        tc=[np.ascontiguousarray(a[selected]) for a in coords]
        mass=h['mass'][selected]*factor;flux=h['flux'][selected]*factor
        tseed=sb.nucleus_index(*tc)
        normalization.append(dict(step=step,nstar=nstar,factor=factor,n_original=ninput,n_retained=int(selected.sum()),
             mass_ratio=float(mass.sum()/h['mass'].sum()),
             flux_ratio_V=float(flux[:,0].sum()/h['flux'][:,0].sum()),flux_ratio_g=float(flux[:,1].sum()/h['flux'][:,1].sum())))
        traw,tpixels=rv.project_raw(*tc,flux,scale,n);tmasks={}
        for ia,los in enumerate(rv.LOS):
            for ib,band in enumerate(rv.BAND):
                base=masks[('primary',los,band)]
                tmasks[('fixed',los,band)]=(base['inside'][selected],base['q'][selected],ib)
                image=rv.convolve_band(traw[ia][ib],kernels[ib])
                region,inside,wide,flags=image_selection(image,sb.THRESHOLDS[ib],tpixels[ia],tseed)
                q=fraction(region,kernels[ib],tpixels[ia],len(mass))
                tmasks[('rebuilt',los,band)]=(inside,q,ib)
                targetrows.append(dict(step=step,nstar=nstar,z=h['z'],variant=f'thin{factor}_rebuilt',los=los,band=band,**flags))
        del traw
        factorrows=[]
        def measure(b):
            roots=sb.fof_largest(*tc,sb.link02(mass,b));keep=roots==roots[tseed];del roots
            print('THIN',factor,float(b),int(keep.sum()),flush=True)
            for (imagekind,los,band),(inside,q,ib) in tmasks.items():
                scores=sums(flux[:,ib],keep,inside)
                row=dict(step=step,nstar=nstar,z=h['z'],factor=factor,image=imagekind,los=los,band=band,b=float(b),**scores)
                row.update({f'psf_{k}':v for k,v in sums(flux[:,ib],keep,q).items()})
                factorrows.append(row)
        for b in ORIGINAL:measure(b)
        # Extend the upper domain only when at least one curve has its maximum there.
        def needs_upper():
            return any(max([r for r in factorrows if r['image']==kind and r['los']==los and r['band']==band],key=lambda r:r[score])['b']==max(r['b'] for r in factorrows)
                       for kind in ['fixed','rebuilt'] for los in rv.LOS for band in rv.BAND for score in ['iou','psf_iou'])
        if needs_upper():
            for b in EXTENDED:measure(b)
        thinrows.extend(factorrows)
        write(dest/f'thin{factor}.tsv',factorrows)
        del tc,mass,flux,tmasks;gc.collect()
    write(dest/'thin_curves.tsv',thinrows)
    write(dest/'sampling_totals.tsv',normalization)
    # Supply a uniform schema despite aperture_flux being specific to primary images.
    allfields=list(dict.fromkeys(k for r in targetrows for k in r))
    targetrows=[{k:r.get(k,np.nan) for k in allfields} for r in targetrows]
    write(dest/'targets.tsv',targetrows)
    (dest/'DONE').write_text('Consistent-distance primary scores and PSF aperture fluxes reproduced.\n')
    print('COMPLETE',step,nstar,flush=True)


if __name__=='__main__':
    import fcntl
    if socket.gethostname().lower()!='lageunha':
        raise SystemExit('Particle calculations are authorized on Lageunha only.')
    step,nstar=sys.argv[1],int(sys.argv[2])
    dest=ROOT/'results'/f'{step}_{nstar}';dest.mkdir(exist_ok=True)
    with (dest/'.calculation.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        run(step,nstar)
