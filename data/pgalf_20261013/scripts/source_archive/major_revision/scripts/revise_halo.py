"""Recompute boundary scores, audit central selection, and rescale distances.

Read-only access to catalogues and archived flux grids. All new results go to
the revision directory. Full baseline curves must reproduce archived values.
"""
import csv
import gc
import json
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent
sys.path.insert(0, str(ROOT/'scripts/source'))
import sb_revise as rv
import sb_vs_b as sb

SUPP = ROOT/'inputs/supplement'


def table(path):
    lines=path.read_text().splitlines()
    body=[s for s in lines if s.strip() and not s.startswith('#')]
    if body[0].split()[0].isdigit():
        header=next(s[2:].split() for s in lines if s.startswith('# step '))
        return [dict(zip(header,s.split(),strict=True)) for s in body]
    return list(csv.DictReader(body,delimiter='\t'))


def write_table(path, rows):
    with path.open('w') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def prepare(step, nstar, dest):
    meta = next(r for r in table(SUPP/'halo_metadata.tsv')
                if r['step']==step and int(r['nstar'])==nstar)
    a, z, gpath, bpath, tag = rv.catalogue_paths(step)
    bgcache = sb.OUT/f"offsets_{tag}_bg.npz"
    if not bgcache.exists():
        raise RuntimeError('Missing precomputed catalogue offsets. No implicit scan permitted.')
    with np.load(bgcache) as offsets:
        bg_offset = int(offsets['off'][int(meta['catalogue_row'])])
    pieces = sb.load_halo(gpath, int(meta['byte_offset']))
    bg = sb.load_halo(bpath, bg_offset, background=True)[0]
    with np.load(SUPP/f'fsps_z{z:.3f}.npz') as grid:
        grids = (*sb.build_time(a), grid['age'], grid['logz'], grid['flux'])
    central = int(np.argmax([p[0] for p in pieces]))
    centre = pieces[central][1]
    count = sum(p[0] for p in pieces)+bg[0]
    if count != nstar:
        raise RuntimeError(f'Parent reconstruction count mismatch {count} != {nstar}')
    gal = []
    for index, p in enumerate(pieces):
        if p[2] is None:
            light = np.zeros(2)
            current_mass = 0.
        else:
            flux = rv.unwrap_flux(p[2], centre, a, grids)[-1]
            light = flux.sum(axis=0)
            current_mass = float(p[2]['mass'].sum())
            del flux
        gal.append(dict(step=step, nstar=nstar, galaxy_index=index, n_stars=p[0],
                        selected=int(index==central), current_mass=current_mass,
                        flux_V=light[0], flux_g=light[1],
                        x=p[1][0], y=p[1][1], z=p[1][2]))
    write_table(dest/'galaxies.tsv', gal)
    cat, centre, nsat, nobj = sb.parent_stars(pieces, bg)
    header_mass=float(meta['mstar_hinv_msun'])
    allmass=sum(r['current_mass'] for r in gal)
    if bg[2] is not None:
        allmass+=float(bg[2]['mass'].sum())
    audit=dict(step=step, nstar=nstar, reconstructed_count=count,
               parent_count=len(cat['mass']), satellite_count=nsat,
               selected_index=central,
               brightest_V=int(np.argmax([r['flux_V'] for r in gal])),
               brightest_g=int(np.argmax([r['flux_g'] for r in gal])),
               most_massive=int(np.argmax([r['current_mass'] for r in gal])),
               catalogue_stellar_mass=header_mass, summed_stellar_mass=allmass)
    (dest/'catalogue_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    del pieces,bg
    dx,dy,dz,mass,flux=rv.unwrap_flux(cat,centre,a,grids)
    del cat
    return dict(a=a,z=z,dx=dx,dy=dy,dz=dz,mass=mass,flux=flux)


def run(step, nstar):
    dest=ROOT/'results'/f'{step}_{nstar}'
    dest.mkdir(parents=True, exist_ok=True)
    if (dest/'DONE').exists():
        print('ALREADY_DONE',step,nstar,flush=True)
        return
    audit=json.loads((ROOT/'results/photometry_audit.json').read_text())
    if not audit['grid_reproduction_passed']:
        raise RuntimeError('Photometric normalization has not been recovered.')
    print('START',step,nstar,flush=True)
    h=prepare(step,nstar,dest)
    factor=next(r['flux_factor'] for r in audit['records']
                if abs(r['z']-h['z'])<1e-6)
    scale,n=rv.angular_scale(h['a'],h['z'])
    coords=(h['dx'],h['dy'],h['dz'])
    raw,pixs=rv.project_raw(*coords,h['flux'],scale,n)
    kernels=[sb.psf_kernel(sb.PIX_AS,sb.FWHM_AS,sb.MOFFAT_BETA,v,sb.KERNEL_AS)[0]
             for v in sb.WING_A]
    regions={}
    areas={}
    masks={}
    for ia in range(3):
        for ib in range(2):
            img=rv.convolve_band(raw[ia][ib],kernels[ib])
            for variant,mult in [('baseline',1.),('consistent_distance',factor)]:
                region,area=sb.central_isophote(img*mult,sb.THRESHOLDS[ib],sb.PIX_AS)
                key=(variant,ia,ib)
                regions[key]=region
                areas[key]=area
                masks[key]=(rv.sample_inside(region,*pixs[ia]),
                            rv.sample_inside(rv.ring_of(region),*pixs[ia]))
            del img
    del raw,kernels
    gc.collect()
    published={(r['los'],r['band'],round(float(r['b']),2)):r for r in
               table(SUPP/'sightline_baseline.tsv')
               if r['step']==step and int(r['nstar'])==nstar}
    cache=dest/'membership.npz'
    old=np.load(cache) if cache.exists() else None
    seed=int(old['seed']) if old is not None else sb.nucleus_index(*coords)
    packed={}
    rows=[]
    profiles=[]
    for b in rv.B_GRID:
        tag=f'b{b:.2f}'
        if old is None:
            root=sb.fof_largest(*coords,sb.link02(h['mass'],b))
            keep=sb.peak_mask(root,*coords,seed)
            del root
            packed[tag]=np.packbits(keep)
        else:
            keep=np.unpackbits(old[tag],count=len(h['mass'])).astype(bool)
        print('GROUP',step,nstar,float(b),int(keep.sum()),flush=True)
        for (variant,ia,ib),(inside,ring) in masks.items():
            w=h['flux'][:,ib] if variant=='baseline' else h['flux'][:,ib]*factor
            iou,completeness,precision=rv.light_scores(w,keep,inside)
            riou,rc,rp=rv.boundary_scores(w,keep,inside,ring)
            if variant=='baseline':
                ref=published[(rv.LOS[ia],rv.BAND[ib],round(float(b),2))]
                for key,got in [('iou',iou),('recall',completeness),('purity',precision)]:
                    if abs(got-float(ref[key]))>5.1e-5:
                        raise RuntimeError(f'Baseline mismatch {step} {nstar} {key} {got} {ref[key]}')
                if int(keep.sum())!=int(ref['n_group']):
                    raise RuntimeError('Baseline membership count differs.')
            rows.append(dict(step=step,nstar=nstar,z=h['z'],variant=variant,
                             los=rv.LOS[ia],band=rv.BAND[ib],b=float(b),
                             area_arcsec2=areas[(variant,ia,ib)],iou=iou,
                             completeness=completeness,precision=precision,
                             ring_iou=riou,ring_completeness=rc,ring_precision=rp,
                             L_group=float(w[keep].sum()),L_inside=float(w[inside].sum()),
                             L_both=float(w[keep&inside].sum()),
                             L_group_ring=float(w[keep&ring].sum()),
                             L_inside_ring=float(w[inside&ring].sum()),
                             L_both_ring=float(w[keep&inside&ring].sum())))
            if variant=='baseline' and (np.isclose(b,.1) or np.isclose(b,.2)):
                j,k=sb.LOS_AXES[ia]
                radius=np.hypot(coords[j],coords[k])*scale
                req=np.sqrt(areas[(variant,ia,ib)]/np.pi)
                annulus=np.searchsorted(np.array([.5,1.,1.5])*req,radius,side='right')
                for band_index in range(4):
                    selected=annulus==band_index
                    profiles.append(dict(step=step,nstar=nstar,los=rv.LOS[ia],band=rv.BAND[ib],
                                         b=float(b),annulus=band_index,req_arcsec=req,
                                         L_group=float(w[selected&keep].sum()),
                                         L_inside=float(w[selected&inside].sum()),
                                         L_both=float(w[selected&keep&inside].sum())))
                del radius,annulus
    if old is None:
        np.savez(cache,seed=seed,n_stars=len(h['mass']),**packed)
    else:
        old.close()
    write_table(dest/'curves.tsv',rows)
    write_table(dest/'radial.tsv',profiles)
    (dest/'DONE').write_text('All baseline scores and membership counts reproduced.\n')
    print('DONE',step,nstar,flush=True)


if __name__=='__main__':
    import fcntl
    step,nstar=sys.argv[1],int(sys.argv[2])
    dest=ROOT/'results'/f'{step}_{nstar}'
    dest.mkdir(parents=True,exist_ok=True)
    # Two finite queues may reach the same halo. Only one calculation may
    # hold the lock. A later caller then observes DONE and reuses the result.
    with (dest/'.calculation.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        run(step,nstar)
