#!/usr/bin/env python3
"""Fixed-b image sensitivity on the locked 49-host CAMELS sample.

The scenarios are one-at-a-time stress tests, not a retuning of the Darwin
linking parameter or a claim to reproduce a particular survey's detector.
"""
from __future__ import annotations
import argparse
import csv
import gc
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy import ndimage
from scipy.signal import fftconvolve
import measure_camels_halo as base

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"results/image_sensitivity_20261030"
CONFIG={"sample":"Locked sample_manifest.csv, 49 groups in CV20-CV26",
        "b_values":[.10,.20],"aperture_projected_pkpc":30.,
        "image_centre":"central subhalo minimum-potential position",
        "baseline_pixel_arcsec":.2,"stamp_half_width_arcsec":180.,
        "thresholds_mag_arcsec2":{"V":26.5,"g":27.0},
        "threshold_offsets_mag":[-.5,.5],"coarse_pixel_arcsec":.4,
        "no_wing":"Set the faint PSF wing amplitude to zero and renormalize the PSF",
        "sky_residual_mag_arcsec2":29.,
        "noise_threshold_to_pixel_rms":[5.,3.],"noise_realisations":3,
        "noise_seed":20261030,
        "noise_definition":"Independent zero-mean Gaussian noise added after PSF convolution, with per-pixel rms equal to threshold pixel flux divided by 5 or 3; not a survey-calibrated noise model",
        "inference":"Views and noise realisations are repeated measurements within host; no linking-parameter fitting on the validation sample"}
CONFIG_SHA=hashlib.sha256(json.dumps(CONFIG,sort_keys=True).encode()).hexdigest()


def scenarios():
    return [("baseline",.2,0.,False,None,None,0),
            ("limit_brighter",.2,-.5,False,None,None,0),
            ("limit_fainter",.2,.5,False,None,None,0),
            ("pixel_coarse",.4,0.,False,None,None,0),
            ("no_psf_wing",.2,0.,True,None,None,0),
            ("sky_residual_29",.2,0.,False,29.,None,0)] + [
            (f"noise_{int(snr)}sigma",.2,0.,False,None,snr,seed)
            for snr in (5.,3.) for seed in range(3)]


def identify(image,threshold,pixel,ix,iy,on):
    per=image/(pixel*pixel)
    positive=per>0
    mu=np.full(image.shape,np.inf)
    mu[positive]=-2.5*np.log10(per[positive])
    mask=mu<=threshold
    labels,nlabels=ndimage.label(mask,structure=base.rv.STRUCT)
    centre=image.shape[0]//2
    component,fallback,clipped=base.rv.component_at(labels,mask,centre,centre)
    region,area=base.rv.region_of(labels,component,pixel)
    target=base.rv.sample_inside(region,ix,iy,on)
    return target,bool(fallback or clipped or not region.any()),area,int(nlabels)


def baseline_audit(rows,box,group):
    reference=ROOT/"results"/f"CV{box:02d}_group{group:05d}_blind20261030_v4_potential_centered/metrics.csv"
    old={}
    for r in csv.DictReader(reference.open()):
        if r["target_variant"]!="potential":continue
        method=r["method"]
        b="aperture" if method=="aperture30pkpc" else float(r["b_star"])
        if b not in ("aperture",.1,.2):continue
        old[(r["line_of_sight"],r["band"],method,b)]=r
    checked=0
    maximum=0.
    for r in rows:
        if r["condition"]!="baseline":continue
        b="aperture" if r["method"]=="aperture30pkpc" else float(r["b_star"])
        ref=old[(r["line_of_sight"],r["band"],r["method"],b)]
        for metric in ("iou","completeness","precision","L_selected","L_contour","L_common","R50_selected_kpc","R50_contour_kpc"):
            actual,expected=float(r[metric]),float(ref[metric])
            if not np.isclose(actual,expected,rtol=1e-9,atol=1e-25,equal_nan=True):
                raise ValueError(f"Baseline failed CV{box} group{group} {metric}: {actual} != {expected}")
            if np.isfinite(actual):maximum=max(maximum,abs(actual-expected))
        checked+=1
    if checked!=18:raise ValueError(f"Expected 18 baseline rows, found {checked}")
    return {"checked_rows":checked,"maximum_absolute_metric_difference":maximum}


def measure(task):
    box,group=task
    directory=OUT/f"CV{box:02d}_group{group:05d}"
    complete=directory/"metadata.json"
    if complete.exists():
        meta=json.loads(complete.read_text())
        if meta["config_sha256"]!=CONFIG_SHA or not meta["baseline_reproduced"]:
            raise ValueError(f"Invalid existing checkpoint {directory}")
        return str(directory)
    directory.mkdir(parents=True,exist_ok=True)
    manifest=base.read_manifest_row(box,group)
    coords,coords_h,mass,formed,birth,zmetal,input_meta=base.read_input_stars(box,group,int(manifest["central_subhalo_index"]),"catalogue_ids")
    cosmo=FlatLambdaCDM(H0=100*base.HUBBLE,Om0=base.OMEGA_M)
    old_cosmo=FlatLambdaCDM(H0=72.,Om0=.27)
    with np.load(base.FSPS_GRID) as grid:
        if not np.isclose(float(grid["z"]),base.OBSERVED_REDSHIFT,rtol=0,atol=1e-12):
            raise ValueError("Unexpected FSPS grid")
        flux,_,_=base.star_fluxes(formed,birth,zmetal,(grid["age"],grid["logz"],grid["flux"]),cosmo,old_cosmo)
    arcsec_per_mpc=206265./cosmo.angular_diameter_distance(base.OBSERVED_REDSHIFT).value
    seed=base.sb.nucleus_index(*[coords_h[:,a] for a in range(3)])
    selections={}
    for b in (.1,.2):
        link=b*(mass/(2.7755e11*base.OMEGA_M))**(1./3.)
        roots=base.sb.fof_largest(coords_h[:,0],coords_h[:,1],coords_h[:,2],link)
        selections[b]=base.sb.peak_mask(roots,*[coords_h[:,a] for a in range(3)],seed)
    rows=[]
    for view_index,(los,axis_x,axis_y) in enumerate(base.LOS):
      radius=np.hypot(coords[:,axis_x],coords[:,axis_y])*1000
      aperture=radius<=30
      for band,band_index,limit in base.BANDS:
        image_cache={}
        for condition,pixel,offset,no_wing,sky,snr,realisation in scenarios():
            geometry=(pixel,no_wing)
            if geometry not in image_cache:
                n=int(round(2*base.STAMP_AS/pixel))
                scale=arcsec_per_mpc/pixel
                ix=np.rint(coords[:,axis_x]*scale).astype(np.int32)+n//2
                iy=np.rint(coords[:,axis_y]*scale).astype(np.int32)+n//2
                on=(ix>=0)&(ix<n)&(iy>=0)&(iy<n)
                kernel=base.sb.psf_kernel(pixel,base.sb.FWHM_AS,base.sb.MOFFAT_BETA,0. if no_wing else base.sb.WING_A[band_index],base.sb.KERNEL_AS)[0]
                image=fftconvolve(base.sb.bin_flux(ix,iy,flux[:,band_index],n),kernel,mode="same")
                image[image<0]=0.
                image_cache[geometry]=(image,ix,iy,on)
            image,ix,iy,on=image_cache[geometry]
            threshold=limit+offset
            rendered=image
            if sky is not None:rendered=image+10**(-.4*sky)*pixel**2
            if snr is not None:
                rng=np.random.default_rng(CONFIG["noise_seed"]+box*1000003+group*9176+view_index*137+band_index*17+realisation)
                noise=rng.standard_normal(image.shape)
                rendered=image+noise*(10**(-.4*threshold)*pixel**2/snr)
            target,failed,area,nlabels=identify(rendered,threshold,pixel,ix,iy,on)
            for method,b,selected in [("FoF",.1,selections[.1]),("FoF",.2,selections[.2]),("aperture30pkpc","aperture",aperture)]:
                metrics=base.score(selected,target,flux[:,band_index],radius)
                if failed:
                    metrics.update(iou=0.,completeness=0.)
                rows.append({"cv_box":box,"fof_group_index":group,"condition":condition,"noise_realisation":realisation,
                             "line_of_sight":los,"band":band,"method":method,"b_star":b,
                             "threshold_mag_arcsec2":threshold,"pixel_arcsec":pixel,
                             "noise_threshold_to_pixel_rms":snr if snr is not None else "",
                             "target_failed":int(failed),"target_area_arcsec2":area,"n_components":nlabels,**metrics})
        del image_cache
    audit=baseline_audit(rows,box,group)
    with (directory/"metrics.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    metadata={"config_sha256":CONFIG_SHA,"baseline_reproduced":True,"baseline_audit":audit,
              "n_rows":len(rows),"n_stars":len(mass),"input_metadata":input_meta}
    complete.write_text(json.dumps(metadata,indent=2)+"\n")
    gc.collect()
    return str(directory)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot",action="store_true")
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    frozen=OUT/"configuration.json"
    if frozen.exists() and json.loads(frozen.read_text())["config_sha256"]!=CONFIG_SHA:
        raise ValueError("Frozen sensitivity configuration differs")
    frozen.write_text(json.dumps({"config_sha256":CONFIG_SHA,**CONFIG},indent=2)+"\n")
    manifest=list(csv.DictReader(base.MANIFEST.open()))
    if len(manifest)!=49:raise ValueError("Locked sample has changed")
    tasks=[(int(r["cv_box"]),int(r["fof_group_index"])) for r in manifest]
    if args.pilot:
        print(measure(tasks[0]),flush=True)
        return
    pilot=OUT/f"CV{tasks[0][0]:02d}_group{tasks[0][1]:05d}/metadata.json"
    if not pilot.exists():raise ValueError("Pilot baseline gate must pass before full sample")
    for box in sorted(set(b for b,_ in tasks)):
        with ProcessPoolExecutor(max_workers=2) as pool:
            futures=[pool.submit(measure,t) for t in tasks if t[0]==box]
            for done in as_completed(futures):print("completed",done.result(),flush=True)
        base.clear_snapshot_id_cache()
    print("Completed locked 49-host image-sensitivity sample",flush=True)


if __name__=="__main__":main()
