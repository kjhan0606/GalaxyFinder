#!/usr/bin/env python3
"""Host-unit summaries of fixed-b image stress tests, without recalibration."""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import measure_image_sensitivity as measurement

ROOT=Path(__file__).resolve().parent
OUT=measurement.OUT
NBOOT=5000
METRICS=("iou","completeness","contamination","logL","abs_logL","logR50","abs_logR50","failure_fraction")
HOST_KEYS=["cv_box","fof_group_index"]


def resampling_weights(hosts):
    rng=np.random.default_rng(20261030)
    n=len(hosts)
    plain=rng.multinomial(n,np.full(n,1/n),size=NBOOT).astype(float)/n
    boxes=sorted(hosts.cv_box.unique())
    by_box={b:np.flatnonzero(hosts.cv_box.to_numpy()==b) for b in boxes}
    nested=np.zeros((NBOOT,n))
    for k in range(NBOOT):
        for b in rng.choice(boxes,len(boxes),replace=True):
            indices=by_box[b]
            np.add.at(nested[k],rng.choice(indices,len(indices),replace=True),1)
    nested/=nested.sum(axis=1,keepdims=True)
    return plain,nested


def paired(left,right,weights,metadata):
    merged=left.merge(right,on=HOST_KEYS,validate="one_to_one",suffixes=("_left","_right")).sort_values(HOST_KEYS)
    if len(merged)!=49:raise ValueError("Paired contrast lost hosts")
    rows=[]
    for metric in METRICS:
        values=(merged[metric+"_left"]-merged[metric+"_right"]).to_numpy(float)
        good=np.isfinite(values)
        row={**metadata,"metric":metric,"n_hosts":49,"n_hosts_with_valid_metric":int(good.sum()),
             "mean_paired_difference":float(np.mean(values[good])) if good.any() else np.nan}
        for label,w in zip(("host","box_then_host"),weights):
            denominator=w@good.astype(float)
            means=np.divide(w@np.nan_to_num(values,nan=0),denominator,
                            out=np.full(NBOOT,np.nan),where=denominator>0)
            finite=np.isfinite(means)
            row[label+"_bootstrap95_low"]=float(np.quantile(means[finite],.025)) if finite.any() else np.nan
            row[label+"_bootstrap95_high"]=float(np.quantile(means[finite],.975)) if finite.any() else np.nan
        rows.append(row)
    return rows


def main():
    if (OUT/"summary_metadata.json").exists():raise FileExistsError("Completed sensitivity summary already exists")
    manifest=list(csv.DictReader(measurement.base.MANIFEST.open()))
    if len(manifest)!=49:raise ValueError("Locked manifest changed")
    frames=[]
    for r in manifest:
        directory=OUT/f"CV{int(r['cv_box']):02d}_group{int(r['fof_group_index']):05d}"
        meta=json.loads((directory/"metadata.json").read_text())
        if not meta["baseline_reproduced"] or meta["config_sha256"]!=measurement.CONFIG_SHA:
            raise ValueError(f"Unvalidated input {directory}")
        frame=pd.read_csv(directory/"metrics.csv")
        if len(frame)!=18*len(measurement.scenarios()):raise ValueError(f"Incomplete stress-test rows {directory}")
        frames.append(frame)
    raw=pd.concat(frames,ignore_index=True)
    raw.b_star=raw.b_star.astype(str)
    if raw.duplicated(HOST_KEYS+["condition","noise_realisation","line_of_sight","band","method","b_star"]).any():
        raise ValueError("Duplicate projection metrics")
    failed=raw.target_failed.astype(bool)
    raw["failure_fraction"]=failed.astype(float)
    raw["contamination"]=(1-raw.precision).where(~failed)
    for source,target in [("L","logL"),("R50","logR50")]:
        if source=="L":num,den=raw.L_selected,raw.L_contour
        else:num,den=raw.R50_selected_kpc,raw.R50_contour_kpc
        valid=(~failed)&num.gt(0)&den.gt(0)&np.isfinite(num)&np.isfinite(den)
        raw[target]=np.log10((num/den).where(valid))
        raw["abs_"+target]=raw[target].abs()
    raw["image_area_arcsec2"]=raw.target_area_arcsec2.where(~failed)
    real_keys=HOST_KEYS+["condition","noise_realisation","band","method","b_star"]
    counts=raw.groupby(real_keys).size()
    if not counts.eq(3).all():raise ValueError("A host-realisation lacks three projections")
    aggregations={k:"median" for k in METRICS if k!="failure_fraction"}
    aggregations.update(failure_fraction="mean",image_area_arcsec2="median")
    real=raw.groupby(real_keys,as_index=False).agg(aggregations)
    host_keys=HOST_KEYS+["condition","band","method","b_star"]
    host=real.groupby(host_keys,as_index=False)[list(METRICS)+["image_area_arcsec2"]].mean()
    host["n_noise_realisations"]=real.groupby(host_keys).size().to_numpy()
    expected=np.where(host.condition.str.startswith("noise_"),3,1)
    if not np.array_equal(host.n_noise_realisations.to_numpy(),expected):
        raise ValueError("Noise repetition counts differ from frozen design")
    archive=pd.read_csv(ROOT/"results/summary_blind20261030_v5_potential_centered/halo_metrics.csv")
    archive.b_star=archive.b_star.astype(str)
    archive.loc[archive.method.eq("aperture30pkpc"),"b_star"]="aperture"
    archive.loc[archive.method.eq("FoF"),"b_star"]=archive.loc[archive.method.eq("FoF"),"b_star"].map(lambda x:str(float(x)))
    audit=host[host.condition.eq("baseline")].merge(archive,on=HOST_KEYS+["band","method","b_star"],validate="one_to_one")
    if len(audit)!=294:raise ValueError("Baseline host aggregation mismatch")
    for now,old in [("iou","median_iou"),("completeness","median_completeness"),
                    ("logL","median_log10_luminosity_ratio"),("abs_logL","median_abs_log10_luminosity_error"),
                    ("logR50","median_log10_R50_ratio"),("abs_logR50","median_abs_log10_R50_error")]:
        if not np.allclose(audit[now],audit[old],rtol=1e-8,atol=1e-12,equal_nan=True):
            raise ValueError(f"Archived host statistic not reproduced for {now}")
    summary=[]
    for (condition,band,method,b),g in host.groupby(["condition","band","method","b_star"]):
        if len(g)!=49:raise ValueError("Condition lacks 49 independent host entries")
        summary.append({"condition":condition,"band":band,"method":method,"b_star":b,"n_hosts":49,
                        "mean_host_iou":g.iou.mean(),"median_host_iou":g.iou.median(),
                        "fraction_hosts_iou_below_0p90":g.iou.lt(.9).mean(),
                        "fraction_hosts_any_failed_view":g.failure_fraction.gt(0).mean(),
                        "mean_failed_projection_fraction":g.failure_fraction.mean(),
                        "median_host_completeness":g.completeness.median(),
                        "median_host_contamination":g.contamination.median(),
                        "mean_signed_logL_bias":g.logL.mean(),"p90_absolute_logL_error":g.abs_logL.quantile(.9),
                        "mean_signed_logR50_bias":g.logR50.mean(),"p90_absolute_logR50_error":g.abs_logR50.quantile(.9),
                        "n_hosts_valid_luminosity_error":int(g.abs_logL.notna().sum()),
                        "n_hosts_valid_radius_error":int(g.abs_logR50.notna().sum())})
    hosts=host[HOST_KEYS].drop_duplicates().sort_values(HOST_KEYS).reset_index(drop=True)
    if len(hosts)!=49:raise ValueError("Host manifest mismatch")
    weights=resampling_weights(hosts)
    changes=[];methods=[]
    for (condition,band,method,b),g in host.groupby(["condition","band","method","b_star"]):
        if condition=="baseline":continue
        baseline=host[host.condition.eq("baseline")&host.band.eq(band)&host.method.eq(method)&host.b_star.eq(b)]
        changes.extend(paired(g,baseline,weights,{"condition":condition,"band":band,"method":method,"b_star":b,
                                                 "difference":"condition minus baseline"}))
    for (condition,band),g in host.groupby(["condition","band"]):
        fixed=g[g.method.eq("FoF")&g.b_star.eq("0.1")]
        for method,b in [("FoF","0.2"),("aperture30pkpc","aperture")]:
            other=g[g.method.eq(method)&g.b_star.eq(b)]
            methods.extend(paired(fixed,other,weights,{"condition":condition,"band":band,
                                                      "baseline_method":method,"baseline_b_star":b,
                                                      "difference":"fixed b=0.10 minus comparison method"}))
    raw.to_csv(OUT/"projection_metrics.csv",index=False)
    real.to_csv(OUT/"host_realisation_metrics.csv",index=False)
    host.to_csv(OUT/"host_metrics.csv",index=False)
    pd.DataFrame(summary).to_csv(OUT/"sensitivity_summary.csv",index=False)
    pd.DataFrame(changes).to_csv(OUT/"paired_condition_changes.csv",index=False)
    pd.DataFrame(methods).to_csv(OUT/"paired_method_contrasts.csv",index=False)
    metadata={"config_sha256":measurement.CONFIG_SHA,"n_hosts":49,"n_boxes":7,"n_projection_metric_rows":len(raw),
              "baseline_host_aggregation_reproduced":True,"bootstrap_draws":NBOOT,
              "aggregation":"Median over three views for each noise realisation, then mean over the three noise realisations within host. Noiseless conditions have one realisation. Hosts are the statistical units.",
              "failed_target_policy":"Failed centre component gives zero overlap/completeness; precision and luminosity/radius errors are undefined. Error quantiles condition on measurable targets, with failures reported separately.",
              "bootstrap":"Shared host and box-then-host resampling weights for paired differences; views and noise seeds are not independent hosts.",
              "psf_parameters":{"fwhm_arcsec":measurement.base.sb.FWHM_AS,"moffat_beta":measurement.base.sb.MOFFAT_BETA,
                                "wing_amplitudes":measurement.base.sb.WING_A,"kernel_halfwidth_arcsec":measurement.base.sb.KERNEL_AS},
              "source_hashes":{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                               [Path(measurement.__file__),Path(measurement.base.__file__),Path(__file__)]},
              "limits":"Post-baseline one-factor sensitivity at z_obs=0.1904 using z=0 structures and a known simulated centre. Noise is an explicit idealized Gaussian scenario, not a survey data product. Linking parameters are fixed, and the original blinded sample is not used to retune them."}
    (OUT/"summary_metadata.json").write_text(json.dumps(metadata,indent=2)+"\n")
    print(pd.DataFrame(summary).to_string(index=False),flush=True)
    print(json.dumps(metadata,indent=2),flush=True)


if __name__=="__main__":main()
