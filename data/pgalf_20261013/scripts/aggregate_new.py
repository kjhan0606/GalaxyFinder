"""Summarize the completed consistent-distance fine/PSF/sampling calculation."""
import os
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
import csv,json,argparse
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1]
SUPP=ROOT/'inputs/supplement'
OUT=ROOT/'results'
GRID=np.array([.02,.03,.05,.06,.07,.08,.09,.10,.11,.12,.13,.14,.15,.18,.20])
BOX=44.171267

def read(path):
    return list(csv.DictReader([s for s in path.read_text().splitlines() if s.strip() and not s.startswith('#')],delimiter='\t'))
def write(name,rows):
    with (OUT/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def q(a):
    a=np.asarray(a,float);a=a[np.isfinite(a)]
    return np.quantile(a,[.16,.5,.84]).tolist() if len(a) else []

def main():
    global OUT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,help='Regenerate from released aggregate tables instead of private per-halo products.')
    parser.add_argument('--out',type=Path,default=ROOT/'results')
    parser.add_argument('--primary-only',action='store_true',help='Summarize all completed fine curves without sampling products.')
    args=parser.parse_args();OUT=args.out;OUT.mkdir(parents=True,exist_ok=True)
    metadata_dir=args.data if args.data else SUPP
    meta=read(metadata_dir/'halo_metadata.tsv');keys=[(r['step'],int(r['nstar'])) for r in meta]
    masslookup={(r['step'],int(r['nstar'])):float(r['group_mass_hinv_msun']) for r in meta}
    dirs=[ROOT/'results'/f'{k[0]}_{k[1]}' for k in keys]
    missing=[d.name for d in dirs if not (d/('fine_curves.tsv' if args.primary_only else 'DONE')).is_file()]
    if missing and not args.data:raise SystemExit('Calculation incomplete: '+str(len(missing))+' missing entries')
    fine=[];thin=[];targets=[];totals=[];mass=[];properties=[]
    if args.data:
        fine=read(args.data/'photometric_curves.tsv');thin=[] if args.primary_only else read(args.data/'conserved_sampling_curves.tsv')
        targets=read(args.data/'photometric_targets.tsv');totals=[] if args.primary_only else read(args.data/'sampling_total_ratios.tsv')
        mass=json.loads((args.data/'stellar_mass_linking_quantiles.json').read_text())
        properties=read(args.data/'stellar_population_metrics.tsv')
    else:
        for d in dirs:
            fine+=read(d/'fine_curves.tsv');targets+=read(d/'targets.tsv')
            if not args.primary_only:
                thin+=read(d/'thin_curves.tsv');totals+=read(d/'sampling_totals.tsv')
            mass.append(json.loads((d/'mass_scale.json').read_text()))
            properties+=read(d/'stellar_properties.tsv')
    assert len(keys)==76 and {(r['step'],int(r['nstar'])) for r in fine}==set(keys)
    (OUT/'stellar_mass_linking_quantiles.json').write_text(json.dumps(mass,indent=2)+'\n')
    write('photometric_curves.tsv',fine)
    if thin:write('conserved_sampling_curves.tsv',thin)
    write('photometric_targets.tsv',targets)
    if totals:write('sampling_total_ratios.tsv',totals)
    for r in properties:
        for kind in ['position','aperture']:
            age=float(r[f'age_{kind}_Gyr'])
            r[f'fractional_age_error_{kind}']=float(r['age_group_Gyr'])/age-1 if age>0 else np.nan
    write('stellar_population_metrics.tsv',properties)
    proplookup={(r['step'],int(r['nstar']),r['band'],r['los'],float(r['b'])):r for r in properties}
    lookup=defaultdict(list)
    for r in fine:lookup[(r['step'],int(r['nstar']),r['variant'],r['band'],r['los'])].append(r)
    halo=[];view=[];curves={};selected={}
    for key,rs in lookup.items():
        rs.sort(key=lambda r:float(r['b']))
        assert np.allclose([float(r['b']) for r in rs],GRID)
        for score in (['iou','psf_iou'] if key[2]=='primary' else ['iou']):
            values=np.array([float(r[score]) for r in rs]);j=int(np.argmax(values))
            prefix='psf_' if score=='psf_iou' else ''
            out=dict(step=key[0],nstar=key[1],variant=key[2],band=key[3],los=key[4],weighting='PSF' if prefix else 'position',
                     z=float(rs[0]['z']),b_best=float(GRID[j]),edge=int(j in [0,len(GRID)-1]),iou=float(values[j]),
                     width=float(np.ptp(GRID[values>=values[j]-.02])))
            for b in [.08,.09,.10,.20]:out[f'iou_b{b:.2f}']=float(values[int(np.argmin(abs(GRID-b)))])
            view.append(out);selected[(*key,score)]=out
    for k in keys:
        for variant in ['primary','threshold_bright','threshold_faint','wing_off','light_peak']:
            for band in ['V','g']:
                rs=[lookup[(*k,variant,band,los)] for los in 'xyz']
                for score in (['iou','psf_iou'] if variant=='primary' else ['iou']):
                    prefix='psf_' if score=='psf_iou' else ''
                    matrix=np.array([[float(r[score]) for r in a] for a in rs])
                    curves[(*k,variant,band,score)]=np.median(matrix,axis=0)
                    vv=[selected[(*k,variant,band,los,score)] for los in 'xyz']
                    best=np.median([r['b_best'] for r in vv])
                    row=dict(step=k[0],nstar=k[1],z=float(rs[0][0]['z']),variant=variant,band=band,
                             weighting='PSF' if prefix else 'position',interior=int(all(not r['edge'] for r in vv)),
                             b=float(best),fitted_view_iou=float(np.median(matrix.max(axis=1))),
                             best_common_b_iou=float(np.median(matrix,axis=0).max()),
                             selected_halo_iou=float(np.median(matrix[:,int(np.argmin(abs(GRID-best)))])),
                             width=float(np.median([r['width'] for r in vv])))
                    for b in [.08,.09,.10,.20]:
                        j=int(np.argmin(abs(GRID-b)))
                        row[f'iou_b{b:.2f}']=float(np.median(matrix[:,j]))
                        row[f'error_b{b:.2f}']=float(np.median([abs(float(a[j][prefix+'luminosity_error'])) for a in rs]))
                    halo.append(row)
    write('photometric_view_selection.tsv',view);write('photometric_halo_metrics.tsv',halo)
    positions={(r['step'],int(r['nstar'])):np.array([float(r[c]) for c in ['x_hinv_cMpc','y_hinv_cMpc','z_hinv_cMpc']]) for r in read(metadata_dir/'halo_positions.tsv')}
    folds=[];pred=[]
    for step in sorted({k[0] for k in keys}):
        held=[k for k in keys if k[0]==step];train=[]
        for k in keys:
            if k[0]==step:continue
            delta=np.array([positions[k]-positions[h] for h in held]);delta-=BOX*np.rint(delta/BOX)
            if np.min(np.linalg.norm(delta,axis=1))>=.5:train.append(k)
        for band in ['V','g']:
            for score in ['iou','psf_iou']:
                objective=np.median([curves[(*k,'primary',band,score)] for k in train],axis=0)
                j=int(np.argmax(objective));weighting='PSF' if score=='psf_iou' else 'position'
                folds.append(dict(step=step,band=band,weighting=weighting,n_train=len(train),n_held=len(held),selected_b=float(GRID[j])))
                for k in held:
                    a=curves[(*k,'primary',band,score)]
                    pp=[proplookup[(*k,band,los,float(GRID[j]))] for los in 'xyz']
                    cc=[lookup[(*k,'primary',band,los)][j] for los in 'xyz']
                    prefix='psf_' if weighting=='PSF' else ''
                    pred.append(dict(step=k[0],nstar=k[1],band=band,weighting=weighting,selected_b=float(GRID[j]),
                                     held_out_iou=float(a[j]),iou_b010=float(a[np.argmin(abs(GRID-.1))]),
                                     iou_b020=float(a[-1]),best_common_b_iou=float(a.max()),shortfall_from_best_common_b=float(a.max()-a[j]),
                                     absolute_luminosity_error=float(np.median([abs(float(r[prefix+'luminosity_error'])) for r in cc])),
                                     absolute_mass_error_position=float(np.median([abs(float(r['mass_error'])) for r in pp])),
                                     absolute_fractional_age_error_position=float(np.median([abs(float(r['fractional_age_error_position'])) for r in pp])),
                                     absolute_fractional_age_error_aperture=float(np.median([abs(float(r['fractional_age_error_aperture'])) for r in pp]))))
    write('fine_validation_folds.tsv',folds);write('fine_validation_predictions.tsv',pred)
    summary={};strata=[];performance=[];controls=[];boundary=[]
    for band in ['V','g']:
        summary[band]={}
        for weighting in ['position','PSF']:
            rs=[r for r in halo if r['variant']=='primary' and r['band']==band and r['weighting']==weighting]
            fit=[r for r in rs if r['interior']]
            summary[band][weighting]=dict(n_all=len(rs),n_fitted=len(fit),b=q([r['b'] for r in fit]),fitted_iou=q([r['fitted_view_iou'] for r in fit]),width=q([r['width'] for r in fit]),
                        selected_halo_iou=q([r['selected_halo_iou'] for r in fit]),view_fit_loss=q([r['fitted_view_iou']-r['selected_halo_iou'] for r in fit]),
                        fixed={str(b):dict(iou=q([r[f'iou_b{b:.2f}'] for r in rs]),error=q([r[f'error_b{b:.2f}'] for r in rs]),n_below090=sum(r[f'iou_b{b:.2f}']<.9 for r in rs)) for b in [.08,.09,.10,.20]})
            pp=[r for r in pred if r['band']==band and r['weighting']==weighting]
            summary[band][weighting]['validation']=dict(held_out_iou=q([r['held_out_iou'] for r in pp]),gain_over_b020=q([r['held_out_iou']-r['iou_b020'] for r in pp]),
                shortfall_from_best_common_b=q([r['shortfall_from_best_common_b'] for r in pp]),n_below090=sum(r['held_out_iou']<.9 for r in pp),fold_b=sorted({r['selected_b'] for r in pp}))
            for field in ['absolute_luminosity_error','absolute_mass_error_position','absolute_fractional_age_error_position','absolute_fractional_age_error_aperture']:
                values=[r[field] for r in pp]
                summary[band][weighting]['validation'][field]=dict(n_defined=int(np.count_nonzero(np.isfinite(values))),quantiles=q(values))
            for step in ['all']+sorted({k[0] for k in keys}):
                pp=[r for r in pred if r['band']==band and r['weighting']==weighting and (step=='all' or r['step']==step)]
                performance.append(dict(step=step,band=band,weighting=weighting,n=len(pp),b_min=min(r['selected_b'] for r in pp),b_max=max(r['selected_b'] for r in pp),
                    held_out_iou=float(np.median([r['held_out_iou'] for r in pp])),best_common_b_iou=float(np.median([r['best_common_b_iou'] for r in pp])),
                    paired_shortfall=float(np.median([r['shortfall_from_best_common_b'] for r in pp])),gain=float(np.median([r['held_out_iou']-r['iou_b020'] for r in pp])),n_below090=sum(r['held_out_iou']<.9 for r in pp)))
            for label,filter_ in [('all',lambda r:True),('low_z',lambda r:r['z']<=.23),('high_z',lambda r:r['z']>=.5),
                                    ('small_region',lambda r:True),('large_region',lambda r:True)]:
                primary_targets={(r['step'],int(r['nstar']),r['band']):[] for r in targets if r['variant']=='primary'}
                for r in targets:
                    if r['variant']=='primary':primary_targets[(r['step'],int(r['nstar']),r['band'])].append(float(r['req_arcsec']))
                ss=[r for r in rs if filter_(r)]
                if label in ['small_region','large_region']:
                    ss=[r for r in ss if (np.median(primary_targets[(r['step'],r['nstar'],r['band'])])<3)==(label=='small_region')]
                ff=[r for r in ss if r['interior']]
                bq=q([r['b'] for r in ff]);strata.append(dict(stratum=label,band=band,weighting=weighting,n=len(ss),n_fitted=len(ff),
                     b16=bq[0] if bq else np.nan,b50=bq[1] if bq else np.nan,b84=bq[2] if bq else np.nan,
                     iou_b010=float(np.median([r['iou_b0.10'] for r in ss])),error_b010=float(np.median([r['error_b0.10'] for r in ss])),n_below090=sum(r['iou_b0.10']<.9 for r in ss)))
        base={(r['step'],r['nstar']):r for r in halo if r['variant']=='primary' and r['band']==band and r['weighting']=='position'}
        for variant in ['threshold_bright','threshold_faint','wing_off','light_peak']:
            rs=[r for r in halo if r['variant']==variant and r['band']==band and r['interior'] and base[(r['step'],r['nstar'])]['interior']]
            controls.append(dict(variant=variant,band=band,n_fitted=len(rs),b_change=q([r['b']-base[(r['step'],r['nstar'])]['b'] for r in rs]),
                                 score_change=q([r['fitted_view_iou']-base[(r['step'],r['nstar'])]['fitted_view_iou'] for r in rs])))
        for width in ['1.0','2.0']:
            for b in [.08,.09,.10,.20]:
                for metric in ['iou','completeness','precision']:
                    vals=[np.median([float(next(r for r in lookup[(*k,'primary',band,los)] if float(r['b'])==b)[f'band{width}_{metric}']) for los in 'xyz']) for k in keys]
                    qq=q(vals);boundary.append(dict(band=band,width_arcsec=float(width),b=b,metric=metric,n_defined=int(np.count_nonzero(np.isfinite(vals))),p16=qq[0] if qq else np.nan,median=qq[1] if qq else np.nan,p84=qq[2] if qq else np.nan))
    write('fine_validation_summary.tsv',performance);write('photometric_strata.tsv',strata);write('wide_boundary_summary.tsv',boundary)
    tgroups=defaultdict(list)
    for r in thin:tgroups[(r['factor'],r['image'],r['step'],int(r['nstar']),r['band'],r['los'])].append(r)
    thin_summary=[]
    for factor in ([] if args.primary_only else ['2','4','8']):
        for image in ['fixed','rebuilt']:
            for band in ['V','g']:
                for weighting,score in [('position','iou'),('PSF','psf_iou')]:
                    points=[];nlow=nupper=0;fixed=[];best=[]
                    for k in keys:
                        rr=[sorted(tgroups[(factor,image,*k,band,los)],key=lambda r:float(r['b'])) for los in 'xyz']
                        js=[int(np.argmax([float(r[score]) for r in a])) for a in rr]
                        lower=any(j==0 for j in js);upper=any(j==len(a)-1 for a,j in zip(rr,js))
                        nlow+=lower;nupper+=upper
                        if not lower and not upper:points.append(float(np.median([float(a[j]['b']) for a,j in zip(rr,js)])))
                        fixed.append(float(np.median([float(next(r for r in a if float(r['b'])==.1)[score]) for a in rr])))
                        best.append(float(np.median([float(a[j][score]) for a,j in zip(rr,js)])))
                    thin_summary.append(dict(factor=int(factor),image=image,band=band,weighting=weighting,n_all=len(keys),n_fitted=len(points),n_lower=nlow,n_upper=nupper,
                        b=q(points),iou_b010=q(fixed),n_below090_b010=sum(v<.9 for v in fixed),best_iou=q(best)))
    # Failures carry descriptive image flags, allowing overlapping categories.
    flags={(r['step'],int(r['nstar']),r['band'],r['los']):r for r in targets if r['variant']=='primary'}
    failed=[];flagcounts=[]
    predictions={(r['step'],r['nstar'],r['band'],r['weighting']):r for r in pred}
    for band in ['V','g']:
        for weighting in ['position','PSF']:
            vv=[v for v in view if v['variant']=='primary' and v['band']==band and v['weighting']==weighting]
            for v in vv:
                p=predictions[(v['step'],v['nstar'],band,weighting)]
                row=next(r for r in lookup[(v['step'],v['nstar'],'primary',band,v['los'])] if float(r['b'])==p['selected_b'])
                v['held_out_b']=p['selected_b'];v['held_out_view_iou']=float(row['psf_iou' if weighting=='PSF' else 'iou'])
            lows=[v for v in vv if v['iou_b0.10']<.9]
            count=dict(band=band,weighting=weighting,n_views=len(vv),n_failed=len(lows),n_failed_held_out=sum(v['held_out_view_iou']<.9 for v in vv),small_region=0,seed_outside=0,fallback=0,offset_gt3=0)
            for v in vv:
                if v['iou_b0.10']>=.9 and v['held_out_view_iou']>=.9:continue
                t=flags[(v['step'],v['nstar'],band,v['los'])];f=dict(v)
                for key in ['area_arcsec2','req_arcsec','fallback','seed_inside','seed_offset_arcsec','n_components','empty','one_pixel']:f[key]=t[key]
                f['group_mass_hinv_msun']=masslookup[(v['step'],v['nstar'])]
                f['failed_common_b010']=int(v['iou_b0.10']<.9);f['failed_held_out']=int(v['held_out_view_iou']<.9)
                f['small_region']=int(float(t['req_arcsec'])<3);f['offset_gt3']=int(float(t['seed_offset_arcsec'])>3)
                if f['failed_common_b010']:
                    count['small_region']+=f['small_region'];count['seed_outside']+=1-int(t['seed_inside']);count['fallback']+=int(t['fallback']);count['offset_gt3']+=f['offset_gt3']
                failed.append(f)
            flagcounts.append(count)
    write('failure_view_catalogue.tsv',failed);write('failure_flags_summary.tsv',flagcounts)
    failedhalos=[]
    for h in halo:
        if h['variant']!='primary':continue
        key=(h['step'],h['nstar'],h['band'],h['weighting']);p=predictions[key]
        if h['iou_b0.10']>=.9 and p['held_out_iou']>=.9:continue
        score='psf_iou' if h['weighting']=='PSF' else 'iou'
        vv=[selected[(h['step'],h['nstar'],'primary',h['band'],los,score)] for los in 'xyz']
        tt=[flags[(h['step'],h['nstar'],h['band'],los)] for los in 'xyz']
        failedhalos.append(dict(step=h['step'],nstar=h['nstar'],z=h['z'],band=h['band'],weighting=h['weighting'],
             group_mass_hinv_msun=masslookup[(h['step'],h['nstar'])],iou_b010=h['iou_b0.10'],held_out_b=p['selected_b'],held_out_iou=p['held_out_iou'],
             n_failed_views_b010=sum(v['iou_b0.10']<.9 for v in vv),n_failed_views_held_out=sum(v['held_out_view_iou']<.9 for v in vv),
             median_req_arcsec=float(np.median([float(t['req_arcsec']) for t in tt])),max_seed_offset_arcsec=max(float(t['seed_offset_arcsec']) for t in tt),
             n_seed_outside=sum(1-int(t['seed_inside']) for t in tt),n_fallback=sum(int(t['fallback']) for t in tt)))
    write('failure_halo_catalogue.tsv',failedhalos)
    correlations=[]
    for band in ['V','g']:
        for step in ['all']+sorted({k[0] for k in keys}):
            rr=[r for r in halo if r['variant']=='primary' and r['weighting']=='position' and r['band']==band and r['interior'] and (step=='all' or r['step']==step)]
            if len(rr)>2:
                rho,p=spearmanr([np.log10(masslookup[(r['step'],r['nstar'])]) for r in rr],[r['b'] for r in rr])
            else:rho=p=np.nan
            correlations.append(dict(step=step,band=band,n=len(rr),spearman_rho=float(rho),two_sided_p=float(p)))
    write('fine_mass_correlations.tsv',correlations)
    matched=[]
    hlookup={(r['step'],r['nstar'],r['band'],r['weighting']):r for r in halo if r['variant']=='primary'}
    for pair in read((args.data if args.data else ROOT/'results')/'known_positional_pairs.tsv'):
        for band in ['V','g']:
            for weighting in ['position','PSF']:
                a=hlookup[(pair['step_a'],int(pair['nstar_a']),band,weighting)]
                b=hlookup[(pair['step_b'],int(pair['nstar_b']),band,weighting)]
                matched.append(dict(**pair,band=band,weighting=weighting,interior_a=a['interior'],interior_b=b['interior'],
                    b_a=a['b'],b_b=b['b'],b_difference=b['b']-a['b'],
                    iou_b010_difference=b['iou_b0.10']-a['iou_b0.10']))
    write('fine_matched_pair_comparison.tsv',matched)
    propsummary={}
    for band in ['V','g']:
        propsummary[band]={}
        for b in [.08,.09,.10,.20]:
            pp=[r for r in properties if r['band']==band and float(r['b'])==b]
            byhalo=defaultdict(list)
            for r in pp:byhalo[(r['step'],r['nstar'])].append(r)
            propsummary[band][str(b)]={}
            for metric in ['mass_error','age_difference_position_Gyr','age_difference_aperture_Gyr','fractional_age_error_position','fractional_age_error_aperture']:
                vals=[np.median([abs(float(r[metric])) for r in rr]) for rr in byhalo.values()]
                propsummary[band][str(b)][metric]=dict(n_defined=int(np.count_nonzero(np.isfinite(vals))),quantiles=q(vals))
    report=dict(grid=GRID.tolist(),primary=summary,sampling=thin_summary,photometric_controls=controls,stellar_properties=propsummary,
                mass_scale=mass,sampling_totals={str(f):{field:q([float(r[field]) for r in totals if int(r['factor'])==f]) for field in ['mass_ratio','flux_ratio_V','flux_ratio_g']} for f in [2,4,8]},
                failure_flags=flagcounts,sampling_scope='Inverse-probability reweighting preserves total mass/light in expectation, one nested realization. Original fixed image and rebuilt image. Not hydrodynamic convergence.')
    (OUT/'new_summary.json').write_text(json.dumps(report,indent=2)+'\n');print('Final numerical summaries written to',OUT)

if __name__=='__main__':main()
