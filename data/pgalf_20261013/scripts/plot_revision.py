"""Draw revised figures from compact completed calculation products."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
from pathlib import Path
import json,argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy import ndimage
from paper_figure_style import save_paper_figure
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data',type=Path,default=ROOT/'results')
parser.add_argument('--out',type=Path,default=ROOT/'figures')
args=parser.parse_args()
DATA=args.data
OUT=args.out
OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':9,'axes.labelsize':9,'axes.titlesize':9,'legend.fontsize':9,
                     'xtick.labelsize':9,'ytick.labelsize':9,'savefig.bbox':'tight'})
COLORS=['#6b3a95','#d55e00','#0072b2','#009e73','#a65f82','#4d4d4d']
STEPS=['00081','00178','00475','00786','00967','01001']
ZLABEL=['3.002','1.990','1.002','0.500','0.227','0.190']
def read(name):return pd.read_csv(DATA/name,sep='\t',dtype={'step':str})
def save(fig,name):
    save_paper_figure(fig,OUT/name);plt.close(fig)
def inside_title(ax,text,x=.04,y=.96,ha='left',va='top',fontsize=9):
    ax.text(x,y,text,transform=ax.transAxes,ha=ha,va=va,fontsize=fontsize,
            bbox=dict(facecolor='white',edgecolor='none',alpha=.88,pad=1.4),zorder=20)
def grouped(frame,metric):
    return frame.groupby(['step','nstar','b'])[metric].agg(lambda a:float(np.median(a.to_numpy()))).reset_index()
def main():
    curves=read('photometric_curves.tsv');curves=curves[curves.variant=='primary']
    views=read('photometric_view_selection.tsv');views=views[(views.variant=='primary')&(views.weighting=='position')]
    halos=read('photometric_halo_metrics.tsv');halos=halos[(halos.variant=='primary')&(halos.weighting=='position')]
    fig,axes=plt.subplots(1,2,figsize=(7.1,2.5),sharey=True)
    for ax,band in zip(axes,['V','g']):
        for i,step in enumerate(STEPS):
            hh=halos[(halos.band==band)&(halos.step==step)&(halos.interior==1)]
            for offset,(_,h) in zip(np.linspace(-.18,.18,max(len(hh),1)),hh.iterrows()):
                vv=views[(views.band==band)&(views.step==step)&(views.nstar==h.nstar)]
                ax.plot([i+offset]*2,[vv.b_best.min(),vv.b_best.max()],color=COLORS[i],lw=.65)
                bad=h.fitted_view_iou<.9
                ax.scatter(i+offset,h.b,s=24 if bad else 12,marker='s' if bad else 'o',
                           facecolor='white' if bad else COLORS[i],edgecolor=COLORS[i],lw=.8,zorder=5 if bad else 2)
        ax.axhline(.1,color='.35',ls='--',lw=.8,label='$b=0.10$')
        ax.axhline(.2,color='.35',ls=':',lw=.8,label='$b=0.20$')
        ax.set(xticks=range(6),xticklabels=ZLABEL,xlabel='Output redshift',ylim=(0,.215))
        inside_title(ax,band)
    axes[0].set_ylabel('Selected linking parameter $b$');axes[1].legend(loc='upper right',ncol=2)
    fig.tight_layout();save(fig,'b_match.pdf')
    # All 76 entries remain in the displayed curves, including endpoint cases.
    fig,axes=plt.subplots(4,2,figsize=(7.1,7.0),sharex=True,sharey='row')
    for j,band in enumerate(['V','g']):
        cc=curves[curves.band==band]
        for row,metric in enumerate(['iou','completeness','precision']):
            ax=axes[row,j];hh=grouped(cc,metric)
            for _,h in hh.groupby(['step','nstar']):ax.plot(h.b,h[metric],color='.82',lw=.4,zorder=0)
            for i,step in enumerate(STEPS):
                h=hh[hh.step==step].groupby('b')[metric].median()
                ax.plot(h.index,h.values,color=COLORS[i],lw=1.15,label='$z='+ZLABEL[i]+'$')
            ax.set_ylim(0,1.02);ax.axvline(.1,color='.45',ls=':',lw=.65)
            if j==0:ax.set_ylabel(metric.capitalize() if metric!='iou' else 'Overlap')
        # Precision labels occupy the empty lower right, away from the curves.
        for row in range(4):
            label=f'({chr(97+2*row+j)}) {band}'
            if row==2:inside_title(axes[row,j],label,x=.96,y=.08,ha='right',va='bottom')
            elif row==1:inside_title(axes[row,j],label,x=.04,y=.08,va='bottom')
            else:inside_title(axes[row,j],label)
        ax=axes[3,j]
        for width,ls in [('1.0','-'),('2.0','--')]:
            for metric,color in [('completeness','#286090'),('precision','#bf6532')]:
                h=grouped(cc,f'band{width}_{metric}').groupby('b')[f'band{width}_{metric}'].median()
                ax.plot(h.index,h.values,ls,color=color,label=f'{metric.capitalize()}, {width}"')
        ax.set_ylim(0,1.02);ax.set_xlabel('Linking parameter $b$');ax.axvline(.1,color='.45',ls=':',lw=.65)
    axes[3,0].set_ylabel('Boundary light fraction')
    axes[0,1].legend(loc='lower right',ncol=2);axes[3,1].legend(loc='lower right',ncol=1)
    for row in range(4):
        for ax in axes[row]:
            ax.set_ylim(0,1)
            ax.set_yticks([0,.2,.4,.6,.8,1])
            ax.set_yticklabels(['0' if row==3 else '', '.2','.4','.6','.8', '1' if row==0 else ''])
            ax.tick_params(axis='both',which='both',direction='in',top=True,bottom=True,left=True,right=True)
            ax.tick_params(axis='x',labelbottom=row==3)
    fig.tight_layout()
    fig.subplots_adjust(hspace=0)
    save(fig,'b_curves.pdf')
    fig,axes=plt.subplots(2,2,figsize=(7.1,5.0))
    for j,band in enumerate(['V','g']):
        cc=curves[(curves.band==band)&(curves.b==.1)]
        for row in range(2):
            ax=axes[row,j]
            for i,step in enumerate(STEPS):
                h=cc[cc.step==step];good=h.iou>=.9
                ax.scatter(h.loc[good,'completeness'],h.loc[good,'precision'],s=15,color=COLORS[i],lw=0)
                ax.scatter(h.loc[~good,'completeness'],h.loc[~good,'precision'],s=20,facecolor='white',edgecolor=COLORS[i],marker='s',lw=.8,zorder=5)
            lo=.88 if row==0 else -.025
            hi=1.005 if row==0 else 1.025
            ax.plot([max(lo,0),1],[max(lo,0),1],color='.55',ls=':',lw=.7)
            ax.set(xlim=(lo,hi),ylim=(lo,hi),xlabel='Completeness',ylabel='Precision',aspect='equal')
        inside_title(axes[1,j],f'{band}, common $b=0.10$')
    fig.tight_layout();save(fig,'b_residual.pdf')
    purpose(curves)
    failures(curves)
    rest_control(curves)
def rest_control(curves):
    from collections import defaultdict
    rest=defaultdict(list);mass={}
    for line in (DATA/'rest_frame.tsv').read_text().splitlines():
        if not line.strip() or line.startswith('#'):continue
        p=line.split()
        if p[1]=='P':mass[(p[0],int(p[3]))]=float(p[7])
        if p[1]=='I':rest[(p[0],int(p[3]),p[5],round(float(p[2]),3))].append((p[4],float(p[6]),float(p[7])))
    fig,axes=plt.subplots(1,2,figsize=(7.1,2.6));styles=['-','--','-.',':']
    for ax,band in zip(axes,['V','g']):
        for i,(key,rows) in enumerate(sorted((k,v) for k,v in rest.items() if k[2]==band)):
            step,nstar,_,z=key;color=COLORS[STEPS.index(step)]
            xx=sorted({r[1] for r in rows});yy=np.array([[r[2] for r in sorted(rows,key=lambda r:r[1]) if r[0]==los] for los in 'xyz'])
            baseline=curves[(curves.step==step)&(curves.nstar==nstar)&(curves.band==band)].groupby('b').iou.median()
            ax.plot(baseline.index,baseline.values,color=color,ls=styles[i],lw=.7,marker='o',ms=2)
            ax.plot(xx,np.median(yy,axis=0),color=color,ls=styles[i],lw=1.2,
                    label=f'$z={z:.3f}$, $M={mass[(step,nstar)]/1e13:.3f}'+r'\times10^{13}$')
            ax.fill_between(xx,yy.min(axis=0),yy.max(axis=0),color=color,alpha=.12)
        ax.set(xlim=(.01,.51),ylim=(0,1.03),xlabel='Linking parameter $b$',ylabel='Overlap')
        inside_title(ax,band,x=.8,y=.72,ha='center')
    axes[1].legend(loc='lower right',fontsize=9)
    fig.tight_layout();save(fig,'b_rest.pdf')
def purpose(curves):
    thin=read('conserved_sampling_curves.tsv')
    fig,axes=plt.subplots(2,2,figsize=(7.1,4.55))
    for j,band in enumerate(['V','g']):
        cc=curves[curves.band==band]
        for prefix,color,label in [('', '#286090','Projected positions'),('psf_','#bf6532','PSF aperture flux')]:
            h=grouped(cc.assign(error=abs(cc[prefix+'luminosity_error'])),'error')
            med=h.groupby('b').error.median();axes[0,j].plot(med.index,100*med.values,color=color,label=label)
            xx=[];yy=[];ll=[];uu=[]
            for factor in [1,2,4,8]:
                for kind,marker,offset in [('fixed','o',-.035),('rebuilt','s',.035)]:
                    if factor==1:values=grouped(cc[cc.b==.1],prefix+'iou')[prefix+'iou'].values
                    else:
                        a=thin[(thin.band==band)&(thin.factor==factor)&(thin.image==kind)&(thin.b==.1)]
                        values=a.groupby(['step','nstar'])[prefix+'iou'].median().values
                    q=np.quantile(values,[.16,.5,.84]);x=np.log2(factor)+offset
                    axes[1,j].errorbar(x,q[1],yerr=[[q[1]-q[0]],[q[2]-q[1]]],fmt=marker,color=color,
                        ms=4,capsize=2,lw=.7,mfc=color if kind=='fixed' else 'white',
                        label=f'{label}, {kind} image' if factor==2 else None)
        axes[0,j].set(xlabel='Linking parameter $b$',ylabel='Median absolute luminosity error (%)')
        inside_title(axes[0,j],band,x=.5,ha='center')
        axes[0,j].axvline(.1,color='.55',ls=':',lw=.7)
        axes[1,j].set(xticks=range(4),xticklabels=['1','1/2','1/4','1/8'],xlabel='Retained particle fraction',
                      ylabel='Overlap at $b=0.10$',ylim=(0,1.02))
    axes[0,1].legend();axes[1,1].legend(fontsize=9)
    fig.tight_layout();save(fig,'photometry_sampling.pdf')
def failures(curves):
    cases=json.loads((DATA/'failure_drawing.json').read_text())
    drawing=np.load(DATA/'failure_drawing.npz')
    fig,axes=plt.subplots(2,3,figsize=(7.1,5.1))
    for ip,(ax,r) in enumerate(zip(axes.flat,cases)):
        band=r['band'];los=r['los'];half=r['half_arcsec'];mu=drawing[f'mu_{ip}'];mask=drawing[f'mask_{ip}']
        # Magnitudes already encode logarithmic flux. Display the same flux
        # with LogNorm and extend the faint end to 30 mag arcsec^-2.
        flux=np.power(10.,-.4*(mu-26.5))
        display=ax.imshow(np.ma.masked_less_equal(flux,0),origin='lower',
            extent=[-half,half,-half,half],cmap='gray',
            norm=LogNorm(vmin=10**(-.4*(30-26.5)),vmax=10**(-.4*(20-26.5))),
            interpolation='nearest')
        x=(np.arange(mu.shape[0])+.5)*.2-half
        ax.contour(x,x,mask,levels=[.5],colors=['#e4685b'],linewidths=.8)
        ax.scatter(drawing[f'x_{ip}'],drawing[f'y_{ip}'],color='#5199bd',s=1.2,lw=0,alpha=.6,rasterized=True)
        seed=r['seed_offsets_arcsec']
        ax.scatter([0],[0],marker='+',color='cyan',s=45,lw=1)
        ax.scatter([seed[0]],[seed[1]],marker='x',color='#edaa3e',s=35,lw=1)
        if max(abs(seed[0]),abs(seed[1]))>half:
            tip=np.asarray(seed)*(.82*half/max(abs(seed[0]),abs(seed[1])))
            ax.annotate('Seed outside',xy=tip,xytext=tip*.45,color='#edaa3e',fontsize=9,
                        arrowprops=dict(arrowstyle='->',color='#edaa3e',lw=.8))
        ax.set(xlim=(-half,half),ylim=(-half,half),xlabel='Offset (arcsec)',ylabel='Offset (arcsec)',
               aspect='equal')
        inside_title(ax,f'({chr(97+ip)}) {band}, {los}, $z={r["z"]:.3f}$\nOverlap {r["iou"]:.3f}, $b=0.10$',
                     x=.97,y=.03 if ip==0 else .97,ha='right',
                     va='bottom' if ip==0 else 'top',fontsize=9)
    fig.tight_layout(rect=[0,.10,1,1])
    cax=fig.add_axes([.28,.055,.44,.018]);bar=fig.colorbar(display,cax=cax,orientation='horizontal')
    levels=np.array([20,22,24,26,28,30])
    bar.set_ticks(10**(-.4*(levels-26.5)));bar.set_ticklabels(levels)
    bar.ax.invert_xaxis()
    bar.ax.minorticks_off()
    bar.set_label(r'Surface brightness (mag arcsec$^{-2}$)',fontsize=9)
    save(fig,'fail_maps.pdf')
if __name__=='__main__':main()
