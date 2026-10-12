#!/usr/bin/env python3
"""Draw publication figures from archived tables without particle calculations."""
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures'
Z=[.190,.227,.500,1.002,1.990,3.002]
COL=dict(zip(Z,['#1F4E79','#148F9A','#5B8C3E','#E0A100','#D4652F','#8E3A73']))
B=np.array([.02,.03,.05,.08,.10,.12,.15,.18,.20])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'pdf.fonttype':42})
import csv
curves=defaultdict(list)
for row in csv.DictReader((ROOT/'supplement/recomputed_curves.tsv').open(),delimiter='\t'):
 if row['variant']=='consistent_distance':
  key=(row['step'],int(row['nstar']),row['los'],row['band']);curves[key].append(row)
views=[]
for key,rows in curves.items():
 rows.sort(key=lambda row:float(row['b']))
 a=np.array([[float(row[k]) for k in ['b','iou','completeness','precision']] for row in rows]);assert np.allclose(a[:,0],B)
 j=int(np.argmax(a[:,1]));best=float(B[j]);edge=int(j in (0,len(B)-1))
 views.append(dict(key=key,z=round(float(rows[0]['z']),3),best=best,edge=edge,a=a,score=a[j,1:],rows=rows))
assert len(views)==456
legend=[Line2D([],[],color=COL[z],label=f'$z={z:.3f}$',lw=1.8) for z in Z]
def base(ax):
 ax.tick_params(direction='out',length=3,width=.6); ax.grid(alpha=.14,lw=.5); ax.set_axisbelow(True)
def save(fig,name):
 fig.savefig(OUT/name,bbox_inches='tight',pad_inches=.035); plt.close(fig)
# Three rows separate the three ratios. Full-range IoU insets retain all failures.
fig,axes=plt.subplots(4,2,figsize=(7.05,7.65),sharex=True)
for c,band in enumerate(['V','g']):
 rows=[v for v in views if v['key'][3]==band and not v['edge']]; assert len(rows)==226
 for r,metric in enumerate(['Overlap','Completeness','Precision']):
  ax=axes[r,c]; base(ax)
  if r==0:
   for v in rows: ax.plot(B,v['a'][:,1],color='.65',lw=.35,alpha=.3)
  for z in Z:
   a=np.array([v['a'][:,r+1] for v in rows if v['z']==z]); ax.plot(B,np.median(a,axis=0),color=COL[z],lw=1.5)
  for v in rows:
   if v['key'][1]==7062464: ax.plot(B,v['a'][:,r+1],color='.1',lw=.7,alpha=.8)
  ax.set_xlim(.015,.205); ax.set_ylim((.85,1.006) if r==0 else (0,1.035)); ax.set_xticks([.02,.05,.10,.15,.20]); ax.set_ylabel(metric)
  ax.text(.03,.97,f'({chr(97+r*2+c)})',transform=ax.transAxes,va='top')
  if r==0:
   ax.set_title('$V$' if band=='V' else 'SDSS $g$'); inset=ax.inset_axes([.10,.12,.37,.37]);base(inset)
   for v in rows: inset.plot(B,v['a'][:,1],color='.7',lw=.25,alpha=.3)
   for z in Z: inset.plot(B,np.median([v['a'][:,1] for v in rows if v['z']==z],axis=0),color=COL[z],lw=.8)
   for v in rows:
    if v['key'][1]==7062464: inset.plot(B,v['a'][:,1],color='.1',lw=.45)
   inset.set(xlim=(.015,.205),ylim=(0,1.03),xticks=[.02,.10,.20],yticks=[0,.5,1]);inset.tick_params(labelsize=6,pad=1,length=2)

for c,band in enumerate(['V','g']):
 ax=axes[3,c];base(ax);groups=defaultdict(list)
 for v in views:
  if v['key'][3]==band:groups[v['key'][:2]].append(v)
 good=[vv for vv in groups.values() if not any(v['edge'] for v in vv)]
 for field,color,label in [('ring_iou','#1F4E79','Overlap'),('ring_completeness','#148F9A','Completeness'),('ring_precision','#D4652F','Precision')]:
  data=np.array([np.median([[float(row[field]) for row in v['rows']] for v in vv],axis=0) for vv in good])
  median=np.nanmedian(data,axis=0)
  if field=='ring_precision':median[np.sum(np.isfinite(data),axis=0)<len(good)]=np.nan
  ax.plot(B,median,color=color,lw=1.5,label=label)
 ax.set(xlim=(.015,.205),ylim=(0,1.035),xlabel='Linking parameter $b$',ylabel='Boundary-band ratios')
 ax.text(.03,.97,f'({chr(103+c)})',transform=ax.transAxes,va='top')
 ax.axvline(.10,color='.5',lw=.7,ls=':')
 if c==1:ax.legend(loc='lower right',fontsize=6.5,frameon=False)
fig.legend(handles=legend,loc='upper center',bbox_to_anchor=(.5,1.01),ncol=6,frameon=False,columnspacing=.9,handlelength=1.4)
fig.tight_layout(rect=(0,0,1,.97),h_pad=.7);save(fig,'b_curves.pdf')
# Halo medians use the published best nodes, including the original tie convention.
fig,axes=plt.subplots(1,2,figsize=(7.05,2.35),sharey=True)
for c,band in enumerate(['V','g']):
 ax=axes[c];base(ax); ax.axhline(.20,color='.45',ls=':',lw=.8)
 for zi,z in enumerate(Z):
  grouped=defaultdict(list)
  for v in views:
   if v['z']==z and v['key'][3]==band: grouped[v['key'][:2]].append(v)
  good=[vv for vv in grouped.values() if len(vv)==3 and not any(v['edge'] for v in vv)]
  assert len(good)==[12,9,33,17,3,1][zi]
  for vv,xoff in zip(good,np.linspace(-.30,.30,len(good)) if len(good)>1 else [0]):
   bb=[v['best'] for v in vv]; xx=zi+xoff; yy=np.median(bb); low=np.median([v['score'][0] for v in vv])<.90
   ax.plot([xx,xx],[min(bb),max(bb)],color=COL[z],lw=.8,alpha=.8)
   ax.plot(xx,yy,marker='s' if low else 'o',ms=5.2 if low else 3.5,mfc='white' if low else COL[z],mec=COL[z],mew=1,zorder=5 if low else 3)
 ax.set(xticks=range(6),xticklabels=[f'{z:.3f}' for z in Z],ylim=(.025,.215),xlabel='Simulation output (redshift)',title='$V$' if band=='V' else 'SDSS $g$')
 ax.text(.025,.95,f'({chr(97+c)})',transform=ax.transAxes,va='top')
axes[0].set_ylabel('Matched linking parameter $b$');fig.tight_layout(w_pad=.6);save(fig,'b_match.pdf')
# Completeness-precision planes, with all views in insets.
fig,axes=plt.subplots(1,2,figsize=(7.05,4.0))
fig.subplots_adjust(left=.085,right=.98,bottom=.46,top=.89,wspace=.28)
for c,band in enumerate(['V','g']):
 ax=axes[c]; pos=ax.get_position(); inset=fig.add_axes([pos.x0+.07,.05,.23,.23]); rows=[v for v in views if v['key'][3]==band]
 for target,small in [(ax,False),(inset,True)]:
  base(target); target.plot([0,1],[0,1],color='.65',lw=.7)
  for v in rows:
   iou,rec,pur=v['score']; marker='x' if v['edge'] else ('s' if iou<.90 else 'o'); color=COL[v['z']]
   target.plot(rec,pur,marker=marker,ms=(3 if small else 5) if marker!='s' else (4 if small else 6),mfc='none' if marker=='s' else color,mec=color,mew=.8,ls='none',zorder=5 if marker!='o' else 2,alpha=.8)
  target.set_xlim((-.02,1.025) if small else (.90,1.002));target.set_ylim((-.02,1.025) if small else (.90,1.002))
  if small:target.set_xticks([0,.5,1]);target.set_yticks([0,.5,1]);target.tick_params(labelsize=6,length=2,pad=1);target.set_title('All views',fontsize=6,pad=2)
 ax.set(xlabel='Completeness',ylabel='Precision',title='$V$' if band=='V' else 'SDSS $g$');ax.text(.03,.96,f'({chr(97+c)})',transform=ax.transAxes,va='top')
fig.legend(handles=legend,loc='upper center',bbox_to_anchor=(.5,1.04),ncol=6,frameon=False,columnspacing=.9,handlelength=1.3)
save(fig,'b_residual.pdf')
# Separate control figure. Stable redshift colours, styles distinguish equal-redshift halos.
rest=defaultdict(list); masses={}
for line in (ROOT/'supplement/rest_frame.tsv').read_text().splitlines():
 if line.startswith('#') or not line.strip():continue
 p=line.split()
 if p[1]=='P':masses[(p[0],int(p[3]))]=float(p[7])
 if p[1]=='I':rest[(p[0],int(p[3]),p[5],round(float(p[2]),3))].append([p[4],float(p[6]),float(p[7])])
fig,axes=plt.subplots(1,2,figsize=(7.05,2.5)); styles=['-','--','-.',':']
for c,band in enumerate(['V','g']):
 ax=axes[c];base(ax)
 for k,((step,ns,bb,z),rs) in enumerate(sorted((k,v) for k,v in rest.items() if k[2]==band)):
  xx=sorted(set(row[1] for row in rs)); y=np.array([[row[2] for row in sorted(rs,key=lambda row:row[1]) if row[0]==los] for los in 'xyz'])
  baseline=[v for v in views if v['key'][:2]==(step,ns) and v['key'][3]==band]
  ax.plot(B,np.median([v['a'][:,1] for v in baseline],axis=0),color=COL[z],ls=styles[k],lw=.8,marker='o',ms=2,alpha=.75)
  mass=masses[(step,ns)]; ax.plot(xx,np.median(y,axis=0),color=COL[z],ls=styles[k],lw=1.4,label=f'$z={z:.3f}$, $M={mass/1e13:.3f}'+r'\times10^{13}$')
  ax.fill_between(xx,y.min(axis=0),y.max(axis=0),color=COL[z],alpha=.09)
 ax.set(xlim=(.01,.51),ylim=(0,1.035),xlabel='Linking parameter $b$',ylabel='Overlap',title='$V$' if band=='V' else 'SDSS $g$');ax.text(.03,.96,f'({chr(97+c)})',transform=ax.transAxes,va='top')
axes[1].legend(loc='lower right',fontsize=7,frameon=False);fig.tight_layout(w_pad=.8);save(fig,'b_rest.pdf')
print('Primary score figures use consistent-distance curves. Rest-frame control retains unredshifted absolute flux.')

# Place total and boundary scores together at the same linking parameters.
fig,axes=plt.subplots(2,2,figsize=(7.05,4.5),sharex=True)
for c,band in enumerate(['V','g']):
 grouped=defaultdict(list)
 for v in views:
  if v['key'][3]==band:grouped[v['key'][:2]].append(v)
 good=[vv for vv in grouped.values() if not any(v['edge'] for v in vv)]
 for r in range(2):
  ax=axes[r,c];base(ax)
  for field,color,label in [('iou','#1F4E79','Overlap'),('completeness','#148F9A','Completeness'),('precision','#D4652F','Precision')]:
   col=field if r==0 else 'ring_'+field
   arr=np.array([np.median([[float(row[col]) for row in v['rows']] for v in vv],axis=0) for vv in good])
   median=np.nanmedian(arr,axis=0)
   if col=='ring_precision':median[np.sum(np.isfinite(arr),axis=0)<len(good)]=np.nan
   ax.plot(B,median,color=color,label=label,lw=1.6)
  ax.axvline(.1,color='.4',lw=.8,ls=':');ax.set_ylim(0,1.03);ax.set_xlim(.015,.205)
  ax.text(.03,.96,f'({chr(97+r*2+c)})',transform=ax.transAxes,va='top')
  ax.set_ylabel('All stellar light' if r==0 else 'Light near the isophote')
  if r==0:ax.set_title('$V$' if band=='V' else 'SDSS $g$')
  else:ax.set_xlabel('Linking parameter $b$')
axes[0,0].legend(loc='lower right',frameon=False,fontsize=8)
fig.tight_layout();save(fig,'b_boundary.pdf')
