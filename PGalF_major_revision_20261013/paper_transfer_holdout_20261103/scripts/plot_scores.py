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
curves=defaultdict(list)
for line in (ROOT/'supplement/sightline_baseline.tsv').read_text().splitlines():
 if line.startswith('#') or not line.strip(): continue
 p=line.split(); key=(p[0],int(p[2]),p[3],p[4]); curves[key].append(p)
views=[]
for key,rows in curves.items():
 rows.sort(key=lambda p:float(p[8])); a=np.array([[float(x) for x in p[8:12]] for p in rows]); assert np.allclose(a[:,0],B)
 best=float(rows[0][6]); edge=int(rows[0][7]); j=np.flatnonzero(np.isclose(B,best))[0]
 views.append(dict(key=key,z=round(float(rows[0][1]),3),best=best,edge=edge,a=a,score=a[j,1:]))
assert len(views)==456
legend=[Line2D([],[],color=COL[z],label=f'$z={z:.3f}$',lw=1.8) for z in Z]
def base(ax):
 ax.tick_params(direction='in',top=True,right=True,length=3,width=.6); ax.grid(alpha=.14,lw=.5); ax.set_axisbelow(True)
def save(fig,name):
 fig.savefig(OUT/name,bbox_inches='tight',pad_inches=.035); plt.close(fig)
# Three rows separate the three ratios. Full-range IoU insets retain all failures.
fig,axes=plt.subplots(3,2,figsize=(7.05,6.15),sharex=True)
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
  ax.set_xlim(.015,.205); ax.set_ylim((.85,1.006) if r==0 else (0,1.035)); ax.set_xticks([.02,.05,.10,.15,.20]); ax.set_ylabel(r'$IoU$' if r==0 else r'$C$' if r==1 else r'$P$')
  ax.text(.03,.97,f'({chr(97+r*2+c)})',transform=ax.transAxes,va='top',fontsize=10)
  if r==0:
   ax.set_title('$V$' if band=='V' else 'SDSS $g$',fontsize=10); inset=ax.inset_axes([.10,.12,.37,.37]);base(inset)
   for v in rows: inset.plot(B,v['a'][:,1],color='.7',lw=.25,alpha=.3)
   for z in Z: inset.plot(B,np.median([v['a'][:,1] for v in rows if v['z']==z],axis=0),color=COL[z],lw=.8)
   for v in rows:
    if v['key'][1]==7062464: inset.plot(B,v['a'][:,1],color='.1',lw=.45)
   inset.set(xlim=(.015,.205),ylim=(0,1.03),xticks=[.02,.10,.20],yticks=[0,.5,1]);inset.tick_params(labelsize=6,pad=1,length=2)
  if r==2: ax.set_xlabel(r'$b_{\mathrm{link}}$')
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
 ax.set(xticks=range(6),xticklabels=[f'{z:.3f}' for z in Z],ylim=(.025,.215),xlabel=r'$z$',title='$V$' if band=='V' else 'SDSS $g$')
 ax.text(.025,.95,f'({chr(97+c)})',transform=ax.transAxes,va='top',fontsize=10)
axes[0].set_ylabel(r'$b_{\mathrm{link}}$');fig.tight_layout(w_pad=.6);save(fig,'b_match.pdf')
# Completeness-precision planes, with all views in insets.
fig,axes=plt.subplots(1,2,figsize=(7.05,3.55))
fig.subplots_adjust(left=.085,right=.98,bottom=.43,top=.89,wspace=.28)
for c,band in enumerate(['V','g']):
 ax=axes[c]; pos=ax.get_position(); inset=fig.add_axes([pos.x0+.12,.05,.17,.16]); rows=[v for v in views if v['key'][3]==band]
 for target,small in [(ax,False),(inset,True)]:
  base(target); target.plot([0,1],[0,1],color='.65',lw=.7)
  for v in rows:
   iou,rec,pur=v['score']; marker='x' if v['edge'] else ('s' if iou<.90 else 'o'); color=COL[v['z']]
   target.plot(rec,pur,marker=marker,ms=(3 if small else 5) if marker!='s' else (4 if small else 6),mfc='none' if marker=='s' else color,mec=color,mew=.8,ls='none',zorder=5 if marker!='o' else 2,alpha=.8)
  target.set_xlim((-.02,1.025) if small else (.90,1.002));target.set_ylim((-.02,1.025) if small else (.90,1.002))
  if small:target.set_xticks([0,.5,1]);target.set_yticks([0,.5,1]);target.tick_params(labelsize=6,length=2,pad=1);target.set_title('All views',fontsize=6,pad=2)
 ax.set(xlabel=r'$C$',ylabel=r'$P$',title='$V$' if band=='V' else 'SDSS $g$');ax.title.set_fontsize(10);ax.text(.03,.96,f'({chr(97+c)})',transform=ax.transAxes,va='top',fontsize=10)
fig.legend(handles=legend,loc='upper center',bbox_to_anchor=(.5,1.04),ncol=6,frameon=False,columnspacing=.9,handlelength=1.3)
save(fig,'b_residual.pdf')
# Separate control figure. Stable redshift colours, styles distinguish equal-redshift halos.
rest=defaultdict(list); masses={}
for line in (ROOT/'supplement/rest_frame.tsv').read_text().splitlines():
 if line.startswith('#') or not line.strip():continue
 p=line.split()
 if p[1]=='P':masses[(p[0],int(p[3]))]=float(p[7])
 if p[1]=='I':rest[(p[0],int(p[3]),p[5],round(float(p[2]),3))].append([p[4],float(p[6]),float(p[7])])
fig,axes=plt.subplots(2,1,figsize=(3.45,4.5),sharex=True,sharey=True); styles=['-','--','-.',':']
for c,band in enumerate(['V','g']):
 ax=axes[c];base(ax)
 for k,((step,ns,bb,z),rs) in enumerate(sorted((k,v) for k,v in rest.items() if k[2]==band)):
  xx=sorted(set(row[1] for row in rs)); y=np.array([[row[2] for row in sorted(rs,key=lambda row:row[1]) if row[0]==los] for los in 'xyz'])
  mass=masses[(step,ns)]; ax.plot(xx,np.median(y,axis=0),color=COL[z],ls=styles[k],lw=1.4,label=f'$z={z:.3f}$, $M={mass/1e13:.3f}'+r'\times10^{13}$')
  ax.fill_between(xx,y.min(axis=0),y.max(axis=0),color=COL[z],alpha=.09)
 ax.set(xlim=(.01,.51),ylim=(0,1.035),ylabel=r'$\mathrm{IoU}$');ax.text(.03,.96,f'({chr(97+c)})',transform=ax.transAxes,va='top',fontsize=10)
 ax.text(.97,.65,'$V$' if band=='V' else 'SDSS $g$',transform=ax.transAxes,ha='right',va='top',fontsize=10)
axes[1].set_xlabel(r'$b_{\mathrm{link}}$')
axes[1].legend(loc='lower right',fontsize=8.5,frameon=False);fig.subplots_adjust(left=.17,right=.99,bottom=.10,top=.99,hspace=0)
for ax in axes:
 ax.tick_params(axis='x',bottom=True,top=True)
 ax.tick_params(axis='y',left=True,right=True)
axes[0].tick_params(axis='x',labelbottom=False)
save(fig,'b_rest.pdf')
print('Four figures redrawn from archived tables. Published endpoint flags retained.')
