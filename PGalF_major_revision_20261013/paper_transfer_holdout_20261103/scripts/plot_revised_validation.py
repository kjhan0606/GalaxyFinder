#!/usr/bin/env python3
"""Regenerate primary figures from the corrected-distance Darwin and centred CAMELS tables."""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[1]; SUP=ROOT/'supplement'; OUT=ROOT/'figures'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'xtick.direction':'in','ytick.direction':'in','pdf.fonttype':42})
colors={.190:'#1F4E79',.227:'#148F9A',.500:'#5B8C3E',1.002:'#E0A100',1.990:'#D4652F',3.002:'#8E3A73'}
d=pd.read_csv(SUP/'recomputed_curves.tsv',sep='\t'); d=d[d.variant=='consistent_distance']; d_all=d.copy()
sel=pd.read_csv(SUP/'selected_parameters.tsv',sep='\t'); sel=sel[sel.variant=='consistent_distance']
meta=pd.read_csv(SUP/'halo_metadata.tsv',sep='\t',comment='#')[['step','nstar','z']]
sel=sel.merge(meta,on=['step','nstar'],how='left')
# Curves: only entries retained by the stated interior-maximum criterion.
valid=sel[sel.edge==0].groupby(['step','nstar','band']).size(); valid=valid[valid==3].index
keys=set((int(a),int(b),c) for a,b,c in valid)
d=d[d.apply(lambda r:(int(r.step),int(r.nstar),r.band) in keys,axis=1)]
fig,axs=plt.subplots(2,2,figsize=(7.05,5.2),sharex=True,sharey=True)
for j,band in enumerate(['V','g']):
 for i,metric in enumerate(['iou','completeness']):
  ax=axs[i,j]
  for z,g in d[d.band==band].groupby('z'):
   # equal-weight entries: first median views per entry at each b, then median entries
   med=g.groupby(['step','nstar','b'])[metric].median().groupby('b').median()
   ax.plot(med.index,med.values,color=colors[round(z,3)],lw=1.4,label=f'{z:.3f}')
  ax.set_ylim(.70,1.00);ax.set_yticks(np.arange(.70,1.001,.05));ax.grid(alpha=.18)
  if j==0: ax.set_ylabel(r'$\mathrm{IoU}$' if metric=='iou' else r'$C$')
  ax.text(.97,.95,'$V$' if band=='V' else 'SDSS $g$',transform=ax.transAxes,ha='right',va='top',fontsize=10,weight='bold');ax.text(.03,.95,f'({chr(97+i*2+j)})',transform=ax.transAxes,va='top',fontsize=10)
  if i==1: ax.set_xlabel(r'$b_{\mathrm{link}}$')
axs[0,0].legend(title='Output $z$',title_fontsize=8.5,fontsize=8.5,
                loc='lower left',ncol=2,frameon=True,framealpha=.88,
                borderpad=.25,labelspacing=.2,handlelength=1.1)
for ax in axs.flat:
 ax.set_xlim(.02,.205);ax.set_xticks(np.arange(.05,.201,.05))
for i in range(2):
 for j in range(2):
  axs[i,j].tick_params(axis='x',bottom=True,top=True,labelbottom=(i==1),labeltop=False)
  axs[i,j].tick_params(axis='y',left=True,right=True,labelleft=(j==0),labelright=(j==1))
fig.canvas.draw()
for ax,edge in ((axs[0,0],.70),(axs[0,1],.70),(axs[1,0],1.0),(axs[1,1],1.0)):
 for tick in ax.yaxis.get_major_ticks():
  if np.isclose(tick.get_loc(),edge):
   tick.label1.set_visible(False);tick.label2.set_visible(False)
fig.subplots_adjust(left=.15,right=.99,bottom=.11,top=.985,wspace=0.0,hspace=0.0);fig.savefig(OUT/'b_curves.pdf',bbox_inches='tight');plt.close(fig)
# Selected b by redshift; open square indicates selected maximum below 0.9.
fig,axs=plt.subplots(2,1,figsize=(3.45,4.3),sharex=True,sharey=True)
zorder=[.190,.227,.500,1.002,1.990,3.002]
for j,band in enumerate(['V','g']):
 ax=axs[j];q=sel[sel.band==band]
 for x,z in enumerate(zorder):
  gg=q[np.isclose(q.z.round(3),z)]
  hs=gg.groupby(['step','nstar']).agg(b=('b_best','median'),iou=('iou_best','median'),edge=('edge','max'))
  for k,(_,r) in enumerate(hs.iterrows()):
   xx=x+np.linspace(-.24,.24,len(hs))[k] if len(hs)>1 else x
   ax.plot(xx,r.b,marker='x' if r.edge else ('s' if r.iou<.9 else 'o'),ms=5 if r.iou<.9 else 3.5,mfc='white' if r.iou<.9 else colors[z],mec=colors[z],color=colors[z],ls='none')
 ax.set_xticks(range(6),[f'{z:.3f}' for z in zorder]);ax.set_ylim(.015,.205);ax.text(.97,.95,'$V$' if band=='V' else 'SDSS $g$',transform=ax.transAxes,ha='right',va='top',fontsize=10,weight='bold');ax.grid(axis='y',alpha=.18);ax.text(.03,.95,f'({chr(97+j)})',transform=ax.transAxes,va='top',fontsize=10)
 ax.tick_params(axis='x',bottom=True,top=True,labelbottom=(j==1),labeltop=False)
 ax.tick_params(axis='y',left=True,right=True)
for ax in axs: ax.set_ylabel(r'$b_{\mathrm{link}}$')
axs[1].set_xlabel(r'$z$')
axs[1].tick_params(axis='x',labelsize=7)
fig.subplots_adjust(left=.18,right=.99,bottom=.11,top=.99,hspace=0.0);fig.savefig(OUT/'b_match.pdf',bbox_inches='tight');plt.close(fig)
# Completeness and precision at each corrected-distance view optimum.
lookup=d_all.set_index(['step','nstar','los','band','b'])[['completeness','precision']]
fig,axs=plt.subplots(1,2,figsize=(7.05,3.0),sharey=True)
for j,band in enumerate(['V','g']):
 ax=axs[j];q=sel[sel.band==band]
 points=[]
 for _,r in q.iterrows():
  try: rec,pre=lookup.loc[(r.step,r.nstar,r.los,band,r.b_best)].values
  except KeyError: continue
  points.append((rec,pre,r))
 for rec,pre,r in points:
  marker='x' if r.edge else ('s' if r.iou_best<.9 else 'o')
  fill='none' if r.iou_best<.9 else colors[round(r.z,3)]
  color=colors[round(r.z,3)]
  if marker=='x': ax.scatter(rec,pre,s=5.5,marker=marker,color=color,linewidths=.65,alpha=.78)
  else: ax.scatter(rec,pre,s=5.5,marker=marker,facecolors=fill,edgecolors=color,linewidths=.6,alpha=.78)
 ax.plot([.6,1],[.6,1],color='.6',lw=.7);ax.set(xlim=(.6,1),ylim=(.6,1),xlabel=r'$C$');ax.text(.12,.96,'$V$' if band=='V' else 'SDSS $g$',transform=ax.transAxes,ha='left',va='top',fontsize=10,weight='bold');ax.grid(alpha=.18);ax.text(.03,.95,f'({chr(97+j)})',transform=ax.transAxes,va='top',fontsize=10)
 if j==0: ax.set_ylabel(r'$P$')
 ax.tick_params(axis='x',bottom=True,top=True)
 ax.tick_params(axis='y',left=True,right=True,labelleft=(j==0),labelright=(j==1))
ticks=np.arange(.6,1.001,.1)
axs[0].set_xticks(ticks,['0.6','0.7','0.8','0.9',''])
axs[1].set_xticks(ticks,['','0.7','0.8','0.9','1.0'])
axs[0].set_yticks(ticks)
z_handles=[Line2D([],[],marker='o',ls='none',color=colors[z],label=f'$z={z:.3f}$')
           for z in sorted(colors)]
axs[0].legend(handles=z_handles,loc='lower left',ncol=2,title='$z$',
              title_fontsize=8,fontsize=7.5,frameon=True,framealpha=.88,
              borderpad=.25,columnspacing=.65,handlelength=1.0,labelspacing=.15)
fig.subplots_adjust(left=.10,right=.995,bottom=.17,top=.98,wspace=0.0);fig.savefig(OUT/'b_residual.pdf',bbox_inches='tight');plt.close(fig)
# CAMELS external transfer distributions, with all hosts shown.
h=pd.read_csv(SUP/'camels_host_metrics.csv')
h=h.copy()
# FoF b values are stored in b_star; aperture has no linking parameter.
methods=[('FoF',0.10),('FoF',0.20),('aperture30pkpc',None)]
fig,axs=plt.subplots(3,2,figsize=(7.05,6.0),sharex='col')
metrics=[('median_iou',r'$\mathrm{IoU}_{\rm med}$'),('median_abs_log10_luminosity_error',r'$|\Delta\log_{10}L|$'),('median_abs_log10_R50_error',r'$|\Delta\log_{10}R_{50}|$')]
for col,band in enumerate(['V','g']):
 for row,(metric,label) in enumerate(metrics):
  ax=axs[row,col]; vals=[]; labs=[]
  for m,bval in methods:
   mask=(h.band==band)&(h.method==m)
   if bval is not None: mask &= np.isclose(pd.to_numeric(h.b_star,errors='coerce'),bval)
   a=h[mask][metric].dropna().values
   vals.append(a);labs.append('30 pkpc' if bval is None else fr'$b_{{\mathrm{{link}}}}={bval:.2f}$')
  ax.boxplot(vals,tick_labels=labs,showfliers=False,widths=.6)
  for i,a in enumerate(vals,1): ax.scatter(np.full(len(a),i)+np.linspace(-.08,.08,len(a)),a,s=7,alpha=.42,color='#1F4E79' if band=='V' else '#D4652F')
  ax.set_ylabel(label);ax.grid(axis='y',alpha=.18);ax.tick_params(axis='x',labelrotation=18,labelbottom=(row==2))
  ax.text(.03,.95,f'({chr(97+row*2+col)})',transform=ax.transAxes,va='top',fontsize=10)
  if row==0:
   ax.set_ylim(.5,1.05);ax.set_yticks(np.arange(.5,1.01,.1))
   ax.text(.97,.95,'$V$' if band=='V' else 'SDSS $g$',transform=ax.transAxes,
           ha='right',va='top',fontsize=9.5,weight='bold',
           bbox=dict(facecolor='white',edgecolor='none',alpha=.88,pad=1.2))
  ax.tick_params(axis='x',bottom=True,top=True,labelbottom=(row==2),labeltop=False)
  ax.tick_params(axis='y',left=True,right=True,labelleft=(col==0),labelright=(col==1))
fig.canvas.draw()
for row in range(2):
 for col in range(2):
  for tick in axs[row,col].yaxis.get_major_ticks():
   if np.isclose(tick.get_loc(),min(axs[row,col].get_yticks())):
    tick.label1.set_visible(False);tick.label2.set_visible(False)
  for tick in axs[row+1,col].yaxis.get_major_ticks():
   if np.isclose(tick.get_loc(),max(axs[row+1,col].get_yticks())):
    tick.label1.set_visible(False);tick.label2.set_visible(False)
for ax in axs[:,1]: ax.yaxis.set_label_position('right')
fig.subplots_adjust(left=.13,right=.87,bottom=.17,top=.95,wspace=0.0,hspace=0.0);fig.savefig(OUT/'validation_transfer.pdf',bbox_inches='tight');plt.close(fig)
# Artificial-redshift overlap across all views and high-z detection accounting.
art=pd.read_csv(SUP/'artificial_redshift_summary.csv');art=art[art.flux_mode=='observer_full'];av=art[art.band=='V'].sort_values('z_obs');ag=art[art.band=='g'].sort_values('z_obs');z=av.z_obs.tolist();v=av.mean_iou_all_views.tolist();g=ag.mean_iou_all_views.tolist()
fig,ax=plt.subplots(figsize=(3.45,2.75))
ax.plot(z,v,'o-',label='$V$',color='#1F4E79')
ax.plot(z,g,'s-',label='SDSS $g$',color='#D4652F')
ax.set(xlabel=r'$z_{\rm obs}$',ylabel=r'$\langle\mathrm{IoU}\rangle_{\rm all}$',ylim=(.60,.98))
ax.grid(alpha=.2);ax.legend(frameon=False,fontsize=9,loc='lower left')
ax.tick_params(direction='in',top=True,right=True)
fig.tight_layout();fig.savefig(OUT/'artificial_redshift.pdf',bbox_inches='tight');plt.close(fig)
