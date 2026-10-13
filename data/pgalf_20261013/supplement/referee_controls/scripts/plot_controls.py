import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
from pathlib import Path
import sys,numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(Path(__file__).parent))
from paper_figure_style import save_paper_figure
n=pd.read_csv(ROOT/'results/outer_cut.tsv',sep='\t').set_index('nstar')
p=pd.read_csv(ROOT/'results/population_controls.tsv',sep='\t')
a=p[p.control=='all_stars'].set_index('nstar');old=p[p.control=='age_ge_1Gyr'].set_index('nstar');joint=p[p.control=='radius_0p05dex_age_0p5Gyr'].set_index('nstar')
plt.rcParams.update({'font.size':9,'axes.labelsize':9,'xtick.labelsize':9,'ytick.labelsize':9})
fig,ax=plt.subplots(1,3,figsize=(7.1,2.8))
ax[0].scatter(100*n.removed_fraction,-100*n.actual_change,s=18,color='#286090',label='Linked stars')
ax[0].scatter(100*n.removed_fraction,-100*n.null_change,s=35,facecolors='none',edgecolors='#bf6532',marker='s',label='Outer truncation')
ax[0].set(xlabel='Stellar mass removed (%)',ylabel='Half-mass radius decrease (%)')
ax[0].legend(loc='upper left',bbox_to_anchor=(.02,.72),fontsize=9,frameon=False)
ax[1].scatter(100*a.age_fractional,100*old.loc[a.index].age_fractional,s=18,color='#286090')
lims=[-3,95];ax[1].plot(lims,lims,':',color='.5',lw=.7);ax[1].set(xlim=lims,ylim=lims,xlabel='Age difference, all stars (%)',ylabel='Age difference, ages >1 Gyr (%)')
x=a.logZ_delta;y=joint.loc[a.index].logZ_delta
lims=[min(x.min(),y.min())-.02,.01]
ax[2].plot(lims,lims,':',color='.5',lw=.7);ax[2].scatter(x,y,s=18,color='#286090');ax[2].set(xlim=lims,ylim=lims,xlabel=r'$\Delta\langle\log_{10} Z\rangle$ (dex)'+'\nRadius matched',ylabel=r'$\Delta\langle\log_{10} Z\rangle$ (dex)'+'\nRadius and age matched')
for panel,title in zip(ax,['(a) Mass\nand size','(b) Age cut','(c) Age matching']):
 panel.text(.04,.12 if panel is ax[2] else .97,title,transform=panel.transAxes,ha='left',va='top',bbox={'facecolor':'white','edgecolor':'none','alpha':.85,'pad':1})
 panel.tick_params(direction='in',top=True,right=True)
fig.tight_layout(w_pad=.8)
(ROOT/'figures').mkdir(exist_ok=True)
save_paper_figure(fig,ROOT/'figures/stellar_selection_controls.pdf',width_tex_pt=508)
