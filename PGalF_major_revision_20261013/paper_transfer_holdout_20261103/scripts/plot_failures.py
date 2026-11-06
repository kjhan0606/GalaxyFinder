#!/usr/bin/env python3
"""Rearrange saved diagnostic images and add readable coordinate markers.
No new surface-brightness measurements. Input raster colours are preserved.
The generating save_figure function uses gray_r normalized from 20 to 28.
"""
from pathlib import Path
import csv
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
ROOT=Path(__file__).resolve().parents[1]
PANELS=[(3.,2.4266,3.6011),(180.,25.3363,-21.6020),(180.,19.8194,25.3363),(8*np.sqrt(3.40/np.pi),.1215,.8936),(8*np.sqrt(2.00/np.pi),.8936,-1.2479),(8*np.sqrt(3.64/np.pi),-1.2479,.1215)]
KEYS=[('00475',10033180,'y','g'),('00967',17494810,'x','V'),('00967',17494810,'z','g'),('00178',7062464,'x','V'),('00178',7062464,'y','V'),('00178',7062464,'z','g')]
with (ROOT/'supplement/boundary_ring.tsv').open() as fh:
 scores={(r['step'],int(r['nstar']),r['los'],r['band']):r for r in csv.DictReader(fh,delimiter='\t')}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42})
fig,axes=plt.subplots(2,3,figsize=(7.05,5.1));fig.subplots_adjust(left=.005,right=.995,bottom=.08,top=.992,wspace=0,hspace=0)
for i,(ax,(half,cx,cy)) in enumerate(zip(axes.ravel(),PANELS)):
 letter=chr(97+i);rgb=np.flipud(np.array(Image.open(ROOT/f'supplement/failure_panels/{letter}.png'))[:,:,:3]);ax.imshow(rgb,origin='lower',extent=[-half,half,-half,half],interpolation='nearest')
 window=60 if i in [1,2] else half
 # Recover only the existing orange seed marker location from the diagnostic raster.
 mask=(rgb[:,:,0]>230)&(rgb[:,:,1]>165)&(rgb[:,:,1]<220)&(rgb[:,:,2]<70)
 ys,xs=np.nonzero(mask)
 if len(xs):
  sx=(np.median(xs)+.5)/rgb.shape[1]*2*half-half;sy=(np.median(ys)+.5)/rgb.shape[0]*2*half-half
  ax.plot(sx,sy,'x',ms=3.5,mew=.8,color='#FFC107',zorder=7)
 ax.plot(0,0,'+',ms=4,mew=.8,color='#43CAFF',zorder=8)
 if abs(cx)<=window and abs(cy)<=window:
  ax.plot(cx,cy,'o',ms=4.5,mfc='none',mec='#D08CFF',mew=.8,zorder=6)
 else:
  scale=.92*window/max(abs(cx),abs(cy));tx,ty=cx*scale,cy*scale
  ax.annotate('',xy=(tx,ty),xytext=(.67*tx,.67*ty),arrowprops=dict(arrowstyle='-|>',color='#D08CFF',lw=1.5),zorder=7)
 ax.set(xlim=(-window,window),ylim=(-window,window));ax.set_axis_off()
 row=scores[KEYS[i]]; overlap=float(row['iou_best'])
 value=f'{overlap:.4f}' if overlap>=1e-4 else r'<10^{-4}'
 ax.text(.035,.97,f"({letter})  ${row['band']}$, ${row['los']}$ view\n$b_{{\\mathrm{{link}}}}={float(row['b_best']):.2f}$, IoU ${value}$",transform=ax.transAxes,va='top',fontsize=10,bbox=dict(fc='white',ec='none',alpha=.85,pad=1.5))
 length=1 if i==0 else (20 if i in [1,2] else 5)
 x0=.94-length/(2*window)
 ax.plot([x0,.94],[.07,.07],transform=ax.transAxes,color='white',lw=2.5,zorder=10)
 ax.text((x0+.94)/2,.09,fr'${length}\,\mathrm{{arcsec}}$',transform=ax.transAxes,ha='center',va='bottom',fontsize=9,bbox=dict(fc='white',ec='none',alpha=.85,pad=1),zorder=10)
 # Thin white separators preserve a visible boundary without panel gaps.
 ax.plot([0,1,1,0,0],[0,0,1,1,0],transform=ax.transAxes,color='white',lw=.7,zorder=20)
handles=[Line2D([],[],marker='+',ms=8,mew=1.4,color='#43CAFF',ls='none',label='Image centre'),Line2D([],[],marker='x',ms=7,mew=1.4,color='#FFC107',ls='none',label='Density seed'),Line2D([],[],marker='o',ms=9,mfc='none',mec='#D08CFF',ls='none',label='Watershed core')]
axes[1,0].legend(handles=handles,loc='lower left',bbox_to_anchor=(.01,.16),fontsize=9,framealpha=.9,handlelength=1.1,borderpad=.25,labelspacing=.25)
cax=fig.add_axes([.30,.035,.40,.018]);cb=fig.colorbar(ScalarMappable(norm=Normalize(20,28),cmap='gray_r'),cax=cax,orientation='horizontal');cb.set_label(r'$\mu\,[\mathrm{mag\,arcsec}^{-2}]$',labelpad=2);cb.ax.tick_params(labelsize=8,length=2,pad=1,direction='in',top=True,bottom=True)
fig.savefig(ROOT/'figures/fail_maps.pdf',bbox_inches='tight',pad_inches=.03)
