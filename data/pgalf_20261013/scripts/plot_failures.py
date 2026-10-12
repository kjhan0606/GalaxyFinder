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
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42})
fig,axes=plt.subplots(2,3,figsize=(7.05,4.7));fig.subplots_adjust(left=.065,right=.99,bottom=.16,top=.88,wspace=.28,hspace=.38)
for i,(ax,(half,cx,cy)) in enumerate(zip(axes.ravel(),PANELS)):
 letter=chr(97+i);rgb=np.flipud(np.array(Image.open(ROOT/f'supplement/failure_panels/{letter}.png'))[:,:,:3]);ax.imshow(rgb,origin='lower',extent=[-half,half,-half,half],interpolation='nearest')
 window=60 if i in [1,2] else half
 # Recover only the existing orange seed marker location from the diagnostic raster.
 mask=(rgb[:,:,0]>230)&(rgb[:,:,1]>165)&(rgb[:,:,1]<220)&(rgb[:,:,2]<70)
 ys,xs=np.nonzero(mask)
 if len(xs):
  sx=(np.median(xs)+.5)/rgb.shape[1]*2*half-half;sy=(np.median(ys)+.5)/rgb.shape[0]*2*half-half
  ax.plot(sx,sy,'x',ms=3.5,mew=.8,color='#FFC107',zorder=7)
 # Locate the saved coral boundary, then outline and label the selected region.
 coral=(rgb[:,:,0]>190)&(rgb[:,:,1]>60)&(rgb[:,:,1]<170)&(rgb[:,:,2]>40)&(rgb[:,:,2]<170)
 yy,xx=np.nonzero(coral)
 if i in [0,1] and len(xx):
  x0=(np.min(xx)+.5)/rgb.shape[1]*2*half-half;x1=(np.max(xx)+.5)/rgb.shape[1]*2*half-half
  y0=(np.min(yy)+.5)/rgb.shape[0]*2*half-half;y1=(np.max(yy)+.5)/rgb.shape[0]*2*half-half
  tx=(x0+x1)/2;ty=(y0+y1)/2
  ax.add_patch(plt.Rectangle((x0-.10,y0-.10),max(x1-x0,.2)+.2,max(y1-y0,.2)+.2,fill=False,ec='#e04126',lw=1.0,zorder=10))
  ax.annotate('Selected region',xy=(tx,ty),xytext=(.06,.12),textcoords='axes fraction',fontsize=6.5,color='#a01b12',bbox=dict(fc='white',ec='none',alpha=.9,pad=1),arrowprops=dict(arrowstyle='->',color='#e04126',lw=.9),zorder=11)
 if i==1:
  ax.annotate('Selected region\nnear centre',xy=(0,0),xytext=(-53,40),fontsize=6.5,color='#a01b12',bbox=dict(fc='white',ec='none',alpha=.9,pad=1),arrowprops=dict(arrowstyle='->',color='#e04126',lw=.9),zorder=11)
 ax.plot(0,0,'+',ms=4,mew=.8,color='#43CAFF',zorder=8)
 if abs(cx)<=window and abs(cy)<=window:
  ax.plot(cx,cy,'o',ms=4.5,mfc='none',mec='#D08CFF',mew=.8,zorder=6)
 else:
  scale=.92*window/max(abs(cx),abs(cy));tx,ty=cx*scale,cy*scale
  ax.annotate('',xy=(tx,ty),xytext=(.67*tx,.67*ty),arrowprops=dict(arrowstyle='-|>',color='#D08CFF',lw=1.5),zorder=7)
 ax.set(xlim=(-window,window),ylim=(-window,window));ax.tick_params(labelsize=8,length=3,pad=2)
 ax.text(.035,.95,f'({letter})',transform=ax.transAxes,va='top',bbox=dict(fc='white',ec='none',alpha=.9,pad=1))
 row=scores[KEYS[i]]; overlap=float(row['iou_best'])
 value=f'{overlap:.4f}' if overlap>=1e-4 else r'<10^{-4}'
 ax.set_title(f"${row['band']}$, ${row['los']}$ view, $b={float(row['b_best']):.2f}$\nOverlap ${value}$",fontsize=7,pad=3)
 if i>=3:ax.set_xlabel('Offset (arcsec)')
 if i%3==0:ax.set_ylabel('Offset (arcsec)')
handles=[Line2D([],[],marker='+',ms=8,mew=1.4,color='#43CAFF',ls='none',label='Image centre'),Line2D([],[],marker='x',ms=7,mew=1.4,color='#FFC107',ls='none',label='Density seed'),Line2D([],[],marker='o',ms=9,mfc='none',mec='#D08CFF',ls='none',label='Watershed core')]
fig.legend(handles=handles,loc='upper center',ncol=3,frameon=False,bbox_to_anchor=(.52,1.005))
cax=fig.add_axes([.33,.05,.42,.025]);cb=fig.colorbar(ScalarMappable(norm=Normalize(20,28),cmap='gray_r'),cax=cax,orientation='horizontal');cb.set_label(r'Surface brightness (mag arcsec$^{-2}$)',labelpad=2);cb.ax.tick_params(labelsize=7,length=2,pad=1)
fig.savefig(ROOT/'figures/fail_maps.pdf',bbox_inches='tight',pad_inches=.03)
