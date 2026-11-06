#!/usr/bin/env python3
"""Show fixed-calibration overlap and luminosity-error sensitivities."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
data=pd.read_csv(ROOT/"sensitivity_summary.csv")
conditions=["baseline","limit_brighter","limit_fainter","pixel_coarse","no_psf_wing",
            "sky_residual_29","noise_5sigma","noise_3sigma","noise_1sigma"]
labels=["Base",r"$-0.5$",r"$+0.5$",r"$0.4''$","No wing","Sky 29",r"$5\sigma$",r"$3\sigma$",r"$1\sigma$"]
styles=[("FoF","0.1",r"$b_{\rm link}=0.10$","#1F4E79","o"),
        ("FoF","0.2",r"$b_{\rm link}=0.20$","#D4652F","s"),
        ("aperture30pkpc","aperture","30 pkpc",".45","^")]
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"pdf.fonttype":42})
fig,axes=plt.subplots(2,2,figsize=(7.05,4.7),sharex=True,sharey="row")
for col,band in enumerate(("V","g")):
 for row,metric in enumerate(("mean_host_iou","p90_absolute_logL_error")):
    ax=axes[row,col]
    for method,b,label,color,marker in styles:
        group=data[data.band.eq(band)&data.method.eq(method)&data.b_star.astype(str).eq(b)].set_index("condition")
        if set(group.index)!=set(conditions):raise ValueError("Incomplete figure conditions")
        ax.plot(range(len(conditions)),group.loc[conditions,metric],marker=marker,ms=3.5,lw=1.,color=color,label=label)
    ax.text(.03,.96,f"({chr(97+2*row+col)})",transform=ax.transAxes,va="top")
    ax.text(.97,.96,"$V$" if band=="V" else "SDSS $g$",transform=ax.transAxes,ha="right",va="top")
    ax.grid(axis="y",alpha=.18,lw=.5)
    ax.tick_params(direction="in",top=True,right=True)
    ax.set_xticks(range(len(conditions)),labels,rotation=42,ha="right")
    ax.set_xlim(-.5,len(conditions)-.5)
    if col==0:ax.set_ylabel(r"$\langle\mathrm{IoU}\rangle_{\rm host}$" if row==0 else r"$P_{90}(|\Delta\log_{10}L|)$")
minimum=float(data.mean_host_iou.min())
axes[0,0].set_ylim(max(0.,np.floor(minimum*10)/10-.05),1.03)
maximum=float(data.p90_absolute_logL_error.max())
axes[1,0].set_ylim(0,np.ceil(maximum/.05)*.05+.02)
handles,legend_labels=axes[0,0].get_legend_handles_labels()
axes[1,0].legend(handles,legend_labels,loc="upper center",bbox_to_anchor=(.58,.995),
                 ncol=2,fontsize=9,framealpha=.92,borderpad=.2,
                 columnspacing=.65,handlelength=1.0,labelspacing=.2)
fig.subplots_adjust(left=.115,right=.99,bottom=.16,top=.99,wspace=0,hspace=0)
out=ROOT/"fixed_b_image_sensitivity.pdf"
fig.savefig(out,bbox_inches="tight")
fig.savefig(out.with_suffix(".png"),dpi=170,bbox_inches="tight")
print(out)
