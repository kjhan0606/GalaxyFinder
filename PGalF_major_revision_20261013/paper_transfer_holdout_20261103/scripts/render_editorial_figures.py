#!/usr/bin/env python3
"""Regenerate current figures from saved measurements with journal typography."""
from pathlib import Path
import runpy
import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
_original_save = Figure.savefig
LABELS = {'Median light-weighted overlap': r'$\widetilde{\mathrm{IoU}}$', 'Median overlap': r'$\widetilde{\mathrm{IoU}}$', r'Linking parameter $b_{\rm link}$': r'$b_{\mathrm{link}}$', r'Distance threshold $\Delta$ (arcsec)': r'$\Delta$ (arcsec)'}
def save_journal_figure(fig, *args, **kwargs):
    for ax in fig.axes:
        for spine in ax.spines.values(): spine.set_linewidth(1.2)
        ax.tick_params(which='both', direction='in', top=True, bottom=True, left=True, right=True, width=1.0)
        ax.set_xlabel(LABELS.get(ax.get_xlabel(), ax.get_xlabel()))
        ax.set_ylabel(LABELS.get(ax.get_ylabel(), ax.get_ylabel()))
        for text in ax.texts: text.set_fontsize(9)
        legend = ax.get_legend()
        if legend:
            for text in legend.get_texts(): text.set_fontsize(9)
    return _original_save(fig, *args, **kwargs)
Figure.savefig = save_journal_figure
JOBS = [
 ('scripts/plot_revised_validation.py', []),
 ('scripts/plot_b_rest.py', []),
 ('supplement/tng300_growth_bias_20261106/plot_growth_r1r4_20261106.py', ['--force']),
 ('supplement/astra_review_response_20261106/plot_branch_growth_recovery.py', ['--force']),
 ('supplement/astra_review_response_20261106/recompute_corrected_psf_and_boundary.py', ['--summaries-only']),
 ('supplement/camels_image_checks/measurement_sources/plot_image_sensitivity.py', []),
 ('supplement/camels_multiz_20261030/plot_camels_multiz_transfer.py', [])]
if __name__ == '__main__':
    for name, args in JOBS:
        sys.argv = [str(ROOT/name), *args]
        runpy.run_path(str(ROOT/name), run_name='__main__')
        plt.close('all')
