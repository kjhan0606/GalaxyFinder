"""Match figure text to the 9 TeX pt body at the MNRAS 508 TeX pt width."""
from matplotlib.text import Text
BODY_TEX_PT = 9.0
TEXT_WIDTH_TEX_PT = 508.0

def save_paper_figure(fig, path, dpi=300, tight=True):
    # Matplotlib uses 72 points/inch, TeX uses 72.27. The ratio below
    # compensates for both units and the includegraphics width scaling.
    box = None
    for _ in range(5):
        fig.canvas.draw()
        if tight:
            box = fig.get_tightbbox(fig.canvas.get_renderer()).padded(.1)
            width_in = box.width
        else:
            width_in = fig.get_figwidth()
        source_size = BODY_TEX_PT * width_in * 72 / TEXT_WIDTH_TEX_PT
        for label in fig.findobj(Text):
            label.set_fontsize(source_size)
    fig.canvas.draw()
    if tight:
        box = fig.get_tightbbox(fig.canvas.get_renderer()).padded(.1)
    fig.savefig(path, dpi=dpi, bbox_inches=box)
