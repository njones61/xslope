"""Regenerate composite_surface.png for docs/lem/overview.md.

Draws the Fredlund & Krahn weak-seam benchmark (verification corpus VP22): the
trial circle bottoms out five feet below the model base, so XSLOPE truncates it
and runs the surface along the base between the two crossings.

    PYTHONPATH=. python3 docs/lem/images/composite_surface.py
"""
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.text import Annotation

from xslope.fileio import load_slope_data
from xslope.plot import plot_inputs
from xslope.slice import build_composite_surface, generate_slices

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
XLSX = os.path.join(HERE, "..", "..", "verification", "files", "rocscience", "vp022a.xlsx")
OUT = os.path.join(HERE, "composite_surface.png")


def main():
    sd = load_slope_data(XLSX)
    circle = sd["circles"][0]
    ok, (df, fs) = generate_slices(sd, circle=circle, num_slices=40, composite=True)
    if not ok:
        raise RuntimeError(fs)

    # Drawn to be shown at 100 px per figure inch (half its 200 dpi pixel width),
    # where 12 pt text renders at about the page's own text size.
    plt.rcParams.update({"font.size": 12, "xtick.labelsize": 11, "ytick.labelsize": 11})
    fig = plt.figure(figsize=(8.4, 4.8))
    plot_inputs(sd, fig=fig, show_title=False, show_mesh=False)
    ax = fig.axes[0]

    # Drop plot_inputs' radius arrow and center marker — the center (120, 90) sits
    # far above the crop, so the arrow would read as a stray diagonal line.
    for t in list(ax.texts):
        if isinstance(t, Annotation) and t.arrow_patch is not None:
            t.remove()
    for ln in list(ax.lines):
        if len(ln.get_ydata()) and max(ln.get_ydata()) > 85:
            ln.remove()

    fx, fy = zip(*fs.coords)
    ax.plot(fx, fy, lw=2.8, color="crimson", zorder=7,
            label="composite failure surface")

    comp = build_composite_surface(sd, circle, min(fx), max(fx))
    for i, xc in enumerate(comp.crossings):
        ax.plot([xc], [comp.y(np.array([xc]))[0]], "o", ms=8, mfc="white",
                mec="crimson", mew=2.0, zorder=8,
                label="arc meets the base" if i == 0 else None)

    ax.annotate("the trial circle passes\nbelow the model base",
                xy=(120, 10.2), xytext=(96, -7), fontsize=12, color="0.3",
                ha="center", arrowprops=dict(arrowstyle="->", color="0.5", lw=1.0,
                                             shrinkA=3, shrinkB=6))
    ax.set_xlim(-5, 190)
    ax.set_ylim(-12, 66)
    ax.tick_params(labelsize=11)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), fontsize=11.5,
              frameon=False, ncol=3)
    ax.set_title("A circle deeper than the model base is truncated: it follows the arc down,\n"
                 "runs along the base between the two crossings, and climbs back out on the arc",
                 fontsize=12)
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
