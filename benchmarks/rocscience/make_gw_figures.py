"""Render figures for the Groundwater corpus (docs/verification/rocscience_groundwater.md):
one two-panel PNG per problem — the inputs (geometry, material zones, seepage boundary
conditions) on the left, and the solved seepage solution (total-head fill contours, the
phreatic surface, flow lines) on the right. Solves live, on the same pipeline as the
type=seep test tags.

Both panels come from XSLOPE's own plotting stack — plot_inputs(mode='seep') and
plot_seep_solution — NOT from hand-rolled matplotlib. That matters: those functions
already handle the material fills, the BC symbols, the phreatic surface, the flow net
and the colorbar, and they keep these figures looking like the ones a user gets from
the package. (An earlier version of this script drew its own triplot/tricontour and
drifted away from the package style.)

Run from the repo root:  python benchmarks/rocscience/make_gw_figures.py [gwNNN ...]
"""

import io
import os
import sys
import contextlib
import warnings

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

from xslope.fileio import load_slope_data
from xslope.mesh import get_material_polygons, build_mesh_from_polygons
from xslope.seep import build_seep_data, run_seepage_analysis
from xslope.plot import plot_inputs
from xslope.plot_seep import plot_seep_solution

SRC = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'verification', 'files', 'rocscience_gw')
OUT = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'verification', 'images')

# (stem, target_size or None for span/120, max_iter) — mesh sizes match the test tags
CASES = [
    ('gw001', 0.20, 400),
    ('gw002', 0.10, 400),
    ('gw003', 0.40, 400),
    ('gw004', 0.06, 2500),
    ('gw005', 0.5, 400),
    ('gw006a', 1.0, 2000),
    ('gw006b', 1.0, 2000),
    ('gw006c', 1.0, 2000),
    ('gw006d', 1.0, 2000),
    ('gw006e', 1.0, 2000),
    ('gw007', 0.04, 1000),
    ('gw008', 0.025, 2000),
    ('gw009a', None, 400),
    ('gw009b', None, 400),
    ('gw010', 0.25, 1500),
    ('gw011', 1.0, 2000),
    ('gw012', 1.0, 1500),
    ('gw013', 1.0, 1500),
]


#: The cushion between the section and the panel frame, as a fraction of the
#: larger domain dimension — one number in DATA units applied to both axes, so
#: under equal aspect the margin reads the same on all four sides. 3.5% is the
#: middle of the project figure-frame spec's 3–4% band, and it is what
#: ``plot_inputs(frame="content")`` uses by default.
CUSHION = 0.035
#: How many passes the figure-fitting loop takes, and how close in inches counts
#: as fitted. The furniture is laid out in points, so each pass removes most of
#: the remaining error and a handful of passes is plenty.
FIT_ROUNDS, FIT_TOL = 6, 0.01


def domain_frame(seep_data, cushion=CUSHION):
    """The axis limits both panels of a pair are framed to.

    Taken from the meshed DOMAIN — the section itself — rather than from each
    panel's own ink, so a water-level line drawn out over the reservoir or a
    hatched base cannot stretch one panel's frame past the other's. One cushion,
    in data units, on both axes.
    """
    nodes = np.asarray(seep_data['nodes'])
    x0, x1 = float(nodes[:, 0].min()), float(nodes[:, 0].max())
    y0, y1 = float(nodes[:, 1].min()), float(nodes[:, 1].max())
    pad = cushion * max(x1 - x0, y1 - y0)
    return (x0 - pad, x1 + pad), (y0 - pad, y1 + pad)


def axes_box_inches(fig, ax):
    """The axes box's size in inches, as drawn."""
    fig.canvas.draw()
    box = ax.get_position()
    fw, fh = fig.get_size_inches()
    return box.width * fw, box.height * fh


def fit_axes_box(fig, ax, width_in=None, height_in=None):
    """Resize the FIGURE until the axes box measures the requested inches.

    The panel's furniture — title, legend, colorbar, tick labels — is laid out in
    points, so its size in inches barely moves when the figure is resized and the
    axes box takes up the difference. Correcting the figure by the measured
    shortfall therefore converges in a few passes, and it needs no per-case
    constant: a wide-thin section ends in a wide-thin figure and a compact one in
    a compact figure. Passing ``None`` for a dimension leaves the figure's own
    extent in that direction alone.
    """
    for _ in range(FIT_ROUNDS):
        have_w, have_h = axes_box_inches(fig, ax)
        fw, fh = fig.get_size_inches()
        dw = 0.0 if width_in is None else width_in - have_w
        dh = 0.0 if height_in is None else height_in - have_h
        if abs(dw) < FIT_TOL and abs(dh) < FIT_TOL:
            break
        fig.set_size_inches(max(fw + dw, 1.0), max(fh + dh, 1.0))
    return fig


#: The gap between the x-axis decoration and the top of the legend, as a figure
#: fraction — the same one ``xslope.plot._legend_below`` leaves.
LEGEND_GAP = 0.035


def repin_legend(fig, ax, gap=LEGEND_GAP):
    """Re-hang the legend under the axes box after the panel has been re-framed.

    The plotting functions pin the legend just below the axes box as they draw
    it. Re-framing the panel afterwards moves that box, so the legend has to
    follow or it strands at the old position and leaves a white band across the
    figure. The rule is the plotting module's own: the legend's top sits one gap
    below the lowest x-axis decoration.
    """
    leg = ax.get_legend()
    if leg is None:
        return
    fig.canvas.draw()
    ax.apply_aspect()
    renderer = fig.canvas.get_renderer()
    spine = ax.get_window_extent().y0
    low = spine
    for artist in list(ax.get_xticklabels()) + [ax.xaxis.get_label()]:
        if not artist.get_text():
            continue
        low = min(low, artist.get_window_extent(renderer).y0)
    leg.set_bbox_to_anchor((0.5, low / fig.bbox.height - gap),
                           transform=fig.transFigure)


def frame_panel(fig, ax, xlim, ylim, width_in=None):
    """Put one panel on the pair's shared frame, at the section's true shape.

    The frame is the pair's; the axes box is fitted to the section's own aspect
    so the panel is no taller than the section needs. ``width_in`` asks for a
    particular drawn width — the second panel of a pair is fitted to the first
    panel's, which is how the two come out at one scale even though only one of
    them carries a colorbar. The first panel is left at the figure width it was
    created with, which is the width its legend was laid out to fit.
    """
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect('equal', adjustable='box')
    aspect = (ylim[1] - ylim[0]) / (xlim[1] - xlim[0])
    want_w = width_in if width_in is not None else axes_box_inches(fig, ax)[0]
    fit_axes_box(fig, ax, width_in, want_w * aspect)
    repin_legend(fig, ax)
    return axes_box_inches(fig, ax)[0]


def trim_border(img, margin_frac=CUSHION):
    """Crop the pasted pair's outer whitespace back to a uniform margin.

    Each panel is fitted to its own content, but the two end up different
    heights and the taller one sets the pasted image's height, so the pair
    carries a band of empty page above and below. This is the LAST step and it
    touches only the outside border — no gutter between the panels, no panel's
    own margins, so the two panels stay in the register the paste put them in.
    """
    grey = img.convert('L')
    ink = np.asarray(grey) < 250
    rows, cols = np.where(ink.any(axis=1))[0], np.where(ink.any(axis=0))[0]
    if not rows.size or not cols.size:
        return img
    margin = int(round(margin_frac * max(img.width, img.height)))
    box = (max(int(cols[0]) - margin, 0), max(int(rows[0]) - margin, 0),
           min(int(cols[-1]) + 1 + margin, img.width),
           min(int(rows[-1]) + 1 + margin, img.height))
    return img.crop(box)


def make_figure(stem, target_size, max_iter, panel_size=(8.0, 5.0), dpi=150):
    sd = load_slope_data(os.path.join(SRC, f'{stem}.xlsx'))
    polygons = get_material_polygons(sd)
    if target_size is None:
        xs = [x for x, _ in sd['ground_surface'].coords]
        target_size = (max(xs) - min(xs)) / 120
    mesh = build_mesh_from_polygons(polygons, target_size, 'tri3')
    seep_data = build_seep_data(mesh, sd)
    with contextlib.redirect_stdout(io.StringIO()):
        solution = run_seepage_analysis(seep_data, tol=1e-5, max_iter=max_iter)

    # One frame per case, derived once and given to both panels: identical axis
    # ranges, identical scale, so the pair reads as one figure.
    xlim, ylim = domain_frame(seep_data)

    q = solution.get('flowrate')
    paths = []
    drawn_width = None          # set by the inputs panel, matched by the solution
    for which in ('inputs', 'solution'):
        fig = plt.figure(figsize=panel_size)
        if which == 'inputs':
            # mode='seep' draws the specified-head and exit-face BC lines
            plot_inputs(sd, fig=fig, mode='seep', mat_table=False,
                        show_title=True, title=f'{stem} — inputs',
                        frame='content', pad_frac=CUSHION)
        else:
            # This script sets its own solution-panel title (below), so suppress the
            # built-in one: the modern steady title is a two-line variant (main line +
            # a smaller flowrate sub-line drawn as a separate annotation), and a bare
            # set_title() override would replace only the main line and leave the
            # flowrate annotation dangling under our title.
            plot_seep_solution(seep_data, solution, fig=fig, show_title=False,
                               fill_contours=True, phreatic=True, flowlines=True,
                               mesh=False, pad_frac=CUSHION)
            head_title = f'{stem} — total head'
            if q is not None:
                head_title += f'  (Q = {q:.3e})'
            fig.axes[0].set_title(head_title)
        drawn_width = frame_panel(fig, fig.axes[0], xlim, ylim, drawn_width)
        p = os.path.join(OUT, f'_{stem}_{which}.png')
        # NO bbox_inches='tight': it would crop each panel by its own ink extent,
        # and the two carry very different ink (the solution has a colorbar and a
        # legend). The frame above already sizes each figure to its content, and
        # both panels' axes boxes are the same size in inches, so the sections
        # paste together at one scale.
        fig.savefig(p, dpi=dpi)
        plt.close(fig)
        paths.append(p)

    imgs = [Image.open(p) for p in paths]
    h = max(im.height for im in imgs)
    combo = Image.new('RGB', (sum(im.width for im in imgs), h), 'white')
    x = 0
    for im in imgs:
        combo.paste(im, (x, (h - im.height) // 2))   # vertically centered
        x += im.width
    out = os.path.join(OUT, f'{stem}.png')
    trim_border(combo).save(out)
    for p in paths:
        os.remove(p)
    return out


if __name__ == '__main__':
    only = set(sys.argv[1:])
    for stem, ts, mi in CASES:
        if only and stem not in only:
            continue
        try:
            make_figure(stem, ts, mi)
            print('ok  ', stem, flush=True)
        except Exception as e:                      # noqa: BLE001
            print('FAIL', stem, ':', repr(e)[:140], flush=True)
