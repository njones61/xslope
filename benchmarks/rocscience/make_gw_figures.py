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
from PIL import Image

from xslope.fileio import load_slope_data
from xslope.mesh import get_material_polygons, build_mesh_from_polygons
from xslope.seep import build_seep_data, run_seepage_analysis
from xslope.plot import plot_inputs
from xslope.plot_seep import plot_seep_solution

SRC = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'files', 'rocscience_gw')
OUT = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'verification', 'images')

# (stem, target_size or None for span/120, max_iter) — mesh sizes match the test tags
CASES = [
    ('gw002', 0.10, 400),
    ('gw003', 0.40, 400),
    ('gw004', 0.147, 400),
    ('gw006a', 1.0, 2000),
    ('gw008', 0.025, 2000),
    ('gw009a', None, 400),
    ('gw010', 0.25, 1500),
    ('gw011', 1.0, 2000),
    ('gw012', 1.0, 1500),
    ('gw013', 1.0, 1500),
]


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

    q = solution.get('flowrate')
    paths = []
    for which in ('inputs', 'solution'):
        fig = plt.figure(figsize=panel_size)
        if which == 'inputs':
            # mode='seep' draws the specified-head and exit-face BC lines
            plot_inputs(sd, fig=fig, mode='seep', mat_table=False,
                        show_title=True, title=f'{stem} — inputs')
        else:
            plot_seep_solution(seep_data, solution, fig=fig, show_title=True,
                               fill_contours=True, phreatic=True, flowlines=True,
                               mesh=False)
            if q is not None:
                fig.axes[0].set_title(f'{stem} — total head  (Q = {q:.3e})')
        p = os.path.join(OUT, f'_{stem}_{which}.png')
        # NO bbox_inches='tight': it crops each panel by its own ink extent, and the
        # two panels have very different ink (the solution carries a colorbar and a
        # legend). Cropping them independently and then scaling to a common height
        # shrinks the inputs panel to a sliver. Saving both at the same figsize/dpi
        # gives two identically-sized images that paste straight together.
        fig.savefig(p, dpi=dpi)
        plt.close(fig)
        paths.append(p)

    imgs = [Image.open(p) for p in paths]
    h = max(im.height for im in imgs)
    combo = Image.new('RGB', (sum(im.width for im in imgs), h), 'white')
    x = 0
    for im in imgs:
        combo.paste(im, (x, (h - im.height) // 2))   # vertically centred
        x += im.width
    out = os.path.join(OUT, f'{stem}.png')
    combo.save(out)
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
