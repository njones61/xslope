"""Render figures for the Groundwater corpus (docs/verification/rocscience_groundwater.md):
one two-panel PNG per problem — the mesh with boundary conditions on the left, solved
total-head contours with the phreatic surface on the right. Solves live (same pipeline
as the type=seep tags).

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
import matplotlib.tri as mtri
import numpy as np
from PIL import Image

from xslope.fileio import load_slope_data
from xslope.mesh import get_material_polygons, build_mesh_from_polygons
from xslope.seep import build_seep_data, run_seepage_analysis

SRC = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'files', 'rocscience_gw')
OUT = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'verification', 'images')

# (stem, target_size or None for span/120, max_iter)
CASES = [
    ('gw002', 0.10, 400),
    ('gw003', 0.40, 400),
    ('gw004', 0.147, 400),
    ('gw009a', None, 400),
    ('gw010', 0.25, 1500),
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
        sol = run_seepage_analysis(seep_data, tol=1e-5, max_iter=max_iter)

    nodes = seep_data['nodes']
    elems = np.asarray(mesh['elements'])[:, :3]
    triang = mtri.Triangulation(nodes[:, 0], nodes[:, 1], elems)
    head = np.asarray(sol['head'])
    u = head - nodes[:, 1]          # pressure head, for the phreatic contour

    paths = []
    for which in ('mesh', 'solution'):
        fig, ax = plt.subplots(figsize=panel_size)
        if which == 'mesh':
            ax.triplot(triang, color='0.6', lw=0.3)
            bc = sd.get('seepage_bc') or {}
            for spec in bc.get('specified_heads', []):
                xs_, ys_ = zip(*spec['coords'])
                ax.plot(xs_, ys_, color='b', lw=2.5,
                        label=f"head = {spec['head']:g}")
            if bc.get('exit_face'):
                xs_, ys_ = zip(*bc['exit_face'])
                ax.plot(xs_, ys_, color='r', lw=2.5, label='exit face')
            ax.set_title(f'{stem} — mesh and boundary conditions')
            ax.legend(loc='best', fontsize=8)
        else:
            levels = np.linspace(head.min(), head.max(), 12)
            cs = ax.tricontour(triang, head, levels=levels, colors='k',
                               linewidths=0.6)
            ax.clabel(cs, fontsize=6, fmt='%.2f')
            try:
                ax.tricontour(triang, u, levels=[0.0], colors='b', linewidths=2.0)
            except ValueError:
                pass                        # fully confined: no phreatic line
            ax.triplot(triang, color='0.85', lw=0.2, zorder=0)
            q = sol.get('flowrate')
            ax.set_title(f'{stem} — total head'
                         + (f'  (Q = {q:.3e})' if q is not None else ''))
        ax.set_aspect('equal')
        fig.tight_layout()
        p = os.path.join(OUT, f'_{stem}_{which}.png')
        fig.savefig(p, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        paths.append(p)

    imgs = [Image.open(p) for p in paths]
    h = min(im.height for im in imgs)
    imgs = [im.resize((int(im.width * h / im.height), h)) for im in imgs]
    combo = Image.new('RGB', (sum(im.width for im in imgs) + 20, h), 'white')
    x = 0
    for im in imgs:
        combo.paste(im, (x, 0))
        x += im.width + 20
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
            print('ok ', stem)
        except Exception as e:                      # noqa: BLE001
            print('FAIL', stem, ':', repr(e)[:100])
