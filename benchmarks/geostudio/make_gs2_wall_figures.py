"""Figures for the SIGMA/W wall benchmark (docs/verification/geostudio.md,
"Slope stabilization with a wall").

Three PNGs into docs/verification/images/:

  gs2_wall_none.png    inputs | the strength-reduction mechanism, no wall
  gs2_wall.png         inputs | the strength-reduction mechanism, wall in place
  gs2_wall_forces.png  the wall's bending moment and shear down its length,
                       drawn against elevation the way the example's own
                       Figures 17 and 18 are, with the published peaks marked

Nothing is solved here. Every field is read back from the companions committed
beside the models by tools/make_gs2_wall_sidecars.py, so the figures, the page's
numbers and the regression locks all describe one run. The two-panel composites
are stitched the same way the rest of the corpus's are -- two equal-height panels
with a white gutter between them -- because the page's figure check reads the
layout back off the PNG and compares it with the caption.

    PYTHONPATH=. python3 benchmarks/geostudio/make_gs2_wall_figures.py
"""
import contextlib
import io
import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, '..', '..'))
sys.path.insert(0, _REPO)

import matplotlib                                                    # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                      # noqa: E402
import numpy as np                                                   # noqa: E402
from PIL import Image                                                # noqa: E402

from xslope.fem import build_fem_data, import_fem_solution           # noqa: E402
from xslope.fileio import load_slope_data                            # noqa: E402
from xslope.plot import plot_inputs                                  # noqa: E402
from xslope.plot_fem import plot_fem_results                         # noqa: E402

FILES = os.path.join(_REPO, 'docs', 'verification', 'files', 'geostudio')
OUT = os.path.join(_REPO, 'docs', 'verification', 'images')

#: The example's published peak wall actions, read off its Figures 17 and 18.
#: They are drawn on the profile figure as reference markers, and the docs entry
#: states them as read-off values.
PUB_M_PEAK = 1650.0      # kN m/m, at the weak layer
PUB_V_NEG = -750.0       # kN/m, above the weak layer
PUB_V_POS = 1300.0       # kN/m, below it

PANEL = (8.0, 5.0)
DPI = 150


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _load(stem):
    """The model, its fem_data, the solution companions beside it and the run
    record. ``import_fem_solution`` returns the converged field with the captured
    mechanism attached under ``failure_solution``."""
    path = os.path.join(FILES, f'{stem}.xlsx')
    sd = _quiet(load_slope_data, path)
    fem_data = _quiet(build_fem_data, sd, sd['mesh'])
    solution = _quiet(import_fem_solution, fem_data, os.path.join(FILES, stem))
    meta_path = os.path.join(FILES, f'{stem}_fem_meta.json')
    meta = {}
    if os.path.exists(meta_path):
        with open(meta_path) as f:
            meta = json.load(f)
    return sd, fem_data, solution, meta


def _stitch(paths, out):
    imgs = [Image.open(p) for p in paths]
    h = min(im.height for im in imgs)
    imgs = [im.resize((int(im.width * h / im.height), h)) for im in imgs]
    combo = Image.new('RGB', (sum(im.width for im in imgs) + 20, h), 'white')
    x = 0
    for im in imgs:
        combo.paste(im, (x, 0))
        x += im.width + 20
    combo.save(out)
    for p in paths:
        os.remove(p)
    return out


def composite(stem):
    """inputs | the at-failure mechanism, side by side."""
    sd, fem_data, solution, meta = _load(stem)
    failure = solution.get('failure_solution')
    fs = meta.get('FS')

    paths = []
    fig = plt.figure(figsize=PANEL)
    _quiet(plot_inputs, sd, fig=fig, mat_table=False, show_title=True,
           title=f'{stem} — inputs')
    p = os.path.join(OUT, f'_{stem}_inputs.png')
    fig.savefig(p, dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    paths.append(p)

    fig = plt.figure(figsize=PANEL)
    _quiet(plot_fem_results, fem_data, solution, plot_type=['shear_strain'],
           fig=fig, fs=fs, failure_solution=failure,
           field_state='failure' if failure is not None else None,
           show_title=True)
    p = os.path.join(OUT, f'_{stem}_solution.png')
    fig.savefig(p, dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    paths.append(p)

    out = _stitch(paths, os.path.join(OUT, f'{stem}.png'))
    print(f'{stem}: {out}  (FS = {fs})')
    return out


def _profile(fem_data, solution):
    """The wall's actions against elevation.

    Each beam element carries a shear and a moment at each of its two ends; the
    moment profile is the end moments in order down the wall, which is continuous
    because adjacent elements share a node.
    """
    nodes = np.asarray(fem_data['nodes'])
    pairs = list(fem_data['pile_node_pairs'])
    M = np.asarray(solution['forces_pile_moment'])
    V = np.asarray(solution['forces_pile_lateral'])

    order = sorted(range(len(pairs)),
                   key=lambda p: -max(nodes[pairs[p][0]][1], nodes[pairs[p][1]][1]))
    y_m, m_vals, y_v, v_vals = [], [], [], []
    for p in order:
        a, b = pairs[p]
        ya, yb = float(nodes[a][1]), float(nodes[b][1])
        # end 1 is the node at (x1, y1) -- the head end, higher up
        if ya < yb:
            ya, yb = yb, ya
            m1, m2 = M[p][1], M[p][0]
        else:
            m1, m2 = M[p][0], M[p][1]
        # the bending moment at an end is the NEGATIVE of the stiffness
        # convention's nodal moment at end 1, and the nodal moment at end 2
        y_m += [ya, yb]
        m_vals += [-m1, m2]
        y_v += [ya, yb]
        v_vals += [V[p], V[p]]
    return (np.array(y_m), np.array(m_vals), np.array(y_v), np.array(v_vals))


def forces(stem='gs2_wall'):
    """The wall's moment and shear profiles, against the published peaks.

    Drawn as two panels stitched with a white gutter rather than one figure with
    two subplots: the page's figure check reads the layout off the PNG by looking
    for a blank column down the middle, and a title spanning both panels would
    ink that column.
    """
    _sd, fem_data, converged, _meta = _load(stem)
    failure = converged.get('failure_solution')
    y_m, m, y_v, v = _profile(fem_data, failure or converged)

    panels = [
        ('moment', y_m, m, 'Bending moment (kN·m/m)',
         [(PUB_M_PEAK, f'SIGMA/W peak ≈ {PUB_M_PEAK:,.0f}')]),
        ('shear', y_v, v, 'Shear (kN/m)',
         [(PUB_V_NEG, f'SIGMA/W peaks ≈ {PUB_V_NEG:,.0f} / {PUB_V_POS:,.0f}'),
          (PUB_V_POS, None)]),
    ]
    paths = []
    for key, yy, vals, xlabel, pubs in panels:
        fig, ax = plt.subplots(figsize=(4.9, 5.6))
        # the weak clay band, which is where both actions turn over
        ax.axhspan(5.0, 6.0, color='#f39c12', alpha=0.18, zorder=0,
                   label='weak clay band')
        ax.plot(vals, yy, '-', color='#1f4e79', lw=1.8, label='XSLOPE')
        ax.axvline(0.0, color='0.6', lw=0.8)
        for xv, lb in pubs:
            ax.axvline(xv, color='#c0392b', ls='--', lw=1.2,
                       label=lb if lb else None)
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Elevation (m)')
        ax.set_title(f'gs2_wall — wall {key} at the mechanism', fontsize=10)
        ax.grid(alpha=0.25)
        ax.legend(loc='lower left', fontsize=8, frameon=False)
        fig.tight_layout()
        p = os.path.join(OUT, f'_gs2_wall_{key}.png')
        fig.savefig(p, dpi=DPI, bbox_inches='tight')
        plt.close(fig)
        paths.append(p)

    out = _stitch(paths, os.path.join(OUT, 'gs2_wall_forces.png'))
    print(f'wall actions: M {m.min():.0f} .. {m.max():.0f} kN·m/m, '
          f'V {v.min():.0f} .. {v.max():.0f} kN/m -> {out}')
    return out


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    only = set(sys.argv[1:])
    if not only or 'gs2_wall_none' in only:
        composite('gs2_wall_none')
    if not only or 'gs2_wall' in only:
        composite('gs2_wall')
        forces('gs2_wall')
