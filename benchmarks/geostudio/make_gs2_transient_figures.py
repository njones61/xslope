"""Render the isochrone comparison figures for the SEEP/W transient corpus
(SEEPW-T01 consolidation, SEEPW-T02 infiltration) in
docs/verification/geostudio.md, "Transient seepage (SEEP/W)".

Both are LINE profiles (the native form of the published figures), so the
field-plot frame spec (equal aspect / colorbars) does not apply; each panel
shares axes and a uniform margin, following the GW transient figures.

  gs2_cons.png   Terzaghi excess PWP vs elevation at t = 150/604/1460 s
                 -- closed form (lines), XSLOPE (open markers), SEEP/W (dots)
  gs2_infil.png  pressure head vs elevation at t = 46 800 s
                 -- SEEP/W node.csv oracle (line), XSLOPE (open markers)

The SEEP/W profiles are the vendor .gsz node.csv results, DISTILLED to plain
arrays here (downsampled); the .gsz itself is never committed.

Run from the repo root:
  PYTHONPATH=. python3 benchmarks/geostudio/make_gs2_transient_figures.py
"""

import io
import os
import sys
import contextlib
import warnings

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from xslope.fileio import load_slope_data
from xslope.mesh import get_material_polygons, build_mesh_from_polygons
from xslope.seep import (build_seep_data, build_tseep_data,
                         run_transient_seepage, transient_frame_index)

import build_gs2_cons as C
import build_gs2_infil as I

SRC = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'verification',
                   'files', 'geostudio')
OUT = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'verification', 'images')

_COLORS = ['#1f77b4', '#d62728', '#2ca02c']

# --- SEEP/W node.csv oracle, distilled (downsampled) from the vendor .gsz -----
# T01 excess PWP (kPa) at y = 0..0.05 (11 pts), t = 150/604/1460 s
_SW_CONS_Y = np.arange(0.0, 0.0501, 0.005)
_SW_CONS = {
    150.0: [0.0, 5.366, 8.235, 9.374, 9.766, 9.859, 9.766, 9.374, 8.235, 5.366, 0.0],
    604.0: [0.0, 2.736, 5.022, 6.635, 7.564, 7.864, 7.564, 6.635, 5.022, 2.736, 0.0],
    1460.0: [0.0, 1.516, 2.860, 3.899, 4.549, 4.770, 4.549, 3.899, 2.860, 1.516, 0.0],
}
# T02 pressure head (m) at y = 0..1 (21 pts), t = 46 800 s
_SW_INF_Y = np.arange(0.0, 1.001, 0.05)
_SW_INF_PSI = [-7.9976, -7.9991, -7.9998, -8.0, -8.0, -8.0, -7.9997, -7.4188,
               -1.6328, -0.7241, -0.4171, -0.2599, -0.1656, -0.1047, -0.0642,
               -0.0372, -0.0197, -0.0090, -0.0033, -0.0011, 0.0]


def _solve(stem, target_size, frac):
    sd = load_slope_data(os.path.join(SRC, f'{stem}.xlsx'))
    ts = build_tseep_data(sd)
    mesh = build_mesh_from_polygons(get_material_polygons(sd), target_size, 'tri3')
    seep = build_seep_data(mesh, sd)
    with contextlib.redirect_stdout(io.StringIO()):
        sol = run_transient_seepage(seep, ts, max_head_change_frac=frac, verbose=False)
    return seep['nodes'], sol


def _sample(nodes, h, xq, yq):
    d2 = (nodes[:, 0] - xq) ** 2 + (nodes[:, 1] - yq) ** 2
    idx = np.argsort(d2)[:4]
    w = 1.0 / np.maximum(d2[idx], 1e-12)
    return float(np.sum(w * h[idx]) / np.sum(w))


def fig_cons():
    cv = C._cv()
    ys = np.linspace(0.0, 0.05, 200)
    ys_s = np.linspace(0.004, 0.046, 11)
    nodes, sol = _solve('gs2_cons', 0.001, 0.02)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for c, t in zip(_COLORS, C._SAVES):
        Tv = cv * t / (C._H * C._H)
        ax.plot(C._U0 * C.terzaghi_ue(ys / C._H, Tv), ys, '-', color=c, lw=1.7,
                label=f't = {t:g} s  ($T_v$ = {Tv:.2f})')
        ax.plot(_SW_CONS[t], _SW_CONS_Y, ':', color=c, lw=1.2, marker='.', ms=5)
        h = np.asarray(sol['frames'][transient_frame_index(sol, t)]['head'])
        ue = np.array([C._GW * (_sample(nodes, h, C._WIDTH / 2, y) - C._H_REF)
                       for y in ys_s])
        ax.plot(ue, ys_s, 'o', color=c, ms=5.5, mfc='white', mew=1.3)
    ax.set_xlabel('excess pore-water pressure  $u_e$  (kPa)')
    ax.set_ylabel('elevation  y  (m)   [drained at y = 0 and y = 0.05]')
    ax.set_title('SEEPW-T01 — Terzaghi 1-D consolidation (double drainage)',
                 fontsize=11)
    ax.grid(alpha=0.3)
    ax.set_xlim(-0.3, 10.3)
    handles = ax.get_legend_handles_labels()[0]
    handles += [Line2D([], [], marker='o', color='0.3', ls='none', mfc='white',
                       mew=1.3, ms=6, label='XSLOPE'),
                Line2D([], [], color='0.3', ls=':', marker='.', label='SEEP/W')]
    ax.legend(handles=handles, loc='upper left', fontsize=8.5,
              title='Terzaghi Eq 17.3 (lines)')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'gs2_cons.png'), dpi=150)
    plt.close(fig)
    return 'gs2_cons.png'


def fig_infil():
    ys_s = np.linspace(0.05, 0.98, 20)
    nodes, sol = _solve('gs2_infil', 0.01, 0.01)
    h = np.asarray(sol['frames'][transient_frame_index(sol, I._SAVE)]['head'])
    psi = np.array([_sample(nodes, h, I._WIDTH / 2, y) - y for y in ys_s])
    fig, ax = plt.subplots(figsize=(6.5, 6))
    ax.plot(_SW_INF_PSI, _SW_INF_Y, '-', color='#1f77b4', lw=1.8,
            label='SEEP/W (node.csv oracle)')
    ax.plot(psi, ys_s, 'o', color='#d62728', ms=5.5, mfc='white', mew=1.3,
            label='XSLOPE')
    ax.axvline(0.0, color='0.6', ls='--', lw=0.9)
    ax.set_xlabel('pressure head  $\\psi$  (m)')
    ax.set_ylabel('elevation  y  (m)   [ponded at y = 1]')
    ax.set_title('SEEPW-T02 — infiltration into dry soil, t = 46 800 s',
                 fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(loc='upper left', fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'gs2_infil.png'), dpi=150)
    plt.close(fig)
    return 'gs2_infil.png'


if __name__ == '__main__':
    for fn in (fig_cons, fig_infil):
        print('ok  ', fn(), flush=True)
