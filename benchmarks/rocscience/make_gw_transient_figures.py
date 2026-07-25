"""Render the isochrone comparison figures for the transient Groundwater corpus
(GW#15, #16, #21) in docs/verification/rocscience_groundwater.md.

Each figure overlays the CLOSED-FORM / recomputed-series target (solid lines) with
XSLOPE's transient solve sampled at the same points (markers), at several times:

  gw015.png  Terzaghi 1-D consolidation — ue/u0 vs depth, double & single drainage
  gw016.png  Pyrah two-layer consolidation — ue/u0 vs depth, uniform / A-B / B-A
  gw021.png  Ferris confined aquifer — head rise vs distance at 600 hr, two ICs

These are LINE plots (a profile vs depth/distance), the native form of the
published figures (Terzaghi Fig 17-3/17-5, Pyrah Fig 18-2/3/4, Ferris Fig 23-6/7),
so the field-plot frame spec (equal aspect / colorbars) does not apply; the panels
share axes, a single legend, and a uniform margin.

Run from the repo root:  python benchmarks/rocscience/make_gw_transient_figures.py
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

_XSLOPE_HANDLE = Line2D([], [], marker='o', color='0.25', ls='none',
                        mfc='white', mew=1.4, ms=6, label='XSLOPE')

from xslope.fileio import load_slope_data
from xslope.mesh import get_material_polygons, build_mesh_from_polygons
from xslope.seep import (build_seep_data, build_tseep_data,
                         run_transient_seepage, transient_frame_index)

import build_groundwater as B

SRC = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'verification',
                   'files', 'rocscience_gw')
OUT = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'verification', 'images')

_COLORS = ['#1f77b4', '#d62728', '#2ca02c', '#9467bd', '#ff7f0e']


def _solve(stem, target_size, frac=None):
    sd = load_slope_data(os.path.join(SRC, f'{stem}.xlsx'))
    ts = build_tseep_data(sd)
    mesh = build_mesh_from_polygons(get_material_polygons(sd), target_size, 'tri3')
    seep = build_seep_data(mesh, sd)
    kw = {'verbose': False}
    if frac is not None:
        kw['max_head_change_frac'] = frac
    with contextlib.redirect_stdout(io.StringIO()):
        sol = run_transient_seepage(seep, ts, **kw)
    return seep['nodes'], sol


def _sample(nodes, h, xq, yq):
    d2 = (nodes[:, 0] - xq) ** 2 + (nodes[:, 1] - yq) ** 2
    idx = np.argsort(d2)[:4]
    w = 1.0 / np.maximum(d2[idx], 1e-12)
    return float(np.sum(w * h[idx]) / np.sum(w))


def fig_gw15():
    cv = B._GW15_K / B._gw15_ss()
    ys = np.linspace(0.02, 0.98, 200)
    ys_s = np.linspace(0.1, 0.9, 9)
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    for ax, (stem, H, saves, Zof, title) in zip(axes, [
            ('gw015a', 0.5, B._GW15_SAVES['a'], lambda y: 2 * (1 - y),
             'Case 1 — double drainage (H = 0.5 m)'),
            ('gw015b', 1.0, B._GW15_SAVES['b'], lambda y: 1 - y,
             'Case 2 — single drainage (H = 1.0 m)')]):
        nodes, sol = _solve(stem, 0.02)
        for c, t in zip(_COLORS, saves):
            Tv = cv * t / (H * H)
            ax.plot(B.terzaghi_ue(Zof(ys), Tv), ys, '-', color=c, lw=1.6,
                    label=f'$T_v$ = {Tv:.2f}')
            h = np.asarray(sol['frames'][transient_frame_index(sol, t)]['head'])
            ue = np.array([(_sample(nodes, h, 0.125, y) - B._H_REF) / B._GW15_U0
                           for y in ys_s])
            ax.plot(ue, ys_s, 'o', color=c, ms=5, mfc='white', mew=1.3)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel('excess pore pressure  $u_e/u_0$')
        ax.grid(alpha=0.3)
        ax.set_xlim(-0.02, 1.02)
    axes[0].set_ylabel('elevation  y  (m)   [drained top at y = 1]')
    axes[0].legend(loc='center left', fontsize=9, title='Terzaghi Eq 17.3')
    axes[1].legend(handles=[_XSLOPE_HANDLE], loc='lower right', fontsize=9)
    fig.suptitle('GW15 — Terzaghi 1-D consolidation: analytical (lines) vs XSLOPE (points)',
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(os.path.join(OUT, 'gw015.png'), dpi=150)
    plt.close(fig)
    return 'gw015.png'


def fig_gw16():
    ys = np.linspace(0.001, 0.999, 400)
    ys_s = np.linspace(0.08, 0.92, 11)
    saves = B._GW16_SAVES
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 5), sharey=True)
    specs = [
        ('gw016a', 'Case 1 — uniform Soil A', None),
        ('gw016b', 'Case 2 — A (top) / B (bottom)', (B._GW16_A, B._GW16_B)),
        ('gw016c', 'Case 3 — B (top) / A (bottom)', (B._GW16_B, B._GW16_A)),
    ]
    for ax, (stem, title, layers) in zip(axes, specs):
        nodes, sol = _solve(stem, 0.02, frac=B._GW16_FRAC)
        if layers is None:
            ana = lambda yy, t: B.terzaghi_ue(1 - yy, t)
        else:
            top, bottom = layers
            betas = B._pyrah_betas(bottom[0], top[0])
            ana = (lambda yy, t, top=top, bottom=bottom, betas=betas:
                   B.pyrah_ue(yy, t, bottom[0], top[0], bottom[1], top[1],
                              betas=betas, u0=1.0))
            ax.axhline(0.5, color='0.6', ls='--', lw=1)
        for c, t in zip(_COLORS, saves):
            ax.plot(ana(ys, t), ys, '-', color=c, lw=1.6, label=f't = {t:g}')
            h = np.asarray(sol['frames'][transient_frame_index(sol, t)]['head'])
            ue = np.array([(_sample(nodes, h, 0.25, y) - B._H_REF) / B._GW16_U0
                           for y in ys_s])
            ax.plot(ue, ys_s, 'o', color=c, ms=5, mfc='white', mew=1.3)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel('excess pore pressure  $u_e/u_0$')
        ax.grid(alpha=0.3)
        ax.set_xlim(-0.02, 1.02)
    axes[0].set_ylabel('elevation  y  (m)   [drained top at y = 1]')
    axes[0].legend(loc='lower right', fontsize=9, title='time (cv=1)')
    axes[2].legend(handles=[_XSLOPE_HANDLE], loc='lower right', fontsize=9)
    fig.suptitle('GW16 — Pyrah two-layer consolidation: recomputed series (lines) vs '
                 'XSLOPE (points)', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(os.path.join(OUT, 'gw016.png'), dpi=150)
    plt.close(fig)
    return 'gw016.png'


def fig_gw21():
    D = B._gw21_D()
    xs = np.linspace(0, 100, 300)
    xs_s = np.array([10, 20, 30, 40, 50, 60, 70.])
    fig, ax = plt.subplots(figsize=(8, 5))
    for c, (stem, ic, lbl) in zip(_COLORS, [
            ('gw021a', 0.0, 'Case 1: IC = 0, step to 5 ft'),
            ('gw021b', 5.0, 'Case 2: IC = 5 ft, step to 10 ft')]):
        nodes, sol = _solve(stem, 0.8)
        ax.plot(xs, ic + B.ferris_dh(xs, B._GW21_T, D, B._GW21_DH), '-',
                color=c, lw=1.8, label=lbl)
        h = np.asarray(sol['frames'][transient_frame_index(sol, B._GW21_T)]['head'])
        hv = np.array([_sample(nodes, h, x, 2.5) - B._H_REF for x in xs_s])
        ax.plot(xs_s, hv, 'o', color=c, ms=6, mfc='white', mew=1.4)
    ax.plot([], [], 'ko', mfc='white', label='XSLOPE')
    ax.set_xlabel('distance from stepped face,  x  (ft)')
    ax.set_ylabel('head above datum  (ft)')
    ax.set_title('GW21 — Ferris confined aquifer at t = 600 hr: erfc (lines) vs XSLOPE (points)',
                 fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'gw021.png'), dpi=150)
    plt.close(fig)
    return 'gw021.png'


# Digitized Ref [1] (RS2/FlexPDE) toe-slope total head from Slide2/RS2 Fig 20.5
# (read off the published chart; used only as the visual overlay, not the lock).
_GW18_FIG205 = {
    0.6: ([28, 30, 32, 35, 38, 40, 42, 45, 48, 50],
          [2.88, 2.78, 2.68, 2.50, 2.36, 2.22, 2.08, 1.83, 1.50, 1.05]),
    19656: ([28, 30, 32, 35, 38, 40, 42, 45, 48, 50, 52],
            [8.35, 7.95, 7.55, 7.00, 6.35, 5.65, 4.85, 3.45, 2.00, 1.20, 0.0]),
}


def _toe_y(x):
    return 12.0 - (x - 28.0) / 2.0          # downstream 2:1 face, x in [28,52]


def fig_gw18():
    """Toe-slope total head vs x: XSLOPE (solid, dense sampling of the downstream
    face) against the digitized Ref [1] Fig 20.5 profile (markers), early 0.6 h and
    near-steady (locked at 1000 h; the manual reports 19656 h — the same steady
    curve)."""
    nodes, sol = _solve('gw018', 1.5, frac=0.25)
    xs = np.linspace(28.0, 52.0, 120)
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    for c, (t_solve, t_pub, lbl) in zip(_COLORS, [
            (0.6, 0.6, 't = 0.6 h (early transient)'),
            (1000.0, 19656, 'near-steady (t = 19656 h)')]):
        h = np.asarray(sol['frames'][transient_frame_index(sol, t_solve)]['head'])
        th = np.array([_sample(nodes, h, x, _toe_y(x)) for x in xs])
        ax.plot(xs, th, '-', color=c, lw=1.8, label=f'XSLOPE — {lbl}')
        px, ph = _GW18_FIG205[t_pub]
        ax.plot(px, ph, 's', color=c, ms=5, mfc='white', mew=1.3)
    ax.plot([], [], 's', color='0.35', mfc='white', mew=1.3,
            label='Ref [1] — digitized Fig 20.5')
    ax.set_xlabel('x coordinate along toe slope  (m)')
    ax.set_ylabel('total head  (m)')
    ax.set_title('GW18 — toe-slope total head: XSLOPE (lines) vs digitized Fig 20.5 (points)',
                 fontsize=11)
    ax.set_xlim(25, 55)
    ax.set_ylim(0, 9)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9, loc='upper right')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'gw018.png'), dpi=150)
    plt.close(fig)
    return 'gw018.png'


def fig_gw17():
    """XSLOPE near-steady total-head field (500 h) for the toe-drain dam, rendered
    through the package's own plot_seep_solution — the visual analog of the
    published Fig 19-5 total-head contours (reservoir 10 drawn down to the toe drain
    at head 0)."""
    from xslope.plot_seep import plot_seep_solution
    sd = load_slope_data(os.path.join(SRC, 'gw017.xlsx'))
    ts = build_tseep_data(sd)
    mesh = build_mesh_from_polygons(get_material_polygons(sd), 1.0, 'tri3')
    seep = build_seep_data(mesh, sd)
    with contextlib.redirect_stdout(io.StringIO()):
        sol = run_transient_seepage(seep, ts, verbose=False, max_head_change_frac=0.25)
    fr = sol['frames'][transient_frame_index(sol, 500.0)]
    frame_solution = {
        'head': np.asarray(fr['head']), 'u': np.asarray(fr['u']),
        # No stream function is stored for a transient frame (no flow net); this
        # figure draws flowlines=False anyway.
        'phi': None, 'flowrate': fr.get('inflow'),
        'inflow': fr.get('inflow'), 'outflow': fr.get('outflow'),
        'unconfined': True,
    }
    fig = plt.figure(figsize=(9.0, 4.0))
    with contextlib.redirect_stdout(io.StringIO()):
        plot_seep_solution(seep, frame_solution, fig=fig, show_title=True,
                           fill_contours=True, phreatic=True, flowlines=False,
                           mesh=False)
    fig.axes[0].set_title('GW17 — near-steady total head (t = 500 h): reservoir 10 '
                          'drawn to the toe drain (cf. Fig 19-5)', fontsize=10)
    fig.savefig(os.path.join(OUT, 'gw017.png'), dpi=150)
    plt.close(fig)
    return 'gw017.png'


def fig_gw19():
    """Pressure head along the top boundary (y = 10) vs x at the four report times
    — the digitizable Fig 21.9 analog: the lagoon-leak pressure mound spreading from
    the centerline (x = 0) toward the far field as the lined pond fills.  The lagoon
    footprint (x in [0,2]) is shaded; the four lock stations are marked."""
    nodes, sol = _solve('gw019', 0.8, frac=0.25)
    times = [73.0, 416.0, 792.0, 11340.0]
    xs = np.linspace(0.0, 19.0, 160)
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    for c, t in zip(_COLORS, times):
        h = np.asarray(sol['frames'][transient_frame_index(sol, t)]['head'])
        ph = np.array([_sample(nodes, h, x, 10.0) - 10.0 for x in xs])
        ax.plot(xs, ph, '-', color=c, lw=1.8, label=f't = {t:g} min')
    ax.axvspan(0.0, 2.0, color='0.85', zorder=0)
    ax.text(1.0, 0.9, 'lagoon', ha='center', va='top', fontsize=8,
            transform=ax.get_xaxis_transform())
    ax.set_xlabel('x along top boundary  (m)   [centerline at x = 0]')
    ax.set_ylabel('pressure head  (m)')
    ax.set_title('GW19 — pressure head along the top boundary as the lagoon fills '
                 '(cf. Fig 21.9)', fontsize=10)
    ax.set_xlim(0, 19)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9, loc='upper right', title='XSLOPE')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'gw019.png'), dpi=150)
    plt.close(fig)
    return 'gw019.png'


def fig_gw20():
    """Total head along a vertical query line at x = 2.0 (through the perched zone
    above the fine lens) vs elevation at the three report times — the Fig 22.7
    analog: as rainfall switches on, the perched mound builds above the low-k fine
    lens (shaded y in [0.6,0.7]) and the water table rises from its initial el 0.3."""
    nodes, sol = _solve('gw020', 0.04, frac=0.25)
    times = [4.6, 31.0, 208.0]
    ys = np.linspace(0.2, 1.0, 160)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for c, t in zip(_COLORS, times):
        h = np.asarray(sol['frames'][transient_frame_index(sol, t)]['head'])
        th = np.array([_sample(nodes, h, 2.0, yy) for yy in ys])
        ax.plot(th, ys, '-', color=c, lw=1.8, label=f't = {t:g} s')
    ax.axhspan(0.6, 0.7, color='0.85', zorder=0)
    ax.text(0.02, 0.65, 'fine lens', ha='left', va='center', fontsize=8,
            transform=ax.get_yaxis_transform())
    ax.set_xlabel('total head  (m)')
    ax.set_ylabel('elevation  y  (m)   [query line at x = 2.0]')
    ax.set_title('GW20 — total head along the query line as rainfall perches on the '
                 'lens (cf. Fig 22.7)', fontsize=9.5)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9, loc='upper right', title='XSLOPE')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'gw020.png'), dpi=150)
    plt.close(fig)
    return 'gw020.png'


if __name__ == '__main__':
    for fn in (fig_gw15, fig_gw16, fig_gw21, fig_gw18, fig_gw17, fig_gw19, fig_gw20):
        print('ok  ', fn(), flush=True)
