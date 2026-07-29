"""Field-diff validation figure for the RS2-67 own-flow transient march (Cases 3c/3d).

Runs XSLOPE's OWN transient drawdown seepage on the RS2-67 dam with RS2's own hydraulics
(k = 1e-7 m/s, mv = 2e-4; '#067_*.slw'): initial condition = steady full pool (el 24.4),
reservoir stepped to the tailwater (el 7.3) at t=0, marched to 90 h. The 90 h own-flow
phreatic surface is overlaid on the phreatic surface of RS2's OWN imported 90 h field
(the committed rs2_67c sidecars) — the direct fidelity check that XSLOPE's transient FLOW
solver reproduces RS2's solved drawdown field. On the upstream face (where drawdown drives
the response) the two phreatic surfaces coincide to <~0.1 m; they part only in the thin
crest/core, and the march is still far from drained at 1500 h, which is why RS2 renders the
"1500 h" Case-4 stage as the drained steady limit (see build_rs2.py rs2_67e/f).

Run:  python benchmarks/rocscience/make_rs2_67_fielddiff.py
Out:  docs/verification/images/rs2_67_fielddiff.png
"""
import os
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation, LinearTriInterpolator

from build_rs2 import _RS2_67_DAM_90H, _poly_slope_data
from xslope.mesh import get_material_polygons, build_mesh_from_polygons
from xslope.seep import (build_seep_data, build_tseep_data, run_transient_seepage,
                         transient_frame_index)

HERE = os.path.dirname(__file__)
IMG = os.path.join(HERE, '..', '..', 'docs', 'verification', 'images')
FILES = os.path.join(HERE, '..', '..', 'docs', 'verification', 'files', 'rocscience')
K_MHR = 1e-7 * 3600.0            # RS2-67 SK = 1e-7 m/s (TimeUnit sec) -> m/hr

# The march runs on the SAME geometry frame as the field it is compared against: the
# 90 h files' re-imported section (benches 6.663 / 6.86), not the authored one Cases 1,
# 2 and 4 are built on. Comparing two phreatic surfaces across a 0.44-0.64 m difference
# in bench elevation would measure the frames, not the flow solver.
_UP_DAYLIGHT = (87.813, 24.4)    # el 24.4 on (33.787, 6.663) -> (99.272, 28.162)
_DN_DAYLIGHT = (156.420, 7.3)    # el 7.3  on (157.459, 6.86) -> (107.161, 28.162)


def own_flow_90h(target_size=3.0):
    """One transient march (IC = steady full pool 24.4, step to 7.3 at t=0), 90 h frame."""
    m = dict(name='rock1', c=13.8, phi=37.0, gamma=18.2, gamma_sat=18.2, E=1.0e5, nu=0.3,
             u='seep', k1=K_MHR, k2=K_MHR, alpha=0.0, unsat='lf', kr0=0.01, h0=-1.0,
             Ss=9.81 * 0.0002, Sy=0.4)
    sd = _poly_slope_data(polygons=[(0, _RS2_67_DAM_90H)], materials=[m],
                          circle={'Xo': 130.0, 'Yo': 34.0, 'Depth': 4.0, 'R': 30.0},
                          max_depth=0.0)
    sd['time_unit'] = 'hr'
    sd['unit_system'] = 'si'
    sd['seepage_bc'] = {
        'specified_heads': [
            # upstream reservoir drawdown series (submerged-only): el 24.4 -> el 7.3
            {'head': 'res', 'kind': 'reservoir',
             'coords': [(0.266, 0.091), (0.256, 6.663),
                        (33.787, 6.663), _UP_DAYLIGHT]},
            {'head': 7.3, 'coords': [(191.582, 0.091), (191.582, 6.86),
                                     (157.459, 6.86), _DN_DAYLIGHT]},
        ],
        'exit_face': [_DN_DAYLIGHT, (107.161, 28.162)],
    }
    sd['tseep'] = {'times': [0.0, 0.0], 'series': {'res': [24.4, 7.3]},
                   'duration': 1500.0 * 1.001, 'save_interval': None,
                   'stage_1': None, 'stage_2': None, 'save_times': [90.0, 1500.0]}
    mesh = build_mesh_from_polygons(get_material_polygons(sd), target_size=target_size,
                                    element_type='tri6')
    sol = run_transient_seepage(build_seep_data(mesh, sd), build_tseep_data(sd), verbose=False)
    f90 = sol['frames'][transient_frame_index(sol, 90.0)]
    f1500 = sol['frames'][transient_frame_index(sol, 1500.0)]
    return mesh, f90, f1500


def _load_vendor(base):
    m = json.load(open(os.path.join(FILES, f'{base}_mesh.json')))
    u = []
    with open(os.path.join(FILES, f'{base}_seep.csv')) as f:
        for row in csv.DictReader(f):
            if row['node_id'].startswith('#'):
                continue
            u.append(float(row['u']))
    return np.asarray(m['nodes']), np.asarray(u), m


def _sub_tris(elements, element_types):
    tris = []
    for el, et in zip(elements, element_types):
        if et == 6:
            a, b, c, d, e, f = el[:6]
            tris += [[a, d, f], [d, b, e], [f, e, c], [d, e, f]]
        else:
            tris.append([el[0], el[1], el[2]])
    return np.asarray(tris)


def _phreatic(nodes, u, tris, xs, ymin=0.0, ymax=28.2):
    interp = LinearTriInterpolator(Triangulation(nodes[:, 0], nodes[:, 1], triangles=tris), u)
    yy = np.linspace(ymin, ymax, 400)
    out = []
    for x in xs:
        uu = np.ma.masked_invalid(np.array([interp(x, y) for y in yy]))
        pos = np.where((uu >= 0) & ~np.ma.getmaskarray(uu))[0]
        out.append(yy[pos.max()] if len(pos) else np.nan)
    return np.asarray(out)


def make_figure():
    mesh, f90, f1500 = own_flow_90h()
    own_nodes = np.asarray(mesh['nodes'])
    own_tris = _sub_tris(mesh['elements'], mesh['element_types'])
    vend_nodes, vend_u, vend_m = _load_vendor('rs2_67c')
    vend_tris = _sub_tris(vend_m['elements'], vend_m['element_types'])

    xs = np.linspace(35, 155, 41)
    own_ph = _phreatic(own_nodes, f90['u'], own_tris, xs)
    vend_ph = _phreatic(vend_nodes, vend_u, vend_tris, xs)
    up = xs <= 87.813          # upstream face (drawdown-driven)
    d = own_ph - vend_ph
    up_valid = up & ~np.isnan(d)
    all_valid = ~np.isnan(d)
    print(f'upstream-face phreatic match (x<=87.8): mean {np.nanmean(d[up_valid]):+.3f} m, '
          f'max|d| {np.nanmax(np.abs(d[up_valid])):.3f} m')
    print(f'whole-section: RMS {np.sqrt(np.nanmean(d[all_valid]**2)):.2f} m, '
          f'max|d| {np.nanmax(np.abs(d[all_valid])):.2f} m')
    print(f'1500 h march crest head = {float(np.max(f1500["head"])):.1f} m '
          f'(full pool 24.4 -> tailwater 7.3; still draining)')

    fig, ax = plt.subplots(figsize=(8.4, 3.9))
    dam = np.array(_RS2_67_DAM_90H + [_RS2_67_DAM_90H[0]])
    ax.fill(dam[:, 0], dam[:, 1], facecolor='#f2efe9', edgecolor='0.35', lw=1.1, zorder=0)
    ax.plot(xs[all_valid], vend_ph[all_valid], '-', color='#1f6feb', lw=1.8, marker='o',
            ms=3.5, label='RS2 imported field (rs2_67c, 90 h)')
    ax.plot(xs[all_valid], own_ph[all_valid], '--', color='#cf222e', lw=1.8, marker='s',
            ms=3.5, label='XSLOPE own transient flow (90 h)')
    ax.axhline(7.3, color='0.55', lw=0.8, ls=':')
    ax.text(2, 7.3, 'tailwater el 7.3', va='bottom', ha='left', fontsize=7, color='0.4')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('elevation (m)')
    ax.set_title('RS2-67 90 h phreatic surface — own transient flow vs RS2 imported field\n'
                 'reservoir el 24.4 $\\to$ 7.3 at $t=0$; k = 1e-7 m/s (Huang & Jia 2009)',
                 fontsize=9.5)
    ax.set_aspect('equal')
    xr = dam[:, 0].max() - dam[:, 0].min()
    yr = 28.2
    ax.set_xlim(dam[:, 0].min() - 0.03 * xr, dam[:, 0].max() + 0.03 * xr)
    ax.set_ylim(-0.04 * yr, yr * 1.06)
    ax.grid(alpha=0.25)
    ax.legend(loc='upper right', fontsize=8, framealpha=0.95)
    fig.tight_layout()
    out = os.path.join(IMG, 'rs2_67_fielddiff.png')
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


if __name__ == '__main__':
    print(make_figure())
