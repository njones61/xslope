"""Builder for Rocscience Slide2 Verification Problem #46, stages 2 and 3 —
Baker, Frydman & Talesnick (1993) three-stage validation dam.

Stage 1 (dry / end-of-construction) is built by ``vp046()`` in build_problems.py
(circular search on the c'=0 upstream skin, FS 2.500). Stages 2 and 3 solve their
own FE seepage / carry a spatially varying undrained strength, so — like VP38 —
they live in a self-contained builder that writes the mesh/seep sidecars.

Source: Slide_SlopeStabilityVerification.pdf §46 (printed pp. 170-176) and its
primary source Baker, Frydman & Talesnick (1993), "Slope stability analysis for
undrained loading conditions", Int. J. Numer. Anal. Methods Geomech. 17:15-43
(the Givat Brener reservoir, §6, pp. 32-42).

GEOMETRY (Slide Fig 46.1/46.2, axis-tick calibrated; same as stage 1): a
compacted-clay embankment crest el 101 sits on deep natural clay; the UPSTREAM
(reservoir, right) face drops on a 4H:1V natural-clay face from (128,95) to the
toe bench (148,90) and runs flat to x=220; the base (regional water table) is
el 0. Ground surface (0,95)-(80,95)-(100,101)-(128,95)-(148,90)-(220,90).

MATERIALS (Table 46.1 / Baker p. 32): compacted clay c'=6.5 kPa, phi'=40, gamma=18;
natural clay c'=0, phi'=32, gamma=18.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/rocscience/build_vp046.py
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx  # noqa: E402
from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,  # noqa: E402
                         export_mesh_to_json)
from xslope.seep import (build_seep_data, run_seepage_analysis,  # noqa: E402
                         export_seep_solution)
sys.path.insert(0, os.path.dirname(__file__))
from vendor_tcut import apply_vendor_t_cut  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'verification', 'files', 'rocscience')
ACADS_1A = os.path.join(os.path.dirname(__file__), '..', '..',
                        'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')

GW = 9.81
# Reservoir surface el 100; waterline on the compacted upstream face
# (100,101)-(128,95) is at y=100, i.e. x = 100 + (101-100)/6 * 28.
XW = 100.0 + (101.0 - 100.0) / 6.0 * 28.0
# Table 46.1 seepage parameters: Ks = 7e-5 both clays (its absolute value is
# FS-irrelevant — the steady head field depends only on the published ratios:
# equal K between clays and K2/K1 = 0.1, i.e. the 10:1 horizontal:vertical
# anisotropy Baker states on p. 32); natural-clay unsaturated fringe Gardner
# a = 0.06, n = 2 (the manual's printed estimate). Compacted fill takes no
# unsaturated reduction (a = n = 0 in the table -> kr = 1, here the linear
# front with kr0 = 1).
KS, RATIO = 7e-5, 0.1
GROUND = [(0.0, 95.0), (80.0, 95.0), (100.0, 101.0), (128.0, 95.0),
          (148.0, 90.0), (220.0, 90.0)]
# Upstream-slope search seeds. Slide keeps stage 1's search limits (Fig 46.2
# note: "Limits are as they were before"), which target the upstream/reservoir
# slope Baker's whole analysis is about; a global grid instead rides an
# unsupported downstream-toe mechanism the published analyses do not report.
STAGE2_SEEDS = [{'Xo': 150.0, 'Yo': 118.0, 'Depth': 55.0, 'R': 38.0},
                {'Xo': 140.0, 'Yo': 122.0, 'Depth': 58.0, 'R': 40.0},
                {'Xo': 158.0, 'Yo': 112.0, 'Depth': 55.0, 'R': 42.0}]


def _stage2_slope_data():
    sd = load_slope_data(ACADS_1A)
    base = dict(sd['materials'][0])
    m0 = dict(base)
    m0.update(name='Compacted clay', c=6.5, phi=40.0, gamma=18.0, gamma_sat=18.0,
              option='mc', u='seep', k1=KS, k2=KS * RATIO, alpha=0.0,
              unsat='lf', kr0=1.0, h0=-1.0)
    m1 = dict(base)
    m1.update(name='Natural clay', c=0.0, phi=32.0, gamma=18.0, gamma_sat=18.0,
              option='mc', u='seep', k1=KS, k2=KS * RATIO, alpha=0.0,
              unsat='gard', vg_a=0.06, vg_n=2.0)
    sd['materials'] = [m0, m1]
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(80.0, 95.0), (100.0, 101.0), (128.0, 95.0)]},
        {'mat_id': 1, 'coords': [(0.0, 95.0), (80.0, 95.0), (128.0, 95.0),
                                 (148.0, 90.0), (220.0, 90.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = GW
    sd['piezo_line'] = []
    # Full reservoir at el 100 as a distributed load on the submerged upstream
    # boundary (Normal = gamma_w * depth of water); same idiom as VP102b.
    sd['dloads'] = [[{'X': XW, 'Y': 100.0, 'Normal': 0.0},
                     {'X': 128.0, 'Y': 95.0, 'Normal': GW * 5.0},
                     {'X': 148.0, 'Y': 90.0, 'Normal': GW * 10.0},
                     {'X': 220.0, 'Y': 90.0, 'Normal': GW * 10.0}]]
    sd['seepage_bc'] = {
        'specified_heads': [
            {'head': 100.0, 'coords': [(XW, 100.0), (128.0, 95.0),
                                       (148.0, 90.0), (220.0, 90.0)]},
            {'head': 0.0, 'coords': [(0.0, 0.0), (220.0, 0.0)]},
        ],
        'exit_face': [(0.0, 95.0), (80.0, 95.0), (100.0, 101.0), (XW, 100.0)],
    }
    sd['circular'] = True
    sd['non_circ'] = []
    sd['circles'] = [dict(c) for c in STAGE2_SEEDS]
    return sd


def _write_seep(stem, sd):
    """Write the xlsx plus the FE-seepage sidecars a plain reload picks up."""
    path = os.path.join(OUT, f'{stem}.xlsx')
    # No vendor cap on these stages; the call clears any t_cut inherited from ACADS_1A.
    apply_vendor_t_cut(sd.get('materials', []), path)
    save_slope_data_to_xlsx(sd, path)
    sd2 = load_slope_data(path)
    polygons = get_material_polygons(sd2)
    mesh = build_mesh_from_polygons(polygons, 220.0 / 70.0, 'tri6')
    p = Path(path)
    export_mesh_to_json(mesh, str(p.parent / f'{p.stem}_mesh.json'))
    seep = build_seep_data(mesh, sd2)
    sol = run_seepage_analysis(seep, tol=1e-5, max_iter=600)
    export_seep_solution(seep, sol, str(p.parent / f'{p.stem}_seep.csv'))
    return f'{stem}.xlsx'


def vp046b():
    """Slide #46 stage 2: steady-state seepage, FULL reservoir (el 100).
    Effective-stress (c'/phi') analysis on the FE seepage field. Upstream-slope
    circular search -> Spencer 7.086 / Bishop 7.093 vs Slide 7.003, Baker 6.98."""
    return _write_seep('vp046b', _stage2_slope_data())


BUILDERS = [vp046b]


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    for b in BUILDERS:
        print('built', b(), '(+ _mesh.json, _seep.csv)')
