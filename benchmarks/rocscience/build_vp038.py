"""Builder for Rocscience Slide2 Verification Problem #38 —
Excavated slope, homogeneous soil, FE groundwater seepage, matric suction.

Source: Slide_SlopeStabilityVerification.pdf §38 (printed pp. 144-146), after
Ng & Shi (1998), "A numerical investigation of the stability of unsaturated soil
slopes subjected to transient seepage". A 28-degree cut slope in Hong Kong: a
24 m soil layer over a 6 m bedrock layer. Steady-state FE groundwater supplies
BOTH the positive AND the negative (matric-suction) pore pressures, and the
Bishop simplified method assesses the cut. The negative pore water pressure
above the water table raises the soil shear strength via the modified
Mohr-Coulomb (Fredlund) criterion:

    tau = c' + (sigma_n - u_a) tan(phi') + (u_a - u_w) tan(phi_b)

i.e. an apparent cohesion (u_a - u_w) tan(phi_b) on the unsaturated slices.

GEOMETRY (from the vendor RS2 mesh's external boundary; Fig 38.1 coords match to
the printed precision (0,19)-(41,41)-(51,41)-(57,49)-(98,71)):
  ground surface (0.13,19.15)-(40.74,40.6)-(50.95,40.6)-(57.06,49.31)-(97.89,70.78);
  the domain closes down the right edge to (97.89,40.49), along an impermeable
  base to the toe (0.0,-11.51), and up the left edge to (0.13,19.15). One
  homogeneous soil fills the domain -- the 6 m "bedrock" band is far below the
  critical surface (which sits entirely in the upper cut), so giving it the soil
  strength cannot affect the specified-surface factor of safety, and the base is
  no-flow either way. Units metric (kPa, kN/m3, m).

MATERIAL (Table 38.1): c'=10 kPa, phi'=38 deg, phi_b=15 deg, gamma=16 kN/m3.
  phi_b is applied as the opt-in suction-strength angle in the LEM tag, NOT baked
  into the file. The permeability function (Fig 38.2 legend / Ng 1998) is a
  Gardner k_r(psi) = 1 / (1 + a |psi|^n) fit in log space to RS2's own .slw
  k-function (ks = 4.19 m/day; a = 7.479, n = 2.908; see scratchpad swcc_fit).
  These SWCC/k parameters are DATA, not tuned.

SEEPAGE BCs (§38.2): constant total head H on the right side of the slope
  (H = 61 / 62 / 63 m, the three cases), constant head 6 m on the left side, the
  ground surface a seepage/exit face, the base no-flow. xslope solves the steady
  unsaturated field itself (u='seep'), writing the sidecars vpNNN_mesh.json /
  vpNNN_seep.csv that a plain reload picks up.

RESULTS (Table 38.2, Bishop simplified):
  H=61m  Slide 1.621  Ng & Shi 1.636
  H=62m  Slide 1.538  Ng & Shi 1.527
  H=63m  Slide 1.407  Ng & Shi 1.436

Slide prints its H=61 critical surface (Fig 38.2): Bishop FS 1.621, center
(47.490, 56.311), radius 16.087, endpoints (50.953, 40.601)-(63.120, 52.500) --
a shallow circle in the upper, unsaturated part of the cut. Only H=61 has a
figure; the critical geometry is head-invariant to that precision, so all three
cases carry Slide's printed H=61 circle and are locked on it as specified
surfaces (single_circle + suction_phi_b). A free search WITH suction localizes
to a slightly deeper circle a few percent below the published shallow-circle FS
(the search-selection difference the specified-surface lock is immune to).

Run from the repo root:  PYTHONPATH=. python3 benchmarks/rocscience/build_vp038.py
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from shapely.geometry import LineString, Polygon  # noqa: E402

from xslope.fileio import load_slope_data  # noqa: E402
from xslope.fileio import save_slope_data_to_xlsx as _write_xlsx  # noqa: E402
from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,  # noqa: E402
                         export_mesh_to_json)
from xslope.seep import (build_seep_data, run_seepage_analysis,  # noqa: E402
                         export_seep_solution)
sys.path.insert(0, os.path.dirname(__file__))
from vendor_tcut import apply_vendor_t_cut  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'verification', 'files', 'rocscience')
ACADS_1A = os.path.join(os.path.dirname(__file__), '..', '..',
                        'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')

SOIL = 'Cut soil'
# Ground surface and closed domain, from the vendor RS2 mesh external boundary.
GS = [(0.13, 19.15), (40.74, 40.6), (50.95, 40.6), (57.06, 49.31), (97.89, 70.78)]
DOMAIN = [(0.0, -11.51), (97.89, 40.49), (97.89, 70.78), (57.06, 49.31),
          (50.95, 40.6), (40.74, 40.6), (0.13, 19.15)]
# Gardner k-function fit (log space) to RS2's .slw curve, ks = 4.19 m/day.
KS, GARD_A, GARD_N = 4.19, 7.479, 2.908
# Slide's printed H=61 critical surface (Fig 38.2); shared across the three heads.
SLIDE_CIRCLE = {'Xo': 47.490, 'Yo': 56.311, 'R': 16.087, 'Depth': 56.311 - 16.087}
# Right-edge total head per case; left edge is always 6 m.
HEADS = {'vp038a': 61.0, 'vp038b': 62.0, 'vp038c': 63.0}
LEFT_HEAD = 6.0


def _slope_data(head):
    sd = load_slope_data(ACADS_1A)
    m = dict(sd['materials'][0])
    m.update(name=SOIL, c=10.0, phi=38.0, gamma=16.0, gamma_sat=16.0,
             option='mc', u='seep', k1=KS, k2=KS, alpha=0.0, kr0=0.001, h0=-1.0,
             unsat='gard', vg_a=GARD_A, vg_n=GARD_N)
    sd['materials'] = [m]
    sd['profile_lines'] = []
    sd['polygons'] = [{'mat_id': 0, 'polygon': Polygon(DOMAIN)}]
    sd['max_depth'] = None
    sd['gamma_water'] = 9.81
    sd['ground_surface'] = LineString(GS)
    sd['domain_polygon'] = Polygon(DOMAIN)
    sd['piezo_line'] = []
    sd['dloads'] = []
    sd['circular'] = True
    sd['non_circ'] = []
    sd['circles'] = [dict(SLIDE_CIRCLE)]
    sd['seepage_bc'] = {
        'specified_heads': [
            {'head': float(head), 'coords': [(97.89, 40.49), (97.89, 70.78)]},
            {'head': LEFT_HEAD, 'coords': [(0.0, -11.51), (0.13, 19.15)]},
        ],
        'exit_face': [(0.13, 19.15), (40.74, 40.6), (50.95, 40.6),
                      (57.06, 49.31), (97.89, 70.78)],
    }
    return sd


def _build_one(stem, head):
    path = os.path.join(OUT, f'{stem}.xlsx')
    sd0 = _slope_data(head)
    # No vendor cap on this row; the call clears any t_cut inherited from ACADS_1A.
    apply_vendor_t_cut(sd0.get('materials', []), path)
    _write_xlsx(sd0, path)
    # Seepage sidecars so a plain reload finds the unsaturated field.
    sd = load_slope_data(path)
    polygons = get_material_polygons(sd)
    target = 97.89 / 70.0
    mesh = build_mesh_from_polygons(polygons, target, 'tri6')
    p = Path(path)
    export_mesh_to_json(mesh, str(p.parent / f'{p.stem}_mesh.json'))
    seep = build_seep_data(mesh, sd)
    sol = run_seepage_analysis(seep, tol=1e-5, max_iter=600)
    export_seep_solution(seep, sol, str(p.parent / f'{p.stem}_seep.csv'))
    return f'{stem}.xlsx'


def vp038a():
    """Slide #38, H=61 m (right-side total head). Published Bishop 1.621."""
    return _build_one('vp038a', HEADS['vp038a'])


def vp038b():
    """Slide #38, H=62 m. Published Bishop 1.538."""
    return _build_one('vp038b', HEADS['vp038b'])


def vp038c():
    """Slide #38, H=63 m. Published Bishop 1.407."""
    return _build_one('vp038c', HEADS['vp038c'])


BUILDERS = [vp038a, vp038b, vp038c]


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    for b in BUILDERS:
        print('built', b(), '(+ _mesh.json, _seep.csv)')
