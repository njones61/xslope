"""Builder for GeoStudio SLOPE/W Verification Manual problem 2.46 —
Eurocode 7: Earth Dam.

Source: GeoStudio Slope Stability Verification Manual (Oct 2022) §2.46
(printed pp. 92-93), after example 5.12 of "Smith's Elements of Soil Mechanics"
(8th ed.). A homogeneous clay dam with a reservoir on the upstream (right)
face; a coupled SEEP/W steady-state analysis supplies the pore pressures and
SLOPE/W checks the downstream (left) slope under Eurocode 7 Design Approach 1,
Combination 2 (DA1-C2).

  Geometry (Fig 136, coords from the .gsz): ground surface
    (0,0)-(4.42,2.6)-(10.37,6.1)-(11.7,6.1)-(13.39,5.1)-(22,0); base at y=0
    (impenetrable). Units metric (kPa, kN/m3, m).
  Material (Table 116, characteristic): Clay c'=12 kPa, phi'=20 deg,
    gamma=19.2 kN/m3.
  Water: reservoir total head 5.1 m against the upstream (right) face; the
    downstream (left) face is a seepage/exit face; phreatic drops to the toe.

PORE PRESSURE — xslope's OWN FE seepage (not SLOPE/W's field)
  The pore pressures come from an xslope steady-state seepage solve (u='seep'),
  written to the sidecars gs2_46_mesh.json / gs2_46_seep.csv and read back on
  load. BCs: specified total head 5.1 on the upstream face (13.39,5.1)-(22,0);
  exit face on the downstream face (0,0)-(4.42,2.6)-(10.37,6.1). The solved
  phreatic reproduces SLOPE/W's SEEP/W phreatic within ~0.1-0.2 m across the
  section (cross-checked against the .gsz field, dev-side only). Homogeneous k,
  so the head field is k-independent.

FACTORED-PARAMETER EMULATION (EC7 DA1-C2, set M2)
  Partial factors read from the .gsz DesignFactors block ("Eurocode 7 - DA1, C2")
  = EN 1997 DA1 combination 2 (set M2 on materials, set A2 = 1.0 on permanent
  actions):
      EffectiveCohesion  gamma_c'    = 1.25  -> c*   = 12 / 1.25       = 9.6 kPa
      EffectiveTanPhi    gamma_phi'  = 1.25  -> phi* = atan(tan20/1.25) = 16.235 deg
      SoilWeight         gamma_gamma = 1.0   -> gamma unchanged (19.2)
  There are no external loads, so the ordinary FS of this factored model equals
  GeoStudio's reported Overdesign Factor (ODF).

Reservoir weight: the reservoir stands only on the upstream (right) face
  (x >= 13.39); the critical downstream circle (center ~(2.56, 9.66), R~9.65)
  spans x in [-7, 12.2] and never reaches it, so the ponded-water weight cannot
  affect this failure. It enters the model only through the seepage field (pore
  pressures), which the FE solve captures. No surface water load is authored.
  Under automatic water loads the engine supplies it anyway, from the head-5.1
  boundary this file already carries -- 250.3 kN/m over the upstream face, which
  is the load the omission left out -- and the three searched factors of safety
  are unchanged to ten digits (M-P 1.073035, Bishop 1.074405, Spencer 1.073572),
  because the critical circle still never reaches x = 13.39. This file is a
  clean candidate for the automatic mode: it gains the physically present
  reservoir without moving a locked number.

Targets (Table 117, printed p.93): ODF M-P 1.091 / Bishop 1.092 / Janbu 1.017;
book (Bishop) 1.07. The .gsz solves M-P ODF = 1.091 on its critical (composite,
base-truncated) surface, and 1.101 on the best circle that xslope can build.
"""

import math
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
# appended, not inserted: this dir is only for the sibling _gs2_donor, and a
# leading entry would let a sibling name shadow a package module.
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from shapely.geometry import LineString, Polygon  # noqa: E402

from xslope.fileio import load_slope_data  # noqa: E402
from xslope.fileio import save_slope_data_to_xlsx as _write_xlsx  # noqa: E402
from xslope.mesh import get_material_polygons, build_mesh_from_polygons, export_mesh_to_json  # noqa: E402
from xslope.seep import build_seep_data, run_seepage_analysis, export_seep_solution  # noqa: E402
from _gs2_donor import donor_material, load_donor  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'verification', 'files', 'geostudio')
ACADS_1A = os.path.join(os.path.dirname(__file__), '..', '..',
                        'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')

# --- DA1-C2 (set M2) partial factors, from the .gsz DesignFactors block ---
GAMMA_C = 1.25
GAMMA_PHI = 1.25

C_CHAR, PHI_CHAR, GAMMA = 12.0, 20.0, 19.2
C_STAR = C_CHAR / GAMMA_C
PHI_STAR = math.degrees(math.atan(math.tan(math.radians(PHI_CHAR)) / GAMMA_PHI))

GS = [(0.0, 0.0), (4.42, 2.6), (10.37, 6.1), (11.7, 6.1), (13.39, 5.1), (22.0, 0.0)]
RES_HEAD = 5.1


def _slope_data():
    sd = load_donor(ACADS_1A)
    m = donor_material(sd)
    m.update(name='Clay (DA1-C2 factored)', c=round(C_STAR, 4), phi=round(PHI_STAR, 4),
             gamma=GAMMA, gamma_sat=GAMMA, option='mc', u='seep',
             k1=1.0e-6, k2=1.0e-6, alpha=0.0, kr0=0.001, h0=-1.0)
    sd['materials'] = [m]
    sd['profile_lines'] = [{'mat_id': 0, 'coords': GS}]
    sd['max_depth'] = 0.0
    sd['gamma_water'] = 9.807
    sd['ground_surface'] = LineString(GS)
    sd['domain_polygon'] = Polygon(GS)
    sd['piezo_line'] = []
    sd['dloads'] = []
    sd['seepage_bc'] = {
        'specified_heads': [{'head': RES_HEAD, 'coords': [(13.39, 5.1), (22.0, 0.0)]}],
        'exit_face': [(0.0, 0.0), (4.42, 2.6), (10.37, 6.1)],
    }
    # downstream (left-slope) critical circle, from SLOPE/W's own solution;
    # the true minimum truncates along the base (composite), so solve with
    # composite=True. Depth = Yo - R = ~0.05 (nearly tangent to the base).
    sd['circular'] = True
    sd['non_circ'] = []
    sd['circles'] = [
        {'Xo': 2.558, 'Yo': 9.659, 'Depth': 0.012, 'R': 9.647},
        {'Xo': 2.0, 'Yo': 8.5, 'Depth': 0.0, 'R': 8.5},
        {'Xo': 3.5, 'Yo': 11.0, 'Depth': 0.2, 'R': 10.8},
    ]
    return sd


def build():
    sd = _slope_data()
    path = os.path.join(OUT, 'gs2_46.xlsx')
    _write_xlsx(sd, path)

    # generate the seepage sidecars so a plain reload finds them
    sd = load_slope_data(path)
    polygons = get_material_polygons(sd)
    xr = [min(x for x, _ in GS), max(x for x, _ in GS)]
    target = (xr[1] - xr[0]) / 80.0
    mesh = build_mesh_from_polygons(polygons, target, 'tri6')
    stem = Path(path)
    export_mesh_to_json(mesh, str(stem.parent / f'{stem.stem}_mesh.json'))
    seep = build_seep_data(mesh, sd)
    sol = run_seepage_analysis(seep, tol=1e-5, max_iter=400)
    export_seep_solution(seep, sol, str(stem.parent / f'{stem.stem}_seep.csv'))
    return path


if __name__ == '__main__':
    p = build()
    print('wrote', p, '(+ _mesh.json, _seep.csv)')
    print('factored: c*=%.4f kPa, phi*=%.4f deg, gamma=%.1f' % (C_STAR, PHI_STAR, GAMMA))
