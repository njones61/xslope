"""Builder for the SIGMA/W example "Slope stabilization with piles using SRS"
-> docs/verification/geostudio.md, "Slope stabilization with a wall (SIGMA/W)".

This is the one benchmark in the corpus that exercises XSLOPE's structural beam
element inside a real slope: a marginally stable slope over a weak clay band,
stabilized by a continuous sheet pile wall driven from the lower bench through
the weak layer into the competent material below. GeoStudio solves it with a
Strength Reduction Stability (SRS) analysis in SIGMA/W, with and without the
wall, and publishes the wall's moment and shear distributions.

Two workbooks, identical but for the wall:

  gs2_wall_none.xlsx   the slope alone      -- SIGMA/W SRS ~1.025
  gs2_wall.xlsx        the same slope + wall -- SIGMA/W SRS ~1.4

MODEL (transcribed from the vendor .gsz and the example write-up; the .gsz is a
read-only oracle and is never committed)

  Extents x in [1, 65], y in [-4, 20]. Ground surface, left to right:
      (1,20) (20,20) (36,10) (38,10) (40,10) (42,9) (44,8) (65,8)
  a flat crest at el. 20, the slope face down to a bench at el. 10 spanning
  x = 36-40 (where the wall head sits), then down to the toe at (44,8) and flat
  ground out to x = 65.

  Three zones, all gamma = 20 kN/m3, nu' = 0.4:
      Sandy Clay   el. 6 -> ground surface   E' = 10,000 kPa  c' = 20  phi' = 30
      Weak Clay    el. 5 -> 6 (1 m band)     E' =  5,000 kPa  c' =  0  phi' = 10
      OC Soil      el. -4 -> 5               E' = 500,000 kPa  LINEAR ELASTIC
  The OC soil is GeoStudio's "Isotropic Elastic" stress model with a "Bedrock"
  slope model -- it carries no strength because it cannot fail. XSLOPE expresses
  exactly that with the mat sheet's option = 'elastic'.

  The wall: a vertical line from (38,10) to (38,1), E = 2e8 kPa (200 GPa),
  A = 0.02 m2/m, I = 0.0005 m4/m. It is CONTINUOUS out of plane -- the .gsz
  carries no spacing parameter at all, because a SIGMA/W beam is a plane-strain
  element whose EA and EI are already per metre of wall. XSLOPE's S column is the
  same idea inverted (it DIVIDES a per-pile stiffness by the spacing), so a
  continuous wall is S = 1 and the transcribed EA/EI ride through unchanged.

PORE PRESSURE -- XSLOPE's OWN seepage solve

  GeoStudio drives this problem from a steady-state SEEP/W analysis whose field
  is then handed to every SIGMA/W analysis. That field is not published as
  numbers, so it is not transcribable; what IS transcribable is the boundary
  value problem that produced it, and the example states it plainly: the water
  table falls from el. 14 on the left edge of the domain to the ground surface on
  the right. The .gsz boundary conditions are

      total head = 14   on the left edge (x = 1), from y = -4 up to y = 14
      pressure head = 0 on the flat ground beyond the toe (x = 44 -> 65, y = 8)
                        -- i.e. total head = 8
      potential seepage face on the slope face (x = 36 -> 44)
      no flow on the right edge, the base, and the crest above el. 14
      k_sat = 1e-5 m/s, the same in all three zones

  so those are authored here and XSLOPE solves them with its own steady FE
  seepage (u = 'seep'), writing the nodal field to the _seep.csv sidecar. The
  conductivity is homogeneous, so the head field does not depend on its value.
  This is the one modelling decision in the transcription and the docs entry
  states it: the comparison is of two programs solving the same boundary value
  problem, not of one reading the other's answer.

MESH

  Both models are meshed here rather than by the tag, because a u='seep' model
  must run its stability solve on the SAME mesh the seepage field was written on
  (the field is nodal). tri6 throughout -- quadratic elements are the only ones
  the SSRM path is verified on. The weak clay is a 1 m band that governs the whole
  problem, so it carries a local size override: a global size alone would leave
  one element across it. GeoStudio does the same thing (1 m globally, 0.5 m in
  the weak clay, 0.253 m along the wall).

Run:  PYTHONPATH=. python3 benchmarks/geostudio/build_gs2_wall.py
"""

import os
import sys
import warnings

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from shapely.geometry import LineString, Polygon  # noqa: E402

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx  # noqa: E402
from xslope.mesh import (build_mesh_from_polygons, export_mesh_to_json,  # noqa: E402
                         extract_constraint_line_geometry,
                         extract_size_regions, get_material_polygons)
from xslope.seep import (build_seep_data, export_seep_solution,  # noqa: E402
                         run_seepage_analysis)

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
# A plain string, and every path below is built with os.path.join: the rebuild
# guard redirects a builder's output by REPLACING this attribute with its scratch
# directory, and it hands over a str.
OUT = os.path.join(REPO, 'docs', 'verification', 'files', 'geostudio')
DONOR = os.path.join(REPO, 'docs', 'fem', 'files', 'xslope_griffiths1.xlsx')

# --- geometry ------------------------------------------------------------- #
GROUND = [(1.0, 20.0), (20.0, 20.0), (36.0, 10.0), (38.0, 10.0), (40.0, 10.0),
          (42.0, 9.0), (44.0, 8.0), (65.0, 8.0)]
WEAK_TOP = [(1.0, 6.0), (65.0, 6.0)]
WEAK_BOT = [(1.0, 5.0), (65.0, 5.0)]
BASE_Y = -4.0
X_LEFT, X_RIGHT = 1.0, 65.0

# --- water ---------------------------------------------------------------- #
GAMMA_W = 9.807          # the .gsz UnitWaterWeight
HEAD_FAR = 14.0          # total head on the left edge
HEAD_TOE = 8.0           # zero pressure head on the flat ground beyond the toe
K_SAT = 1.0e-5           # m/s, all three zones

# --- wall ----------------------------------------------------------------- #
WALL = dict(x1=38.0, y1=10.0, x2=38.0, y2=1.0, E=2.0e8, area=0.02, I=0.0005, S=1.0)

# --- mesh ----------------------------------------------------------------- #
TARGET = 1.5             # global tri6 size
WEAK_SIZE = 0.3          # local size across the 1 m weak clay band


def _slope_data(with_wall):
    sd = load_slope_data(str(DONOR))
    sd['unit_system'] = 'metric'
    sd['gamma_water'] = GAMMA_W

    donor_mat = dict(sd['materials'][0])

    def mat(name, **kw):
        m = dict(donor_mat)
        m.update(dict(name=name, gamma=20.0, gamma_sat=20.0, nu=0.4, u='seep',
                      k1=K_SAT, k2=K_SAT, alpha=0.0, kr0=0.001, h0=-1.0,
                      t_cut=None, ru=0.0, phi_b=None, s_cap=None))
        m.update(kw)
        return m

    # t_cut = 0: the tensile strength the example's own material legends print
    # for both Mohr-Coulomb zones. Without it a c'-phi' material can carry
    # tension all the way to the Mohr-Coulomb apex (c'/tan phi' = 34.6 kPa in
    # the sandy clay), which lets a back scarp stand that the vendor's material
    # cannot, and the factor of safety comes out high.
    sd['materials'] = [
        mat('Sandy Clay', option='mc', c=20.0, phi=30.0, E=10000.0, t_cut=0.0),
        mat('Weak Clay', option='mc', c=0.0, phi=10.0, E=5000.0, t_cut=0.0),
        # Isotropic elastic / "Bedrock": no strength, cannot fail.
        mat('OC Soil', option='elastic', c=0.0, phi=0.0, E=500000.0),
    ]

    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': GROUND, 'size': None},
        {'mat_id': 1, 'coords': WEAK_TOP, 'size': WEAK_SIZE},
        {'mat_id': 2, 'coords': WEAK_BOT, 'size': None},
    ]
    sd['max_depth'] = BASE_Y
    sd['ground_surface'] = LineString(GROUND)
    sd['domain_polygon'] = Polygon(GROUND + [(X_RIGHT, BASE_Y), (X_LEFT, BASE_Y)])
    sd['polygons'] = []
    sd['piezo_line'] = []
    sd['piezo_phreatic'] = False
    sd['dloads'] = []
    sd['dloads2'] = []
    sd['circles'] = []
    sd['non_circ'] = []
    sd['reinforce_lines'] = []
    sd['reinforcement_lines'] = []
    sd['k_seismic'] = 0.0
    sd['tcrack_depth'] = 0.0

    sd['seepage_bc'] = {
        'specified_heads': [
            {'head': HEAD_FAR, 'coords': [(X_LEFT, BASE_Y), (X_LEFT, HEAD_FAR)]},
            {'head': HEAD_TOE, 'coords': [(44.0, 8.0), (X_RIGHT, 8.0)]},
        ],
        'exit_face': [(36.0, 10.0), (38.0, 10.0), (40.0, 10.0), (42.0, 9.0),
                      (44.0, 8.0)],
    }

    sd['pile_lines'] = [dict(WALL, H=None, theta_p=None, D_pile=None,
                             V_cap=None, M_cap=None, fixity='free',
                             appl='passive', label='sheet pile wall')
                        ] if with_wall else []
    return sd


def _build_one(name, with_wall):
    path = os.path.join(OUT, f'{name}.xlsx')
    save_slope_data_to_xlsx(_slope_data(with_wall), path)

    # Reload through the writer's own output, then mesh and solve the seepage on
    # that mesh so the stability run reads a field written on its own nodes.
    sd = load_slope_data(path)
    lines, _n_reinf, _n_pile = extract_constraint_line_geometry(sd)
    polygons = get_material_polygons(sd, reinf_lines=lines)
    mesh = build_mesh_from_polygons(polygons, target_size=TARGET,
                                    element_type='tri6', lines=lines,
                                    size_regions=extract_size_regions(sd))
    export_mesh_to_json(mesh, os.path.join(OUT, f'{name}_mesh.json'))
    seep = build_seep_data(mesh, sd)
    sol = run_seepage_analysis(seep, tol=1e-5, max_iter=400)
    export_seep_solution(seep, sol, os.path.join(OUT, f'{name}_seep.csv'))
    print(f'  {name}: {len(mesh["nodes"])} nodes, {len(mesh["elements"])} elements, '
          f'{len(mesh.get("elements_1d", []))} beam elements')
    return f'{name}.xlsx'


def gs2_wall_none():
    """The slope alone. SIGMA/W SRS ~1.025; its FE stability run 1.035."""
    return _build_one('gs2_wall_none', with_wall=False)


def gs2_wall():
    """The same slope with the continuous sheet pile wall. SIGMA/W SRS ~1.4."""
    return _build_one('gs2_wall', with_wall=True)


BUILDERS = [gs2_wall_none, gs2_wall]


if __name__ == '__main__':
    for b in BUILDERS:
        print('wrote', b())
