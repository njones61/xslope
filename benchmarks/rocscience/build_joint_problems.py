"""Builders for the RS2 JOINT-analysis verification corpus (docs/verification/rs2_joints.md).

Rocscience's *RS2 Joint Verification* manual carries 23 problems on jointed rock:
block and flexural toppling, plane failure, Alejano's sliding and ploughing slabs,
step-path failure, a Voronoi-tessellated mass, a jointed tunnel, and two problems
that exercise the joint's own constitutive law rather than a slope. Twenty-one of
them report a factor of safety (or, for problem 16, a tilt angle); problems 22 and
23 report neither, being shear-box tests of the joint model itself.

**Naming.** One file per problem, ``rjNNN.xlsx`` with the manual's own problem
number — ``rj018.xlsx`` is problem 18 — and a letter suffix where the manual
carries lettered cases (``rj001a``…``rj001d``). Variants the manual does not
letter take a descriptive suffix (``rj016_3deg``). The ``rj`` prefix separates
this corpus from the ``vpNNN`` Slide2 one and the ``rs2_NN`` RS2 slope-stability
one, both of which are numbered by a different manual.

**Where the inputs come from.** Every geometry, material, joint property and
restraint is read from the vendor's own ``.fez`` model, not from the manual's
tables: the tables have known errata (problem 7's slope angle, problem 19's joint
inclination), and the model is what RS2 actually solved. The manual supplies the
referee — the closed form or the UDEC run each problem is scored against — and
RS2's own two reported factors. The vendor files are NOT redistributed; the page
records what was read from them.

**Units.** The vendor models are Metric MPa with unit weight in MN/m^3. These
files are written in metric kPa / kN/m^3, which is the unit system the rest of the
corpus uses, so every stress is multiplied by 1000 and every unit weight by 1000.

Run this script to regenerate every built problem.
"""

import os
import sys
import warnings

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.dirname(__file__))

from shapely.geometry import LineString, Polygon                    # noqa: E402

from xslope.fileio import load_slope_data                           # noqa: E402
from xslope.fileio import save_slope_data_to_xlsx as _write_xlsx    # noqa: E402
from xslope.fileio import build_ground_surface_from_polygons        # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', '..',
                   'docs', 'verification', 'files', 'rocscience', 'joints')

#: A shipped FEM model, read only for the boilerplate every slope_data carries —
#: the unit system, the water unit weight, the solver options. Every builder here
#: replaces the geometry, the materials and the lines outright.
DONOR = os.path.join(os.path.dirname(__file__), '..', '..',
                     'docs', 'fem', 'files', 'xslope_griffiths1.xlsx')

#: The vendor's joint stiffnesses, in the files' own MPa/m, and what they are in
#: the kPa/m these files carry. Kn = 100 000 MPa/m and Ks = 10 000 MPa/m are the
#: pair every problem in the set uses except 8, 15, 16 and 21.
KN_STD = 1.0e8
KS_STD = 1.0e7


def _base():
    """A loaded model stripped to its boilerplate."""
    sd = load_slope_data(DONOR)
    sd['unit_system'] = 'metric'
    sd['profile_lines'] = []
    sd['circles'] = []
    sd['non_circ'] = []
    sd['piezo_line'] = []
    sd['piezo_phreatic'] = False
    sd['dloads'] = []
    sd['dloads2'] = []
    sd['reinforcement_lines'] = []
    sd['reinforce_lines'] = []
    sd['pile_lines'] = []
    sd['joint_lines'] = []
    sd['max_depth'] = None
    sd['k_seismic'] = 0.0
    return sd


def _rock(name, gamma, E, nu, c, phi, t_cut=0.0):
    """One Mohr-Coulomb rock, in kPa and kN/m^3."""
    return {'name': name, 'gamma': gamma, 'gamma_sat': gamma,
            'option': 'mc', 'c': c, 'phi': phi, 'psi': 0.0, 'r_elev': 0.0,
            'u': 'none', 'ru': 0.0, 'E': E, 'nu': nu, 't_cut': t_cut,
            'sigma_gamma': 0.0, 'sigma_c': 0.0, 'sigma_phi': 0.0,
            'sigma_cp': 0.0, 'sigma_d': 0.0, 'sigma_psi': 0.0}


def _joint(label, p1, p2, c, phi, kn=KN_STD, ks=KS_STD, t_cut=0.0,
           c_res=None, phi_res=None, dil=None, jred=''):
    """One row of the joints sheet, in kPa."""
    nan = float('nan')
    return {'label': label,
            'x1': float(p1[0]), 'y1': float(p1[1]),
            'x2': float(p2[0]), 'y2': float(p2[1]),
            'c': c, 'phi': phi,
            'c_res': nan if c_res is None else c_res,
            'phi_res': nan if phi_res is None else phi_res,
            'dil': nan if dil is None else dil,
            't_cut': t_cut, 'kn': kn, 'ks': ks, 'jred': jred}


def _finish(sd, rings_and_ids, materials):
    """Install the geometry and derive the surface the loader would."""
    sd['materials'] = materials
    sd['polygons'] = [{'polygon': Polygon(r), 'mat_id': i}
                      for r, i in rings_and_ids]
    gs, dom = build_ground_surface_from_polygons(sd['polygons'])
    sd['ground_surface'], sd['domain_polygon'] = gs, dom
    return sd


def _write(sd, name):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name)
    _write_xlsx(sd, path)
    return name


# ---------------------------------------------------------------------------
# Problem 18 — step-path failure through continuous joints
# ---------------------------------------------------------------------------

def rj018():
    """RJ-18 — step-path failure, continuous joints (vendor `joint #018.fez`).

    A 45 x 20 m section with a slope face rising from (17, 8.2) to (26.9, 20),
    cut by three parallel joints at 36.1 degrees running from the face to the
    crest at a perpendicular spacing of 0.883 m. The rock is one Mohr-Coulomb
    material (gamma 19.62 kN/m^3, E 20 GPa, nu 0.3, c 25 kPa, phi 25, no tensile
    capacity); the joints carry c = 1 kPa, phi = 35, Kn = 1e8 kPa/m,
    Ks = 1e7 kPa/m, and are reduced with the rock in the strength reduction
    (RS2 `CoupledSSR`). Referee: UDEC 1.01. RS2 reports 1.01 without joint
    improvement and 1.00 with it.

    The external boundary carries the three joints' upper and lower endpoints as
    vertices of its own, which is how the vendor file states them and what lets
    each joint end ON the boundary rather than a hair inside it.
    """
    sd = _base()
    ring = [(45.0, 0.0), (45.0, 20.0), (33.2, 20.0), (31.7, 20.0),
            (30.2, 20.0), (30.0, 20.0), (26.9, 20.0),
            (21.7143, 13.819), (19.3571, 11.0095), (17.0, 8.2),
            (15.5, 8.2), (14.0, 8.2), (0.0, 8.2), (0.0, 0.0)]
    mats = [_rock('Rock', 19.62, 2.0e7, 0.3, 25.0, 25.0, t_cut=0.0)]
    _finish(sd, [(ring, 0)], mats)
    sd['joint_lines'] = [
        _joint('joint-1', (17.0, 8.2), (33.2, 20.0), 1.0, 35.0),
        _joint('joint-2', (19.3571, 11.0095), (31.7, 20.0), 1.0, 35.0),
        _joint('joint-3', (21.7143, 13.819), (30.2, 20.0), 1.0, 35.0),
    ]
    # Inert for a strength reduction, and there so the file opens as a complete
    # model: the toe circle the starting-circle rule gives on this face.
    sd['circles'] = [{'Xo': 21.95, 'Yo': 31.8, 'R': 24.11, 'Depth': 31.8 - 24.11}]
    return _write(sd, 'rj018.xlsx')


# ---------------------------------------------------------------------------
# Problem 19 — bi-planar step-path failure
# ---------------------------------------------------------------------------

def rj019():
    """RJ-19 — bi-planar step-path failure (vendor `joint #019.fez`).

    A 120 x 70 m section with a slope face from (30, 20) to (60, 70), cut by two
    discontinuous joints with a rock bridge between them: a basal joint from
    (39.0149, 35.0248) to (63, 48) at 28.4 degrees, and an upper joint from
    (62, 49) to (76, 70) at 56.3 degrees. The rock is Mohr-Coulomb
    (gamma 27 kN/m^3, E 20 GPa, nu 0.3, c 10 500 kPa, phi 35, tensile capacity
    200 kPa); the joints carry c = 0, phi = 40 and the standard stiffness pair.
    Referee: UDEC 1.46. RS2 reports 1.50 without joint improvement and 1.41 with
    it.

    The manual's own table for this problem states one joint inclination as 59
    degrees; the figure dimensions 56 and 28, and the vendor model's endpoints
    give 56.3 and 28.4. The model is what is transcribed.
    """
    sd = _base()
    ring = [(120.0, 0.0), (120.0, 70.0), (76.0, 70.0), (60.0, 70.0),
            (39.0149, 35.0248), (30.0, 20.0), (0.0, 20.0), (0.0, 0.0)]
    mats = [_rock('Rock', 27.0, 2.0e7, 0.3, 10500.0, 35.0, t_cut=200.0)]
    _finish(sd, [(ring, 0)], mats)
    sd['joint_lines'] = [
        _joint('basal', (39.0149, 35.0248), (63.0, 48.0), 0.0, 40.0),
        _joint('upper', (62.0, 49.0), (76.0, 70.0), 0.0, 40.0),
    ]
    sd['circles'] = [{'Xo': 45.0, 'Yo': 120.0, 'R': 101.12,
                      'Depth': 120.0 - 101.12}]
    return _write(sd, 'rj019.xlsx')


#: Every builder in this module, in problem-number order. ``verify_rebuild.py``'s
#: ``joints`` group is this list, so a builder missing here is a corpus file
#: nothing guards.
BUILDERS = [rj018, rj019]


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    for b in BUILDERS:
        print('built', b())
