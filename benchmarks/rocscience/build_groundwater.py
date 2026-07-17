"""Builders for the Rocscience Slide2 Groundwater verification corpus
(docs/verification/rocscience_groundwater.md). Steady-state problems only;
the transient tier (GW#15-21) waits on a transient solver.

Unlike the LEM corpus, these tags run the seepage solver live (run_tests.py
type=seep / seep_head mesh and solve per run), so the committed artifact is
the xlsx alone - no mesh/solution sidecars.

Run from the repo root:  python benchmarks/rocscience/build_groundwater.py
"""

import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), '..', '..',
                   'docs', 'files', 'rocscience_gw')
ACADS_1A = os.path.join(os.path.dirname(__file__), '..', '..',
                        'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')


def _base_sd(k1=1e-5, kr0=1e-3, h0=-0.4):
    """Seepage-only base model: one material with u='seep'; the LEM circle is
    a placeholder (these problems are never solved for a factor of safety)."""
    sd = load_slope_data(ACADS_1A)
    m = dict(sd['materials'][0])
    m.update(name='Soil', c=1.0, phi=30.0, gamma=20.0, gamma_sat=20.0,
             option='mc', u='seep', k1=k1, k2=k1, alpha=0.0, kr0=kr0, h0=h0)
    sd['materials'] = [m]
    sd['gamma_water'] = 9.81
    sd['dloads'] = []
    sd['piezo_line'] = []
    sd['circular'] = True
    sd['non_circ'] = []
    sd['circles'] = [{'Xo': 5.0, 'Yo': 10.0, 'Depth': 0.0, 'R': 10.0}]
    return sd


def gw001():
    """GW#1: shallow unconfined flow with rainfall (Haar 1990; Dupuit). Flow
    between two parallel rivers 10 m apart with a rainfall recharge on the land
    surface, the manual's first specified-flux problem after the flux BC landed.
    Domain 10 x 5 m (vendor RS2 groundwater #001_01: base y=0 impermeable, ground
    surface y=5); river heads h1=3.75 on the left edge (0,0)-(0,3.75) and h2=3.0
    on the right (10,0)-(10,3.0), exactly as the RS2 nodal 'total head' cards;
    uniform rainfall P=2.5e-6 m/s on the top ('vert infilt' 2.5e-6), k=1e-5 m/s.
    The whole recharge mounds an INTERNAL free surface between the rivers - there
    is no daylighting seepage face, so the top edge is declared BOTH the flux and
    an exit face: it stays inactive (dry, the mound crest el ~4.6 sits below the
    y=5 surface) but switches the solver onto the unsaturated free-surface path,
    which the internal mound needs. Targets (Table 1.1): the free-surface crest
    (x_a, h_max) - Slide 4.06 / 4.49, Haar eqs 1.2-1.3 3.98 / 4.25. xslope reads
    x_a~=4.1, h_max~=4.61 (a touch above Slide, the same free-surface-family bias
    the SEEP2D cross-check documents); Q=P*L=2.5e-5 m3/s per m is exact by
    construction. The free surface is k- and unsat-model independent (mass balance
    sets it); the flowrate lock plus a mound-guarding head regression are taken."""
    sd = _base_sd(k1=1e-5)
    from shapely.geometry import Polygon
    sd['profile_lines'] = []
    sd['polygons'] = [{'mat_id': 0, 'polygon': Polygon(
        [(0.0, 0.0), (10.0, 0.0), (10.0, 5.0), (0.0, 5.0)])}]
    sd['max_depth'] = None
    sd['seepage_bc'] = {
        'specified_heads': [
            {'head': 3.75, 'coords': [(0.0, 0.0), (0.0, 3.75)]},
            {'head': 3.0, 'coords': [(10.0, 0.0), (10.0, 3.0)]},
        ],
        'specified_fluxes': [
            {'flux': 2.5e-6, 'coords': [(0.0, 5.0), (10.0, 5.0)]},
        ],
        # Same polyline as the flux: an internal recharge mound has no daylight
        # face, so the top doubles as an (inactive) exit face to select the
        # unsaturated solver. It never activates - the mound crest is below y=5.
        'exit_face': [(0.0, 5.0), (10.0, 5.0)],
    }
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'gw001.xlsx'))
    return 'gw001.xlsx'


def gw007():
    """GW#7: seepage within a layered slope (Rulon & Freeze sandbox, after
    Fredlund & Rahardjo 1993). A medium-sand slope with a thin fine-sand lens;
    the fine sand's ks (5.5e-5 m/s) is BELOW the applied rainfall (2.1e-4 m/s),
    so recharge perches on the lens and daylights as a seepage face - the classic
    perched-water-table result. Geometry from vendor RS2 groundwater #007 (a 2.4 x
    1.0 m box, 2:1 downstream slope from crest (1.6,1.0) to toe (0,0.2)); the fine
    lens is the y=0.6-0.7 band from the slope face to the right wall. Conductivity
    functions are the RS2 'Custom' curves (Fig 7.2), fit by Mualem-vG here: medium
    ks=0.0014 (vg_a=1.774, vg_n=2.328), fine ks=5.5e-5 (vg_a=1.672, vg_n=2.197).
    BCs from the RS2 cards: rainfall 2.1e-4 'vert infilt' on the crest top
    (1.6-2.4, y=1.0); tailwater 'total head' 0.3 on the submerged toe (0,0.2)-
    (0.2,0.3); the rest of the slope an exit (seepage) face; base and right wall
    no-flow. All published targets are chart curves (Fig 7.4 water table, 7.7/7.8
    head profiles) with no tabulated value, so - as the methodology note allows for
    GW6/GW7 - only the flowrate is locked (Q=q*L=1.68e-4, exact by construction),
    with a head regression guarding the field. xslope reproduces the stated water
    table (daylights at el 0.30 at the toe) and the perched zone above the lens."""
    sd = _base_sd()
    med = dict(sd['materials'][0])
    med.update(name='Medium sand', k1=0.0014, k2=0.0014, alpha=0.0,
               kr0=1e-3, h0=-0.4, unsat='vg', vg_a=1.7745, vg_n=2.3276)
    fin = dict(med)
    fin.update(name='Fine sand', k1=5.5e-5, k2=5.5e-5,
               unsat='vg', vg_a=1.6722, vg_n=2.1965)
    sd['materials'] = [med, fin]                # mat 0 = medium, mat 1 = fine lens
    from shapely.geometry import Polygon
    sd['profile_lines'] = []
    sd['polygons'] = [
        {'mat_id': 0, 'polygon': Polygon(
            [(0.0, 0.0), (2.4, 0.0), (2.4, 0.6), (0.8, 0.6), (0.0, 0.2)])},
        {'mat_id': 1, 'polygon': Polygon(
            [(0.8, 0.6), (2.4, 0.6), (2.4, 0.7), (1.0, 0.7)])},
        {'mat_id': 0, 'polygon': Polygon(
            [(1.0, 0.7), (2.4, 0.7), (2.4, 1.0), (1.6, 1.0)])},
    ]
    sd['max_depth'] = None
    sd['seepage_bc'] = {
        'specified_heads': [{'head': 0.3, 'coords': [(0.0, 0.2), (0.2, 0.3)]}],
        'specified_fluxes': [{'flux': 2.1e-4, 'coords': [(1.6, 1.0), (2.4, 1.0)]}],
        'exit_face': [(0.2, 0.3), (0.8, 0.6), (1.0, 0.7), (1.6, 1.0)],
    }
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'gw007.xlsx'))
    return 'gw007.xlsx'


def gw002():
    """GW#2: confined potential flow around a cylinder (Streeter analytical;
    Desai & Kundu FE). Half-domain 8 x 4 m with a 24-segment semicircular
    notch (r=1 at (4,0)); heads 1.0 / 0.0 on the left/right edges (Fig 2.2,
    all dimensions printed); fully saturated -> confined solve. Targets:
    heads at the 5 printed points - xslope matches Slide within 0.0013 m
    everywhere and the closed form within its own idealization error.
    Homogeneous confined flow: the head field is k-independent."""
    sd = _base_sd()
    from shapely.geometry import Polygon
    arc = [(4.0 - math.cos(t), math.sin(t))
           for t in [math.pi * i / 24 for i in range(25)]]
    coords = ([(0.0, 0.0), (3.0, 0.0)] + arc[1:-1] + [(5.0, 0.0), (8.0, 0.0)])
    sd['profile_lines'] = []
    sd['polygons'] = [{'mat_id': 0,
                       'polygon': Polygon(coords + [(8.0, 4.0), (0.0, 4.0)])}]
    sd['max_depth'] = None
    sd['seepage_bc'] = {
        'specified_heads': [
            {'head': 1.0, 'coords': [(0.0, 0.0), (0.0, 4.0)]},
            {'head': 0.0, 'coords': [(8.0, 0.0), (8.0, 4.0)]},
        ],
        'exit_face': [],
    }
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'gw002.xlsx'))
    return 'gw002.xlsx'


def gw003():
    """GW#3: confined flow under a dam foundation (Rushton & Redshaw). Soil
    block 40 x 10 m; head 5 on the ground x 0-8, head 0 on x 20-40, dam base
    x 8-20 impervious (no BC). Targets: head profiles along line 1-1 (y=-4)
    and line 2-2 (x=20) from Fig 3.5/3.6 - xslope within 0.08 m of the
    chart (Slide's markers coincide with Rushton & Redshaw's)."""
    sd = _base_sd()
    from shapely.geometry import Polygon
    sd['profile_lines'] = []
    sd['polygons'] = [{'mat_id': 0, 'polygon': Polygon(
        [(0.0, -10.0), (40.0, -10.0), (40.0, 0.0), (0.0, 0.0)])}]
    sd['max_depth'] = None
    sd['seepage_bc'] = {
        'specified_heads': [
            {'head': 5.0, 'coords': [(0.0, 0.0), (8.0, 0.0)]},
            {'head': 0.0, 'coords': [(20.0, 0.0), (40.0, 0.0)]},
        ],
        'exit_face': [],
    }
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'gw003.xlsx'))
    return 'gw003.xlsx'


def gw004():
    """GW#4: steady unconfined flow through an earth dam with a toe drain
    (Kozeny basic parabola). Dam (0,0)-(10,5)-(12.5,5)-(22,1.5)-(22,0)
    (Fig 4.2 dimensions), upstream head 4 below el 4, exit face on the
    downstream face and drain end; k arbitrary (the free-surface shape is
    k-independent). The solved phreatic matches the parabola within 1-2%
    over the dam body; at the drain tip y1 = 0.50 vs Slide's measured 0.442
    and the parabola's 0.480 - the published pair itself spreads 9%, and
    the parabola is an idealization exact only at the drain."""
    sd = _base_sd(h0=-0.25)
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 0.0), (10.0, 5.0), (12.5, 5.0),
                                 (22.0, 1.5), (22.0, 0.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['seepage_bc'] = {
        'specified_heads': [
            {'head': 4.0, 'coords': [(0.0, 0.0), (8.0, 4.0)]},
        ],
        'exit_face': [(12.5, 5.0), (22.0, 1.5), (22.0, 0.0)],
    }
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'gw004.xlsx'))
    return 'gw004.xlsx'


def gw009a():
    """GW#9 dam 1: Bowles' homogeneous dam via Chapuis et al. (2001). Base
    100 m, crest 10 m at el 20 (upstream 2.5:1, downstream 2:1), reservoir
    head 18.5; ks = 6.67e-6 m/s with the manual's printed 8-point k-function
    fit by Mualem-van Genuchten (vg_a=0.2835, vg_n=2.765, log-space, max
    deviation ~1.5x mid-range). Q = 1.379e-3 m3/(min*m) vs Slide 1.378e-3,
    SEEP/W fine-mesh 1.37e-3, Bowles' flow nets 1.10-1.28e-3. Dam 2 (with
    drain) is not built: its k-function and reservoir level are chart-only
    and the published Q implies a k two decades below the chart - needs
    Chapuis et al. (2001), CGJ 38:1113."""
    sd = _base_sd(k1=6.67e-6)
    m = sd['materials'][0]
    m.update(name='Dam fill', c=10.0, phi=30.0, kr0=0.0, h0=0.0,
             unsat='vg', vg_a=0.2835, vg_n=2.765)
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 0.0), (50.0, 20.0), (60.0, 20.0),
                                 (100.0, 0.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['seepage_bc'] = {
        'specified_heads': [
            {'head': 18.5, 'coords': [(0.0, 0.0), (46.25, 18.5)]},
        ],
        'exit_face': [(60.0, 20.0), (100.0, 0.0)],
    }
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'gw009a.xlsx'))
    return 'gw009a.xlsx'




def gw010():
    """GW#10: Clement, Wise, Molz & Wen (1996) unconfined square domain with
    van Genuchten conductivity (vg_a=0.64, vg_n=4.65, ks=1.1574e-5 m/s):
    10 x 10 m block, head 10 on the left edge, tailwater 2 on the right,
    exit face above the tailwater, no-flow top and base. Q = 6.070e-5 vs
    Slide 6.066e-5 (+0.07%) / Clement 6.076e-5 (-0.1%); phreatic exit el
    4.87 vs Clement 4.8 / Slide 5.0 (the manual's "seepage face" column is
    the exit ELEVATION, not a face length). Only the tailwater-2 case has
    published numbers."""
    sd = _base_sd()
    sd['materials'][0].update(u='none', k1=1.1574e-5, k2=1.1574e-5, alpha=0.0,
                              unsat='vg', vg_a=0.64, vg_n=4.65,
                              kr0=0.001, h0=-1.0)
    sd['profile_lines'] = [{'mat_id': 0, 'coords': [(0.0, 10.0), (10.0, 10.0)]}]
    sd['max_depth'] = 0.0
    sd['seepage_bc'] = {
        'specified_heads': [
            {'head': 10.0, 'coords': [(0.0, 0.0), (0.0, 10.0)]},
            {'head': 2.0, 'coords': [(10.0, 0.0), (10.0, 2.0)]},
        ],
        'exit_face': [(10.0, 2.0), (10.0, 10.0)],
    }
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'gw010.xlsx'))
    return 'gw010.xlsx'


def gw012():
    """GW#12: Vedernikov's trapezoidal ditch seeping into a deep drainage
    layer, half-model by symmetry: ground (0,40)-(15,40)-(25,50)-(100,50),
    head 50 on the ditch perimeter, head 0 on the base (the deep drain),
    symmetry/right edges no-flow. Q_half = 4.137e-4 vs Slide 4.093e-4
    (+1.1%) / Vedernikov k(B+AH)/2 = 4.0e-4 (+3.4%); flow-bulb half-width
    ~42 vs Slide 41 / theory 40. The detached bulb converges at
    max_iter=1500 (tag key), not the default 400."""
    sd = _base_sd()
    sd['materials'][0].update(u='none', k1=1e-5, k2=1e-5, alpha=0.0,
                              unsat='lf', kr0=0.001, h0=-1.0)
    sd['profile_lines'] = [{'mat_id': 0, 'coords': [(0.0, 40.0), (15.0, 40.0),
                                                    (25.0, 50.0), (100.0, 50.0)]}]
    sd['max_depth'] = 0.0
    sd['seepage_bc'] = {
        'specified_heads': [
            {'head': 50.0, 'coords': [(0.0, 40.0), (15.0, 40.0), (25.0, 50.0)]},
            {'head': 0.0, 'coords': [(0.0, 0.0), (100.0, 0.0)]},
        ],
        'exit_face': [(25.0, 50.0), (100.0, 50.0)],
    }
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'gw012.xlsx'))
    return 'gw012.xlsx'


def gw013():
    """GW#13: Vedernikov's triangular ditch into a deep drainage layer,
    half-model: apex (0,40) to (10,50), k=1e-3. Q_half = 2.087e-2 vs Slide
    2.050e-2 (+1.8%) / Vedernikov 2.0e-2 (+4.3%). max_iter=1500 as GW#12."""
    sd = _base_sd()
    sd['materials'][0].update(u='none', k1=1e-3, k2=1e-3, alpha=0.0,
                              unsat='lf', kr0=0.001, h0=-1.0)
    sd['profile_lines'] = [{'mat_id': 0, 'coords': [(0.0, 40.0), (10.0, 50.0),
                                                    (100.0, 50.0)]}]
    sd['max_depth'] = 0.0
    sd['seepage_bc'] = {
        'specified_heads': [
            {'head': 50.0, 'coords': [(0.0, 40.0), (10.0, 50.0)]},
            {'head': 0.0, 'coords': [(0.0, 0.0), (100.0, 0.0)]},
        ],
        'exit_face': [(10.0, 50.0), (100.0, 50.0)],
    }
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'gw013.xlsx'))
    return 'gw013.xlsx'




def gw011():
    """GW#11 case 1: uniform earth/rock-fill dam with the Gardner
    kr = 1/(1 + a*psi^n) conductivity law, after Zhang, Xu & Chen (2001).
    Geometry from Fig 11.1 (p.45, printed dimensions): 45 m high, 17 m crest,
    upstream run 89.10 (1:1.98), downstream run 76.90 -> base 183.0 m. (The
    text's "1:1.171" downstream is a typo: 76.90/45 = 1.709, and the printed
    76.90 governs.) Reservoir at el 40 (upstream elevation head 40, downstream
    0); exit face on the whole downstream slope; impermeable base.
    Gardner a = 0.15, n = 6 exactly as printed on p.45. ks is not printed for
    the uniform dam and does not matter: the dam is homogeneous, so the free
    surface (the published target) is k-independent - 1e-7 m/s is the dam value
    from Table 11.1 (p.46) for the companion non-homogeneous case.
    Target: release point on the downstream face - Slide 19.397 m,
    ABAQUS/Zhang 19.64 m (p.46)."""
    sd = _base_sd(k1=1e-7)
    sd['materials'][0].update(name='Dam fill', kr0=0.0, h0=0.0,
                              unsat='gard', vg_a=0.15, vg_n=6.0)
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 0.0), (89.10, 45.0), (106.10, 45.0),
                                 (183.0, 0.0)]},
    ]
    sd['max_depth'] = 0.0
    sd['seepage_bc'] = {
        'specified_heads': [
            {'head': 40.0, 'coords': [(0.0, 0.0), (79.20, 40.0)]},
        ],
        'exit_face': [(106.10, 45.0), (183.0, 0.0)],
    }
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'gw011.xlsx'))
    return 'gw011.xlsx'


def gw006a():
    """GW#6 case 1: Fredlund & Rahardjo (1993) isotropic earth dam with a
    12 m horizontal drain (Fig 6.1: 12 m high, 4 m crest, symmetric 2:1
    faces, base 52 m; reservoir 10 m). k-function digitized from the Fig 6.2
    chart (ks=1e-7 m/s, air entry ~9 kPa, 3 decades down by ~97 kPa) and fit
    by Mualem-vG (vg_a=0.4029, vg_n=1.9156; max deviation 0.36 decades in
    the controlling 1-6 m suction band - three fit variants move the answer
    <0.03 m). Target: pressure head along line 1-1 (the crest centerline,
    x=26) from Fig 6.6, where Slide and F&R coincide. xslope reads a
    k-fit- and mesh-insensitive +0.25-0.5 m above the published curves -
    the same free-surface family as task #30; locked at xslope's own values
    with the caveat stated. Case 4 (infiltration) is blocked on the flux BC
    (task #28); cases 2 (9:1 anisotropy), 3 (core), 5 (seepage-face variant)
    are buildable-deferred with chart targets like this one."""
    sd = _base_sd(k1=1e-7)
    m = sd['materials'][0]
    m.update(name='Dam fill', c=10.0, phi=30.0, kr0=0.0, h0=0.0,
             unsat='vg', vg_a=0.4029, vg_n=1.9156)
    sd['profile_lines'] = [{'mat_id': 0, 'coords': [(0.0, 0.0), (24.0, 12.0),
                                                    (28.0, 12.0), (52.0, 0.0)]}]
    sd['max_depth'] = 0.0
    sd['seepage_bc'] = {
        'specified_heads': [{'head': 10.0, 'coords': [(0.0, 0.0), (20.0, 10.0)]}],
        'exit_face': [(40.0, 0.0), (52.0, 0.0)],
    }
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'gw006a.xlsx'))
    return 'gw006a.xlsx'


def gw008():
    """GW#8: flow through ditch-drained soils (Gureghian 1981), the corpus'
    first specified-flux (Neumann) problem. Fig 8.1 (p.34): half-drain
    spacing 1.0 m wide, 0.5 m deep to the impermeable base, two layers -
    Soil A is the LOWER 0.1 m (coarse, k = 1.11e-3 m/s, Gardner a = 1000,
    n = 4.5) and Soil B the upper 0.4 m (fine, k = 1.11e-4 m/s, a = 2777.7,
    n = 4.2), both from Table 8.1 (p.34). Rainfall infiltration 4.4e-6 m/s
    on the top boundary as a specified flux (positive = inflow, Fig 8.2
    p.35). The water-free ditch is the left wall: exit face over its full
    height, so the invert node carries zero head when active (Slide draws
    the same thing as "seepage face" + "zero head" at the invert). Base and
    right-hand symmetry edge are no-flow, i.e. simply unspecified.
    Targets are chart-only (the manual prints no point value and no
    discharge): pressure head -0.10 to -0.20 m in the unsaturated zone
    (Fig 8.3) and total head 0.05 to 0.29 m (Fig 8.4)."""
    sd = _base_sd()
    soil_b = dict(sd['materials'][0])
    soil_b.update(name='Soil B', k1=1.11e-4, k2=1.11e-4, alpha=0.0,
                  kr0=0.0, h0=0.0, unsat='gard', vg_a=2777.7, vg_n=4.2)
    soil_a = dict(soil_b)
    soil_a.update(name='Soil A', k1=1.11e-3, k2=1.11e-3,
                  unsat='gard', vg_a=1000.0, vg_n=4.5)
    sd['materials'] = [soil_b, soil_a]          # mat 0 = upper, mat 1 = lower
    sd['profile_lines'] = [
        {'mat_id': 0, 'coords': [(0.0, 0.5), (1.0, 0.5)]},
        {'mat_id': 1, 'coords': [(0.0, 0.1), (1.0, 0.1)]},
    ]
    sd['max_depth'] = 0.0
    sd['seepage_bc'] = {
        'specified_heads': [],
        'specified_fluxes': [
            {'flux': 4.4e-6, 'coords': [(0.0, 0.5), (1.0, 0.5)]},
        ],
        'exit_face': [(0.0, 0.0), (0.0, 0.5)],
    }
    save_slope_data_to_xlsx(sd, os.path.join(OUT, 'gw008.xlsx'))
    return 'gw008.xlsx'


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    for fn in (gw001, gw002, gw003, gw004, gw006a, gw007, gw008, gw009a,
               gw010, gw011, gw012, gw013):
        print(fn())
