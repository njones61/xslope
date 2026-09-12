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

import math
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
from xslope.joints import cross_jointed, parallel_set               # noqa: E402
from benchmarks.tag_k0 import apply_tag_k0                          # noqa: E402

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
    # The donor carries the imperial unit weight of water; these files are metric
    # and none of them has water, so the value is inert but it must not be wrong.
    sd['gamma_water'] = 9.81
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
    # RS2 restrains the SIDES of every model in this manual in both directions,
    # not just in x: the restraint list of each .fez names every node on the base
    # and on both sides with tx = ty = 0, and mapping those node numbers back to
    # coordinates puts all of them on those three faces and nowhere else. XSLOPE's
    # default side restraint is a roller, so the files declare the vendor's.
    # (Problem 21, the tunnel, is the one model in the set that uses rollers.)
    sd['side_bc'] = 'fixed'
    return sd


def _rock(name, gamma, E, nu, c, phi, t_cut=0.0, option='mc'):
    """One rock, in kPa and kN/m^3.

    ``option='elastic'`` is a rock the FEM holds out of the strength reduction
    entirely — it cannot yield, so only the joints can fail. That is what the
    vendor's ``Plasticity Specifications: Non`` says, and several of these
    problems carry it (see :func:`rj002`).
    """
    # gamma_sat is left unset: none of these models has water, and a saturated
    # unit weight on a dry model is a value nothing reads.
    return {'name': name, 'gamma': gamma, 'gamma_sat': float('nan'),
            'option': option, 'c': c, 'phi': phi, 'psi': 0.0, 'r_elev': 0.0,
            'u': 'none', 'ru': 0.0, 'E': E, 'nu': nu, 't_cut': t_cut,
            'sigma_gamma': 0.0, 'sigma_c': 0.0, 'sigma_phi': 0.0,
            'sigma_cp': 0.0, 'sigma_d': 0.0, 'sigma_psi': 0.0}


def _offset_through(dip_deg, x, y):
    """The ``parallel_set`` offset that puts one trace through ``(x, y)``.

    ``offset`` is measured along the set's own normal, and ``offset=0`` is the
    trace through the origin, so the offset a stated location asks for is that
    location projected onto the normal. The vendor files state a set as a dip, a
    spacing and an ``init_joint_loc``, which is exactly this triple.
    """
    th = math.radians(float(dip_deg))
    nx, ny = -math.sin(th), math.cos(th)
    return float(x) * nx + float(y) * ny


def _toe_circle(toe, crest):
    """A starting circle for a section whose mechanism is not a single plane.

    Centred above mid-slope at the toe elevation plus twice the slope height and
    passing through the toe, which is the starting circle a user would draw on
    this section. It is what makes the file a complete slope-stability model;
    a strength reduction never reads it.
    """
    (tx, ty), (cx, cy) = toe, crest
    xo = 0.5 * (tx + cx)
    yo = ty + 2.0 * (cy - ty)
    r = ((xo - tx) ** 2 + (yo - ty) ** 2) ** 0.5
    # Depth is the elevation of the circle's lowest point, and the loader reads
    # it in preference to R, so the two have to agree or the stated radius is
    # replaced by one that reaches elevation zero.
    return [{'Xo': float(xo), 'Yo': float(yo), 'Depth': float(yo - r),
             'R': float(r)}]


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


def _surface(points):
    """A non-circular failure surface through the stated points.

    The two ends are Free — they are on the ground surface and the slicer finds
    where — and every point between them is Fixed. Each carries an explicit Y,
    which is what a Free end needs.
    """
    n = len(points)
    return [{'X': float(x), 'Y': float(y),
             'Movement': 'Free' if i in (0, n - 1) else 'Fixed'}
            for i, (x, y) in enumerate(points)]


def _finish(sd, rings_and_ids, materials):
    """Install the geometry and derive the surface the loader would."""
    sd['materials'] = materials
    sd['polygons'] = [{'polygon': Polygon(r), 'mat_id': i}
                      for r, i in rings_and_ids]
    gs, dom = build_ground_surface_from_polygons(sd['polygons'])
    sd['ground_surface'], sd['domain_polygon'] = gs, dom
    return sd


def _write(sd, name):
    """Write one corpus file, declaring the K0 its own test tag names.

    The initial stress has to live in the FILE, not only in the tag, or the
    locked factor is reproducible from the suite and from nothing a user would
    open. ``apply_tag_k0`` reads the page's own tag and clears K0 where no tag
    names one, which is what stops a donor's value riding into a problem that
    never asked for it.
    """
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name)
    apply_tag_k0(sd, path)
    _write_xlsx(sd, path)
    return name


#: The section problems 3 to 7 share: a 700 x 400 m block with a 260 m face at
#: 55 degrees running from the crest at (377.786, 400) down to the toe at
#: (560, 140). The four problems differ only in the joint network cut into it.
LV_RING = [(0.0, 0.0), (700.0, 0.0), (700.0, 140.0), (560.0, 140.0),
           (377.785587105239, 400.0), (0.0, 400.0)]

#: The joint the same four problems carry: c = 0.1 MPa, phi = 40 degrees. The
#: manual's tables for problems 3 to 7 state the friction angle and no cohesion
#: at all; 0.1 MPa is the models' own value and it is what is built.
LV_JOINT = {'c': 100.0, 'phi': 40.0, 't_cut': 0.0, 'kn': KN_STD, 'ks': KS_STD}



# ---------------------------------------------------------------------------
# Problem 1 — Goodman & Bray block toppling
#
# NOT in BUILDERS, and no file is shipped. The section transcribes cleanly — the
# sixteen column outlines below are the vendor model's own, and the joints are
# derived from where the columns touch — but the mesh split refuses the result
# at every target size from 1.0 to 4.0 m:
#
#     The mesh edge (0, N) on jointed line 1 is carried by 4 two-dimensional
#     element(s), not two.
#
# The cause is the topology, not the geometry. A stepped base puts each column's
# basal contact at a point PARTWAY ALONG its downslope neighbour's side joint,
# so the section has fifteen T-junctions: one joint ending in the interior of
# another. The split copies a shared node once per wedge of material around it,
# which is right for a crossing (two joints, four wedges, two elements per edge)
# and wrong for a termination (three wedges), and the edge at the T comes out
# carrying four elements. Every crossing in the rest of this corpus is an X.
#
# The code is kept because it is the finished transcription: the round that
# teaches the split about terminations can build these four rows by registering
# ``rj001a`` and ``rj001c`` and reading the toe force off the manual for b and d.
# ---------------------------------------------------------------------------

#: Goodman & Bray's sixteen columns, each a closed outline read verbatim
#: from the vendor model's own joint boundaries. The columns sit on a
#: STEPPED base: every column's base is 10 m long at 30 degrees and one
#: metre, measured perpendicular, above its downslope neighbour's.
GB_BLOCKS = [
    [(-0.5, 0.866025403784439), (8.16025403784439, 5.86602540378444), (7.66025403784439, 6.73205080756888), (6.16025403784439, 9.33012701892219), (-2.5, 4.33012701892219)],
    [(7.66025403784439, 6.73205080756888), (16.3205080756888, 11.7320508075689), (15.8205080756888, 12.5980762113533), (12.3205080756888, 18.6602540378444), (3.66025403784439, 13.6602540378444), (6.16025403784439, 9.33012701892219)],
    [(15.8205080756888, 12.5980762113533), (24.4807621135332, 17.5980762113533), (23.9807621135332, 18.4641016151378), (18.4807621135332, 27.9903810567666), (9.82050807568877, 22.9903810567666), (12.3205080756888, 18.6602540378444)],
    [(23.9807621135332, 18.4641016151378), (32.6410161513775, 23.4641016151378), (32.1410161513775, 24.3301270189222), (24.6410161513775, 37.3205080756888), (15.9807621135332, 32.3205080756888), (18.4807621135332, 27.9903810567666)],
    [(32.1410161513775, 24.3301270189222), (40.8012701892219, 29.3301270189222), (40.3012701892219, 30.1961524227066), (30.8012701892219, 46.650635094611), (22.1410161513775, 41.650635094611), (24.6410161513775, 37.3205080756888)],
    [(40.3012701892219, 30.1961524227066), (48.9615242270663, 35.1961524227066), (48.4615242270663, 36.0621778264911), (36.9615242270663, 55.9807621135332), (28.3012701892219, 50.9807621135332), (30.8012701892219, 46.650635094611)],
    [(48.4615242270663, 36.0621778264911), (57.1217782649107, 41.0621778264911), (56.6217782649107, 41.9282032302755), (43.1217782649107, 65.3108891324553), (34.4615242270663, 60.3108891324553), (36.9615242270663, 55.9807621135332)],
    [(56.6217782649107, 41.9282032302755), (65.2820323027551, 46.9282032302755), (64.7820323027551, 47.7942286340599), (49.2820323027551, 74.6410161513775), (40.6217782649107, 69.6410161513775), (43.1217782649107, 65.3108891324553)],
    [(64.7820323027551, 47.7942286340599), (73.4422863405995, 52.7942286340599), (72.9422863405995, 53.6602540378444), (55.4422863405995, 83.9711431702997), (46.7820323027551, 78.9711431702997), (49.2820323027551, 74.6410161513775)],
    [(72.9422863405995, 53.6602540378444), (81.6025403784439, 58.6602540378444), (81.1025403784439, 59.5262794416288), (64.1025403784439, 88.9711431702997), (61.6025403784439, 93.3012701892219), (52.9422863405995, 88.3012701892219), (55.4422863405995, 83.9711431702997)],
    [(81.1025403784439, 59.5262794416288), (89.7627944162882, 64.5262794416288), (89.2627944162882, 65.3923048454133), (75.2627944162882, 89.6410161513776), (72.7627944162882, 93.9711431702997), (64.1025403784439, 88.9711431702997)],
    [(89.2627944162882, 65.3923048454133), (97.9230484541326, 70.3923048454133), (97.4230484541326, 71.2583302491977), (86.4230484541326, 90.3108891324553), (83.9230484541326, 94.6410161513776), (75.2627944162882, 89.6410161513776)],
    [(97.4230484541326, 71.2583302491977), (106.083302491977, 76.2583302491977), (105.583302491977, 77.1243556529821), (97.583302491977, 90.9807621135331), (95.083302491977, 95.3108891324553), (86.4230484541326, 90.3108891324553)],
    [(105.583302491977, 77.1243556529821), (114.243556529821, 82.1243556529821), (113.743556529821, 82.9903810567666), (108.743556529821, 91.6506350946109), (106.243556529821, 95.9807621135331), (97.583302491977, 90.9807621135331)],
    [(113.743556529821, 82.9903810567666), (122.403810567666, 87.9903810567666), (121.903810567666, 88.856406460551), (119.903810567666, 92.3205080756888), (117.403810567666, 96.6506350946109), (108.743556529821, 91.6506350946109)],
    [(121.903810567666, 88.856406460551), (130.56406460551, 93.856406460551), (128.56406460551, 97.3205080756888), (119.903810567666, 92.3205080756888)],
]


#: The rectangle the vendor model is cut from. Everything above the column
#: stack and the crest shelf is a deleted region in the file — 2 100 elements of
#: air — so the domain built here is the rectangle bounded above by the ground.
GB_EXT = (-65.2820323027551, -46.650635094611,
          195.846096908265, 93.856406460551)


def _gb_section():
    """The Goodman & Bray section: its domain, and the joints inside it.

    The vendor file draws sixteen closed column outlines and deletes the air
    above them. Rebuilt here as seventeen polygons — the sixteen columns and the
    rock they stand on — whose SHARED edges are the joints: sixteen stepped
    basal contacts and fifteen column-to-column ones. Deriving the joints from
    the contacts rather than listing them keeps them exactly where the outlines
    put them, and a column that touches its neighbour over three metres gets a
    three metre joint rather than a full-height one.

    Returns ``(domain_ring, joint_segments)``.
    """
    from shapely.geometry import Polygon
    from shapely.ops import linemerge, unary_union

    x0, y0, x1, y1 = GB_EXT
    cols = [Polygon(b) for b in GB_BLOCKS]
    stack = unary_union(cols)
    ring = list(stack.exterior.coords)
    toe = ring.index((-0.5, 0.866025403784439))
    crest = ring.index((130.56406460551, 93.856406460551))
    # The stepped base is the stack's own lower chain, crest back down to toe.
    chain = ring[crest:toe + 1]
    base = Polygon([(x0, y0), (x1, y0), (x1, y1)] + list(chain)
                   + [(0.0, 0.0), (x0, 0.0)])

    polys = [base] + cols
    segs = []
    for a in range(len(polys)):
        for b in range(a + 1, len(polys)):
            shared = polys[a].boundary.intersection(polys[b].boundary)
            if shared.is_empty:
                continue
            merged = (shared if shared.geom_type == 'LineString'
                      else linemerge(shared))
            parts = ([merged] if merged.geom_type == 'LineString'
                     else list(getattr(merged, 'geoms', [])))
            for part in parts:
                if part.geom_type != 'LineString':
                    continue
                pts = list(part.coords)
                # One contact can be a two-segment polyline — a column's base
                # and the step up to the next one. The sheet holds straight
                # lines, so it is written as the segments it is made of.
                for p, q in zip(pts, pts[1:]):
                    if (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 > 1e-12:
                        segs.append((p, q))
    dom = unary_union(polys)
    return list(dom.exterior.coords), _straight_runs(segs)


def _straight_runs(segs, tol=1.0e-9):
    """Contiguous collinear segments joined into one line each.

    The contacts come out of the union one polygon pair at a time, so the step
    between two columns' bases and the contact between those columns above it
    arrive as two segments of one straight plane. Handing the mesher two
    collinear joint lines that meet end to end puts four elements on the edge at
    their shared node and the split refuses it; they are one plane and are
    written as one line.
    """
    import math

    def ang(p, q):
        return math.atan2(q[1] - p[1], q[0] - p[0]) % math.pi

    def key(pt):
        return (round(pt[0], 9), round(pt[1], 9))

    todo = list(segs)
    runs = []
    while todo:
        p, q = todo.pop(0)
        a = ang(p, q)
        grew = True
        while grew:
            grew = False
            for i, (r, t) in enumerate(todo):
                if abs(((ang(r, t) - a + math.pi / 2) % math.pi)
                       - math.pi / 2) > 1.0e-6:
                    continue
                if key(t) == key(p):
                    p = r
                elif key(r) == key(q):
                    q = t
                elif key(r) == key(p):
                    p = t
                elif key(t) == key(q):
                    q = r
                else:
                    continue
                todo.pop(i)
                grew = True
                break
        runs.append((p, q))
    return runs


def _gb_case(name, phi_joint):
    """One of Goodman & Bray's four cases: the shared section at its own joint
    friction angle."""
    sd = _base()
    ring, segs = _gb_section()
    mats = [_rock('Rock', 25.0, 2.0e7, 0.3, 0.0, 0.0, option='elastic')]
    _finish(sd, [(ring, 0)], mats)
    sd['joint_lines'] = [
        _joint(f'jnt-{i + 1:02d}', p, q, 0.0, phi_joint)
        for i, (p, q) in enumerate(segs)]
    # The seed surface is the basal plane the stack stands on, toe to crest.
    sd['non_circ'] = _surface([(-0.5, 0.866025403784439),
                               (130.56406460551, 93.856406460551)])
    return _write(sd, name)


def rj001a():
    """RJ-1a — Goodman & Bray block toppling, case a (vendor `joint #001_a.fez`).

    Sixteen rock columns on a stepped 30 degree base, rising to a crest at
    (130.564, 93.856) and cut off by the slope face above column ten. The rock
    is ELASTIC (gamma 25 kN/m^3, E 20 GPa, nu 0.3), so the columns cannot break
    and the whole mechanism is sliding and rotation on the joints, which is the
    idealization Goodman & Bray's limit equilibrium makes.

    The joints carry no cohesion and phi = 38.15 degrees — the angle this case
    is posed at. Referees: Goodman & Bray 1.0 and UDEC 0.99. RS2 reports 0.99
    without joint improvement and 0.97 with it.

    The manual describes a stabilizing force at the toe of the lowest column. In
    this case the model carries 0.5 kN, which is 0.05% of that column's own
    weight and does nothing, so the file is built without it; see the page's
    departures table.
    """
    return _gb_case('rj001a.xlsx', 38.15)


def rj001c():
    """RJ-1c — Goodman & Bray block toppling, case c (vendor `joint #001_c.fez`).

    Case a's section and rock at phi = 38.6598 degrees on the joints, and the
    same 0.5 kN toe force that does nothing. Referees: Goodman & Bray 1.02 and
    UDEC 1.01. RS2 reports 1.01 without joint improvement and 0.99 with it.
    """
    return _gb_case('rj001c.xlsx', 38.6598)

# ---------------------------------------------------------------------------
# Problem 2 — Alejano & Alonso block toppling
# ---------------------------------------------------------------------------

def rj002():
    """RJ-2 — Alejano & Alonso block toppling (vendor `joint #002.fez`).

    A 30 x 19.85 m section whose 9.85 m face rises at 58.65 degrees from
    (20, 10) to (14, 19.85). A basal joint runs from the face's toe at (20, 10)
    up to the crest at (2.9393, 19.85) at 30 degrees — the stepped surface the
    columns topple over — and a set of columns at 64 degrees, 1.6 m apart,
    passes through that toe. Both carry phi = 31 degrees and no cohesion.

    The rock is ELASTIC (the vendor's ``Plasticity Specifications: Non``), so
    the columns cannot yield and every mechanism the model has is a joint one.
    gamma = 25 kN/m^3, E = 20 GPa, nu = 0.3.

    Referees: Goodman & Bray's limit-equilibrium 0.76 and UDEC 0.87. RS2 reports
    0.86 without joint improvement and 0.82 with it.

    The set's sign is the file's own: its 24 stored segments run at +64 degrees
    counter-clockwise (ascending to the right, into the slope), which is what
    makes them columns over a basal plane rather than a second sliding set.
    """
    sd = _base()
    ring = [(0.0, 0.0), (30.0, 0.0), (30.0, 10.0), (20.0, 10.0),
            (14.0, 19.85), (2.9393, 19.85), (0.0, 19.85)]
    mats = [_rock('Rock', 25.0, 2.0e7, 0.3, 0.0, 0.0, option='elastic')]
    _finish(sd, [(ring, 0)], mats)
    props = {'c': 0.0, 'phi': 31.0, 't_cut': 0.0, 'kn': KN_STD, 'ks': KS_STD}
    columns = parallel_set(sd, 64.0, 1.6, offset=_offset_through(64.0, 20.0, 10.0),
                           label='col', props=props)
    sd['joint_lines'] = [
        _joint('basal', (20.0, 10.0), (2.9393, 19.85), 0.0, 31.0),
    ] + columns
    # The seed surface is the basal joint: the plane the toppling column stack
    # stands on. Inert for a strength reduction, and it makes the file complete.
    sd['non_circ'] = _surface([(20.0, 10.0), (2.9393, 19.85)])
    return _write(sd, 'rj002.xlsx')


# ---------------------------------------------------------------------------
# Problem 3 — Lorig & Varona forward block toppling
# ---------------------------------------------------------------------------

def rj003():
    """RJ-3 — Lorig & Varona forward block toppling (vendor `joint #003.fez`).

    The shared 260 m / 55 degree section cut by two sets: columns at 70 degrees
    at 20 m spacing, and a cross set at -20 degrees at 30 m, both through the
    origin. The manual states the two as "70 and 160" degrees, which is the
    same pair measured the other way round the half circle.

    The rock is ELASTIC, so only the joints can fail. gamma = 26.0946 kN/m^3,
    E = 9.072 GPa, nu = 0.26. Referee: UDEC 1.13. RS2 reports 1.12 without joint
    improvement and 1.09 with it.
    """
    sd = _base()
    mats = [_rock('Rock', 26.0946, 9.072e6, 0.26, 0.0, 0.0, option='elastic')]
    _finish(sd, [(LV_RING, 0)], mats)
    sd['joint_lines'] = cross_jointed(
        parallel_set(sd, 70.0, 20.0, label='col', props=LV_JOINT),
        parallel_set(sd, -20.0, 30.0, label='cross', props=LV_JOINT))
    sd['circles'] = _toe_circle((560.0, 140.0), (377.785587105239, 400.0))
    return _write(sd, 'rj003.xlsx')


# ---------------------------------------------------------------------------
# Problem 4 — Lorig & Varona flexural toppling
# ---------------------------------------------------------------------------

def rj004():
    """RJ-4 — Lorig & Varona flexural toppling (vendor `joint #004.fez`).

    The shared section cut by ONE set of columns at 70 degrees at 20 m spacing
    through the origin — problem 3's first set without its cross joints, so the
    columns bend rather than topple as blocks.

    The rock is Mohr-Coulomb and carries a tensile cutoff of zero, which is what
    lets a column break in flexure: gamma = 26.1 kN/m^3, E = 9.072 GPa,
    nu = 0.26, c = 675 kPa, phi = 43 degrees, T = 0. Referee: UDEC 1.3. RS2
    reports 1.19 without joint improvement and 1.27 with it.
    """
    sd = _base()
    mats = [_rock('Rock', 26.1, 9.072e6, 0.26, 675.0, 43.0, t_cut=0.0)]
    _finish(sd, [(LV_RING, 0)], mats)
    sd['joint_lines'] = parallel_set(sd, 70.0, 20.0, label='col', props=LV_JOINT)
    sd['circles'] = _toe_circle((560.0, 140.0), (377.785587105239, 400.0))
    return _write(sd, 'rj004.xlsx')


# ---------------------------------------------------------------------------
# Problem 5 — Lorig & Varona backward block toppling
# ---------------------------------------------------------------------------

def rj005():
    """RJ-5 — Lorig & Varona backward block toppling (vendor `joint #005.fez`).

    The shared section cut by a set at -55 degrees at 10 m spacing through the
    toe at (560, 140) — dipping out of the face, so the blocks topple backward —
    and a horizontal set at 40 m spacing through (0, 400).

    The rock is ELASTIC. gamma = 26.1 kN/m^3, E = 9.072 GPa, nu = 0.26.
    Referee: UDEC 1.7. RS2 reports 1.65 without joint improvement and 1.86 with
    it — the widest "joint improvement" spread in the manual.
    """
    sd = _base()
    mats = [_rock('Rock', 26.1, 9.072e6, 0.26, 0.0, 0.0, option='elastic')]
    _finish(sd, [(LV_RING, 0)], mats)
    sd['joint_lines'] = cross_jointed(
        parallel_set(sd, -55.0, 10.0, offset=_offset_through(-55.0, 560.0, 140.0),
                     label='dip', props=LV_JOINT),
        parallel_set(sd, 0.0, 40.0, offset=_offset_through(0.0, 0.0, 400.0),
                     label='bed', props=LV_JOINT))
    sd['circles'] = _toe_circle((560.0, 140.0), (377.785587105239, 400.0))
    return _write(sd, 'rj005.xlsx')


# ---------------------------------------------------------------------------
# Problem 6 — plane failure, daylighting discontinuities
# ---------------------------------------------------------------------------

def rj006():
    """RJ-6 — plane failure with daylighting discontinuities (`joint #006.fez`).

    The shared section cut by one set at -35 degrees at 10 m spacing through the
    origin. The joints dip out of the 55 degree face at a shallower angle than
    the face itself, so every one of them daylights and the slabs between them
    are free to slide out.

    The rock is Mohr-Coulomb: gamma = 26.1 kN/m^3, E = 9.072 GPa, nu = 0.26,
    c = 675 kPa, phi = 43 degrees, no tensile capacity. Referee: UDEC 1.27. RS2
    reports 1.25 without joint improvement and 1.31 with it.
    """
    sd = _base()
    mats = [_rock('Rock', 26.1, 9.072e6, 0.26, 675.0, 43.0, t_cut=0.0)]
    _finish(sd, [(LV_RING, 0)], mats)
    sd['joint_lines'] = parallel_set(sd, -35.0, 10.0, label='jnt', props=LV_JOINT)
    # The seed surface is the daylighting plane itself: the 35 degree joint
    # through the toe of the face, back to where it reaches the crest plateau.
    # Inert for a strength reduction, and it makes the file a complete model.
    sd['non_circ'] = _surface([(188.6478, 400.0), (560.0, 140.0)])
    return _write(sd, 'rj006.xlsx')


# ---------------------------------------------------------------------------
# Problem 7 — plane failure, non-daylighting discontinuities
# ---------------------------------------------------------------------------

def rj007():
    """RJ-7 — plane failure with non-daylighting discontinuities (`joint #007.fez`).

    The same section and the same rock as problem 6, cut by one set at
    -70 degrees at 20 m spacing through the origin. The joints now dip out of
    the face STEEPER than the 55 degree face, so none of them daylights: a slab
    cannot slide out along one without shearing rock, and the slope stands
    higher than problem 6's.

    Referee: UDEC 1.5. RS2 reports 1.57 without joint improvement and 1.59 with
    it. The manual's own table prints the slope angle as 5 degrees; the figure
    and the model are the same 55 degrees problem 6 uses.
    """
    sd = _base()
    mats = [_rock('Rock', 26.1, 9.072e6, 0.26, 675.0, 43.0, t_cut=0.0)]
    _finish(sd, [(LV_RING, 0)], mats)
    sd['joint_lines'] = parallel_set(sd, -70.0, 20.0, label='jnt', props=LV_JOINT)
    sd['circles'] = _toe_circle((560.0, 140.0), (377.785587105239, 400.0))
    return _write(sd, 'rj007.xlsx')


# ---------------------------------------------------------------------------
# Problem 8 — flexural toppling in a base friction model
# ---------------------------------------------------------------------------

def rj008():
    """RJ-8 — flexural toppling, base friction model (vendor `joint #008.fez`).

    Pritchard & Savigny's base-friction table model, scaled up a hundred times:
    a 72.407 x 36.5 m section whose 30.5 m face rises at 78 degrees from (15, 6)
    to (21.48, 36.5). Columns at -60 degrees at 5.08 m spacing pass through the
    crest at (21.48, 36.5); a horizontal joint at y = 6 runs the width of the
    model under them, and a vertical joint at x = 68.4 closes the back of the
    column stack.

    The rock is Mohr-Coulomb with a tensile cutoff of 75 kPa — the strength a
    column has to break in flexure: gamma = 25.506 kN/m^3, E = 22.771 GPa,
    nu = 0.139, c = 60 kPa, phi = 39 degrees. The joints carry no cohesion,
    phi = 39 degrees, and a NORMAL stiffness of 1.5e7 kPa/m rather than the
    set's usual 1e8 — this is the one problem in the corpus whose joints are
    softer than the standard pair. Referee: UDEC 0.76. RS2 reports 0.75 both
    with and without joint improvement.
    """
    sd = _base()
    ring = [(0.0, 0.0), (72.407, 0.0), (72.407, 6.0), (72.407, 36.5),
            (68.4, 36.5), (21.48, 36.5), (15.0, 6.0), (0.0, 6.0)]
    mats = [_rock('Rock', 25.506, 2.2771e7, 0.139, 60.0, 39.0, t_cut=75.0)]
    _finish(sd, [(ring, 0)], mats)
    kn, ks = 1.5e7, 1.0e7
    props = {'c': 0.0, 'phi': 39.0, 't_cut': 0.0, 'kn': kn, 'ks': ks}
    # The column set exists only ABOVE the basal joint, which is how the vendor
    # file states it: its 13 stored segments are clipped to the block bounded
    # below by y = 6 and on the left by the face, not to the whole section. The
    # same region generated over the whole section would put columns under the
    # basal joint, where the model has none.
    stack = [(15.0, 6.0), (72.407, 6.0), (72.407, 36.5), (21.48, 36.5)]
    columns = parallel_set(sd, -60.0, 5.08,
                           offset=_offset_through(-60.0, 21.48, 36.5),
                           region=stack, label='col', props=props)
    sd['joint_lines'] = [
        _joint('base-1', (15.0, 6.0), (68.4, 6.0), 0.0, 39.0, kn=kn, ks=ks),
        _joint('base-2', (68.4, 6.0), (72.407, 6.0), 0.0, 39.0, kn=kn, ks=ks),
        _joint('back', (68.4, 36.5), (68.4, 6.0), 0.0, 39.0, kn=kn, ks=ks),
    ] + columns
    # The seed surface is the basal joint the column stack stands on.
    sd['non_circ'] = _surface([(15.0, 6.0), (72.407, 6.0)])
    return _write(sd, 'rj008.xlsx')



# ---------------------------------------------------------------------------
# Problems 9 to 14 — Alejano et al. sliding and ploughing slabs
#
# NOT in BUILDERS, and no file is shipped, for the reason problem 1 is not:
# a joint that ends ON another joint. Here the termination is a NEAR one. Each
# release trace is meant to run from the crest down to a bedding plane, and the
# vendor states its tip to six decimals, so the tip lands 2e-7 to 2e-6 from the
# bedding trace it belongs on — a part in 10^7 of the section. That leaves a
# sliver between the two lines, and the mesher reports it two ways:
#
#     rj009, rj010   the two elements on mesh edge (N, N+1) of jointed line 48
#                    stand on the same side of it (3.3e-08, -0.437)
#     rj011          the mesh edge (125, 1706) on jointed line 45 is carried by
#                    4 two-dimensional element(s), not two
#     rj012          gmsh never returns: "Impossible to recover edge 113 113",
#                    "2 intersections in the 1D mesh (curves 113 168)", split
#                    and retry, level after level (recorded by a 300 s alarm)
#
# Snapping the tip onto the bedding line would turn the near termination into an
# exact one, which is problem 1's refusal, so both families want the same thing:
# a split that knows what to do where one joint ENDS on another. The
# transcription is kept whole so that round can register these seven functions.
# ---------------------------------------------------------------------------

#: Every one of these six models carries TWO joint strengths. The bedding
#: network runs at one friction angle and the one or two short explicit traces
#: that release the slab at the crest run at another — the vendor states it as a
#: per-segment ``segment joint property`` index into its own joint list, and a
#: single quoted joint friction angle for these problems is incomplete.
#:
#: (ring, bedding dip, bedding spacing, phi_bedding, phi_release, release traces)
ALEJANO = {
    'rj009': ([(58.045, -50.0), (58.045, 0.0), (0.0, 0.0), (-1.4649, 1.74595),
               (-6.76319540951412, 8.060062428408), (-41.9549, 49.9999),
               (-141.955, 50.0), (-141.955, -50.0)],
              -50.0, 3.0, 30.0, 40.0,
              [((-6.76319540951412, 8.060062428408),
                (-9.06132845221374, 6.13169959257877)),
               ((-9.06132845221374, 6.13169959257877), (-1.4649, 1.74595))]),
    'rj010': ([(58.045, -50.0), (58.045, 0.0), (0.0, 0.0), (-1.4649, 1.74595),
               (-9.97696, 11.8902), (-41.9549, 49.9999),
               (-141.955, 50.0), (-141.955, -50.0)],
              -50.0, 3.0, 30.0, 40.0,
              [((-9.06132845221374, 6.13169959257877), (-1.4649, 1.74595)),
               ((-9.97696, 11.8902), (-12.2753, 9.96194))]),
    'rj011': ([(30.0, -30.0), (30.0, 0.0), (0.0, 0.0),
               (-2.43625294577519, 2.90341192441329), (-20.9775, 25.0),
               (-60.0, 25.0), (-60.0, -30.0)],
              -50.0, 1.5, 30.0, 20.0,
              [((0.0, 0.0), (-1.14907, -0.964181)),
               ((-2.43625294577519, 2.90341192441329), (-2.62114, 0.790165))]),
    'rj012': ([(30.0, -30.0), (30.0, 0.0), (0.0, 0.0),
               (-2.57111, 4.45329), (-14.4338, 25.0),
               (-60.0, 25.0), (-60.0, -30.0)],
              -60.0, 1.5, 30.0, 40.0,
              [((0.0, 0.0), (-1.29904, -0.75)),
               ((-2.57111, 4.45329), (-2.79904, 1.84808))]),
    'rj013': ([(30.0, -30.0), (30.0, 0.0), (0.0, 0.0),
               (-2.22963, 3.18425), (-17.5052, 25.0),
               (-60.0, 25.0), (-60.0, -30.0)],
              -55.0, 1.5, 25.0, 20.0,
              [((0.0, 0.0), (-1.22852, -0.860663)),
               ((-2.22963, 3.18425), (-2.42975, 0.854877))]),
    'rj014': ([(30.0, -30.0), (30.0, 0.0), (0.0, 0.0),
               (-2.86498, 4.96228), (-14.4338, 25.0),
               (-60.0, 25.0), (-60.0, -30.0)],
              -60.0, 1.5, 20.0, 30.0,
              [((0.0, 0.0), (-1.29928, -0.749584)),
               ((-2.86498, 4.96228), (-3.09385, 2.35871))]),
}


def _alejano(name):
    """One of the six Alejano slab models.

    The rock is ELASTIC at E = 2 x 10^8 MPa, which is not a rock modulus but the
    manual's own device: the UDEC model these are scored against uses RIGID
    blocks, and the manual says outright that "an artificially high modulus of
    2x10^8 MPa was given to the material" to reproduce that, with a tightened
    convergence tolerance to go with it. gamma = 25 kN/m^3, nu = 0.3.
    """
    ring, dip, spacing, phi_bed, phi_rel, release = ALEJANO[name]
    sd = _base()
    mats = [_rock('Rock', 25.0, 2.0e11, 0.3, 0.0, 0.0, option='elastic')]
    _finish(sd, [(ring, 0)], mats)
    bed = {'c': 0.0, 'phi': phi_bed, 't_cut': 0.0, 'kn': KN_STD, 'ks': KS_STD}
    sd['joint_lines'] = [
        _joint(f'rel-{i + 1:02d}', p, q, 0.0, phi_rel)
        for i, (p, q) in enumerate(release)
    ] + parallel_set(sd, dip, spacing, label='bed', props=bed)
    sd['circles'] = _toe_circle((0.0, 0.0), (ring[-3][0], ring[-3][1]))
    return _write(sd, name + '.xlsx')


def rj009():
    """RJ-9 — Alejano et al. bilinear slab failure, example 1a (`joint #009.fez`).

    A 50 m slope at 50 degrees with bedding at -50 degrees at 3 m spacing
    (phi = 30) and a release trace under the crest at phi = 40. Sliding on the
    basal plane combines with sliding on a shallow joint undercut by the face.
    Referee: UDEC 1.03; the paper's limit equilibrium spans 0.40 to 1.45, which
    is a range rather than an answer. RS2 reports 1.01 without joint improvement
    and 1.09 with it.
    """
    return _alejano('rj009')


def rj010():
    """RJ-10 — bilinear slab failure, example 1b (vendor `joint #010.fez`).

    Example 1a with the release joint moved upslope, which is the whole
    difference between the two. Referee: UDEC 1.03. RS2 reports 0.92 without
    joint improvement and 1.08 with it.
    """
    return _alejano('rj010')


def rj011():
    """RJ-11 — ploughing sliding slab failure (vendor `joint #011.fez`).

    A 25 m slope with bedding at -50 degrees at 1.5 m (phi = 30) and two release
    traces at phi = 20, one of which ends inside the rock at the toe. Referee:
    UDEC 1.21. RS2 reports 1.22 without joint improvement and 1.30 with it.
    """
    return _alejano('rj011')


def rj012():
    """RJ-12 — ploughing toppling slab failure (vendor `joint #012.fez`).

    Bedding at -60 degrees at 1.5 m (phi = 30) with releases at phi = 40.
    Referee: UDEC 1.78. RS2 reports 1.39 without joint improvement and 1.75 with
    it — the family's widest spread between the two RS2 runs.
    """
    return _alejano('rj012')


def rj013():
    """RJ-13 — ploughing sliding slab, example 4 (vendor `joint #013.fez`).

    Bedding at -55 degrees at 1.5 m (phi = 25) with releases at phi = 20.
    Referee: UDEC 1.0. RS2 reports 1.0 without joint improvement and 1.05 with
    it.
    """
    return _alejano('rj013')


def rj014():
    """RJ-14 — ploughing sliding slab, example 5 (vendor `joint #014.fez`).

    Example 4's section with the release joint moved and the two strengths the
    other way round: bedding at -60 degrees at 1.5 m at phi = 20, releases at
    phi = 30. Referee: UDEC 0.9. RS2 reports 0.89 without joint improvement and
    1.09 with it.
    """
    return _alejano('rj014')


# ---------------------------------------------------------------------------
# Problem 15 — partially joint-controlled footwall slope
# ---------------------------------------------------------------------------

def rj015():
    """RJ-15 — partially joint-controlled footwall slope (`joint #015.fez`).

    A 25 m footwall at 40 degrees with bedding dipping the same way at the same
    angle, 2 m apart, so the slabs are parallel to the face and the failure has
    to break rock at the toe to get out. That is the one problem in this family
    whose rock can yield: Mohr-Coulomb, c = 200 kPa, phi = 35 degrees,
    gamma = 28 kN/m^3, E = 1 GPa, nu = 0.3.

    The stated tensile strength of 1000 kPa is above the Mohr-Coulomb apex
    c/tan(phi) = 285.6 kPa, so it never binds and the apex governs.

    The joints are also the corpus's softest: k_n = 5 x 10^6 kPa/m and
    k_s = 5 x 10^5 kPa/m, twenty times below the standard pair, with no cohesion
    and phi = 25 degrees. Referees: Slide's limit equilibrium 1.25 and UDEC 1.6.
    RS2 reports 1.28 without joint improvement and 1.42 with it.
    """
    sd = _base()
    ring = [(40.0, -40.0), (40.0, 0.0), (0.0, 0.0), (-47.6701, 40.0),
            (-100.0, 40.0), (-100.0, -40.0)]
    mats = [_rock('Rock', 28.0, 1.0e6, 0.3, 200.0, 35.0, t_cut=1000.0)]
    _finish(sd, [(ring, 0)], mats)
    props = {'c': 0.0, 'phi': 25.0, 't_cut': 0.0, 'kn': 5.0e6, 'ks': 5.0e5}
    sd['joint_lines'] = parallel_set(sd, -40.0, 2.0, label='bed', props=props)
    sd['circles'] = _toe_circle((0.0, 0.0), (-47.6701, 40.0))
    return _write(sd, 'rj015.xlsx')

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
    # The seed failure surface. A rock slope cut by through-going joints does not
    # fail on a circle, and the surface worth shipping with the file is the one
    # the joints draw: the lowest joint, from the toe of the face to the crest.
    # Inert for a strength reduction, and it makes the file a complete model.
    sd['non_circ'] = _surface([(17.0, 8.2), (33.2, 20.0)])
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
    # The seed surface is the bi-planar path itself: up the basal joint, across
    # the rock bridge, and out along the upper one.
    sd['non_circ'] = _surface([(39.0149, 35.0248), (63.0, 48.0),
                               (76.0, 70.0)])
    return _write(sd, 'rj019.xlsx')


#: Every builder in this module, in problem-number order. ``verify_rebuild.py``'s
#: ``joints`` group is this list, so a builder missing here is a corpus file
#: nothing guards.
BUILDERS = [rj002, rj015, rj003, rj004, rj005, rj006, rj007, rj008, rj018, rj019]


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    for b in BUILDERS:
        print('built', b())
