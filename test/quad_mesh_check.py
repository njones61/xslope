"""Standing checks on the mesh builder.

Five things are guarded here, and each one has a specific way of going wrong that
was seen in practice:

(a) FREE-STYLE QUALITY on a spread of corpus sections. The quad path used to
    inflate the requested element size by 1.4, subdivide everything into quads
    afterwards, and hand gmsh hand-written node-count constraints on every long
    polygon edge. The net effect was elements about 0.72x the size the user asked
    for and, because the node-count constraint's rounding error depends only on how
    long a particular edge happens to be, zones meshed at visibly different sizes
    stitched together — up to 7.7:1 on the worst element. This checks the worst
    element, the fraction of stretched elements, and the delivered element size.

(b) DETERMINISM. Verification factors of safety are locked to specific meshes, so
    the same input has to produce the same mesh, run to run, in both styles.

(c) STRUCTURED SWEEPS on the levee section: the three rectangular foundation
    blocks come out as exact grids, the 5:1 trapezoidal embankment is declined and
    free-meshed, and the interface between the two conforms. The off-by-one that
    this whole exercise exists to fix lives here — gmsh's setTransfiniteCurve takes
    a NODE count, not an element count — so a swept rectangle's element count is
    asserted against rows x columns exactly.

(d) STRUCTURED IS NEVER WORSE THAN FREE. A zone the classifier wrongly accepts
    produces a worse mesh than free meshing (force-sweeping the levee's 5:1
    embankment drives its elements to 0.33 against a 0.58 target). The fallback
    guarantee is that a zone in doubt is simply not swept, so the structured style
    can at worst equal the free style.

(e) THE ELEMENT TYPE DEFAULT. A caller that states no element type — a script, or
    a model whose main!D18 is blank, which reaches the builder as None — must get
    tri6. Silence must not select a linear element: the stress solvers lock on
    them and read the factor of safety high. The explicit 'tri3' request, which is
    the seepage choice, has to keep working alongside it.

Run standalone with ``python3 test/quad_mesh_check.py``.
"""
import contextlib
import hashlib
import io
import os
import sys
import warnings

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

# Sections meshed by (a), (b) and (d). The first five are the range of shapes the
# corpus actually contains — one zone, a zoned dam, a layered section, a benched
# slope with thin inclined seams, a thin weak seam — and the last two are the
# sections the structured sweeps were prototyped on.
GEOMETRIES = [
    ('simple1',  'docs/inputs/slope/xslope_simple1.xlsx'),
    ('johnson',  'docs/seep/files/xslope_johnson_res.xlsx'),
    ('vp042',    'docs/verification/files/rocscience/vp042.xlsx'),
    ('vp005',    'docs/verification/files/rocscience/vp005.xlsx'),
    ('rs2_9',    'docs/verification/files/rocscience/rs2_9.xlsx'),
    ('levee',    'docs/seep/files/xslope_levee_poly.xlsx'),
    ('earthdam', 'docs/seep/files/xslope_earth_dam1.xlsx'),
]

# (a) thresholds. Measured on these seven sections: worst element 3.06, worst
# stretched fraction 0.07 %, delivered size 0.93-1.00 x requested.
AR_MAX_LIMIT = 3.5
AR_GT3_LIMIT_PCT = 0.5
MED_OVER_REQ_RANGE = (0.85, 1.10)

# (d) the structured style may not degrade the shape statistics. AR95 is the
# stable measure of "the mesh got worse"; ARmax is a single element and is
# additionally held under the same absolute ceiling as (a).
AR95_TOLERANCE = 1.15
ARMAX_TOLERANCE = 1.15


def _quiet(fn, *args, **kwargs):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        return fn(*args, **kwargs)


def _load(rel_path):
    """(polygons, target_size) for a corpus section, at the Studio default size."""
    from xslope.fileio import load_slope_data
    from xslope.mesh import get_material_polygons
    data = _quiet(load_slope_data, os.path.join(_REPO, rel_path))
    polys = _quiet(get_material_polygons, data)
    gx = [x for x, _ in data['ground_surface'].coords]
    return polys, (max(gx) - min(gx)) / 100.0


def _corners(mesh, i):
    n = 4 if mesh['element_types'][i] in (4, 8, 9) else 3
    return mesh['nodes'][mesh['elements'][i][:n].astype(int)]


def _metrics(mesh, target_size):
    """Aspect ratio statistics and the delivered element size."""
    ars, sizes = [], []
    for i in range(len(mesh['elements'])):
        v = _corners(mesh, i)
        lengths = np.linalg.norm(np.roll(v, -1, axis=0) - v, axis=1)
        if lengths.min() <= 0:
            continue
        ars.append(lengths.max() / lengths.min())
        x, y = v[:, 0], v[:, 1]
        area = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
        sizes.append(np.sqrt(area))
    ars = np.array(ars)
    return {
        'n': len(mesh['elements']),
        'ar95': float(np.percentile(ars, 95)),
        'armax': float(ars.max()),
        'gt3_pct': float(100.0 * np.mean(ars > 3.0)),
        'med_over_req': float(np.median(sizes) / target_size),
    }


def _mesh_hash(mesh):
    h = hashlib.sha256()
    for key in ('nodes', 'elements', 'element_types', 'element_materials'):
        h.update(np.ascontiguousarray(np.asarray(mesh[key])).tobytes())
    return h.hexdigest()


def _edge_owners(mesh):
    """{(node_a, node_b): [element indices]} over every element edge."""
    owners = {}
    for i in range(len(mesh['elements'])):
        n = 4 if mesh['element_types'][i] in (4, 8, 9) else 3
        ring = [int(v) for v in mesh['elements'][i][:n]]
        for k in range(n):
            a, b = ring[k], ring[(k + 1) % n]
            owners.setdefault((min(a, b), max(a, b)), []).append(i)
    return owners


# --------------------------------------------------------------------- (a)
def check_free_quality(build, failures):
    """Free-style meshes honour the requested size and stay well shaped."""
    for name, path in GEOMETRIES:
        polys, ts = _load(path)
        mesh = _quiet(build, polys, ts, 'quad4', quad_style='free')
        m = _metrics(mesh, ts)
        if m['armax'] > AR_MAX_LIMIT:
            failures.append(f"free/{name}: worst aspect ratio {m['armax']:.2f} "
                            f"> {AR_MAX_LIMIT}")
        if m['gt3_pct'] > AR_GT3_LIMIT_PCT:
            failures.append(f"free/{name}: {m['gt3_pct']:.2f} % of elements beyond "
                            f"3:1 > {AR_GT3_LIMIT_PCT} %")
        lo, hi = MED_OVER_REQ_RANGE
        if not (lo <= m['med_over_req'] <= hi):
            failures.append(f"free/{name}: median element {m['med_over_req']:.3f} x "
                            f"the requested size, outside [{lo}, {hi}]")


# --------------------------------------------------------------------- (b)
def check_determinism(build, failures):
    """The same input produces the same mesh, run to run, in both styles."""
    for name, path in GEOMETRIES[:3] + GEOMETRIES[-2:]:
        polys, ts = _load(path)
        for style in ('free', 'structured'):
            h1 = _mesh_hash(_quiet(build, polys, ts, 'quad4', quad_style=style))
            h2 = _mesh_hash(_quiet(build, polys, ts, 'quad4', quad_style=style))
            if h1 != h2:
                failures.append(f"{style}/{name}: two builds of the same input "
                                f"differ ({h1[:12]} vs {h2[:12]})")


# --------------------------------------------------------------------- (c)
def check_structured_levee(build, failures):
    """The levee: foundation swept as exact grids, embankment free, interface conforms."""
    from xslope.mesh import detect_sweepable_regions

    n_before = len(failures)
    polys, ts = _load('docs/seep/files/xslope_levee_poly.xlsx')
    coords = [p['coords'] for p in polys]
    swept, divisions = detect_sweepable_regions(
        coords, ts, [p.get('size') for p in polys])

    # The trapezoidal embankment is the zone that must NOT be swept: its base is
    # 18 m against a 3.6 m crest, and a transfinite surface would put the same
    # number of divisions on both.
    embankment = [i for i, p in enumerate(polys) if p.get('mat_id') == 0]
    for i in embankment:
        if i in swept:
            failures.append(f"structured/levee: the trapezoidal embankment (zone {i}) "
                            f"was swept; it is 5:1 and must fall back to free meshing")
    # The three rectangular foundation blocks must be.
    blocks = [i for i, p in enumerate(polys) if p.get('mat_id') in (1, 2)]
    for i in blocks:
        if i not in swept:
            failures.append(f"structured/levee: foundation zone {i} was not swept")
    if len(failures) > n_before:
        return          # the rest of this check assumes the classifier's answer

    # The exact grid each block must come out as, from its own edge lengths and
    # the requested element size: 26 m and 6 m wide, 10 m deep, at 0.58.
    expect = {}
    for i in blocks:
        ring = coords[i]
        xs = [p[0] for p in ring]
        n_u = divisions[frozenset([(round(min(xs), 9), 0.0),
                                   (round(max(xs), 9), 0.0)])]
        n_v = divisions[frozenset([(round(min(xs), 9), 0.0),
                                   (round(min(xs), 9), 10.0)])]
        expect[i] = (n_u, n_v)
    if sorted(expect.values()) != [(10, 17), (45, 17), (45, 17)]:
        failures.append(f"structured/levee: division counts {sorted(expect.values())}, "
                        f"expected two 45 x 17 blocks and one 10 x 17 curtain")
        return

    mesh = _quiet(build, polys, ts, 'quad4', quad_style='structured')
    mats = np.asarray(mesh['element_materials'])
    types = np.asarray(mesh['element_types'])
    # element_materials is 1-based on the way out of the builder.
    per_material = {}
    for i in blocks:
        per_material.setdefault(polys[i]['mat_id'] + 1, 0)
        per_material[polys[i]['mat_id'] + 1] += expect[i][0] * expect[i][1]
    for mat, want in per_material.items():
        got = int(np.sum(mats == mat))
        if got != want:
            failures.append(
                f"structured/levee: material {mat} has {got} elements, expected "
                f"exactly {want} (rows x columns) — a swept block is not a full grid. "
                f"setTransfiniteCurve takes a NODE count (divisions + 1); passing the "
                f"element count is the classic way to land here")
        quads = int(np.sum((mats == mat) & np.isin(types, [4, 8, 9])))
        if quads != got:
            failures.append(f"structured/levee: material {mat} has {got - quads} "
                            f"triangles inside a swept block")

    # Every element in a swept block is near-square.
    for i in range(len(mesh['elements'])):
        if mats[i] not in per_material:
            continue
        v = _corners(mesh, i)
        lengths = np.linalg.norm(np.roll(v, -1, axis=0) - v, axis=1)
        ar = lengths.max() / max(lengths.min(), 1e-12)
        if ar > 1.25:
            failures.append(f"structured/levee: a swept element has aspect ratio "
                            f"{ar:.2f}; a grid at the requested size should be square")
            break

    # Interior nodes of a swept grid have exactly four elements around them.
    node_elems = {}
    for i in range(len(mesh['elements'])):
        n = 4 if types[i] in (4, 8, 9) else 3
        for v in mesh['elements'][i][:n]:
            node_elems.setdefault(int(v), []).append(i)
    interior_bad = 0
    for node, elems in node_elems.items():
        x, y = mesh['nodes'][node]
        if not (0.05 < x < 25.95 and 0.05 < y < 9.95):
            continue                       # keep to the interior of the left block
        if len(elems) != 4:
            interior_bad += 1
    if interior_bad:
        failures.append(f"structured/levee: {interior_bad} interior nodes of the swept "
                        f"foundation block are not shared by exactly 4 elements")

    # Conformity across the swept/free interface: no edge may be owned by more
    # than two elements, and every node on the foundation top under the
    # embankment must be shared by the embankment and the foundation alike.
    owners = _edge_owners(mesh)
    over = [e for e, o in owners.items() if len(o) > 2]
    if over:
        failures.append(f"structured/levee: {len(over)} mesh edges are shared by more "
                        f"than two elements")
    interface = [n for n, e in node_elems.items()
                 if abs(mesh['nodes'][n][1] - 10.0) < 1e-6
                 and 20.0 + 1e-6 < mesh['nodes'][n][0] < 38.0 - 1e-6]
    if len(interface) < 20:
        failures.append(f"structured/levee: only {len(interface)} nodes on the "
                        f"swept/free interface — expected about 30")
    for node in interface:
        around = set(int(mats[i]) for i in node_elems[node])
        if 1 not in around or not around & set(per_material):
            failures.append(
                f"structured/levee: node at {tuple(np.round(mesh['nodes'][node], 3))} "
                f"on the swept/free interface is not shared by both sides "
                f"(materials {sorted(around)}) — a hanging node")
            break
    # and the spacing across that interface is the requested one on both sides.
    for (a, b), elems in owners.items():
        ya, yb = mesh['nodes'][a][1], mesh['nodes'][b][1]
        xa, xb = mesh['nodes'][a][0], mesh['nodes'][b][0]
        if abs(ya - 10.0) > 1e-6 or abs(yb - 10.0) > 1e-6:
            continue
        if not (20.0 - 1e-6 <= min(xa, xb) and max(xa, xb) <= 38.0 + 1e-6):
            continue
        length = abs(xa - xb)
        if not (0.7 * ts <= length <= 1.4 * ts):
            failures.append(f"structured/levee: an interface edge is {length:.3f} long "
                            f"against a {ts:.3f} requested size")
            break


def check_structured_single_block(build, failures):
    """The node-count guard, on a geometry with one right answer.

    A 20 x 10 rectangle at a requested size of 1.0 sweeps to exactly 20 x 10 = 200
    square elements. gmsh's setTransfiniteCurve takes the NUMBER OF NODES, so those
    20 divisions are declared as 21 nodes; declaring 20 gives 19 divisions and
    19 x 9 = 171 elements.
    """
    ring = [(0.0, 0.0), (20.0, 0.0), (20.0, 10.0), (0.0, 10.0)]
    mesh = _quiet(build, [{'coords': ring, 'mat_id': 0}], 1.0, 'quad4',
                  quad_style='structured')
    n = len(mesh['elements'])
    if n != 200:
        failures.append(
            f"structured/20x10 block: {n} elements at a requested size of 1.0, "
            f"expected exactly 200 (20 divisions x 10 divisions). "
            f"{n} == 171 means setTransfiniteCurve was handed the element count "
            f"where it wants the node count (divisions + 1)")
    if not np.all(np.isin(mesh['element_types'], [4, 8, 9])):
        failures.append("structured/20x10 block: a swept rectangle left triangles")


# --------------------------------------------------------------------- (d)
def check_structured_never_worse(build, failures):
    """Sweeping may not make any section's mesh worse than free meshing does."""
    for name, path in GEOMETRIES:
        polys, ts = _load(path)
        free = _metrics(_quiet(build, polys, ts, 'quad4', quad_style='free'), ts)
        stru = _metrics(_quiet(build, polys, ts, 'quad4', quad_style='structured'), ts)
        if stru['ar95'] > free['ar95'] * AR95_TOLERANCE:
            failures.append(f"structured/{name}: AR95 {stru['ar95']:.2f} against "
                            f"{free['ar95']:.2f} free — sweeping made the mesh worse")
        ceiling = max(free['armax'] * ARMAX_TOLERANCE, AR_MAX_LIMIT)
        if stru['armax'] > ceiling:
            failures.append(f"structured/{name}: worst element {stru['armax']:.2f} "
                            f"against {free['armax']:.2f} free, over the {ceiling:.2f} "
                            f"ceiling")
        if stru['gt3_pct'] > free['gt3_pct'] + AR_GT3_LIMIT_PCT:
            failures.append(f"structured/{name}: {stru['gt3_pct']:.2f} % of elements "
                            f"beyond 3:1 against {free['gt3_pct']:.2f} % free")
        lo, hi = MED_OVER_REQ_RANGE
        if not (lo <= stru['med_over_req'] <= hi):
            failures.append(f"structured/{name}: median element "
                            f"{stru['med_over_req']:.3f} x the requested size, "
                            f"outside [{lo}, {hi}]")


# --------------------------------------------------------------------- (e)
def check_element_type_default(build, failures):
    """A caller that states no element type must get quadratic triangles.

    The default is what a scripted run gets, and it is the only place the choice
    is made when nothing states one -- a script that never mentions element_type,
    and a model whose main!D18 is blank, which reaches the builder as None. Both
    have to land on tri6: the stress solvers lock on linear elements and read the
    factor of safety high, so silence must not select them.
    """
    square = [{'coords': [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
               'mat_id': 0}]
    ref = _quiet(build, square, 2.5, 'tri6')
    n6 = len(ref['nodes'])
    for label, mesh in (('omitted', _quiet(build, square, 2.5)),
                        ('None', _quiet(build, square, 2.5, None))):
        codes = sorted({int(v) for v in mesh['element_types']})
        if codes != [6]:
            failures.append(f"element_type {label}: built element type code(s) "
                            f"{codes}, not tri6 (6) — the default must be quadratic")
        elif len(mesh['nodes']) != n6:
            failures.append(f"element_type {label}: {len(mesh['nodes'])} nodes "
                            f"against {n6} for an explicit tri6 request")
    # The explicit linear request still works: it is the seepage choice, and the
    # element-type comparisons depend on being able to ask for it.
    tri3 = _quiet(build, square, 2.5, 'tri3')
    if sorted({int(v) for v in tri3['element_types']}) != [3]:
        failures.append("an explicit element_type='tri3' no longer builds "
                        "3-node triangles")


def run():
    """Returns a list of failure messages; empty means every guard passed."""
    warnings.filterwarnings('ignore')
    from xslope.mesh import build_mesh_from_polygons as build

    failures = []
    check_free_quality(build, failures)
    check_determinism(build, failures)
    check_structured_single_block(build, failures)
    check_structured_levee(build, failures)
    check_structured_never_worse(build, failures)
    check_element_type_default(build, failures)
    return failures


if __name__ == '__main__':
    problems = run()
    if problems:
        print(f"FAIL ({len(problems)})")
        for p in problems:
            print(" -", p)
        sys.exit(1)
    print("quad mesh checks: OK")
