"""The mesh split along a jointed reinforcement line.

A reinforcement line flagged as a joint is not a bar sharing the soil's nodes:
it is a slip surface, and the mesh is torn along it. Every node on the line
becomes three at the same point — a copy carried by the soil above, a copy
carried by the bar, and the original, which the soil below keeps — and a pair of
joint elements spans each bar element, upper soil to bar and bar to lower soil.
The two soil faces rejoin at the line's ends, so a tip is one shared soil node
and a crack tip; the bar still gets its own node there, which is what makes a
sheet end free to pull out.

Nine legs, each on a mesh built for it:

  a. a horizontal sheet inside one material polygon, on tri3 and on tri6. The
     counts: n stations along the curve, 2n - 2 nodes added (every station
     triples except the two buried tips, which only double), and two joint elements per
     bar element — 2(n-1) on tri3, where a station pair is a bar element, and
     n-1 on tri6, where a bar element spans two stations and its joints carry
     three node pairs. The classification: every 2D element above the line holds
     only upper copies, every element below only originals, no element straddles,
     and no 2D element stands on a bar node.
  b. the same sheet crossing a material boundary. The crossing node is a station
     like any other and is tripled like any other, and both sides of it read the
     same line — one interface law per line in phase 1.
  c. a sheet ending on the domain boundary, on a real model. That end is a
     crack reaching the surface: the fan of elements around it is cut in two, so
     it carries two soil faces rather than the one a buried tip carries. The
     boundary condition is built from node coordinates, so every copy at the tip,
     which sits at the same point as the original, must come back with the same
     restraint. This leg measures that through build_fem_data rather than
     assuming it.
  d. an inclined sheet, so nothing in the split depends on the line being
     horizontal.
  e. the reinforced tutorial slope with one line flagged in memory: the mesh
     builds, and with the flag off it carries none of the joint keys and its
     node, element and 1D arrays are the ones the mesher always produced. The
     jointed mesh's first N nodes are the unflagged mesh's nodes unchanged, so
     the split is purely additive.
  f. what the split refuses and what it does not: a jointed line met by a BONDED
     constraint line (the pile case) and a line load's application point sitting
     on a jointed line both raise; two jointed lines that MEET are built, and
     what the junction becomes is test/joint_junction_check.py's subject.
  g. a sheet lying ALONG a material boundary, on tri3 and tri6 — the base
     geotextile at a fill/foundation contact, and the wall sheet on a
     fill/facing-block contact. The line coincides with a zone edge, so the
     mesher must reuse that edge's curve rather than embed a second one on top
     of it; the stations are then the nodes the two zones already share and the
     split classifies the faces by which zone's elements hold them. The leg
     asserts the counts and that every upper copy is held only by elements of
     the zone above and every original only by the zone below.
  h. a sheet lying on a material boundary over part of its length and running on
     into the zone above over the rest, the wall geometry where a sheet leaves
     the facing contact and continues into the fill. The bar and both joints run
     the whole length, and the station at the transition is a station like any
     other.
  i. a sheet ending exactly on a material boundary rather than crossing it, and
     the fifteen sheets of vp088 all flagged at once: the mesh builds, and the
     unflagged mesh is the one the corpus ships.

Every leg that builds a jointed mesh also round-trips it through JSON and
compares every array and record.

Run directly:  PYTHONPATH=. python3 test/joint_mesh_check.py
"""

import contextlib
import io
import os
import sys
import warnings

warnings.filterwarnings('ignore')

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np

from xslope.mesh import (JOINT_MESH_KEYS, add_intersection_points_to_polygons,
                         build_mesh_from_polygons, export_mesh_to_json,
                         extract_constraint_line_geometry,
                         extract_point_constraints, extract_size_regions,
                         get_material_polygons, import_mesh_from_json)

#: The shipped reinforced slope, used for legs c and e. Six horizontal
#: geotextile layers inside a two-material section, all well clear of the
#: domain boundary.
MODEL = os.path.join(_ROOT, 'docs/tutorials/files/xslope_reinforced_slope.xlsx')

#: The shipped RS2 geotextile wall, used for leg i. Fifteen horizontal sheets
#: whose front ends run into the facing-block columns, meshed at a target size of
#: 1.0 — the size that reproduces the corpus's own vp088_mesh.json.
WALL_MODEL = os.path.join(_ROOT,
                          'docs/verification/files/rocscience/vp088.xlsx')

#: The keys the mesher wrote before joints existed. Leg e pins that an unflagged
#: model still writes exactly these and nothing else.
BASE_MESH_KEYS = ('element_materials', 'element_materials_1d', 'element_types',
                  'element_types_1d', 'elements', 'elements_1d', 'nodes')

TOL = 1e-9


def _build(polygons, lines, element_type, **kwargs):
    """A mesh, with the mesher's chatter swallowed."""
    with contextlib.redirect_stdout(io.StringIO()):
        return build_mesh_from_polygons(polygons, 2.0, element_type, lines=lines,
                                        element_size_1d=1.0, **kwargs)


def _roundtrip(mesh, tag, failures):
    """Export the mesh to JSON, read it back, and compare every entry."""
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, 'mesh.json')
        with contextlib.redirect_stdout(io.StringIO()):
            export_mesh_to_json(mesh, path)
            back = import_mesh_from_json(path)
    if set(back) != set(mesh):
        failures.append(f"{tag}: JSON round trip changed the key set, "
                        f"{sorted(set(mesh) ^ set(back))} differ")
        return
    for key, value in mesh.items():
        if isinstance(value, np.ndarray):
            got = np.asarray(back[key])
            if got.shape != value.shape or not np.array_equal(got, value):
                failures.append(f"{tag}: JSON round trip did not reproduce "
                                f"'{key}'")
        elif back[key] != value:
            failures.append(f"{tag}: JSON round trip did not reproduce "
                            f"'{key}'")


def _line_frame(line):
    p1 = np.array(line[0], dtype=float)
    p2 = np.array(line[-1], dtype=float)
    d = p2 - p1
    d = d / np.hypot(d[0], d[1])
    return p1, np.array([-d[1], d[0]])


def _check_split(base, mesh, line, tag, failures, expect_line=1, open_ends=0):
    """The counts and the classification, on one jointed mesh.

    ``open_ends`` is how many of the line's two ends stand on the model's
    external boundary. A buried end is a crack tip: the two soil faces rejoin
    there and the station carries one soil node. An end on the boundary is a
    crack that reaches the surface, and the fan of elements around it is cut in
    two, so it carries two — one node more than a buried end, and the leg states
    which of its ends are which rather than reading it off the split.
    """
    joints = mesh.get('joints')
    if not joints:
        failures.append(f"{tag}: the mesh carries no 'joints' record")
        return
    if len(joints) != 1 or joints[0]['line'] != expect_line:
        failures.append(f"{tag}: expected one jointed line {expect_line}, got "
                        f"{[j['line'] for j in joints]}")
        return
    stations = joints[0]['stations']
    n = len(stations)

    # --- counts -----------------------------------------------------------
    added = len(mesh['nodes']) - len(base['nodes'])
    expect_added = 2 * n - 2 + open_ends
    if added != expect_added:
        failures.append(f"{tag}: {n} stations with {open_ends} end(s) on the "
                        f"boundary added {added} nodes, not {expect_added} "
                        f"(three copies per station; a buried tip's soil faces "
                        f"rejoin, an end on the boundary opens)")
    mats_1d = np.asarray(mesh['element_materials_1d'], dtype=int)
    own_bars = np.where(mats_1d == expect_line)[0]
    n_bar = len(own_bars)
    if len(mesh['elements_joint']) != 2 * n_bar:
        failures.append(f"{tag}: {len(mesh['elements_joint'])} joint elements "
                        f"for {n_bar} bar elements, not {2 * n_bar} (an upper "
                        f"and a lower joint on each)")
    pairs = int(mesh['element_types_joint'][0])
    expected_joints = 2 * (n - 1) if pairs == 2 else n - 1
    if len(mesh['elements_joint']) != expected_joints:
        failures.append(f"{tag}: {len(mesh['elements_joint'])} joint elements "
                        f"on {n} stations at {pairs} node pairs each, not "
                        f"{expected_joints}")
    if set(np.asarray(mesh['element_types_joint']).tolist()) != {pairs}:
        failures.append(f"{tag}: mixed joint element types "
                        f"{sorted(set(np.asarray(mesh['element_types_joint']).tolist()))}")
    sides = np.asarray(mesh['element_side_joint'])
    if int((sides == 1).sum()) != n_bar or int((sides == 0).sum()) != n_bar:
        failures.append(f"{tag}: {int((sides == 1).sum())} upper and "
                        f"{int((sides == 0).sum())} lower joints, not {n_bar} "
                        f"of each")
    if set(np.asarray(mesh['element_materials_joint']).tolist()) != {expect_line}:
        failures.append(f"{tag}: the joint elements do not all read line "
                        f"{expect_line}")

    # --- the three copies stand at one point ------------------------------
    nodes = np.asarray(mesh['nodes'], dtype=float)
    for k, (up, bar, low) in enumerate(stations):
        for name, nid in (('upper', up), ('bar', bar)):
            if np.max(np.abs(nodes[nid] - nodes[low])) > TOL:
                failures.append(f"{tag}: station {k}'s {name} copy is not at "
                                f"the same point as the original")
    tips = (0, len(stations) - 1)
    n_open = 0
    for k in tips:
        up, bar, low = stations[k]
        if up != low:
            n_open += 1
        if bar == low:
            failures.append(f"{tag}: the bar shares the soil's node at tip "
                            f"station {k}; the bar keeps its own end node")
    if n_open != open_ends:
        failures.append(f"{tag}: {n_open} of the line's ends carry two soil "
                        f"faces, not {open_ends}; a buried tip must rejoin and "
                        f"an end on the external boundary must not")
    for k in range(1, len(stations) - 1):
        up, bar, low = stations[k]
        if len({up, bar, low}) != 3:
            failures.append(f"{tag}: interior station {k} is not tripled "
                            f"({up}, {bar}, {low})")

    # --- classification ---------------------------------------------------
    upper_ids = {s[0] for s in stations[1:-1]}
    bar_ids = {s[1] for s in stations}
    lower_ids = {s[2] for s in stations[1:-1]}
    p1, normal = _line_frame(line)
    elements = np.asarray(mesh['elements'], dtype=int)
    types = np.asarray(mesh['element_types'], dtype=int)
    n_straddle = n_wrong = n_bar_in_2d = 0
    for i in range(len(elements)):
        et = int(types[i])
        enodes = [int(elements[i, k]) for k in range(et)]
        if any(nid in bar_ids for nid in enodes):
            n_bar_in_2d += 1
        on_line = [nid for nid in enodes
                   if nid in upper_ids or nid in lower_ids]
        if not on_line:
            continue
        n_corner = 3 if et in (3, 6) else 4
        centroid = nodes[enodes[:n_corner], :2].mean(axis=0)
        side = float((centroid - p1) @ normal)
        signs = set()
        for nid in enodes:
            if nid in upper_ids or nid in lower_ids or nid in bar_ids:
                continue
            s = float((nodes[nid, :2] - p1) @ normal)
            if abs(s) > 1e-9:
                signs.add(1 if s > 0 else -1)
        if len(signs) > 1:
            n_straddle += 1
        wanted = upper_ids if side > 0 else lower_ids
        if any(nid not in wanted for nid in on_line):
            n_wrong += 1
    if n_straddle:
        failures.append(f"{tag}: {n_straddle} 2D elements straddle the joint")
    if n_wrong:
        failures.append(f"{tag}: {n_wrong} 2D elements hold a node from the "
                        f"wrong face of the joint")
    if n_bar_in_2d:
        failures.append(f"{tag}: {n_bar_in_2d} 2D elements stand on a bar node; "
                        f"the bar carries its own node set")

    # --- the bar stands on its own nodes ----------------------------------
    e1d = np.asarray(mesh['elements_1d'], dtype=int)
    t1d = np.asarray(mesh['element_types_1d'], dtype=int)
    for i in own_bars:
        for k in range(int(t1d[i])):
            if int(e1d[i, k]) not in bar_ids:
                failures.append(f"{tag}: bar element {i} stands on node "
                                f"{int(e1d[i, k])}, which is not a bar node")
                break

    _roundtrip(mesh, tag, failures)
    return n


def _leg_a(failures):
    """A horizontal sheet inside one polygon, on tri3 and tri6."""
    ring = [(0, 0), (20, 0), (20, 10), (0, 10)]
    line = [(5.0, 5.0), (15.0, 5.0)]
    polys = [{'coords': list(ring), 'mat_id': 0}]
    counts = {}
    for et in ('tri3', 'tri6'):
        base = _build([{'coords': list(ring), 'mat_id': 0}], [line], et)
        mesh = _build([dict(p) for p in polys], [line], et, joint_lines=[0])
        counts[et] = _check_split(base, mesh, line, f"a/{et}", failures)
    return counts


def _leg_b(failures):
    """The same sheet crossing a material boundary at x = 10."""
    left = [(0, 0), (10, 0), (10, 10), (0, 10)]
    right = [(10, 0), (20, 0), (20, 10), (10, 10)]
    line = [(5.0, 5.0), (15.0, 5.0)]
    rings = add_intersection_points_to_polygons([left, right], [line])
    polys = [{'coords': rings[0], 'mat_id': 0}, {'coords': rings[1], 'mat_id': 1}]
    counts = {}
    for et in ('tri3', 'tri6'):
        base = _build([dict(p) for p in polys], [line], et)
        mesh = _build([dict(p) for p in polys], [line], et, joint_lines=[0])
        counts[et] = _check_split(base, mesh, line, f"b/{et}", failures)
        stations = mesh['joints'][0]['stations']
        nodes = np.asarray(mesh['nodes'], dtype=float)
        crossing = [k for k, s in enumerate(stations)
                    if abs(nodes[s[2], 0] - 10.0) < 1e-9]
        if len(crossing) != 1:
            failures.append(f"b/{et}: the material crossing at x = 10 is not a "
                            f"single station ({len(crossing)} found)")
        elif crossing[0] in (0, len(stations) - 1):
            failures.append(f"b/{et}: the material crossing landed on a tip")
        else:
            up, bar, low = stations[crossing[0]]
            if len({up, bar, low}) != 3:
                failures.append(f"b/{et}: the material crossing node is not "
                                f"tripled like any other station")
        mats = set(np.asarray(mesh['element_materials']).tolist())
        if len(mats) != 2:
            failures.append(f"b/{et}: the split lost a material zone; "
                            f"materials are {sorted(mats)}")
    return counts


def _leg_c(failures):
    """A sheet ending on the domain boundary, and the fixity of its copies."""
    from xslope.fem import build_fem_data
    from xslope.fileio import build_reinforce_lines, load_slope_data
    with contextlib.redirect_stdout(io.StringIO()):
        sd = load_slope_data(MODEL)
    x_left = min(x for x, _ in sd['domain_polygon'].exterior.coords)
    row = dict(sd['reinforcement_lines'][0])
    row.update({'x1': x_left, 'y1': -5.0, 'x2': x_left + 20.0, 'y2': -5.0,
                'label': 'sheet on the boundary', 'tend1': 0.0, 'tend2': 0.0})
    sd['reinforcement_lines'].append(row)
    sd['reinforce_lines'] = build_reinforce_lines(sd['reinforcement_lines'])
    lines, _n_r, _n_p = extract_constraint_line_geometry(sd)
    polys = get_material_polygons(sd, reinf_lines=lines)
    kw = dict(target_size=sd['target_size'], element_type=sd['element_type'],
              lines=lines, element_size_1d=sd.get('element_size_1d'),
              point_constraints=extract_point_constraints(sd),
              size_regions=extract_size_regions(sd))
    with contextlib.redirect_stdout(io.StringIO()):
        base = build_mesh_from_polygons(polys, **kw)
        mesh = build_mesh_from_polygons(polys, joint_lines=[len(lines) - 1], **kw)
    n = _check_split(base, mesh, lines[-1], 'c', failures,
                     expect_line=len(lines), open_ends=1)
    stations = mesh['joints'][0]['stations']
    nodes = np.asarray(mesh['nodes'], dtype=float)
    if abs(nodes[stations[0][2], 0] - x_left) > 1e-9:
        failures.append("c: the sheet's end 1 is not on the domain boundary")
    with contextlib.redirect_stdout(io.StringIO()):
        fem_data = build_fem_data(sd, mesh)
    bc = np.asarray(fem_data['bc_type'])
    if int(bc[stations[0][2]]) == 0:
        failures.append("c: the tip on the domain boundary carries no boundary "
                        "condition, so the leg cannot see whether the copies "
                        "inherit it")
    disagreeing = [(k, int(bc[u]), int(bc[b]), int(bc[l]))
                   for k, (u, b, l) in enumerate(stations)
                   if not (bc[u] == bc[b] == bc[l])]
    if disagreeing:
        failures.append(f"c: {len(disagreeing)} stations carry different "
                        f"boundary conditions across their copies, e.g. "
                        f"{disagreeing[0]}")
    return n, int(bc[stations[0][2]])


def _leg_d(failures):
    """An inclined sheet."""
    ring = [(0, 0), (20, 0), (20, 10), (0, 10)]
    line = [(4.0, 3.0), (16.0, 8.0)]
    counts = {}
    for et in ('tri3', 'tri6'):
        base = _build([{'coords': list(ring), 'mat_id': 0}], [line], et)
        mesh = _build([{'coords': list(ring), 'mat_id': 0}], [line], et,
                      joint_lines=[0])
        counts[et] = _check_split(base, mesh, line, f"d/{et}", failures)
    return counts


def _leg_e(failures):
    """The reinforced tutorial slope, one line flagged in memory."""
    from xslope.fileio import load_slope_data
    with contextlib.redirect_stdout(io.StringIO()):
        sd = load_slope_data(MODEL)
    lines, _n_r, _n_p = extract_constraint_line_geometry(sd)
    polys = get_material_polygons(sd, reinf_lines=lines)
    kw = dict(target_size=sd['target_size'], element_type=sd['element_type'],
              lines=lines, element_size_1d=sd.get('element_size_1d'),
              point_constraints=extract_point_constraints(sd),
              size_regions=extract_size_regions(sd))
    flagged = 2                                  # the third layer, mid-height
    with contextlib.redirect_stdout(io.StringIO()):
        base = build_mesh_from_polygons(polys, **kw)
        mesh = build_mesh_from_polygons(polys, joint_lines={flagged: {
            'tend1': 12.0}}, **kw)

    if tuple(sorted(base)) != BASE_MESH_KEYS:
        failures.append(f"e: an unflagged model wrote {sorted(base)}, not the "
                        f"{len(BASE_MESH_KEYS)} keys the mesher has always "
                        f"written — the joint construction is not inert")
    for key in JOINT_MESH_KEYS:
        if key in base:
            failures.append(f"e: an unflagged model wrote the joint key "
                            f"'{key}'")

    # the tutorial's geotextiles start ON the slope face, so end 1 is a crack
    # that reaches the surface and end 2 is buried in the fill.
    n = _check_split(base, mesh, lines[flagged], 'e', failures,
                     expect_line=flagged + 1, open_ends=1)

    # The split only appends: the unflagged mesh's nodes come back unchanged.
    n_base = len(base['nodes'])
    if not np.array_equal(np.asarray(mesh['nodes'])[:n_base],
                          np.asarray(base['nodes'])):
        failures.append("e: the split moved a node of the unflagged mesh; it "
                        "must only append copies")
    if not np.array_equal(np.asarray(mesh['element_types']),
                          np.asarray(base['element_types'])):
        failures.append("e: the split changed the 2D element types")
    if not np.array_equal(np.asarray(mesh['element_materials']),
                          np.asarray(base['element_materials'])):
        failures.append("e: the split changed the 2D element materials")

    ties = mesh['ties']
    if len(ties) != 1 or ties[0]['end'] != 1 or ties[0]['capacity'] != 12.0:
        failures.append(f"e: the tied end 1 (Tend1 = 12) produced {ties}")
    else:
        stations = mesh['joints'][0]['stations']
        if (ties[0]['bar_node'] != stations[0][1]
                or ties[0]['soil_node'] != stations[0][2]):
            failures.append("e: the tie does not join the bar's end node to "
                            "the soil node at the same point")
    return n, len(base['nodes']), len(mesh['nodes'])


def _leg_f(failures):
    """The geometries the split still refuses — and the one it no longer does.

    Two jointed lines that MEET are built, not refused: the split gives the
    shared node one copy per wedge of material around it (test/joint_junction_
    check.py). What stays refused is a jointed line met by a BONDED member, whose
    element would keep one wedge's copy and lose the material on the other side,
    and a line load standing on a jointed line, which has no single node to act
    on.
    """
    ring = [(0, 0), (20, 0), (20, 10), (0, 10)]
    sheet = [(5.0, 5.0), (15.0, 5.0)]
    crossing = [(10.0, 2.0), (10.0, 8.0)]
    refusals = [
        ("a jointed line meeting another constraint line (a pile)",
         'touches another constraint line',
         dict(lines=[sheet, crossing], joint_lines=[0])),
        ("a line load's application point on a jointed line",
         'carries a point constraint',
         dict(lines=[sheet], joint_lines=[0],
              point_constraints=[(10.0, 5.0)])),
    ]
    raised = 0
    for what, phrase, kwargs in refusals:
        lines = kwargs.pop('lines')
        try:
            _build([{'coords': list(ring), 'mat_id': 0}], lines, 'tri6',
                   **kwargs)
        except ValueError as exc:
            if phrase not in str(exc):
                failures.append(f"f: {what} raised, but for another reason: "
                                f"{exc}")
            else:
                raised += 1
        else:
            failures.append(f"f: the mesher accepted {what}; it must refuse it")
    try:
        _build([{'coords': list(ring), 'mat_id': 0}], [sheet, crossing], 'tri6',
               joint_lines=[0, 1])
    except ValueError as exc:
        failures.append(f"f: two jointed lines that meet were refused: {exc}")
    else:
        raised += 1
    return raised


def _material_sides(mesh, line, tag, failures):
    """On a line that lies along a zone edge, the two faces belong to different
    materials: every upper copy must be held only by elements of the zone above
    and every original only by the zone below."""
    stations = mesh['joints'][0]['stations']
    upper_ids = {s[0] for s in stations[1:-1]}
    lower_ids = {s[2] for s in stations[1:-1]}
    nodes = np.asarray(mesh['nodes'], dtype=float)
    elements = np.asarray(mesh['elements'], dtype=int)
    types = np.asarray(mesh['element_types'], dtype=int)
    mats = np.asarray(mesh['element_materials'], dtype=int)
    p1, normal = _line_frame(line)
    seen_upper, seen_lower = set(), set()
    for i in range(len(elements)):
        et = int(types[i])
        enodes = [int(elements[i, k]) for k in range(et)]
        n_corner = 3 if et in (3, 6) else 4
        side = float((nodes[enodes[:n_corner], :2].mean(axis=0) - p1) @ normal)
        if any(nid in upper_ids for nid in enodes):
            seen_upper.add(int(mats[i]))
            if side < 0:
                failures.append(f"{tag}: an element below the line holds an "
                                f"upper copy")
        if any(nid in lower_ids for nid in enodes):
            seen_lower.add(int(mats[i]))
            if side > 0:
                failures.append(f"{tag}: an element above the line holds the "
                                f"original (lower) node")
    if len(seen_upper) != 1 or len(seen_lower) != 1:
        failures.append(f"{tag}: the two faces are not one material each — "
                        f"above {sorted(seen_upper)}, below {sorted(seen_lower)}")
    elif seen_upper == seen_lower:
        failures.append(f"{tag}: both faces read material {seen_upper}; the "
                        f"line is supposed to lie on the boundary between two")
    return sorted(seen_lower)[0], sorted(seen_upper)[0]


def _leg_g(failures):
    """A sheet lying along the boundary between two materials."""
    lower = [(0, 0), (20, 0), (20, 5), (0, 5)]
    upper = [(0, 5), (20, 5), (20, 10), (0, 10)]
    line = [(5.0, 5.0), (15.0, 5.0)]
    counts = {}
    for et in ('tri3', 'tri6'):
        polys = [{'coords': list(lower), 'mat_id': 0},
                 {'coords': list(upper), 'mat_id': 1}]
        base = _build([dict(p) for p in polys], [line], et)
        mesh = _build([dict(p) for p in polys], [line], et, joint_lines=[0])
        n = _check_split(base, mesh, line, f"g/{et}", failures)
        counts[et] = n
        if n and mesh.get('joints'):
            below, above = _material_sides(mesh, line, f"g/{et}", failures)
            if (below, above) != (1, 2):
                failures.append(f"g/{et}: the faces read materials {below} "
                                f"below and {above} above, not 1 and 2")
        # The bar runs the whole line at the stated 1D size, which is what the
        # doubled-curve failure used to cost: the fallback mesher rebuilt the
        # section at the global size instead.
        mats_1d = np.asarray(mesh['element_materials_1d'], dtype=int)
        nodes = np.asarray(mesh['nodes'], dtype=float)
        e1d = np.asarray(mesh['elements_1d'], dtype=int)
        own = np.where(mats_1d == 1)[0]
        xs = sorted(float(nodes[int(e1d[i, k]), 0])
                    for i in own for k in range(2))
        if abs(min(xs) - 5.0) > TOL or abs(max(xs) - 15.0) > TOL:
            failures.append(f"g/{et}: the bar runs from x = {min(xs):.3f} to "
                            f"{max(xs):.3f}, not the stated 5 to 15")
        if len(own) != 10:
            failures.append(f"g/{et}: {len(own)} bar elements over 10 m at a 1D "
                            f"size of 1.0; the line was meshed at the wrong size")
    return counts


def _leg_h(failures):
    """A sheet on a material boundary over part of its length, interior to the
    zone above over the rest."""
    lower = [(0, 0), (20, 0), (20, 3), (12, 3), (12, 5), (0, 5)]
    upper = [(0, 5), (12, 5), (12, 3), (20, 3), (20, 10), (0, 10)]
    line = [(4.0, 5.0), (18.0, 5.0)]
    counts = {}
    for et in ('tri3', 'tri6'):
        polys = [{'coords': list(lower), 'mat_id': 0},
                 {'coords': list(upper), 'mat_id': 1}]
        base = _build([dict(p) for p in polys], [line], et)
        mesh = _build([dict(p) for p in polys], [line], et, joint_lines=[0])
        n = _check_split(base, mesh, line, f"h/{et}", failures)
        counts[et] = n
        if not n or not mesh.get('joints'):
            continue
        nodes = np.asarray(mesh['nodes'], dtype=float)
        stations = mesh['joints'][0]['stations']
        xs = [float(nodes[s[2], 0]) for s in stations]
        transition = [k for k, x in enumerate(xs) if abs(x - 12.0) < 1e-9]
        if len(transition) != 1:
            failures.append(f"h/{et}: the point where the sheet leaves the "
                            f"boundary, x = 12, is not a single station "
                            f"({len(transition)} found)")
        elif transition[0] in (0, len(stations) - 1):
            failures.append(f"h/{et}: the transition landed on a tip")
        else:
            up, bar, low = stations[transition[0]]
            if len({up, bar, low}) != 3:
                failures.append(f"h/{et}: the transition station is not tripled")
        # Both stretches carry bar elements, so the sheet is continuous across
        # the point where it leaves the boundary.
        mats_1d = np.asarray(mesh['element_materials_1d'], dtype=int)
        e1d = np.asarray(mesh['elements_1d'], dtype=int)
        own = np.where(mats_1d == 1)[0]
        on_edge = sum(1 for i in own
                      if max(float(nodes[int(e1d[i, k]), 0]) for k in range(2)) <= 12.0 + 1e-9)
        interior = len(own) - on_edge
        if on_edge == 0 or interior == 0:
            failures.append(f"h/{et}: {on_edge} bar elements on the boundary "
                            f"and {interior} inside the zone above; the sheet "
                            f"must be continuous across x = 12")
        if len(own) != 14:
            failures.append(f"h/{et}: {len(own)} bar elements over 14 m at a 1D "
                            f"size of 1.0; the line was meshed at the wrong size")
    return counts


def _leg_i(failures):
    """A sheet ending on a material boundary, and vp088's fifteen sheets."""
    # --- the end that lands on a zone edge without crossing it ------------
    lower = [(0, 0), (20, 0), (20, 5), (0, 5)]
    upper = [(0, 5), (20, 5), (20, 10), (0, 10)]
    line = [(5.0, 8.0), (10.0, 5.0)]
    ends = {}
    for et in ('tri3', 'tri6'):
        polys = [{'coords': list(lower), 'mat_id': 0},
                 {'coords': list(upper), 'mat_id': 1}]
        base = _build([dict(p) for p in polys], [line], et)
        mesh = _build([dict(p) for p in polys], [line], et, joint_lines=[0])
        ends[et] = _check_split(base, mesh, line, f"i-end/{et}", failures)
        if not mesh.get('joints'):
            continue
        nodes = np.asarray(mesh['nodes'], dtype=float)
        tip = mesh['joints'][0]['stations'][-1]
        if abs(float(nodes[tip[2], 1]) - 5.0) > TOL:
            failures.append(f"i-end/{et}: the sheet's end 2 is not on the "
                            f"material boundary at y = 5")
        n_bar = int((np.asarray(mesh['element_materials_1d'], dtype=int) == 1).sum())
        if n_bar != 6:
            failures.append(f"i-end/{et}: {n_bar} bar elements over 5.83 m at a "
                            f"1D size of 1.0; the line was meshed at the wrong "
                            f"size")

    # --- vp088: fifteen wall sheets, all flagged --------------------------
    from xslope.fileio import load_slope_data
    if not os.path.exists(WALL_MODEL):
        failures.append(f"i: missing model {WALL_MODEL}")
        return ends, None
    with contextlib.redirect_stdout(io.StringIO()):
        sd = load_slope_data(WALL_MODEL)
    lines, n_r, _n_p = extract_constraint_line_geometry(sd)
    polys = get_material_polygons(sd, reinf_lines=lines)
    kw = dict(target_size=sd['target_size'] or 1.0,
              element_type=sd['element_type'], lines=lines,
              element_size_1d=sd.get('element_size_1d'),
              point_constraints=extract_point_constraints(sd),
              size_regions=extract_size_regions(sd))
    with contextlib.redirect_stdout(io.StringIO()):
        base = build_mesh_from_polygons(polys, **kw)
        mesh = build_mesh_from_polygons(polys, joint_lines=list(range(n_r)), **kw)

    joints = mesh.get('joints') or []
    if len(joints) != n_r:
        failures.append(f"i: {len(joints)} of {n_r} sheets were split")
        return ends, None
    added = len(mesh['nodes']) - len(base['nodes'])
    want = sum(2 * len(j['stations']) - 2 for j in joints)
    if added != want:
        failures.append(f"i: vp088's {n_r} sheets added {added} nodes, not "
                        f"{want}")
    mats_1d = np.asarray(mesh['element_materials_1d'], dtype=int)
    n_joint = 0
    for j in joints:
        n_bar = int((mats_1d == j['line']).sum())
        own = int((np.asarray(mesh['element_materials_joint']) == j['line']).sum())
        n_joint += own
        if own != 2 * n_bar:
            failures.append(f"i: line {j['line']} has {own} joint elements for "
                            f"{n_bar} bar elements, not {2 * n_bar}")
    if len(mesh['elements_joint']) != n_joint:
        failures.append(f"i: {len(mesh['elements_joint'])} joint elements, "
                        f"{n_joint} accounted for by line")
    # The unflagged wall mesh is the one the corpus ships.
    shipped = os.path.join(os.path.dirname(WALL_MODEL), 'vp088_mesh.json')
    if os.path.exists(shipped):
        with contextlib.redirect_stdout(io.StringIO()):
            ship = import_mesh_from_json(shipped)
        for key in ('nodes', 'elements', 'element_types', 'element_materials',
                    'elements_1d', 'element_types_1d', 'element_materials_1d'):
            if not np.array_equal(np.asarray(base[key]), np.asarray(ship[key])):
                failures.append(f"i: the unflagged vp088 mesh no longer "
                                f"reproduces the shipped '{key}'")
    _roundtrip(mesh, 'i', failures)
    return ends, (len(base['nodes']), len(mesh['nodes']),
                  sum(len(j['stations']) for j in joints), n_joint)


def main():
    failures = []
    if not os.path.exists(MODEL):
        return [f"missing model {MODEL}"]

    a = _leg_a(failures)
    b = _leg_b(failures)
    c_n, c_bc = _leg_c(failures)
    d = _leg_d(failures)
    e_n, e_base, e_split = _leg_e(failures)
    n_refused = _leg_f(failures)
    g = _leg_g(failures)
    h = _leg_h(failures)
    i_ends, i_wall = _leg_i(failures)

    if not failures:
        print(f"a: horizontal sheet, {a['tri3']} stations on tri3 / "
              f"{a['tri6']} on tri6")
        print(f"b: crossing a material boundary, {b['tri3']} / {b['tri6']} "
              f"stations")
        print(f"c: sheet on the domain boundary, {c_n} stations; the copies at "
              f"the open end all carry boundary condition {c_bc}")
        print(f"d: inclined sheet, {d['tri3']} / {d['tri6']} stations")
        print(f"e: {os.path.basename(MODEL)} with one line flagged, {e_n} "
              f"stations, {e_base} -> {e_split} nodes; unflagged it writes no "
              f"joint key")
        print(f"f: {n_refused - 1} of 2 refused geometries raised; two jointed "
              f"lines that meet were built")
        print(f"g: sheet along a material boundary, {g['tri3']} / {g['tri6']} "
              f"stations")
        print(f"h: sheet part on the boundary, part inside the zone above, "
              f"{h['tri3']} / {h['tri6']} stations")
        print(f"i: sheet ending on a material boundary, {i_ends['tri3']} / "
              f"{i_ends['tri6']} stations")
        if i_wall:
            print(f"i: {os.path.basename(WALL_MODEL)} with all fifteen sheets "
                  f"flagged, {i_wall[2]} stations, {i_wall[0]} -> {i_wall[1]} "
                  f"nodes, {i_wall[3]} joint elements; unflagged it is the mesh "
                  f"the corpus ships")
    return failures


def run():
    """Failures as a list, for run_tests.py."""
    return main()


def _cli():
    failures = main()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nThe mesh splits along a jointed line into three node sets, the "
          "soil faces rejoin at the ends, and nothing changes when no line is "
          "flagged.")


if __name__ == "__main__":
    _cli()
