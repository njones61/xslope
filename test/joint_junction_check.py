"""Joints that meet: the mesh split at a junction, and three closed forms.

A jointed line no longer has to stand alone. Where jointed lines meet — one
ending on another (a T), two crossing (an X), two sharing an endpoint (an L),
two running end to end (a chain) — and where a jointed line reaches the model's
external boundary (a crack to the surface), the split copies the shared node
ONCE PER WEDGE of material around it. A wedge is a connected group of the 2D
elements standing on that node, two of them belonging to the same wedge when
they share an edge at the node and that edge is not part of a jointed line. Away
from junctions nothing changes, so every mesh that carried no junction is the
mesh it was.

The wedge counts this check pins, each one countable by hand from the drawing:

  ordinary station     2   the material above the line and the material below
  buried crack tip     1   the ring of elements is cut in one place and holds
  end on the boundary  2   the fan of elements is cut in two
  chain / L corner     2   two rays leave the node
  T                    3
  X                    4

Six mesh fixtures (a-f), each on tri3 and tri6, assert for every junction that
the number of soil copies at the point equals the wedge count; that the elements
in one angular sector all hold the same copy and no two sectors share one; that
each bar keeps its own node there; that every joint element's two sides stand at
one point and are held by the elements on opposite sides of its own line; and
that the mesh JSON round trip reproduces every array and record. Leg g is the
refusals that survive.

Five more fixtures (h-m) are TERMINATIONS: a joint that ends on another partway
along its segment, at a point that is a vertex of the ending line and of nothing
else. Whether that point is one point is a question about arithmetic rather than
about topology — a tip stated to six decimals lands a part in 10^7 off the trace
it belongs on, and a tip computed from the same angle and spacing as the line it
stops on still misses it by a rounding — so the mesher pulls an end that close
onto the line it stops on and gives that line the same point. Each fixture states
what the snap does to it before gmsh runs: how many ends move, and that none
moves further than the tolerance. The exact fixtures (a-f) must have no end
moved at all, which is the guard that the rule is inert on everything that
already worked.

  h. T, the tip a rounding off the line      3 wedges
  i. T, the tip 1e-7 off the line            3, and the mesh is fixture a's
  j. two tips on one segment                 3 at each
  k. the stepped base, three columns         3 at each step, 2 at each corner
  m. L, the tip 1e-7 off the through line's end   2, and the mesh is fixture c's
  l. a tip on an existing vertex             2, nothing snapped
  n. a joint ending on a zone edge AT a twelve-decimal vertex of it

No 2D element in any of them may name one node twice. A collapsed element is what
a pair of geometry points a rounding apart becomes once the mesher has merged
them as duplicate nodes, and it carries the same edge twice — which is how an
edge on a joint comes back with four elements standing on it rather than two.

Three closed-form rows then run the split through the finite element engine:

  i.   R2's block on a plane, CUT by a vertical joint into two blocks — a T on
       the base joint and a crack to the ground surface. Each block rides its
       own stretch of the same plane, so the strength reduction must still
       return tan phi_j / tan beta.
  ii.  the same slab as a two-course STACK: a joint parallel to the plane at
       mid-thickness, crossed by the vertical cut (an X), so four blocks ride
       three interfaces of the same friction. Same closed form.
  iii. an L-shaped joint around a rectangular ELASTIC block: a horizontal joint
       under it and a vertical joint up its back face, meeting at the corner.
       The block is driven by a horizontal seismic coefficient k -- a known
       lateral load -- with weightless space behind the back face and in front,
       so the sliding-block calculation is

           driving  k W        resisting  (W) tan phi_base

       and the block stands while k < tan phi_base. Reducing the interface
       strength by F, the block stands while tan phi_base / F > k, so

           FS = tan phi_base / k.

       The back-face joint carries the corner: down-slope of it the block leaves
       the corner, which it cannot do while the corner node is one node. Its own
       normal traction goes to zero — the face opens — and that is measured, not
       assumed.

Run directly:  PYTHONPATH=. python3 test/joint_junction_check.py
"""

import contextlib
import copy
import io
import math
import os
import sys
import time
import warnings

warnings.filterwarnings('ignore')

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
from shapely.geometry import LineString, Point, Polygon

from xslope.fem import build_fem_data, solve_fem, solve_ssrm
from xslope.fileio import build_reinforce_lines, load_slope_data
from xslope.mesh import (_joint_line_tol, _snap_joint_line_ends,
                         add_intersection_points_to_polygons,
                         build_mesh_from_polygons, export_mesh_to_json,
                         import_mesh_from_json)

#: A shipped FEM model, for the boilerplate every slope_data carries (units,
#: gamma_water, solver options). Every row replaces geometry and materials.
BASE_FILE = os.path.join(_ROOT, 'docs', 'fem', 'files', 'xslope_griffiths1.xlsx')

GAMMA = 20.0
E_SOIL = 30000.0
NU_SOIL = 0.3
TOL = 1e-9

_base_sd = None


# ---------------------------------------------------------------------------
# mesh fixtures
# ---------------------------------------------------------------------------

def _build(lines, joint_lines, element_type, ring=None, mats=None,
           target_size=2.0, s1d=1.0):
    """A mesh on a rectangular block, with the mesher's chatter swallowed."""
    if mats is None:
        ring = ring or [(0.0, 0.0), (20.0, 0.0), (20.0, 10.0), (0.0, 10.0)]
        polys = [{'coords': list(ring), 'mat_id': 0}]
    else:
        polys = [{'coords': list(r), 'mat_id': i} for r, i in mats]
    with contextlib.redirect_stdout(io.StringIO()):
        return build_mesh_from_polygons(
            polys, target_size, element_type,
            lines=[[tuple(p) for p in ln] for ln in lines],
            element_size_1d=s1d, joint_lines=joint_lines)


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
            failures.append(f"{tag}: JSON round trip did not reproduce '{key}'")


def _at(nodes, p):
    """Every node id standing at point ``p``."""
    q = np.asarray(nodes, dtype=float)[:, :2]
    return [int(i) for i in np.flatnonzero(
        (np.abs(q[:, 0] - p[0]) < 1e-8) & (np.abs(q[:, 1] - p[1]) < 1e-8))]


def _bar_nodes(mesh):
    e1d = np.asarray(mesh['elements_1d'], dtype=int)
    t1d = np.asarray(mesh['element_types_1d'], dtype=int)
    out = set()
    for i in range(len(e1d)):
        for k in range(int(t1d[i])):
            out.add(int(e1d[i, k]))
    return out


def _sector_of(nodes, elements, element_types, ei, p, rays):
    """Which angular sector between the joint rays at ``p`` element ``ei`` is in.

    The sectors are read from the LINE DIRECTIONS at the junction, not from the
    mesh's own connectivity, so this is an independent statement about the split:
    the elements the drawing puts in one wedge must hold one copy of the node.
    """
    et = int(element_types[ei])
    n_corner = 3 if et in (3, 6) else 4
    cen = np.mean([np.asarray(nodes[int(elements[ei, k])][:2], dtype=float)
                   for k in range(n_corner)], axis=0)
    a = math.atan2(cen[1] - p[1], cen[0] - p[0]) % (2 * math.pi)
    ordered = sorted(r % (2 * math.pi) for r in rays)
    for k in range(len(ordered)):
        lo = ordered[k]
        hi = ordered[(k + 1) % len(ordered)]
        span = (hi - lo) % (2 * math.pi) or 2 * math.pi
        if ((a - lo) % (2 * math.pi)) < span - 1e-12:
            return k
    return len(ordered) - 1


def _check_junction(mesh, p, rays, n_wedges, n_bars, tag, failures):
    """One junction: the copies, the wedges the elements land in, the bars."""
    nodes = np.asarray(mesh['nodes'], dtype=float)
    elements = np.asarray(mesh['elements'], dtype=int)
    types = np.asarray(mesh['element_types'], dtype=int)
    here = _at(nodes, p)
    bars = _bar_nodes(mesh)
    soil_copies = [n for n in here if n not in bars]
    bar_copies = [n for n in here if n in bars]
    if len(soil_copies) != n_wedges:
        failures.append(f"{tag}: {len(soil_copies)} soil copies at {p}, not "
                        f"{n_wedges} — one per wedge of material there")
    if len(bar_copies) != n_bars:
        failures.append(f"{tag}: {len(bar_copies)} bar nodes at {p}, not "
                        f"{n_bars} — each jointed line keeps its own")

    # every element standing at the point, by the sector the drawing puts it in
    by_sector = {}
    for ei in range(len(elements)):
        et = int(types[ei])
        hit = [int(elements[ei, k]) for k in range(et)
               if int(elements[ei, k]) in here]
        if not hit:
            continue
        if len(set(hit)) != 1:
            failures.append(f"{tag}: 2D element {ei} holds {len(set(hit))} "
                            f"different copies of the node at {p}")
            continue
        s = _sector_of(nodes, elements, types, ei, p, rays)
        by_sector.setdefault(s, set()).add(hit[0])
    if len(by_sector) != n_wedges:
        failures.append(f"{tag}: the elements at {p} fall in {len(by_sector)} "
                        f"angular sectors, not the {n_wedges} the joint lines "
                        f"cut there")
    used = []
    for s, ids in sorted(by_sector.items()):
        if len(ids) != 1:
            failures.append(f"{tag}: the elements in sector {s} at {p} hold "
                            f"{len(ids)} different copies; a wedge holds one")
        used.extend(ids)
    if len(set(used)) != len(by_sector):
        failures.append(f"{tag}: two wedges at {p} share a copy "
                        f"({sorted(used)}) — the split did not separate them")


def _check_pairs(mesh, tag, failures, lines):
    """Every joint element: its node pairs stand at one point, its soil edge is
    held by the 2D element on its own side of its own line, and the other side
    of every pair is the bar's node.

    The soil edge of a joint element is the pair of soil nodes its first two
    columns name (the upper joint) or its last two (the lower one). Exactly one
    2D element stands on that edge, and it must be on that side of the line — the
    one exception is a buried crack tip, where the two soil faces are one node
    and the edge is shared by the elements on both sides.
    """
    nodes = np.asarray(mesh['nodes'], dtype=float)
    elements = np.asarray(mesh['elements'], dtype=int)
    types = np.asarray(mesh['element_types'], dtype=int)
    conn = np.asarray(mesh['elements_joint'], dtype=int)
    pairs = np.asarray(mesh['element_types_joint'], dtype=int)
    mats = np.asarray(mesh['element_materials_joint'], dtype=int)
    side = np.asarray(mesh['element_side_joint'], dtype=int)
    bars = _bar_nodes(mesh)

    tips = set()
    for rec in mesh.get('joints') or []:
        for st in rec['stations']:
            if int(st[0]) == int(st[2]):
                tips.add(int(st[0]))

    on_edge = {}
    for ei in range(len(elements)):
        et = int(types[ei])
        n_corner = 3 if et in (3, 6) else 4
        for k in range(n_corner):
            a = int(elements[ei, k])
            b = int(elements[ei, (k + 1) % n_corner])
            on_edge.setdefault((min(a, b), max(a, b)), []).append(ei)

    n_off = n_side = n_bar = n_edge = 0
    for i in range(len(conn)):
        li = int(mats[i]) - 1
        p1 = np.asarray(lines[li][0], dtype=float)
        p2 = np.asarray(lines[li][-1], dtype=float)
        d = (p2 - p1) / np.hypot(*(p2 - p1))
        normal = np.array([-d[1], d[0]])
        for k in range(int(pairs[i])):
            a, b = int(conn[i, k]), int(conn[i, 3 + k])
            if np.max(np.abs(nodes[a][:2] - nodes[b][:2])) > TOL:
                n_off += 1
            bar_node = b if int(side[i]) == 1 else a
            if bar_node not in bars:
                n_bar += 1
        cols = (0, 1) if int(side[i]) == 1 else (3, 4)
        want = 1.0 if int(side[i]) == 1 else -1.0
        s0, s1 = int(conn[i, cols[0]]), int(conn[i, cols[1]])
        members = on_edge.get((min(s0, s1), max(s0, s1)), [])
        if not members:
            n_edge += 1
            continue
        if len(members) > 1 and not (tips & {s0, s1}):
            n_edge += 1
        ok = False
        for ei in members:
            et = int(types[ei])
            n_corner = 3 if et in (3, 6) else 4
            cen = np.mean([nodes[int(elements[ei, q])][:2]
                           for q in range(n_corner)], axis=0)
            if float((cen - p1) @ normal) * want > 0.0:
                ok = True
        if not ok:
            n_side += 1
    if n_off:
        failures.append(f"{tag}: {n_off} joint node pairs do not stand at one "
                        f"point")
    if n_side:
        failures.append(f"{tag}: {n_side} joint elements take their soil edge "
                        f"from the wrong side of their own line")
    if n_edge:
        failures.append(f"{tag}: {n_edge} joint elements do not stand on one 2D "
                        f"element's edge")
    if n_bar:
        failures.append(f"{tag}: {n_bar} joint node pairs do not stand on the "
                        f"bar of their line")


def _counts(mesh, tag, failures, lines, joint_idx):
    """The joint element count, two per bar element on every jointed line."""
    mats_1d = np.asarray(mesh['element_materials_1d'], dtype=int)
    mats_j = np.asarray(mesh['element_materials_joint'], dtype=int)
    side = np.asarray(mesh['element_side_joint'], dtype=int)
    for li in joint_idx:
        n_bar = int((mats_1d == li + 1).sum())
        n_j = int((mats_j == li + 1).sum())
        if n_j != 2 * n_bar:
            failures.append(f"{tag}: line {li + 1} has {n_j} joint elements for "
                            f"{n_bar} bar elements, not {2 * n_bar}")
        n_up = int(((mats_j == li + 1) & (side == 1)).sum())
        if n_up != n_bar:
            failures.append(f"{tag}: line {li + 1} has {n_up} upper joints for "
                            f"{n_bar} bar elements")


#: The five junction fixtures: the lines, the junction points with the ray
#: directions the joint lines leave them by, the wedge count and how many bars
#: meet there. Every wedge count is countable from the drawing.
FIXTURES = {
    'a. T': dict(
        lines=[[(2.0, 5.0), (18.0, 5.0)], [(10.0, 5.0), (10.0, 9.0)]],
        joints=[0, 1],
        junctions=[((10.0, 5.0), [0.0, math.pi, math.pi / 2], 3, 2)]),
    'b. X': dict(
        lines=[[(2.0, 5.0), (18.0, 5.0)], [(10.0, 1.0), (10.0, 9.0)]],
        joints=[0, 1],
        junctions=[((10.0, 5.0),
                    [0.0, math.pi, math.pi / 2, 3 * math.pi / 2], 4, 2)]),
    'c. L': dict(
        lines=[[(6.0, 5.0), (14.0, 5.0)], [(6.0, 5.0), (6.0, 9.0)]],
        joints=[0, 1],
        junctions=[((6.0, 5.0), [0.0, math.pi / 2], 2, 2)]),
    'd. chain': dict(
        lines=[[(4.0, 5.0), (10.0, 5.0)], [(10.0, 5.0), (16.0, 5.0)]],
        joints=[0, 1],
        junctions=[((10.0, 5.0), [0.0, math.pi], 2, 2)]),
    'e. crack to the boundary': dict(
        lines=[[(4.0, 5.0), (20.0, 5.0)]],
        joints=[0],
        junctions=[((20.0, 5.0), [math.pi, math.pi / 2, 3 * math.pi / 2], 2, 1)]),
}

# ---------------------------------------------------------------------------
# terminations: a joint that ENDS on another, at a point the through line does
# not carry as a vertex of its own
# ---------------------------------------------------------------------------

#: The frame the stepped fixtures are drawn in: a base at 30 degrees and columns
#: standing on it at 120, which is Goodman & Bray's section. Nothing in it is
#: exactly representable, so a point computed ON a line of this frame misses that
#: line by a rounding — which is what the shipped stepped base does.
_TH = math.radians(30.0)
_D = (math.cos(_TH), math.sin(_TH))
_N = (-math.sin(_TH), math.cos(_TH))


def _P(t, n, origin=(0.0, 0.0)):
    """The point ``t`` along the base and ``n`` perpendicular to it."""
    return (origin[0] + t * _D[0] + n * _N[0], origin[1] + t * _D[1] + n * _N[1])


def _stepped(ncol=3, w=5.0, step=0.7, h=3.0, origin=(0.0, 0.0)):
    """The joint lines of a stepped base, the pattern problem 1 is made of.

    Each column stands on its own stretch of base, one ``step`` perpendicular
    above its downslope neighbour's, so the contact between two columns runs from
    the lower base plane up: the piece below the upper column's base is the step,
    the piece above it is the column-to-column contact, and the two are one
    straight plane, which is how the vendor's outlines give it. The upper
    column's base therefore BEGINS partway along that plane — a T at a point the
    plane does not carry as a vertex of its own, once per column.
    """
    lines = [[_P(i * w, i * step, origin), _P((i + 1) * w, i * step, origin)]
             for i in range(ncol)]
    lines += [[_P(i * w, (i - 1) * step, origin), _P(i * w, i * step + h, origin)]
              for i in range(1, ncol)]
    return lines


def _stepped_spec(ncol=3, w=5.0, step=0.7, h=3.0):
    """The stepped base and the junctions to count on it."""
    lines = _stepped(ncol, w, step, h)
    up, down, along = (math.atan2(_N[1], _N[0]),
                       math.atan2(-_N[1], -_N[0]),
                       math.atan2(_D[1], _D[0]))
    junctions = [(_P(i * w, i * step), [up, down, along], 3, 2)
                 for i in range(1, ncol)]                 # the T at each step
    junctions += [(_P(i * w, (i - 1) * step),
                   [up, math.atan2(-_D[1], -_D[0])], 2, 2)
                  for i in range(1, ncol)]                # the L below each step
    return dict(lines=lines, joints=list(range(len(lines))), junctions=junctions,
                ring=[(-6.0, -5.0), (16.0, -5.0), (16.0, 14.0), (-6.0, 14.0)])


#: The termination fixtures. Each is a joint ENDING on another at a point in the
#: interior of the through line's own segment, which is the case a crossing rule
#: alone cannot make: the two lines have to be given the same point before gmsh
#: sees them, or the mesher is handed a sliver and returns an edge carrying four
#: elements or none.
TERMINATIONS = {
    'h. T, tip a rounding off': dict(
        # the through line and the tip are computed in the same 30 degree frame,
        # so the tip is ON the line in intent and a rounding off it in doubles
        lines=[[_P(0.0, 0.0, (3.0, 2.0)), _P(12.0, 0.0, (3.0, 2.0))],
               [_P(6.0, 0.0, (3.0, 2.0)), _P(6.0, -3.0, (3.0, 2.0))]],
        joints=[0, 1],
        junctions=[(_P(6.0, 0.0, (3.0, 2.0)),
                    [math.atan2(_D[1], _D[0]), math.atan2(-_D[1], -_D[0]),
                     math.atan2(-_N[1], -_N[0])], 3, 2)],
        not_incident=True),
    'i. T, tip 1e-7 off': dict(
        lines=[[(2.0, 5.0), (18.0, 5.0)], [(10.0, 5.0 + 1e-7), (10.0, 9.0)]],
        joints=[0, 1],
        junctions=[((10.0, 5.0), [0.0, math.pi, math.pi / 2], 3, 2)],
        snapped=1, same_as='a. T'),
    'j. two tips on one segment': dict(
        lines=[[(2.0, 5.0), (18.0, 5.0)], [(7.0, 5.0), (7.0, 9.0)],
               [(13.0, 5.0), (13.0, 1.0)]],
        joints=[0, 1, 2],
        junctions=[((7.0, 5.0), [0.0, math.pi, math.pi / 2], 3, 2),
                   ((13.0, 5.0), [0.0, math.pi, 3 * math.pi / 2], 3, 2)]),
    # The stepped base: two terminations and two corners in one section, which is
    # the topology of a column stack on a stepped base. Both tips are computed
    # from the same base angle as the plane they stand on and land 3 x 10^-16 off
    # it, so both snap; and the computed intersection of the two segments lands
    # 4 x 10^-16 from the tip, so the fixture also pins that the point inserted
    # into the through line is the TIP's rather than the computed one.
    'k. stepped base, three columns': {**_stepped_spec(), 'snapped': 2},
    'm. L, tip 1e-7 off the end': dict(
        # the nearest point on the through line is its own END, so the tip takes
        # that vertex rather than a point a hair inside the segment
        lines=[[(6.0, 5.0), (14.0, 5.0)],
               [(6.0 + 1e-7, 5.0 + 1e-7), (6.0, 9.0)]],
        joints=[0, 1],
        junctions=[((6.0, 5.0), [0.0, math.pi / 2], 2, 2)],
        snapped=1, same_as='c. L'),
    'l. tip on an existing vertex': dict(
        # the tip lands on the through line's own END: an L corner, which the
        # split has always made, and which the snap must leave alone
        lines=[[(4.0, 5.0), (14.0, 5.0)], [(14.0, 5.0), (14.0, 9.0)]],
        joints=[0, 1],
        junctions=[((14.0, 5.0), [math.pi, math.pi / 2], 2, 2)],
        snapped=0),
}


def _no_collapsed(mesh, tag, failures):
    """No 2D element may name one node twice.

    A collapsed element is what a pair of geometry points a rounding apart turns
    into: the mesher gives each exact coordinate its own point, merges the two as
    duplicate nodes, and what is left is a triangle with two corners at one node.
    It has no area, it carries the same edge twice — which is how an edge on a
    joint comes back with four elements on it rather than two — and nothing
    downstream can tell it from a real element.
    """
    elements = np.asarray(mesh['elements'], dtype=int)
    types = np.asarray(mesh['element_types'], dtype=int)
    bad = 0
    for ei in range(len(elements)):
        n_corner = 3 if int(types[ei]) in (3, 6) else 4
        corners = [int(elements[ei, k]) for k in range(n_corner)]
        if len(set(corners)) != n_corner:
            bad += 1
    if bad:
        failures.append(f"{tag}: {bad} two-dimensional element(s) name one node "
                        f"twice — a collapsed element, so two geometry points a "
                        f"rounding apart reached the mesher")


def _snapped(spec):
    """What :func:`_snap_joint_line_ends` does to a fixture, and how far.

    Returns ``(ends moved, the largest move)``, on a copy — so a fixture can say
    that nothing of its geometry is touched.
    """
    lines = [[tuple(map(float, p[:2])) for p in ln] for ln in spec['lines']]
    before = [list(ln) for ln in lines]
    tol = _joint_line_tol(lines)
    moved = _snap_joint_line_ends(lines, {i: {} for i in spec['joints']}, tol)
    far = 0.0
    for a, b in zip(before, lines):
        for p, q in ((a[0], b[0]), (a[-1], b[-1])):
            far = max(far, math.hypot(p[0] - q[0], p[1] - q[1]))
    return moved, far, tol


def _leg_fixtures(failures, results):
    for tag, spec in FIXTURES.items():
        moved, _far, _tol = _snapped(spec)
        if moved:
            failures.append(f"{tag}: the end snap moved {moved} end(s) of a "
                            f"fixture whose lines already meet exactly")
        for et in ('tri3', 'tri6'):
            mesh = _build(spec['lines'], spec['joints'], et)
            name = f"{tag} / {et}"
            for (p, rays, n_wedges, n_bars) in spec['junctions']:
                _check_junction(mesh, p, rays, n_wedges, n_bars, name, failures)
            _counts(mesh, name, failures, spec['lines'], spec['joints'])
            _check_pairs(mesh, name, failures, spec['lines'])
            _roundtrip(mesh, name, failures)
            nodes = np.asarray(mesh['nodes'])
            here = _at(nodes, spec['junctions'][0][0])
            results.append(
                f"{tag:26s} {et}: {len(mesh['nodes']):5d} nodes, "
                f"{len(mesh['elements_joint']):3d} joint elements, "
                f"{len(here)} nodes at the junction")


def _leg_terminations(failures, results):
    """h to m. A joint that ENDS on another, partway along its segment.

    The point is a vertex of the ending line and of nothing else, so it is one
    point only if the mesher is told to make it one. Each fixture states what the
    end snap does to it before gmsh runs — how many ends it moves and how far —
    and then the mesh is held to the same wedge counts, joint pairing and round
    trip as every other junction.
    """
    for tag, spec in TERMINATIONS.items():
        moved, far, tol = _snapped(spec)
        want = spec.get('snapped')
        if want is not None and moved != want:
            failures.append(f"{tag}: the end snap moved {moved} end(s), not "
                            f"{want}")
        if far > tol:
            failures.append(f"{tag}: the end snap moved an end {far:.3g}, past "
                            f"its own tolerance of {tol:.3g}")
        if spec.get('not_incident'):
            thr = LineString([tuple(spec['lines'][0][0]),
                              tuple(spec['lines'][0][-1])])
            tip = Point(tuple(spec['lines'][1][0]))
            if thr.distance(tip) == 0.0:
                failures.append(f"{tag}: the tip lies exactly on the through "
                                "line, so the fixture no longer stands for a "
                                "termination the arithmetic misses")
        for et in ('tri3', 'tri6'):
            mesh = _build(spec['lines'], spec['joints'], et,
                          ring=spec.get('ring'))
            name = f"{tag} / {et}"
            for (p, rays, n_wedges, n_bars) in spec['junctions']:
                _check_junction(mesh, p, rays, n_wedges, n_bars, name, failures)
            _counts(mesh, name, failures, spec['lines'], spec['joints'])
            _check_pairs(mesh, name, failures, spec['lines'])
            _no_collapsed(mesh, name, failures)
            _roundtrip(mesh, name, failures)
            twin = spec.get('same_as')
            if twin:
                other = _build(FIXTURES[twin]['lines'], FIXTURES[twin]['joints'],
                               et)
                for key in ('nodes', 'elements', 'elements_joint',
                            'element_side_joint'):
                    if not np.array_equal(np.asarray(mesh[key]),
                                          np.asarray(other[key])):
                        failures.append(
                            f"{name}: the near miss does not mesh as '{twin}' "
                            f"does — '{key}' differs, so the snap did not put "
                            "the tip where the exact fixture has it")
                        break
            nodes = np.asarray(mesh['nodes'])
            here = _at(nodes, spec['junctions'][0][0])
            results.append(
                f"{tag:30s} {et}: {len(mesh['nodes']):5d} nodes, "
                f"{len(mesh['elements_joint']):3d} joint elements, "
                f"{len(here)} nodes at the first junction"
                + (f", {moved} end(s) snapped by {far:.3g}" if moved else ""))


#: A zone-boundary vertex stated to twelve decimals. Nothing rounds it, so a
#: crossing computed against it comes back a rounding away.
_FINE_X = 6.123456789012


def _leg_zone_vertex(failures, results):
    """n. A joint ending on a zone edge AT a vertex stated to twelve decimals.

    The crossing of a constraint line with a polygon edge is inserted as a
    polygon vertex, and the point that is inserted is rounded to six decimals. So
    a joint that ends on a zone edge at a vertex the section states more
    precisely than that meets a polygon which already carries the point — and the
    insertion has to recognise it, or the zone gains a second vertex 2 x 10^-7
    from the first and the sliver between them meshes as a collapsed element.
    """
    left = [(0.0, 0.0), (_FINE_X, 0.0), (_FINE_X, 5.0), (_FINE_X, 10.0),
            (0.0, 10.0)]
    right = [(_FINE_X, 0.0), (14.0, 0.0), (14.0, 10.0), (_FINE_X, 10.0),
             (_FINE_X, 5.0)]
    lines = [[(1.0, 5.0), (_FINE_X, 5.0)]]
    # the crossings of the line with the zone edges become polygon vertices
    # first, which is what get_material_polygons does for a real model
    rings = add_intersection_points_to_polygons([left, right], lines)
    for et in ('tri3', 'tri6'):
        mesh = _build(lines, [0], et, mats=[(rings[0], 0), (rings[1], 1)])
        tag = f"n. joint on a zone vertex / {et}"
        _no_collapsed(mesh, tag, failures)
        _counts(mesh, tag, failures, lines, [0])
        _check_pairs(mesh, tag, failures, lines)
        _roundtrip(mesh, tag, failures)
        mats = np.asarray(mesh['element_materials'], dtype=int)
        if len(set(mats.tolist())) != 2:
            failures.append(f"{tag}: the split lost a material zone")
        nodes = np.asarray(mesh['nodes'], dtype=float)
        here = _at(nodes, (_FINE_X, 5.0))
        # the joint's end is a buried tip standing ON the zone edge: one soil
        # copy per wedge, and the two zones meet there, so it is not one node
        results.append(f"n. joint on a zone vertex {et}: "
                       f"{len(mesh['nodes']):5d} nodes, "
                       f"{len(mesh['elements_joint']):3d} joint elements, "
                       f"{len(here)} nodes at the zone vertex")


def _leg_material_boundary(failures, results):
    """f. A joint ON a material boundary, met by a second joint at a T.

    The vertical line is the zone edge between two materials — the mesher reuses
    that edge rather than laying a curve on it (R2b) — and the horizontal joint
    ends on it. The junction is therefore a T whose three wedges are two zones on
    one side and one on the other, which is the wall's block column met by a
    course joint.
    """
    left = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
    right = [(10.0, 0.0), (20.0, 0.0), (20.0, 10.0), (10.0, 10.0)]
    lines = [[(10.0, 1.0), (10.0, 9.0)], [(10.0, 5.0), (18.0, 5.0)]]
    rings = add_intersection_points_to_polygons([left, right], lines)
    for et in ('tri3', 'tri6'):
        mesh = _build(lines, [0, 1], et,
                      mats=[(rings[0], 0), (rings[1], 1)])
        tag = f"f. joint on a material boundary / {et}"
        _check_junction(mesh, (10.0, 5.0),
                        [math.pi / 2, 3 * math.pi / 2, 0.0], 3, 2, tag, failures)
        _counts(mesh, tag, failures, lines, [0, 1])
        _check_pairs(mesh, tag, failures, lines)
        _roundtrip(mesh, tag, failures)
        # the two zones survive, and the wedge on the left of the vertical joint
        # is the left zone's alone
        mats = np.asarray(mesh['element_materials'], dtype=int)
        if len(set(mats.tolist())) != 2:
            failures.append(f"{tag}: the split lost a material zone")
        nodes = np.asarray(mesh['nodes'])
        elements = np.asarray(mesh['elements'], dtype=int)
        types = np.asarray(mesh['element_types'], dtype=int)
        here = set(_at(nodes, (10.0, 5.0)))
        for ei in range(len(elements)):
            et_ = int(types[ei])
            hit = [int(elements[ei, k]) for k in range(et_)
                   if int(elements[ei, k]) in here]
            if not hit:
                continue
            cen = np.mean([nodes[int(elements[ei, k])][:2]
                           for k in range(3 if et_ in (3, 6) else 4)], axis=0)
            want = 1 if cen[0] < 10.0 else 2
            if int(mats[ei]) != want:
                failures.append(f"{tag}: element {ei} at x = {cen[0]:.2f} reads "
                                f"material {int(mats[ei])}, not {want}")
                break
        results.append(f"f. joint on a zone edge   {et}: "
                       f"{len(mesh['nodes']):5d} nodes, "
                       f"{len(mesh['elements_joint']):3d} joint elements, "
                       f"{len(here)} nodes at the T")


def _leg_refusals(failures, results):
    """g. What is still refused, and what is not."""
    ring = [(0.0, 0.0), (20.0, 0.0), (20.0, 10.0), (0.0, 10.0)]
    cases = [
        ("a jointed line meeting a bonded reinforcement line",
         [[(2.0, 5.0), (18.0, 5.0)], [(10.0, 1.0), (10.0, 9.0)]], [0],
         "touches another constraint line"),
        ("a line load on a jointed line",
         [[(2.0, 5.0), (18.0, 5.0)]], [0], "carries a point constraint"),
        ("two jointed lines lying on one another",
         [[(2.0, 5.0), (14.0, 5.0)], [(8.0, 5.0), (18.0, 5.0)]], [0, 1],
         "lie on one another"),
    ]
    for name, lines, joints, want in cases:
        kw = {}
        if 'line load' in name:
            kw['point_constraints'] = [(10.0, 5.0)]
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                build_mesh_from_polygons(
                    [{'coords': list(ring), 'mat_id': 0}], 2.0, 'tri3',
                    lines=[[tuple(p) for p in ln] for ln in lines],
                    element_size_1d=1.0, joint_lines=joints, **kw)
        except ValueError as e:
            if want not in str(e):
                failures.append(f"g: {name} raised the wrong message: {e}")
            else:
                results.append(f"g. refused: {name}")
            continue
        failures.append(f"g: the mesher accepted {name}; it must refuse it")

    # a jointed line ON the domain boundary has material on one side only
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            build_mesh_from_polygons(
                [{'coords': list(ring), 'mat_id': 0}], 2.0, 'tri3',
                lines=[[(2.0, 0.0), (18.0, 0.0)]], element_size_1d=1.0,
                joint_lines=[0])
    except (ValueError, Exception) as e:          # noqa: BLE001
        results.append("g. refused: a jointed line on the domain boundary")
        if 'element' not in str(e) and 'boundary' not in str(e):
            results[-1] += f" ({type(e).__name__})"
    else:
        failures.append("g: a jointed line lying on the domain boundary was "
                        "accepted; it has material on one side only")


# ---------------------------------------------------------------------------
# the closed-form rows
# ---------------------------------------------------------------------------

def _base():
    global _base_sd
    if _base_sd is None:
        with contextlib.redirect_stdout(io.StringIO()):
            _base_sd = load_slope_data(BASE_FILE)
    return copy.deepcopy(_base_sd)


def _material(name, **kw):
    m = dict(name=name, gamma=GAMMA, gamma_sat=None, option='mc', c=100.0,
             phi=45.0, E=E_SOIL, nu=NU_SOIL, t_cut=None, u='none', ru=0.0)
    m.update(kw)
    return m


def _reinf_row(line, cj, phi_j, label):
    return dict(label=label, x1=line[0][0], y1=line[0][1],
                x2=line[-1][0], y2=line[-1][1],
                t_max=1.0e6, t_res=float('nan'), lp1=0.0, lp2=0.0,
                tend1=0.0, tend2=0.0, E=2.0e5, area=1.0e-5, spacing=None,
                adhesion=cj, delta=phi_j, kn=None, ks=None, jred='yes')


def _model(polys_and_ids, domain, ground, y_bottom, materials, joint_lines,
           cj, phi_j, ts, s1d, k_seismic=0.0, element_type='tri6'):
    """Install the geometry and every jointed line, mesh, build fem_data."""
    d = _base()
    d['unit_system'] = 'metric'
    d['gamma_water'] = 9.81
    d['profile_lines'] = []
    d['polygons'] = [{'polygon': Polygon(r), 'mat_id': i}
                     for r, i in polys_and_ids]
    d['domain_polygon'] = domain
    d['ground_surface'] = ground
    d['max_depth'] = y_bottom
    d['circles'] = []
    d['non_circ'] = []
    d['piezo_line'] = []
    d['piezo_phreatic'] = False
    d['materials'] = materials
    d['reinforcement_lines'] = [
        _reinf_row(ln, cj, phi_j, f'joint {i + 1}')
        for i, ln in enumerate(joint_lines)]
    d['reinforce_lines'] = build_reinforce_lines(d['reinforcement_lines'])
    d['pile_lines'] = []
    lines = [[(float(ln[0][0]), float(ln[0][1])),
              (float(ln[-1][0]), float(ln[-1][1]))] for ln in joint_lines]
    polys = [{'coords': list(r), 'mat_id': i} for r, i in polys_and_ids]
    with contextlib.redirect_stdout(io.StringIO()):
        mesh = build_mesh_from_polygons(polys, target_size=ts,
                                        element_type=element_type, lines=lines,
                                        element_size_1d=s1d,
                                        joint_lines=list(range(len(lines))))
        fem_data = build_fem_data(d, mesh)
    fem_data['k_seismic'] = k_seismic
    return d, mesh, fem_data


def _solve(fem_data, F=1.0, max_iterations=8000):
    with contextlib.redirect_stdout(io.StringIO()):
        return solve_fem(fem_data, F=F, max_iterations=max_iterations,
                         fast_kernel=False)


def _ssrm(fem_data, F_min=1.0, F_max=2.5, tolerance=0.005, **kw):
    kw.setdefault('max_iterations', 3000)
    kw.setdefault('max_iterations_ceiling', 6000)
    with contextlib.redirect_stdout(io.StringIO()):
        return solve_ssrm(fem_data, F_min=F_min, F_max=F_max,
                          tolerance=tolerance, **kw)


#: Rows i and ii: R2's slab on a plane at beta, with the cuts phase 2 adds.
SLAB = dict(beta=20.0, HV=3.0, X1=36.0, XA=6.0, XB=30.0, YB=-6.0, VD=1.0,
            ts=1.5, s1d=1.0, phi_j=30.0, E_void=1.0)


def _slab(cut=False, course=False):
    """The slab on its plane, optionally cut by a vertical joint and stacked.

    ``cut`` adds a vertical joint at mid-length from the base joint to the ground
    surface: a T where it meets the plane and a crack to the surface at the top.
    ``course`` adds a joint parallel to the plane at mid-thickness, which the cut
    crosses at an X.
    """
    p = SLAB
    T = math.tan(math.radians(p['beta']))
    y = lambda x: x * T
    HV, X1, XA, XB, YB, VD = (p['HV'], p['X1'], p['XA'], p['XB'], p['YB'],
                              p['VD'])
    soil = [(0.0, YB), (X1, YB), (X1, y(X1) - VD), (XB, y(XB) - VD),
            (XB, y(XB) + HV), (XA, y(XA) + HV), (XA, y(XA) - VD), (0.0, -VD)]
    vl = [(0.0, -VD), (XA, y(XA) - VD), (XA, y(XA) + HV), (0.0, HV)]
    vr = [(XB, y(XB) - VD), (X1, y(X1) - VD), (X1, y(X1) + HV), (XB, y(XB) + HV)]
    lines = [[(0.0, 0.0), (X1, y(X1))]]
    if course:
        lines.append([(0.0, 0.5 * HV), (X1, y(X1) + 0.5 * HV)])
    if cut:
        xc = 0.5 * (XA + XB)
        lines.append([(xc, y(xc)), (xc, y(xc) + HV)])
    rings = add_intersection_points_to_polygons([soil, vl, vr], lines)
    mats = [_material('soil'),
            _material('void', option='elastic', c=0.0, phi=0.0,
                      E=p['E_void'], gamma=0.001, nu=0.2)]
    ids = [(rings[0], 0), (rings[1], 1), (rings[2], 1)]
    domain = Polygon([(0.0, HV), (X1, y(X1) + HV), (X1, YB), (0.0, YB)])
    ground = LineString([(0.0, HV), (X1, y(X1) + HV)])
    d, mesh, fem_data = _model(ids, domain, ground, YB, mats, lines,
                               0.0, p['phi_j'], p['ts'], p['s1d'])
    return d, mesh, fem_data, lines


def _leg_row_i(failures, results):
    """i. The block on a plane, cut in two by a vertical joint."""
    p = SLAB
    expected = math.tan(math.radians(p['phi_j'])) / math.tan(
        math.radians(p['beta']))
    _d, mesh, fem_data, lines = _slab(cut=True)
    xc = 0.5 * (p['XA'] + p['XB'])
    yc = xc * math.tan(math.radians(p['beta']))
    _check_junction(mesh, (xc, yc),
                    [math.radians(p['beta']), math.radians(p['beta']) + math.pi,
                     math.pi / 2], 3, 2, 'row i (T on the plane)', failures)
    _check_junction(mesh, (xc, yc + p['HV']),
                    [3 * math.pi / 2, math.radians(p['beta']),
                     math.radians(p['beta']) + math.pi], 2, 1,
                    'row i (crack to the surface)', failures)
    _check_pairs(mesh, 'row i', failures, lines)
    res = _ssrm(fem_data, F_min=1.0, F_max=2.5, tolerance=0.005)
    FS = res.get('FS')
    if FS is None:
        failures.append("row i: the strength reduction returned no factor")
        return
    err = (FS - expected) / expected
    if abs(err) > 0.03:
        failures.append(f"row i: the cut slab returns FS = {FS:.4f} against "
                        f"tan phi_j / tan beta = {expected:.4f} "
                        f"({100 * err:+.2f}%)")
    results.append(f"row i   cut into two blocks: SSRM FS = {FS:.4f} vs "
                   f"{expected:.4f}  ({100 * err:+.2f}%), "
                   f"{len(mesh['elements_joint'])} joint elements")


def _leg_row_ii(failures, results):
    """ii. The same slab as a two-course stack, the cut crossing the course."""
    p = SLAB
    expected = math.tan(math.radians(p['phi_j'])) / math.tan(
        math.radians(p['beta']))
    _d, mesh, fem_data, lines = _slab(cut=True, course=True)
    xc = 0.5 * (p['XA'] + p['XB'])
    yc = xc * math.tan(math.radians(p['beta']))
    _check_junction(mesh, (xc, yc + 0.5 * p['HV']),
                    [math.radians(p['beta']), math.radians(p['beta']) + math.pi,
                     math.pi / 2, 3 * math.pi / 2], 4, 2,
                    'row ii (X on the course joint)', failures)
    _check_pairs(mesh, 'row ii', failures, lines)
    res = _ssrm(fem_data, F_min=1.0, F_max=2.5, tolerance=0.005)
    FS = res.get('FS')
    if FS is None:
        failures.append("row ii: the strength reduction returned no factor")
        return
    err = (FS - expected) / expected
    if abs(err) > 0.03:
        failures.append(f"row ii: the stacked slab returns FS = {FS:.4f} "
                        f"against tan phi_j / tan beta = {expected:.4f} "
                        f"({100 * err:+.2f}%)")
    results.append(f"row ii  two courses, four blocks: SSRM FS = {FS:.4f} vs "
                   f"{expected:.4f}  ({100 * err:+.2f}%), "
                   f"{len(mesh['elements_joint'])} joint elements")


#: Row iii: the L-shaped joint around an elastic block.
BLOCK = dict(W=30.0, H=12.0, XA=12.0, XB=18.0, YB=6.0, phi_base=30.0,
             k=0.30, ts=1.0, s1d=0.75, E_void=1.0)


def _block_model(phi_base=None, k=None):
    """An elastic block in a notch: base joint under it, back joint up its back.

    The block occupies XA..XB, YB..H. The soil fills 0..H below YB; the space
    behind the block (XA..0 above YB) and in front of it (XB..W above YB) is the
    near-weightless elastic 'void' the R2 rows use for a free face, so the only
    horizontal load on the block is the seismic one, k W.

    The base joint runs from the corner to the model's right-hand boundary rather
    than stopping under the block's front face. A joint that stops there leaves a
    buried crack tip at the block's own bottom corner — one shared node, which is
    a rigid link between the block and the foundation — and the block cannot
    slide out past it: at that geometry the reduction returns 6.57 rather than
    1.92, which is not the interface's answer but the pin's.
    """
    p = dict(BLOCK)
    if phi_base is not None:
        p['phi_base'] = phi_base
    if k is not None:
        p['k'] = k
    W, H, XA, XB, YB = p['W'], p['H'], p['XA'], p['XB'], p['YB']
    soil = [(0.0, 0.0), (W, 0.0), (W, YB), (0.0, YB)]
    block = [(XA, YB), (XB, YB), (XB, H), (XA, H)]
    behind = [(0.0, YB), (XA, YB), (XA, H), (0.0, H)]
    front = [(XB, YB), (W, YB), (W, H), (XB, H)]
    base = [(XA, YB), (W, YB)]
    back = [(XA, YB), (XA, H)]
    rings = add_intersection_points_to_polygons([soil, block, behind, front],
                                                [base, back])
    mats = [_material('soil', c=1000.0, phi=45.0),
            _material('block', option='elastic', c=0.0, phi=0.0, E=E_SOIL,
                      gamma=GAMMA, nu=0.2),
            _material('void', option='elastic', c=0.0, phi=0.0, E=p['E_void'],
                      gamma=0.001, nu=0.2)]
    ids = [(rings[0], 0), (rings[1], 1), (rings[2], 2), (rings[3], 2)]
    domain = Polygon([(0.0, 0.0), (W, 0.0), (W, H), (0.0, H)])
    ground = LineString([(0.0, H), (W, H)])
    d, mesh, fem_data = _model(ids, domain, ground, 0.0, mats, [base, back],
                               0.0, p['phi_base'], p['ts'], p['s1d'],
                               k_seismic=p['k'])
    return d, mesh, fem_data, [base, back], p


def _leg_row_iii(failures, results):
    """iii. The block slides out of its corner at the base joint's friction."""
    _d, mesh, fem_data, lines, p = _block_model()
    expected = math.tan(math.radians(p['phi_base'])) / p['k']
    _check_junction(mesh, (p['XA'], p['YB']),
                    [0.0, math.pi / 2], 2, 2, 'row iii (the L corner)',
                    failures)
    _check_pairs(mesh, 'row iii', failures, lines)

    # the corner is what lets the block leave: with the two lines a single
    # bonded pair of nodes there would be no mechanism at all.
    res = _ssrm(fem_data, F_min=0.5, F_max=3.0, tolerance=0.005)
    FS = res.get('FS')
    if FS is None:
        failures.append("row iii: the strength reduction returned no factor")
        return
    err = (FS - expected) / expected
    if abs(err) > 0.05:
        failures.append(f"row iii: the block returns FS = {FS:.4f} against "
                        f"tan phi_base / k = {expected:.4f} "
                        f"({100 * err:+.2f}%)")

    # the back face opens: its normal traction is zero where the block has left
    sol = _solve(fem_data, F=max(1.0, expected * 0.9))
    jd = fem_data['joint_data']
    back = np.asarray(jd['line_id']) == 2
    tn = np.asarray(sol['joint_tn'])
    w = np.asarray(jd['w'])
    thrust = float(np.sum(np.maximum(tn[back], 0.0) * w[back]))
    weight = GAMMA * (p['XB'] - p['XA']) * (p['H'] - p['YB'])
    if thrust > 0.02 * weight:
        failures.append(f"row iii: the back face carries {thrust:.2f} kN/m of "
                        f"thrust against a block weight of {weight:.1f}; it "
                        f"must open as the block slides out")
    results.append(f"row iii the block in an L: SSRM FS = {FS:.4f} vs "
                   f"tan phi_base / k = {expected:.4f}  ({100 * err:+.2f}%), "
                   f"back-face thrust {thrust:.3f} kN/m on a "
                   f"{weight:.0f} kN/m block")


def run_mesh_legs():
    """The mesh fixtures and the refusals, without the rows that solve.

    What the --mesh scope runs: every leg here is a mesh build and a count.
    """
    failures, results = [], []
    t0 = time.time()
    _leg_fixtures(failures, results)
    _leg_terminations(failures, results)
    _leg_zone_vertex(failures, results)
    _leg_material_boundary(failures, results)
    _leg_refusals(failures, results)
    print(f"Joint junction check, mesh legs ({time.time() - t0:.0f} s):")
    for line in results:
        print("  " + line)
    return failures


def run():
    """Returns a list of failure strings (empty = pass)."""
    failures, results = [], []
    t0 = time.time()
    _leg_fixtures(failures, results)
    _leg_terminations(failures, results)
    _leg_zone_vertex(failures, results)
    _leg_material_boundary(failures, results)
    _leg_refusals(failures, results)
    _leg_row_i(failures, results)
    _leg_row_ii(failures, results)
    _leg_row_iii(failures, results)
    print(f"Joint junction check ({time.time() - t0:.0f} s):")
    for line in results:
        print("  " + line)
    return failures


def main():
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nJointed lines meet: the split copies a junction node once per wedge "
          "of material, and the closed forms hold through the junctions.")


if __name__ == '__main__':
    main()
