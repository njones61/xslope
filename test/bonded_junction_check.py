# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Bonded reinforcement and pile lines that meet share a mesh node.

A bonded bar and a pile both stand on the soil's nodes. Where two of them meet —
a tieback ending partway down a pile, a bar crossing a pile, two bars crossing, a
bar ending on another bar — they must stand on ONE node at the meeting point, and
every element of both must lie along an edge of the soil elements. gmsh gives a
shared node only where the geometry carries a shared point, so the mesher inserts
each meeting point as a vertex of both lines before gmsh sees them
(``mesh._insert_bonded_junction_points``). Without that, gmsh can honor only one
of the two curves at the meeting point: the other comes back with an element that
cuts across the soil elements, or skips the soil node at the bar's end, and the
two members are joined nowhere.

What is checked, on the FEM-4 sheet-pile wall section (pile at x = 10 from
y = 10 to -10), on tri6 and quad8 meshes:

  A. FOUR MEETINGS — a bar ending partway down the pile, a bar crossing the
     pile, two bars crossing, and a T (a bar ending on another bar's interior),
     plus the partway end stated a rounding (1e-9) off the pile. In each: a node
     stands at the meeting point and belongs to elements of BOTH lines; no 1D
     element lies off a soil edge; every 1D element carries its midside node;
     and each line's elements cover its whole length.
  B. LINES THAT DO NOT MEET ARE NOT TOUCHED — the insertion leaves a pair of
     lines clear of each other exactly as given, and returns no junction.
  C. THE REFUSAL — a mesh with a 1D element that is not an edge of the soil
     elements raises MeshInputError naming the line, and a sound mesh passes.
  D. AN END A ROUNDING ERROR OFF A ZONE CORNER — a bar whose end is stated
     1e-15 off a zone corner, and vp032a's geotextile (end x =
     -1.1069999999999993 against the face corner -1.107): the end is snapped
     onto the corner, so exactly one node stands there and it is on the bar,
     and the mesh passes the soil-edge check. Without the snap gmsh is handed
     two points 1e-15 apart and either cannot recover the edge between them
     (the last bar element then cuts across the soil and the mesh is refused)
     or leaves a second node beside the corner.

Run directly:  PYTHONPATH=. python3 test/bonded_junction_check.py
"""

import contextlib
import copy
import io
import os
import sys
import warnings

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

import numpy as np                                              # noqa: E402

#: The FEM-4 wall section. Its one pile is the sheet-pile wall at x = 10.
WALL = os.path.join(_REPO, "docs", "tutorials", "files", "xslope_pile_wall.xlsx")
PILE = [(10.0, 10.0), (10.0, -10.0)]
SIZE = 2.0
ELEMENT_TYPES = ("tri6", "quad8")

#: name -> (lines, the meeting point, the two line indices that meet there)
CASES = {
    "bar ending partway down the pile": (
        [[(10.0, 5.3), (30.0, 0.3)], PILE], (10.0, 5.3), (0, 1)),
    "bar ending a rounding off the pile": (
        [[(10.0 + 1e-9, 5.3), (30.0, 0.3)], PILE], (10.0, 5.3), (0, 1)),
    "bar crossing the pile": (
        [[(6.0, 4.0), (30.0, -2.0)], PILE], (10.0, 3.0), (0, 1)),
    "two bars crossing": (
        [[(12.0, 8.0), (40.0, 2.0)], [(12.0, 2.0), (40.0, 8.0)], PILE],
        (26.0, 5.0), (0, 1)),
    "a bar ending on another bar (T)": (
        [[(12.0, 8.0), (40.0, 2.0)], [(26.0, 5.0), (40.0, 12.0)], PILE],
        (26.0, 5.0), (0, 1)),
}

#: Leg D. A bar from inside the backfill to the toe corner (0, 0), its end
#: stated 1e-15 off the corner, and the LAST point of the line: the order in
#: which the line's own end, not the corner, survives when the two are merged
#: into the line's point list.
SNAP_LINE = [(20.0, -5.0), (1e-15, 0.0)]
#: Leg D. Where it was found: vp032a's geotextile end, recomputed from the LEM
#: point list, is x = -1.1069999999999993 against the face corner -1.107.
VP032A = os.path.join(_REPO, "docs", "verification", "files", "rocscience",
                      "vp032a.xlsx")

_CORNERS = {3: 3, 6: 3, 4: 4, 8: 4, 9: 4}


def _wall():
    from xslope.fileio import load_slope_data
    with contextlib.redirect_stdout(io.StringIO()):
        return load_slope_data(WALL)


def _mesh(sd, lines, element_type):
    from xslope.mesh import build_mesh_from_polygons, get_material_polygons
    lines = [list(ln) for ln in lines]
    with contextlib.redirect_stdout(io.StringIO()):
        polys = get_material_polygons(sd, reinf_lines=lines)
        return build_mesh_from_polygons(polys, target_size=SIZE,
                                        element_type=element_type, lines=lines)


def _soil_edges(mesh):
    edges = set()
    for e, t in zip(mesh["elements"], mesh["element_types"]):
        c = [int(v) for v in e[:_CORNERS[int(t)]]]
        for k in range(len(c)):
            a, b = c[k], c[(k + 1) % len(c)]
            edges.add((min(a, b), max(a, b)))
    return edges


def _leg_meetings(failures):
    sd = _wall()
    rows = []
    for et in ELEMENT_TYPES:
        for name, (lines, meet, (la, lb)) in CASES.items():
            tag = f"{name} ({et})"
            try:
                mesh = _mesh(sd, lines, et)
            except Exception as exc:
                failures.append(f"{tag}: the mesher raised {type(exc).__name__}: "
                                f"{str(exc)[:160]}")
                continue
            nodes = np.asarray(mesh["nodes"], dtype=float)
            e1 = np.asarray(mesh["elements_1d"], dtype=int)
            t1 = np.asarray(mesh["element_types_1d"], dtype=int)
            m1 = np.asarray(mesh["element_materials_1d"], dtype=int)
            edges = _soil_edges(mesh)

            def line_nodes(li):
                out = set()
                for e, t in zip(e1[m1 == li + 1], t1[m1 == li + 1]):
                    out.update(int(v) for v in e[:t])
                return out

            # 1. one node at the meeting point, on both lines
            d = np.hypot(nodes[:, 0] - meet[0], nodes[:, 1] - meet[1])
            j = int(np.argmin(d))
            if d[j] > 1e-6:
                failures.append(f"{tag}: no node at the meeting point {meet} "
                                f"(nearest {d[j]:.3g} away)")
            elif j not in line_nodes(la) or j not in line_nodes(lb):
                failures.append(f"{tag}: the node at {meet} is not shared by "
                                f"lines {la + 1} and {lb + 1}")
            # 2. no 1D element off a soil edge
            off = [(int(e[0]), int(e[1])) for e in e1
                   if (min(int(e[0]), int(e[1])), max(int(e[0]), int(e[1])))
                   not in edges]
            if off:
                a, b = off[0]
                failures.append(f"{tag}: {len(off)} 1D element(s) lie off every "
                                f"soil edge, e.g. {tuple(np.round(nodes[a], 4))} -> "
                                f"{tuple(np.round(nodes[b], 4))}")
            # 3. every element carries its midside node
            if int(t1.min()) != 3:
                failures.append(f"{tag}: {int((t1 != 3).sum())} 1D element(s) "
                                f"are two-node on a quadratic mesh")
            # 4. each line's elements cover its whole length
            for li, ln in enumerate(lines):
                want = float(np.hypot(ln[-1][0] - ln[0][0], ln[-1][1] - ln[0][1]))
                sel = e1[m1 == li + 1]
                got = float(sum(np.hypot(*(nodes[e[1]] - nodes[e[0]])) for e in sel))
                if abs(got - want) > 1e-6 * max(1.0, want):
                    failures.append(f"{tag}: line {li + 1}'s elements cover "
                                    f"{got:.6g} of its {want:.6g}")
            rows.append(f"{tag:46s} node {j} at {meet} on lines {la + 1} and "
                        f"{lb + 1}; {len(e1)} 1D elements, {len(off)} off a soil "
                        f"edge")
    return rows


def _leg_untouched(failures):
    from xslope.mesh import _insert_bonded_junction_points
    lines = [[(6.0, 4.0), (9.0, 2.0)], [tuple(p) for p in PILE]]
    before = copy.deepcopy(lines)
    ids = [id(ln) for ln in lines]
    got = _insert_bonded_junction_points(lines, None, None)
    if got:
        failures.append(f"lines that do not meet returned junctions {got}")
    if lines != before or [id(ln) for ln in lines] != ids:
        failures.append("lines that do not meet were changed by the junction "
                        "insertion")
    # ... and a jointed pair is not the bonded rule's to touch
    lines = [[(6.0, 4.0), (30.0, -2.0)], [tuple(p) for p in PILE]]
    before = copy.deepcopy(lines)
    _insert_bonded_junction_points(lines, {0: {}}, None)
    if lines != before:
        failures.append("a jointed line was changed by the bonded junction "
                        "insertion")


def _leg_refusal(failures):
    from xslope.mesh import MeshInputError, _check_1d_elements_on_soil_edges
    # Two triangles on a unit square, split along the diagonal 0-2.
    nodes = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    mesh = {"nodes": nodes,
            "elements": np.array([[0, 1, 2], [0, 2, 3]]),
            "element_types": np.array([3, 3]),
            "elements_1d": np.array([[0, 2, 0]]),
            "element_types_1d": np.array([2]),
            "element_materials_1d": np.array([4])}
    try:
        _check_1d_elements_on_soil_edges(mesh)
    except MeshInputError as exc:
        failures.append(f"a 1D element on the soil diagonal was refused: {exc}")
    bad = dict(mesh, elements_1d=np.array([[1, 3, 0]]))      # the other diagonal
    try:
        _check_1d_elements_on_soil_edges(bad)
    except MeshInputError as exc:
        if "line 4" not in str(exc):
            failures.append(f"the refusal does not name the line: {exc}")
    else:
        failures.append("a 1D element across the soil elements was not refused")


def _leg_snap(failures):
    from xslope.fileio import load_slope_data
    from xslope.mesh import extract_constraint_line_geometry, get_material_polygons
    with contextlib.redirect_stdout(io.StringIO()):
        vp = load_slope_data(VP032A)
    cases = (("a bar end 1e-15 off the toe corner", _wall(), [list(SNAP_LINE)]),
             ("vp032a's geotextile end", vp,
              extract_constraint_line_geometry(vp)[0]))
    rows = []
    for name, sd, lines in cases:
        # The premise: the end is off a zone corner by more than nothing and by
        # no more than the 1e-12 under which the snap used to leave it alone.
        with contextlib.redirect_stdout(io.StringIO()):
            polys = get_material_polygons(sd, reinf_lines=[list(ln) for ln in lines])
        corners = np.array([p for poly in polys for p in poly["coords"]], dtype=float)
        end = np.asarray(lines[0][-1], dtype=float)
        dc = np.hypot(corners[:, 0] - end[0], corners[:, 1] - end[1])
        corner = corners[int(np.argmin(dc))]
        if not 0.0 < dc.min() <= 1e-12:
            failures.append(f"{name}: the end is {dc.min():.3g} off the nearest "
                            f"zone corner, not a rounding error off it")
            continue
        for et in ELEMENT_TYPES:
            tag = f"{name} ({et})"
            try:
                mesh = _mesh(sd, lines, et)
            except Exception as exc:
                failures.append(f"{tag}: the mesher raised {type(exc).__name__}: "
                                f"{str(exc)[:160]}")
                continue
            nodes = np.asarray(mesh["nodes"], dtype=float)
            e1 = np.asarray(mesh["elements_1d"], dtype=int)
            t1 = np.asarray(mesh["element_types_1d"], dtype=int)
            m1 = np.asarray(mesh["element_materials_1d"], dtype=int)
            bar = set()
            for e, t in zip(e1[m1 == 1], t1[m1 == 1]):
                bar.update(int(v) for v in e[:t])
            # 1. the end was snapped onto the corner: one node there, on the bar
            at = f"({corner[0]:g}, {corner[1]:g})"
            d = np.hypot(nodes[:, 0] - corner[0], nodes[:, 1] - corner[1])
            near = np.flatnonzero(d <= 1e-9)
            on_bar = len(near) == 1 and int(near[0]) in bar
            if len(near) != 1 or d[near[0]] != 0.0:
                failures.append(f"{tag}: {len(near)} nodes within 1e-9 of the "
                                f"corner {at}; the end was not snapped onto it")
            elif not on_bar:
                failures.append(f"{tag}: the node at the corner {at} is not on "
                                f"the bar")
            # 2. no 1D element off a soil edge (the mesher's own check passed
            # if we got here; this states it independently)
            edges = _soil_edges(mesh)
            off = [e for e in e1 if (min(int(e[0]), int(e[1])),
                                     max(int(e[0]), int(e[1]))) not in edges]
            if off:
                failures.append(f"{tag}: {len(off)} 1D element(s) lie off every "
                                f"soil edge")
            rows.append(f"{tag:46s} end {dc.min():.1e} off {at}: "
                        f"{len(near)} node(s) within 1e-9 of it"
                        f"{', on the bar' if on_bar else ''}; {len(e1)} 1D "
                        f"elements, {len(off)} off a soil edge")
    return rows


LEGS = (
    ("A. four meetings share a node", _leg_meetings),
    ("B. lines that do not meet are not touched", _leg_untouched),
    ("C. a 1D element off the soil edges is refused", _leg_refusal),
    ("D. an end a rounding error off a corner snaps", _leg_snap),
)


def run():
    """Failures as a list, for run_tests.py."""
    failures = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for name, leg in LEGS:
            before = len(failures)
            try:
                rows = leg(failures)
            except Exception as exc:             # a leg that raises is a failure
                failures.append(f"{name}: {type(exc).__name__}: {exc}")
                rows = None
            print(f"  {name:48s} " + ("ok" if len(failures) == before else "FAILED"))
            for r in rows or ():
                print(f"    {r}")
    return failures


def _cli():
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nBonded lines that meet share a node, every 1D element lies along a "
          "soil edge, a mesh where one does not is refused, and a line end a "
          "rounding error off a zone corner is snapped onto it.")


if __name__ == "__main__":
    _cli()
