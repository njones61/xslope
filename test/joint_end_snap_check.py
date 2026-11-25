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

"""Tutorial FEM-5 typed to three decimals meshes, to the completed file's mesh.

Tutorial FEM-5's toppling slope, typed to three decimals: the joint region and the
base plane run from (9, 6.887) to (20.928, 0), and the column set is built over
that region. The toe corner of the section is (20.9282032303, 0), 0.2 mm from the
base plane's end, and the region's face edge runs 0.2 mm inside the slope face
there, so the first column, clipped to the region, stops 0.07 mm under the face.

Two things went wrong on that model, and both are pinned here.

The mesher builds the base plane in two steps. The columns' feet are placed on it
as vertices first, where they stand; then its end is moved onto the toe corner,
since a line end that close to a zone corner would otherwise leave gmsh a sliver
it cannot mesh. The end moves and the feet do not, so the curve gmsh meshes bends
at the feet, by up to 0.1 mm off its new chord. The joint split measured every
node against that chord at a millionth of the section and stopped there, with an
unhandled ValueError the Studio Log printed as a traceback
(``Node 9 of jointed line 6 lies -8.54e-05 off the line``).

And once that meshes, the first column is a joint that stops 0.07 mm inside the
rock: a buried tip, which bonds the column to the face and gives FS 1.152 where
the completed file gives 1.105. A jointed line end that close to the outer
boundary -- within a twentieth of an element, the distance within which the
mesher already moves a line end onto a zone corner -- is now moved onto the
boundary, and the Log says so.

  A. **the reader's model is the completed file's mesh.** The three-decimal
     model, its columns built by the network generator Studio's Build
     network... dialog runs, meshes at tri6 / 1.5 m to the completed file's 947
     nodes, 402 elements and 41 joint elements on 6 lines; the base plane's
     meshed end is the toe corner (so the bent curve is what was split); the
     first column's meshed end lies on the slope face, and the mesher printed
     that it moved it;
  B. **the completed file still meshes** to the counts the tutorial prints;
  C. **a split refusal is a message.** Every refusal of the joint split, and of
     the joint checks before it, is a ``MeshInputError`` -- the error Studio's
     Build Mesh reports as a sentence -- not a bare ValueError it prints as a
     traceback;
  D. **which ends move.** On a square section meshed at 1.0: a jointed line
     ending 0.03 inside the outer boundary is moved onto it; one ending 0.1
     inside is not (beyond a twentieth of an element); one ending 0.03 inside
     near a corner is left for the corner snap; one ending 0.03 short of an
     interior material boundary is not moved (that is a joint stopping in rock).

Run directly:  PYTHONPATH=. python3 test/joint_end_snap_check.py
"""

import contextlib
import io
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FILES = ROOT / 'docs' / 'tutorials' / 'files'
START = FILES / 'xslope_rock_joints_start.xlsx'
DONE = FILES / 'xslope_rock_toppling.xlsx'

#: Tutorial FEM-5 Part 2 at three decimals: the region and base plane as a reader
#: types the page's 9.000 / 6.887 / 20.928 / 0.000, and the column set's
#: Build network... settings (Parallel set, dip 110, spacing 1.5, offset 1.25,
#: c 0, phi 35).
REGION_3DP = [(9.0, 6.887), (20.928, 0.0), (14.0, 12.0), (9.0, 12.0)]
BASE_3DP = (9.0, 6.887, 20.928, 0.0)
TOE = (20.9282032303, 0.0)
#: What the tutorial prints for the completed file's mesh.
DONE_COUNTS = (947, 402, 41)
#: What the mesher printed on the last build (the lines Studio's Log shows).
_LOG = io.StringIO()


def _mesh(sd):
    """Build the model's mesh the way Studio's Build Mesh does at the tutorial's
    settings (tri6, target 1.5, no auto size). Returns ``(mesh, lines)``."""
    from xslope.mesh import (build_mesh_from_polygons,
                             extract_constraint_line_geometry,
                             extract_joint_options, extract_point_constraints,
                             extract_size_regions, get_material_polygons)
    lines, _nr, _np = extract_constraint_line_geometry(sd)
    polys = get_material_polygons(sd, reinf_lines=lines)
    with contextlib.redirect_stdout(_LOG):
        mesh = build_mesh_from_polygons(
            polys, target_size=1.5, element_type='tri6', lines=lines or None,
            element_size_1d=sd.get('element_size_1d'),
            size_regions=extract_size_regions(sd),
            joint_lines=extract_joint_options(sd),
            point_constraints=extract_point_constraints(sd) or None)
    return mesh, lines


def _reader_model():
    """The starter with the three-decimal region, the column set built over it,
    and the three-decimal base plane, as the page's steps leave the model."""
    from xslope.fileio import load_slope_data
    from xslope.joints import parallel_set

    done = load_slope_data(str(DONE))
    sd = load_slope_data(str(START))
    zone = dict(done['joint_zones'][0], polygon=list(REGION_3DP))
    sd['joint_zones'] = [zone]
    base = [r for r in done['joint_lines'] if r['label'] == 'base plane'][0]
    props = {k: base.get(k) for k in ('c', 'phi')}
    cols = parallel_set(sd, 110.0, 1.5, offset=1.25, region=list(REGION_3DP),
                        label='col', props=props)
    x1, y1, x2, y2 = BASE_3DP
    sd['joint_lines'] = list(cols) + [dict(base, x1=x1, y1=y1, x2=x2, y2=y2)]
    return sd


def check_the_reader_model_meshes(fails):
    """A: the three-decimal model meshes to the completed file's mesh."""
    import numpy as np
    try:
        sd = _reader_model()
    except Exception as exc:                                 # pragma: no cover
        fails.append(f"A: could not build the reader's model: {exc!r}")
        return
    n_cols = sum(1 for r in sd['joint_lines'] if str(r['label']).startswith('col'))
    if n_cols != 5:
        fails.append(f"A: the column set over the three-decimal region has "
                     f"{n_cols} rows, not the page's five")
    _LOG.seek(0)
    _LOG.truncate()
    try:
        mesh, lines = _mesh(sd)
    except Exception as exc:
        fails.append(f"A: the three-decimal model does not mesh: "
                     f"{type(exc).__name__}: {exc}")
        return
    got = (len(mesh['nodes']), len(mesh['elements']),
           len(mesh.get('elements_joint', [])))
    if got != DONE_COUNTS:
        fails.append(f"A: the three-decimal model meshes to {got} (nodes, "
                     f"elements, joint elements), not the completed file's "
                     f"{DONE_COUNTS}")
    joints = mesh.get('joints') or []
    if len(joints) != 6:
        fails.append(f"A: {len(joints)} jointed lines split, not 6")
    end = tuple(map(float, lines[-1][-1]))
    if (round(end[0], 9), round(end[1], 9)) != (round(TOE[0], 9), TOE[1]):
        fails.append(f"A: the base plane's meshed end is {end}, not the toe "
                     f"corner {TOE}: the bent curve this check is for was not "
                     f"exercised")
    # The first column's top end, on the slope face from (14, 12) to the toe.
    face_a, face_b = np.array([14.0, 12.0]), np.array(TOE)
    top = max((tuple(map(float, lines[0][0])), tuple(map(float, lines[0][-1]))),
              key=lambda p: p[1])
    q = np.array(top) - face_a
    t = face_b - face_a
    off = abs(q[0] * t[1] - q[1] * t[0]) / float(np.hypot(*t))
    if off > 1e-9:
        fails.append(f"A: column 1's meshed top end {top} lies {off:.3g} off the "
                     f"slope face; it stops inside the rock")
    moved = [ln for ln in _LOG.getvalue().splitlines()
             if ln.startswith("Jointed line 1 ends") and "moved onto the boundary" in ln]
    if not moved:
        fails.append("A: the mesher did not say it moved column 1's end onto the "
                     "boundary")


def check_the_completed_file_meshes(fails):
    """B: the completed file meshes to the counts the tutorial prints."""
    from xslope.fileio import load_slope_data
    try:
        mesh, _lines = _mesh(load_slope_data(str(DONE)))
    except Exception as exc:
        fails.append(f"B: the completed file does not mesh: "
                     f"{type(exc).__name__}: {exc}")
        return
    got = (len(mesh['nodes']), len(mesh['elements']),
           len(mesh.get('elements_joint', [])))
    if got != DONE_COUNTS:
        fails.append(f"B: the completed file meshes to {got} "
                     f"(nodes, elements, joint elements), not {DONE_COUNTS}")


def check_a_split_refusal_is_a_message(fails):
    """C: the joint split's and the joint checks' refusals are MeshInputError."""
    import numpy as np
    from xslope.mesh import (MeshInputError, _validate_joint_lines,
                             split_mesh_along_joints)

    mesh = {'nodes': np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
            'elements': np.array([[0, 1, 2]]), 'element_types': np.array([3]),
            'elements_1d': np.zeros((0, 3), dtype=int)}
    try:
        split_mesh_along_joints(mesh, [[(0.0, 0.0), (1.0, 0.0)]], [0])
        fails.append("C: a joint over a mesh with no 1D elements was not refused")
    except MeshInputError:
        pass
    except Exception as exc:
        fails.append(f"C: the split refused with {type(exc).__name__}, not "
                     f"MeshInputError: {exc}")

    # A jointed line touching a bonded one (a reinforcement line or a pile).
    lines = [[(0.0, 0.0), (10.0, 0.0)], [(5.0, -2.0), (5.0, 2.0)]]
    try:
        _validate_joint_lines(lines, {0: {}})
        fails.append("C: a jointed line crossed by a bonded one was not refused")
    except MeshInputError:
        pass
    except Exception as exc:
        fails.append(f"C: the joint check refused with {type(exc).__name__}, "
                     f"not MeshInputError: {exc}")


def check_which_ends_move(fails):
    """D: only an end within a twentieth of an element of the OUTER boundary,
    and clear of a corner, is moved onto it."""
    try:
        from xslope.mesh import _move_joint_ends_onto_boundary
    except ImportError as exc:
        fails.append(f"D: {exc}")
        return
    left = [(0.0, 0.0), (5.0, 0.0), (5.0, 10.0), (0.0, 10.0)]
    right = [(5.0, 0.0), (10.0, 0.0), (10.0, 10.0), (5.0, 10.0)]
    lines = [
        [(7.0, 5.0), (9.97, 5.0)],     # 0.03 short of the right side: moved
        [(7.0, 3.0), (9.9, 3.0)],      # 0.1 short: a joint stopping in rock
        [(7.0, 9.0), (9.98, 9.97)],    # 0.03 short, beside the corner (10, 10)
        [(1.0, 5.0), (4.97, 5.0)],     # 0.03 short of the zone boundary x = 5
    ]
    want = {0: (10.0, 5.0)}
    with contextlib.redirect_stdout(io.StringIO()):
        moved = _move_joint_ends_onto_boundary(
            lines, {0: {}, 1: {}, 2: {}, 3: {}}, [left, right], 1.0)
    got = {li: tuple(round(c, 9) for c in q) for li, _p, q in moved}
    if got != want:
        fails.append(f"D: the ends moved are {got}, not {want}")
    for li, end in ((1, (9.9, 3.0)), (2, (9.98, 9.97)), (3, (4.97, 5.0))):
        if tuple(lines[li][-1]) != end:
            fails.append(f"D: line {li + 1}'s end moved to {lines[li][-1]}")


def run():
    fails = []
    for check in (check_the_reader_model_meshes, check_the_completed_file_meshes,
                  check_a_split_refusal_is_a_message, check_which_ends_move):
        try:
            check(fails)
        except Exception as exc:                             # pragma: no cover
            fails.append(f"{check.__name__} raised {type(exc).__name__}: {exc}")
    return fails


if __name__ == '__main__':
    failures = run()
    for f in failures:
        print("FAIL", f)
    print("joint_end_snap:", "FAIL" if failures else "PASS")
    sys.exit(1 if failures else 0)
