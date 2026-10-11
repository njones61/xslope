"""Solve Tutorial FEM-5's toppling model and write its companions.

FEM-5 part 2 is the slow half of the page: a base plane with a generated set of
steep columns on it, six jointed lines between them, and a strength reduction
whose standing trials each grow slip over tens of thousands of viscoplastic
sweeps. So the toppling model is the one the tutorial ships **solved** — a reader
opens the workbook, Studio attaches the companions beside it, and the finished
strength reduction is on the screen without a solve.

What is written, beside ``docs/tutorials/files/xslope_rock_toppling.xlsx``:

    xslope_rock_toppling_mesh.json        the mesh the fields were solved on
    xslope_rock_toppling_fem_nodes.csv    the last standing trial's field
    xslope_rock_toppling_fem_elements.csv
    xslope_rock_toppling_fem_joints.csv   every node pair's tractions, slip and
                                          state — the solve's own history, which
                                          nothing can recover from displacements
    xslope_rock_toppling_fem_meta.json    the run record: bracket, tolerance,
                                          criterion, budget, final interval, the
                                          factor of safety, and every trial's
                                          verdict

The run record is not optional. ``export_fem_solution`` writes the CSVs whether
or not it is given one, and a set written without it reloads with the factor of
safety blank: the field is there, the joints draw, and the number the page is
about is missing. ``ssrm_run_record`` is what carries it, so it is built and
passed here, and the check at the end reloads the set and reads the factor of
safety back out of it.

The at-failure snapshot is not written. It is a second solve past the bracket,
and what this page shows is the slope standing at its critical state; the
``capture_failure_state=False`` below is what keeps it out of both the run and
the shipped file set.

Deterministic: same mesh, same bracket, same trials. Solved on the pure-NumPy
reference kernel, the path the page's numbers are measured on.

    PYTHONPATH=. python3 tools/make_rock_joints_sidecars.py
"""
import contextlib
import io
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import run_tests as RT                                                  # noqa: E402
import xslope.fem as _fem                                               # noqa: E402
from xslope.fem import (build_fem_data, export_fem_solution,            # noqa: E402
                        import_fem_solution, solve_ssrm, ssrm_run_record)
from xslope.fileio import load_slope_data                               # noqa: E402
from xslope.mesh import (build_mesh_from_polygons, export_mesh_to_json,  # noqa: E402
                         extract_constraint_line_geometry,
                         extract_joint_options, extract_point_constraints,
                         extract_size_regions, get_material_polygons)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL = os.path.join(REPO_ROOT, "docs", "tutorials", "files",
                     "xslope_rock_toppling.xlsx")

#: Build Mesh and Run FEM exactly as the page tells the reader to set them:
#: quadratic triangles at 1.5 m, a strength reduction over [1.0, 2.0] to a 0.01
#: bracket, and 100,000 sweeps a trial. ``hybrid`` is the criterion the Run FEM
#: dialog opens on for a jointed model, and the one the answer below is read
#: under: a standing trial on a jointed model never settles on force alone,
#: because the pairs at their limit keep slipping and re-sticking, so
#: ``non_convergence`` would count every near-critical trial as a failure.
ELEMENT_TYPE = "tri6"
TARGET_SIZE = 1.5
F_MIN, F_MAX = 1.0, 2.0
TOLERANCE = 0.01
CRITERION = "hybrid"
MAX_ITERATIONS = 100000
#: What the page prints. A run that lands anywhere else is reported and NOT
#: written: a companion that disagrees with the sentence beside it on the page is
#: the defect this script exists to prevent.
EXPECTED_FS = 1.1055


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def build_mesh(slope_data):
    """The mesh Studio's Build Mesh dialog builds for this model.

    The joint lines go in twice: as constraint lines, so every trace lies on
    element edges, and through ``extract_joint_options``, which is what tells the
    mesher to split the mesh along them. Without the second the traces are meshed
    and bonded and the model has no joints in it at all.
    """
    lines, _n_reinf, _n_pile = extract_constraint_line_geometry(slope_data)
    return _quiet(build_mesh_from_polygons,
                  get_material_polygons(slope_data, reinf_lines=lines),
                  TARGET_SIZE, ELEMENT_TYPE, lines=lines or None,
                  element_size_1d=slope_data.get("element_size_1d"),
                  point_constraints=extract_point_constraints(slope_data),
                  size_regions=extract_size_regions(slope_data),
                  joint_lines=extract_joint_options(slope_data))


def verify(fem_data, stem, fs):
    """Read the set back the way Studio reads it when the workbook is opened.

    Three things have to come back, and each of them is a way the set has failed
    before: the joint state (without ``_fem_joints.csv`` a reloaded jointed field
    draws no joints), no ``sidecar_notes`` (a note is a companion the importer
    refused rather than grafted), and the factor of safety out of the meta
    sidecar (which is what a set written with no ``meta=`` loses).
    """
    import json

    back = _quiet(import_fem_solution, fem_data, stem)
    notes = back.get("sidecar_notes")
    if notes:
        raise RuntimeError(f"the reloaded solution carries sidecar notes: {notes}")
    n_pairs = int((fem_data.get("joint_data") or {}).get("n", 0)) * 3
    slip = back.get("joint_slip")
    if slip is None or getattr(slip, "size", 0) != n_pairs:
        raise RuntimeError("the reloaded solution carries no joint state — "
                           f"expected {n_pairs} node pairs, got "
                           f"{None if slip is None else getattr(slip, 'size', None)}")
    with open(f"{stem}_fem_meta.json") as fh:
        meta = json.load(fh)
    if meta.get("FS") is None:
        raise RuntimeError("the meta sidecar records no factor of safety")
    if abs(float(meta["FS"]) - fs) > 1e-9:
        raise RuntimeError(f"the meta sidecar reads FS = {meta['FS']}, the run "
                           f"gave {fs}")
    print(f"  reloaded: {int(slip.size)} node pairs of joint state, no sidecar "
          f"notes, meta FS = {float(meta['FS']):.4f}, criterion "
          f"{meta.get('failure_criterion')!r}, bracket "
          f"[{meta.get('F_min')}, {meta.get('F_max')}]")


def main():
    t0 = time.time()
    slope_data = load_slope_data(MODEL)
    mesh = build_mesh(slope_data)
    fem_data = _quiet(build_fem_data, slope_data, mesh)
    jd = fem_data.get("joint_data") or {}
    print(f"{os.path.basename(MODEL)}: {len(mesh['nodes'])} nodes, "
          f"{len(mesh['elements'])} elements, {int(jd.get('n', 0))} joint "
          f"elements on {len(jd.get('jointed_lines', []) or [])} jointed lines")

    options = {"F_min": F_MIN, "F_max": F_MAX, "tolerance": TOLERANCE,
               "failure_criterion": CRITERION}
    with RT._force_fast_kernel(_fem, False):
        result = _quiet(solve_ssrm, fem_data, debug_level=0,
                        capture_failure_state=False,
                        max_iterations=MAX_ITERATIONS, **options)
    if not result.get("converged"):
        raise RuntimeError("the strength reduction did not resolve a bracket")
    fs = result["FS"]
    if abs(fs - EXPECTED_FS) > TOLERANCE:
        raise RuntimeError(f"solved FS = {fs:.6f}, the page states "
                           f"{EXPECTED_FS:.4f} +/- {TOLERANCE:g} — "
                           f"nothing written")
    field = result["last_solution"]
    meta = ssrm_run_record(result, fem_data, options)
    meta.update({"analysis_type": "ssrm", "FS": fs,
                 "element_type": ELEMENT_TYPE, "target_size": TARGET_SIZE,
                 "max_iter": MAX_ITERATIONS,
                 # The reduction factor the WRITTEN FIELD was solved at — the last
                 # standing trial, not the bracket's midpoint. Neither CSV carries
                 # it, and Studio's own writer records it.
                 "F": field.get("F")})
    stem = os.path.splitext(MODEL)[0]
    _quiet(export_fem_solution, fem_data, field, stem, meta=meta)
    export_mesh_to_json(mesh, f"{stem}_mesh.json")
    print(f"  FS = {fs:.4f} (page {EXPECTED_FS:.4f}) from "
          f"[{result['final_interval'][0]:.6f}, {result['final_interval'][1]:.6f}], "
          f"{len(result['trials'])} trials, field at F = {field.get('F'):.6f}, "
          f"{time.time() - t0:.0f}s")
    verify(fem_data, stem, fs)
    for name in sorted(os.listdir(os.path.dirname(stem))):
        if name.startswith(os.path.basename(stem) + "_"):
            print(f"  {name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
