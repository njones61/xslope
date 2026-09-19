"""Regenerate Tutorial FEM-3's plot figures in ``docs/tutorials/images/``.

FEM-3 puts slip joints into a structure: a segmental block wall standing on the
contacts between its own courses, the same wall with three geogrid layers tied
into it, and two bonded-against-jointed pairs that show what the Joint column
does to a sheet under an embankment.  Every figure here is one standard plot call
on one of the seven models ``tools/build_block_wall.py`` writes, at the run
settings the page tells the reader to use:

    the wall models   tri6, global target size 0.8 m (each block polygon carries
                      its own 0.3 m Size, which the mesher reads from the file)
    part 3's models   tri6, target size 1.2 m
    both              SSRM, F in [1.0, 2.0], tolerance 0.01, failure criterion
                      Hybrid, 100,000 sweeps a trial

Nothing is restyled: the plots are called the way Studio calls them, with
Studio's own titles and Studio's own autoscaling, so the reader's screen and the
committed figure are the same picture.  What cannot be drawn is printed instead —
the mesh and interface counts, the bracket walk trial by trial, and how far each
jointed line has slipped at the last standing trial, which is the measurement the
page's slip table is read from.

The Studio dialog captures are a separate producer, because they need Qt:
``tools/capture_block_wall_screenshots.py``.  The other tutorials' figures are
``tools/make_tutorial_figures.py``, whose ``capture`` helper this module reuses
rather than copies.

Run:  PYTHONPATH=. python3 tools/make_block_wall_figures.py             # everything
      PYTHONPATH=. python3 tools/make_block_wall_figures.py inputs      # one group
      PYTHONPATH=. python3 tools/make_block_wall_figures.py wall grid   # by name
"""

from __future__ import annotations

import contextlib
import io
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import matplotlib
matplotlib.use("Agg")

from xslope.fileio import load_slope_data                            # noqa: E402
from xslope.plot import declared_unit_labels, plot_inputs            # noqa: E402

# The capture idiom — swap ``plt.show`` for a save, one PNG per state — belongs to
# the tutorial figure producer and is used here unchanged, so a change to how a
# tutorial figure is written reaches FEM-3's figures too.
from tools.make_tutorial_figures import OUT_DIR, capture              # noqa: E402,F401

FILES = os.path.join(REPO_ROOT, "docs", "tutorials", "files")
#: Where the per-member profiles behind the detail figures are written.  A
#: profile is the solved field read along one line, and re-deriving one costs the
#: whole strength reduction again — a quarter of an hour on the wall — so every
#: profile a figure is drawn from is also written out as a CSV.  A later question
#: about what the figure shows is then answered from the file rather than from a
#: second run, and never by eye off the PNG.  Outside the repo: these are working
#: measurements of a round, not content of the package.
OUT_DATA = os.path.join(os.path.expanduser("~"), "python_projects", "xslope_private",
                        "reports", "campaign_joints_2026-09", "r34_data",
                        "fem03_figures")
FEM03_START = os.path.join(FILES, "xslope_block_wall_start.xlsx")
FEM03_WALL = os.path.join(FILES, "xslope_block_wall.xlsx")
FEM03_GRID = os.path.join(FILES, "xslope_block_wall_grid.xlsx")
FEM03_SHEET_BONDED = os.path.join(FILES, "xslope_base_geotextile_bonded.xlsx")
FEM03_SHEET_JOINTED = os.path.join(FILES, "xslope_base_geotextile_jointed.xlsx")
FEM03_LINER_BONDED = os.path.join(FILES, "xslope_liner_bonded.xlsx")
FEM03_LINER_JOINTED = os.path.join(FILES, "xslope_liner_jointed.xlsx")

#: Quadratic triangles, because a joint is read on the relative displacement of
#: two faces and a quadratic interface element carries a mid-side pair as well as
#: its ends.  The size is the GLOBAL target; each of the six block polygons
#: carries its own 0.3 m Size in the file, which a 0.6 m course needs and which
#: costs far fewer nodes than refining the whole 24 m section.
FEM03_ELEMENT_TYPE = "tri6"
FEM03_TARGET_SIZE = 0.8
#: Part 3's embankment is a 60 m section with one sheet in it and no thin zone, so
#: it is meshed coarser than the wall.
FEM03_SHEET_TARGET_SIZE = 1.2
#: Run FEM: the dialog's own bracket and bisection tolerance, and ``hybrid`` —
#: which is the criterion the dialog opens on for a jointed model.  On a jointed
#: model almost all of the out-of-balance force sits on the joints, and a pair of
#: faces at its limit alternates between slipping and sticking, so a standing
#: trial never converges on force alone (docs/fem/joints.md, "How a trial is
#: decided").  Part 3's four models are run at the same criterion so the bonded
#: half and the jointed half of a pair are decided the same way.
FEM03_CRITERION = "hybrid"
FEM03_F_MIN, FEM03_F_MAX = 1.0, 2.0
FEM03_TOLERANCE = 0.01
#: 100,000 sweeps a trial — the model checks' floor for a jointed model and the
#: Run FEM box's own ceiling.  A joint reaches equilibrium by growing slip a
#: little per sweep, so a jointed trial settles over tens of thousands of sweeps
#: where a bonded one settles over hundreds, and a trial that runs out of budget
#: is recorded undecided and read as not standing — which reports the budget
#: rather than the wall.
FEM03_MAX_ITERATIONS = 100000
#: The layer the page reads the 1D Details panel on: the middle of the three, the
#: one the block column's own movement develops tension in from both ends.
FEM03_DETAIL_LINE = "grid-02"


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _mesh(model, target_size):
    """The mesh Studio's Build Mesh dialog builds for this model.

    The constraint lines are carried in twice, and both are needed: as lines, so
    every joint trace and every sheet lies on element edges, and through
    ``extract_joint_options``, which is what tells the mesher to SPLIT the mesh
    along the jointed ones and give every node on such a trace one copy per wedge
    of material around it.  Without the second the traces are meshed and bonded
    and the model has no joints in it at all.
    """
    from xslope.mesh import (build_mesh_from_polygons,
                             extract_constraint_line_geometry,
                             extract_joint_options,
                             extract_point_constraints,
                             extract_size_regions, get_material_polygons)

    lines, _n_reinf, _n_pile = extract_constraint_line_geometry(model)
    with contextlib.redirect_stdout(io.StringIO()):
        return build_mesh_from_polygons(
            get_material_polygons(model, reinf_lines=lines),
            target_size, FEM03_ELEMENT_TYPE, lines=lines or None,
            element_size_1d=model.get("element_size_1d"),
            point_constraints=extract_point_constraints(model),
            size_regions=extract_size_regions(model),
            joint_lines=extract_joint_options(model))


def _solve(model, mesh):
    """One strength reduction at the page's settings, on the reference kernel.

    ``fast_kernel=False`` is pinned the way every producer of a committed number
    pins it: the compiled Mohr-Coulomb kernel is built on some machines and not
    others, and on this wall the two paths do not agree — so the figure has to
    carry the answer the reference path gives.  Returns
    ``(fem_data, result, seconds)``.
    """
    import time

    import run_tests as RT
    import xslope.fem as _fem
    from xslope.fem import build_fem_data, solve_ssrm

    fem_data = build_fem_data(model, mesh)
    t0 = time.time()
    with RT._force_fast_kernel(_fem, False):
        with contextlib.redirect_stdout(io.StringIO()):
            result = solve_ssrm(fem_data, F_min=FEM03_F_MIN, F_max=FEM03_F_MAX,
                                tolerance=FEM03_TOLERANCE, debug_level=0,
                                failure_criterion=FEM03_CRITERION,
                                max_iterations=FEM03_MAX_ITERATIONS)
    return fem_data, result, time.time() - t0


def _dump(tag, prof):
    """Write one member's profile out as a CSV, station by station.

    ``fem_details.write_profile_csv`` is the package's own exporter, so the file
    carries the same columns the detail view reads and nothing is transcribed.
    """
    from xslope.fem_details import write_profile_csv

    os.makedirs(OUT_DATA, exist_ok=True)
    path = os.path.join(OUT_DATA, "%s.csv" % tag)
    write_profile_csv(prof, path)
    print("   wrote       %s" % path)
    return path


def _stations(prof, units):
    """Every station of one jointed line, printed.

    The figure shows a shape; this is what the shape is made of — where along the
    line each station sits, what the bar is carrying there against its capacity,
    what the interface is carrying against its own limit, and how far the two
    faces have moved past each other.  The page quotes these rather than reading
    the panel by eye.
    """
    import numpy as np

    s = np.asarray(prof["s"], dtype=float)
    ts = np.asarray(prof["ts"], dtype=float)
    tlim = np.asarray(prof["tlim"], dtype=float)
    slip = np.asarray(prof["slip"], dtype=float)
    slipping = np.asarray(prof["slipping"], dtype=bool)
    opened = np.asarray(prof["open"], dtype=bool)
    bar_s = np.asarray(prof["bar_s"], dtype=float)
    bar_T = np.asarray(prof["bar_T"], dtype=float)
    bar_cap = np.asarray(prof["bar_cap"], dtype=float)
    print("   stations    %s — %d station(s) over %.4f %s"
          % (prof["label"], len(s), prof["length"], units))
    for k in range(len(s)):
        print("        s %6.3f  ts %10.4f  limit %10.4f  |ts|/limit %6.3f  "
              "slip %10.6f  %s"
              % (s[k], ts[k], tlim[k],
                 abs(ts[k]) / tlim[k] if tlim[k] > 1e-12 else float("nan"),
                 slip[k],
                 "OPEN" if opened[k] else ("slipping" if slipping[k] else "")))
    # The bar is sampled at its own element centroids, not at the interface's
    # stations, so it is printed as its own short list rather than folded into
    # the rows above.
    for k in range(len(bar_s)):
        print("        bar s %6.3f  T %10.4f  capacity %10.4f  T/cap %6.3f"
              % (bar_s[k], bar_T[k], bar_cap[k],
                 bar_T[k] / bar_cap[k] if bar_cap[k] > 1e-12 else float("nan")))
    # How far back from end 1 the interface is actually at its limit: the span
    # the flagged stations cover, which is the number the page needs rather than
    # a count of stations.
    flagged = np.flatnonzero(slipping | opened)
    if len(flagged):
        print("        at limit from s %.3f to s %.3f (%.3f %s of %.3f %s), "
              "largest slip %.6f %s"
              % (s[flagged[0]], s[flagged[-1]],
                 s[flagged[-1]] - s[flagged[0]], units, prof["length"], units,
                 float(slip.max()), units))
    else:
        print("        no station at its limit · largest slip %.6f %s"
              % (float(slip.max()) if len(slip) else float("nan"), units))


def _counts(label, model, mesh, fem_data):
    """What the meshed model came to: the counts the page prints."""
    jd = fem_data.get("joint_data") or {}
    print("   %-11s %d nodes · %d elements · %d joint elements on %d jointed "
          "line(s) · %d row(s) on the joints sheet · %d reinforcement line(s) · "
          "%s at %g %s"
          % (label, len(mesh["nodes"]), len(mesh["elements"]),
             int(jd.get("n", 0)), len(jd.get("jointed_lines", []) or []),
             len(model.get("joint_lines") or []),
             len(model.get("reinforcement_lines") or []),
             model["element_type"], model["target_size"],
             declared_unit_labels(model)["length"]))


def _report(label, result, seconds):
    """The bracket walk: what each trial was asked, and what it answered."""
    print("   %-11s FS %.4f from [%.6f, %.6f] (width %.6f) after %d bisection "
          "step(s) · %.0f s"
          % (label, result["FS"], result["final_interval"][0],
             result["final_interval"][1], result["interval_width"],
             result["iterations_ssrm"], seconds))
    for tr in result["trials"]:
        print("        F %.4f  %-6s  %-13s  %s sweeps"
              % (tr["F"], tr.get("role"), tr.get("verdict"), tr.get("iterations")))
    last = result["last_solution"]
    print("        last standing trial F %.4f · %s · %s sweeps (%s)"
          % (last.get("F", float("nan")),
             "equilibrium" if last.get("converged") else "no equilibrium",
             last.get("iterations"), last.get("exit_reason")))


def _joints(model, fem_data, result, dump=None):
    """Which jointed line slipped at the last standing trial, and how far.

    Read station by station rather than element by element: a quadratic interface
    element carries three node pairs and shares its end pairs with the element
    next along it, so counting slots would count the shared pairs twice.
    ``fem_details.joint_profile`` keys its record by the split node itself, which
    makes exactly one record per node pair on the line, and names the line with
    the label the joints sheet gave it.

    The raw arrays are read as well, and the two are printed side by side: the
    per-slot ``solution["joint_slip"]`` / ``["joint_slipping"]`` against
    ``fem_data["joint_data"]["line_id"]`` is the measurement, the profile is what
    the reader sees, and they have to agree.
    """
    import numpy as np

    from xslope import fem_details

    last = result["last_solution"]
    units = declared_unit_labels(model)["length"]
    jd = fem_data["joint_data"]
    line_of = np.asarray(jd["line_id"], dtype=int)
    # Slip is SIGNED — which way the two faces went past each other — so every
    # reading of it is taken in magnitude.  The unsigned maximum of a line whose
    # faces all slid the same, negative way is zero, which reads as a joint that
    # never moved.
    slip_raw = np.abs(np.asarray(last["joint_slip"], dtype=float)).reshape(-1)
    slipping_raw = np.asarray(last["joint_slipping"], dtype=bool).reshape(-1)
    # Every interface element carries the same number of slots, so the per-element
    # line id repeats across its slots.
    per_el = slip_raw.size // int(jd["n"])
    line_slot = np.repeat(line_of, per_el)
    for line_id in fem_details.joint_line_ids(fem_data, last):
        prof = fem_details.joint_profile(fem_data, last, line_id,
                                         slope_data=model,
                                         field_state="converged")
        slipping = np.asarray(prof["slipping"], dtype=bool)
        slip = np.asarray(prof["slip"], dtype=float)
        sel = line_slot == line_id
        print("   joint %-2d    %-11s %d of %d node pair(s) slipping · largest "
              "slip %.6f %s · %d open · peak |ts|/limit %s · %s   "
              "[raw: %d of %d slot(s) slipping, max slip %.6f]"
              % (line_id, prof["label"], int(slipping.sum()), len(slip),
                 float(slip.max()) if len(slip) else float("nan"), units,
                 int(np.asarray(prof["open"], dtype=bool).sum()),
                 "n/a" if prof["peak_utilization"] is None
                 else "%.3f" % prof["peak_utilization"], prof["status"],
                 int(slipping_raw[sel].sum()), int(sel.sum()),
                 float(slip_raw[sel].max()) if sel.any() else float("nan")))
        if dump:
            _dump("%s_joint_%02d_%s" % (dump, line_id,
                                        str(prof["label"]).replace(" ", "_")),
                  prof)


# --------------------------------------------------------------------------- #
# A — the figures that need no solve
# --------------------------------------------------------------------------- #
def fem03_inputs():
    """The section as the reader opens it, the seven joint lines they type into
    it, the three geogrid layers of part 2, and the split mesh part 1 runs on.

    Printed rather than drawn: the four materials, every joint line's own
    geometry and strength, and the mesh and interface counts of both wall models —
    which is what says the 0.3 m block Size in the file reached the mesher.
    """
    from xslope.fem import build_fem_data
    from xslope.plot_fem import plot_fem_data

    # ---- the section as it opens -------------------------------------------- #
    start = load_slope_data(FEM03_START)
    _u = declared_unit_labels(start)
    for m in start["materials"]:
        print("   material    %-16s γ %g %s · %s · c %g %s · φ %g° · E %g %s · ν %g"
              % (m["name"], m["gamma"], _u["unit_weight"], m["option"], m["c"],
                 _u["stress"], m["phi"], m["E"], _u["stress"], m["nu"]))
    print("   starter     %d polygon(s) · %d row(s) on the joints sheet"
          % (len(start["polygons"]), len(start["joint_lines"] or [])))
    capture("fem03_inputs_start.png", plot_inputs, start, mode="fem",
            title="Slope Geometry and Inputs")

    # ---- part 1: the seven contacts ----------------------------------------- #
    wall = load_slope_data(FEM03_WALL)
    for row in wall["joint_lines"]:
        print("   joint line  %-11s (%.4f, %.4f) to (%.4f, %.4f) · c %s %s · "
              "φ %g° · t_cut %s · kn %s · ks %s"
              % (row["label"], row["x1"], row["y1"], row["x2"], row["y2"],
                 row["c"], _u["stress"], row["phi"], row["t_cut"], row["kn"],
                 row["ks"]))
    capture("fem03_inputs_joints.png", plot_inputs, wall, mode="fem",
            title="Slope Geometry and Inputs")

    # ---- part 2: the three layers ------------------------------------------- #
    grid = load_slope_data(FEM03_GRID)
    for row in grid["reinforcement_lines"]:
        print("   geogrid     %-11s (%.4f, %.4f) to (%.4f, %.4f) · Tmax %g %s · "
              "Tend1 %g %s · adhesion %g %s · delta %g° · Joint %r"
              % (row["label"], row["x1"], row["y1"], row["x2"], row["y2"],
                 row["t_max"], _u["force_per_len"],
                 row["tend1"], _u["force_per_len"],
                 row["adhesion"], _u["stress"], row["delta"], row["joint"]))
    capture("fem03_inputs_grid.png", plot_inputs, grid, mode="fem",
            title="Slope Geometry and Inputs")

    # ---- the mesh part 1 runs on -------------------------------------------- #
    mesh = _mesh(wall, FEM03_TARGET_SIZE)
    fem_data = build_fem_data(wall, mesh)
    _counts("wall", wall, mesh, fem_data)
    capture("fem03_mesh.png", plot_fem_data, fem_data)

    # The grid model's mesh is not drawn — the page shows one mesh — but its
    # counts are what say the three layers cost what the page claims they do.
    mesh_g = _mesh(grid, FEM03_TARGET_SIZE)
    fem_data_g = build_fem_data(grid, mesh_g)
    _counts("wall + grid", grid, mesh_g, fem_data_g)
    capture("fem03_mesh_grid.png", plot_fem_data, fem_data_g)


# --------------------------------------------------------------------------- #
# C1 — part 1: the wall on its seven contacts
# --------------------------------------------------------------------------- #
def fem03_wall():
    """The strength reduction on the wall as built, and the two panels the page
    reads it from: the blocks, and the shear strain in the soil around them.

    On a jointed model the displacement panel IS the block picture — the scaled
    deformed mesh with both faces of every joint drawn, because a jointed model's
    mechanism is wedges moving as bodies on their contacts and an arrow field
    sampled at nodes says nothing about the contacts.
    """
    from xslope.plot_fem import plot_fem_results

    wall = load_slope_data(FEM03_WALL)
    mesh = _mesh(wall, FEM03_TARGET_SIZE)
    fem_data, result, seconds = _solve(wall, mesh)
    _counts("wall", wall, mesh, fem_data)
    _report("wall", result, seconds)
    _joints(wall, fem_data, result, dump="wall")
    _both_states("fem03_fem_blocks.png", fem_data, result, "displace_vector")
    _both_states("fem03_fem_shear.png", fem_data, result, "shear_strain")
    _keep_solution("wall", fem_data, result)

    # The back face's 1D details: the panel the page's slip table is read from,
    # shown for the contact that carries the failure so the reader can see the
    # slip curve and the slipping stations the table counts.
    from xslope import fem_details
    from xslope.plot_fem_details import plot_detail

    last = result["last_solution"]
    for line_id in fem_details.joint_line_ids(fem_data, last):
        prof = fem_details.joint_profile(fem_data, last, line_id,
                                         slope_data=wall, field_state="converged")
        if str(prof.get("label", "")).strip().lower() == "back face":
            capture("fem03_1d_details_back_face.png", plot_detail, prof)
            break
    else:
        raise SystemExit("no joint line labelled 'back face' on the wall")


# --------------------------------------------------------------------------- #
# C2 — part 2: the same wall with three geogrid layers
# --------------------------------------------------------------------------- #
def fem03_grid():
    """The same wall with the geogrid tied in: the blocks again, and the 1D
    Details panel for the middle layer.

    A layer declared ``Joint = Yes`` is two members on one chord — the bar, whose
    tension is read against its capacity, and the interface the mesh was split
    along, whose shear traction is read against its Mohr-Coulomb limit — and the
    joint profile carries both, which is why the panel is drawn from it rather
    than from the reinforcement profile.
    """
    from xslope import fem_details
    from xslope.plot_fem_details import plot_detail
    from xslope.plot_fem import plot_fem_results

    grid = load_slope_data(FEM03_GRID)
    mesh = _mesh(grid, FEM03_TARGET_SIZE)
    fem_data, result, seconds = _solve(grid, mesh)
    _counts("wall + grid", grid, mesh, fem_data)
    _report("wall + grid", result, seconds)
    _joints(grid, fem_data, result, dump="grid")
    _both_states("fem03_fem_blocks_grid.png", fem_data, result, "displace_vector")
    _both_states("fem03_fem_shear_grid.png", fem_data, result, "shear_strain")
    _keep_solution("wall_grid", fem_data, result)

    last = result["last_solution"]
    rows = fem_details.list_lines(fem_data, last, grid)
    for row in rows:
        print("   member      %-14s %-11s %d station(s) · utilization %s · %s"
              % (row["kind"], row["label"], row["n_elements"],
                 "n/a" if row["utilization"] is None
                 else "%.3f" % row["utilization"], row["status"]))
    want = [r for r in rows
            if r["kind"] == "joint" and r["label"] == FEM03_DETAIL_LINE]
    if not want:
        raise SystemExit("no joint member named %r among %s"
                         % (FEM03_DETAIL_LINE,
                            [r["label"] for r in rows]))
    # All three layers are read station by station, not just the one the figure
    # is drawn from: the page says whether the middle layer is representative,
    # and that is a claim about the other two.
    units = declared_unit_labels(grid)["length"]
    for row in rows:
        if row["kind"] == "joint" and str(row["label"]).startswith("grid-"):
            _stations(fem_details.joint_profile(fem_data, last, row["index"],
                                                slope_data=grid,
                                                field_state="converged"), units)

    prof = fem_details.joint_profile(fem_data, last, want[0]["index"],
                                     slope_data=grid, field_state="converged")
    print("   detail      %s · %d station(s) · peak bar tension %.4f of "
          "capacity %.4f · peak |ts| %.4f against limit %.4f · largest slip %.6f"
          % (prof["label"], len(prof["s"]),
             float(max(abs(prof["bar_T"]))) if len(prof["bar_T"]) else float("nan"),
             float(max(prof["bar_cap"])) if len(prof["bar_cap"]) else float("nan"),
             float(max(abs(prof["ts"]))) if len(prof["ts"]) else float("nan"),
             float(max(prof["tlim"])) if len(prof["tlim"]) else float("nan"),
             float(max(prof["slip"])) if len(prof["slip"]) else float("nan")))
    capture("fem03_1d_details.png", plot_detail, prof)


# --------------------------------------------------------------------------- #
# C3 — part 3: the two bonded-against-jointed pairs
# --------------------------------------------------------------------------- #
#: The four runs, in the order the page reads them: a sheet the surface has to
#: cross, then a sheet the mass can slide along, each bonded and jointed.
FEM03_PAIRS = (
    ("sheet bonded", FEM03_SHEET_BONDED, "fem03_shear_sheet_bonded.png"),
    ("sheet jointed", FEM03_SHEET_JOINTED, "fem03_shear_sheet_jointed.png"),
    ("liner bonded", FEM03_LINER_BONDED, "fem03_shear_liner_bonded.png"),
    ("liner jointed", FEM03_LINER_JOINTED, "fem03_shear_liner_jointed.png"),
)


def fem03_sheets():
    """Part 3's four runs, each drawn as its own viscoplastic shear strain panel.

    The four panels are autoscaled to their own ranges, the way Studio draws
    every field state.  What the page reads off them is WHERE the band runs —
    across the sheet in the geotextile pair, along it in the liner pair — and a
    shared color scale would be a fifth thing on the page rather than the picture
    the reader's own run produces.
    """
    from xslope.plot_fem import plot_fem_results

    for label, path, name in FEM03_PAIRS:
        model = load_slope_data(path)
        mesh = _mesh(model, FEM03_SHEET_TARGET_SIZE)
        fem_data, result, seconds = _solve(model, mesh)
        _counts(label, model, mesh, fem_data)
        _report(label, result, seconds)
        if (fem_data.get("joint_data") or {}).get("n"):
            _joints(model, fem_data, result)
        # Both field states, so the page can choose: the last converged trial
        # (the slope still standing) and the captured at-failure state (the
        # mechanism). The solve itself is kept beside the report so either can
        # be redrawn without re-solving.
        # The deformation panel first (the blocks with their faces on a jointed
        # model, the deformed mesh on a bonded one), then the strain panel,
        # each in both field states.
        _both_states(name.replace("fem03_shear_", "fem03_deform_"), fem_data,
                     result, _deform_panel(fem_data))
        _both_states(name, fem_data, result, "shear_strain")
        _keep_solution(label, fem_data, result)


def _deform_panel(fem_data):
    """The deformation panel a model gets: the blocks (displace_vector routes
    there on a jointed model) or the deformed mesh on a model with no joint."""
    return ("displace_vector" if (fem_data.get("joint_data") or {}).get("n")
            else "deformation")


def _both_states(name, fem_data, result, plot_type):
    """One results panel in both field states: the last converged trial under
    ``name`` and the captured at-failure state under ``name`` + ``_failure``."""
    from xslope.plot_fem import plot_fem_results

    fail = result.get("failure_solution")
    capture(name, plot_fem_results, fem_data, result["last_solution"],
            plot_type=plot_type, fs=result["FS"], failure_solution=fail,
            field_state="converged")
    if fail is not None:
        capture(name.replace(".png", "_failure.png"), plot_fem_results, fem_data,
                result["last_solution"], plot_type=plot_type, fs=result["FS"],
                failure_solution=fail, field_state="failure")


def _keep_solution(label, fem_data, result):
    """Pickle a part 3 solve under the private campaign directory so both field
    states can be redrawn and compared without the solve."""
    import pickle

    d = os.path.join(os.path.expanduser("~"), "python_projects", "xslope_private",
                     "reports", "campaign_joints_2026-09", "fem03_part3_solutions")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "%s.pkl" % label.replace(" ", "_")), "wb") as fh:
        pickle.dump({"fem_data": fem_data, "result": result}, fh)


GROUPS = {
    "fem03_inputs": fem03_inputs,
    "fem03_wall": fem03_wall,
    "fem03_grid": fem03_grid,
    "fem03_sheets": fem03_sheets,
}


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    os.makedirs(OUT_DIR, exist_ok=True)
    names = [n for n in GROUPS if not argv or any(a in n for a in argv)]
    if not names:
        print("no figure group matching %s; known groups: %s"
              % (argv, ", ".join(sorted(GROUPS))))
        return 1
    for name in names:
        print("== %s" % name)
        GROUPS[name]()
    print("\nwrote %d group(s) to docs/tutorials/images/" % len(names))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
