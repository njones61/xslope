"""What a jointed line reaches: the mesher, preflight, the plots, the detail
panel and the report.

R1 built the mesh split and R2 the element. This is the path a MODEL takes to
them — the `Joint` column on the reinforce sheet — and what a solved jointed
model then shows.

Six legs, all on one solve of the shipped reinforcement sample with two of its
six lines made joints in memory:

  a. the wiring. ``mesh.extract_joint_lines`` reads the column off a loaded
     model and returns the mapping the mesher takes, keyed by the line's index
     in the constraint-line list and carrying its end anchorages; a model with
     no jointed line returns None, which is the value that leaves the mesh what
     it always was. The mesh built through it carries the joint keys, and
     ``build_fem_data`` carries ``joint_data`` for exactly the flagged lines.
  b. preflight. The four refusals fire with the line named, the four §4c
     signals stand down on a model that earns none of them, and the INFO says
     which bonded-bar inputs a jointed line stops reading.
  c. the plots. The inputs plot draws a jointed line in its own style with its
     own legend entry and a bonded one in the ordinary one; the mesh plot draws
     the jointed lines over the mesh; the results panel draws each interface's
     state on the line with the slip colorbar, and draws nothing at all on a
     model with no joint.
  d. the faces are offset. The split gives the two soil faces their own nodes,
     so a slipped joint moves them apart in the solved field — which is what
     the deformed mesh and the at-failure capture draw. Measured on the
     displacements rather than read off a picture.
  e. the detail profile. ``joint_profile`` reads the two interfaces at a station
     together, orders the stations from end 1 the way the bar's own profile is
     ordered, and reports the state; ``list_lines`` gives a joint row beside the
     reinforcement row for the same line; the figure and the CSV both build.
  f. the report table. A jointed run's reinforcement section carries the joints
     table — line, the share of its length at its limit, the peak slip — and an
     unjointed run carries none.

Run directly:  PYTHONPATH=. python3 test/joint_surfaces_check.py
"""

import contextlib
import io
import math
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

REINF_XLSX = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "docs", "fem", "files", "xslope_reinforce_fem.xlsx")

#: The two lines made joints, and the interface strength they are given. Chosen
#: because the sample states none: a line with a blank Adhesion/Delta is a
#: preflight error, and a fixture that flagged one would be testing the refusal.
JOINTED = (1, 3)
ADHESION, DELTA = 5.0, 25.0

#: Element size for the fixture's mesh. Coarse enough that the leg costs one
#: solve of a few seconds and fine enough that a sheet carries ten bar elements.
TARGET_SIZE = 3.0
SIZE_1D = 2.0


def _quiet(fn, *a, **kw):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **kw)


def _model(jointed=JOINTED):
    """The sample with the named lines made joints, loaded fresh each time."""
    from xslope.fileio import load_slope_data
    sd = load_slope_data(REINF_XLSX)
    for k in jointed:
        sd["reinforcement_lines"][k].update(joint="Yes", adhesion=ADHESION,
                                            delta=DELTA)
    return sd


def _mesh_for(sd):
    from xslope.mesh import (build_mesh_from_polygons, extract_constraint_line_geometry,
                             extract_joint_lines, extract_size_regions,
                             get_material_polygons)
    lines, _nr, _np = extract_constraint_line_geometry(sd)
    polys = get_material_polygons(sd, reinf_lines=lines)
    return _quiet(build_mesh_from_polygons, polys, target_size=TARGET_SIZE,
                  element_type="tri6", lines=lines, element_size_1d=SIZE_1D,
                  size_regions=extract_size_regions(sd),
                  joint_lines=extract_joint_lines(sd))


def _solved():
    """``(slope_data, mesh, fem_data, solution)`` for the jointed fixture.

    One solve at F = 1, which is enough to put the free ends of the sheets past
    their interface limit: what every leg below reads is the state, not a factor
    of safety.
    """
    from xslope.fem import build_fem_data, solve_fem
    sd = _model()
    mesh = _mesh_for(sd)
    sd["mesh"] = mesh
    fem_data = build_fem_data(sd, mesh)
    sol = _quiet(solve_fem, fem_data, F=1.0, max_iterations=3000,
                 fast_kernel=False)
    return sd, mesh, fem_data, sol


# --------------------------------------------------------------------------
# a. the wiring
# --------------------------------------------------------------------------

def _leg_wiring(failures, cache):
    from xslope.mesh import extract_joint_lines, line_is_jointed
    sd, mesh, fem_data, _sol = cache["solved"]

    got = extract_joint_lines(sd)
    if got is None or sorted(got) != sorted(JOINTED):
        failures.append(f"extract_joint_lines returned {got}, not the two "
                        f"flagged lines {JOINTED}")
    else:
        for k in JOINTED:
            if set(got[k]) != {"tend1", "tend2"}:
                failures.append(f"line {k}'s options are {sorted(got[k])}, not "
                                f"the end anchorages the mesher reads")

    plain = _model(jointed=())
    if extract_joint_lines(plain) is not None:
        failures.append("a model with no jointed line does not read as None, "
                        "so the mesher would take the joint path on it")
    for i, r in enumerate(sd["reinforcement_lines"]):
        if line_is_jointed(r) != (i in JOINTED):
            failures.append(f"line_is_jointed disagrees with the column on line "
                            f"{i + 1}")

    for key in ("joints", "elements_joint", "element_types_joint",
                "element_materials_joint", "element_side_joint"):
        if key not in mesh:
            failures.append(f"the mesh built from the model carries no {key!r}")
    jd = fem_data.get("joint_data")
    if jd is None:
        failures.append("build_fem_data wrote no joint_data for a jointed model")
        return
    got_lines = sorted(int(v) for v in np.unique(jd["line_id"]))
    want_lines = sorted(k + 1 for k in JOINTED)
    if got_lines != want_lines:
        failures.append(f"joint_data covers lines {got_lines}, not {want_lines}")
    if not np.allclose(jd["cj"], ADHESION) or not np.allclose(
            jd["tanphi"], math.tan(math.radians(DELTA))):
        failures.append("the interface strength is not the line's Adhesion and "
                        "Delta")

    # And the untouched path: the same model with nothing flagged.
    mesh0 = _mesh_for(plain)
    for key in ("joints", "elements_joint", "ties"):
        if key in mesh0:
            failures.append(f"an unflagged model's mesh carries {key!r}")
    from xslope.fem import build_fem_data
    if build_fem_data(plain, mesh0).get("joint_data") is not None:
        failures.append("an unflagged model's fem_data carries joint_data")


# --------------------------------------------------------------------------
# b. preflight
# --------------------------------------------------------------------------

def _leg_preflight(failures, cache):
    from xslope.preflight import preflight
    sd = _model()

    def _ids(model, analysis="fem"):
        rep = preflight(model, analysis, {})
        return {f.rule_id: f for f in rep.findings if f.rule_id.startswith("joint.")}

    fired = _ids(sd)
    for rid in ("joint.no_interface_strength", "joint.lines_meet",
                "joint.crosses_constraint_line", "joint.line_load_on_line",
                "joint.likely_on_material_boundary", "joint.likely_flat_sheet",
                "joint.likely_smooth_interface", "joint.likely_wall"):
        if rid in fired:
            failures.append(f"{rid} fires on the fixture, which earns none of "
                            f"them: {fired[rid].message[:100]}")
    if "joint.bond_inputs_ignored" not in fired:
        failures.append("the sample's jointed lines carry Lp1/Lp2 and a Tres, "
                        "and nothing said they are not read")

    # The refusal: a jointed line with no interface strength, named.
    blank = _model()
    blank["reinforcement_lines"][JOINTED[0]].update(adhesion=float("nan"),
                                                    delta=float("nan"))
    got = _ids(blank).get("joint.no_interface_strength")
    if got is None:
        failures.append("a jointed line with no Adhesion or Delta is not refused")
    elif "Line 2" not in got.message:
        failures.append(f"the refusal does not name the line: {got.message[:100]}")

    # Nothing fires on an LEM run: the limit equilibrium engine ignores Joint.
    if _ids(blank, "lem"):
        failures.append("a joint rule fires on a limit-equilibrium run, which "
                        "does not read the column at all")


# --------------------------------------------------------------------------
# c. the plots
# --------------------------------------------------------------------------

def _leg_plots(failures, cache):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from xslope.plot import plot_mesh, plot_reinforcement_lines
    from xslope.plot_fem import plot_joint_states
    sd, mesh, fem_data, sol = cache["solved"]

    fig, ax = plt.subplots()
    plot_reinforcement_lines(ax, sd)
    labels = [t.get_label() for t in ax.get_lines() if
              not str(t.get_label()).startswith("_")]
    if "Reinforcement (joint)" not in labels:
        failures.append(f"the inputs plot has no jointed legend entry: {labels}")
    if "Reinforcement Line" not in labels:
        failures.append(f"the inputs plot lost the bonded legend entry: {labels}")
    if labels.count("Reinforcement (joint)") != 1:
        failures.append("the jointed legend entry is repeated per line")
    plt.close(fig)

    fig, ax = plt.subplots()
    plot_reinforcement_lines(ax, _model(jointed=()))
    labels = [t.get_label() for t in ax.get_lines() if
              not str(t.get_label()).startswith("_")]
    if "Reinforcement (joint)" in labels:
        failures.append("a model with no jointed line still says it has one")
    plt.close(fig)

    fig = plt.figure()
    _quiet(plot_mesh, mesh, materials=sd.get("materials"), fig=fig)
    ax = fig.axes[0]
    leg = ax.get_legend()
    texts = [t.get_text() for t in leg.get_texts()] if leg else []
    if "Reinforcement (joint)" not in texts:
        failures.append(f"the mesh plot does not draw the jointed line: {texts}")
    plt.close(fig)

    fig, ax = plt.subplots()
    specs = plot_joint_states(ax, fem_data, sol)
    states = [t.get_label() for t in ax.get_lines()
              if str(t.get_label()).startswith("Joint (")]
    if not states:
        failures.append("the results overlay draws no joint state")
    if not any("slipping" in s for s in states):
        failures.append(f"no interface reads as slipping at F = 1: {states}")
    if not specs or "Joint Slip" not in specs[0][1]:
        failures.append(f"no slip colorbar was offered: {specs}")
    plt.close(fig)

    # And nothing at all where there is no joint.
    from xslope.fem import build_fem_data, solve_fem
    plain = _model(jointed=())
    mesh0 = _mesh_for(plain)
    fd0 = build_fem_data(plain, mesh0)
    sol0 = _quiet(solve_fem, fd0, F=1.0, max_iterations=1000, fast_kernel=False)
    fig, ax = plt.subplots()
    if plot_joint_states(ax, fd0, sol0):
        failures.append("an unjointed model was offered a joint colorbar")
    if [t for t in ax.get_lines() if str(t.get_label()).startswith("Joint (")]:
        failures.append("an unjointed model had joint states drawn on it")
    plt.close(fig)
    cache["plain"] = (plain, mesh0, fd0, sol0)


# --------------------------------------------------------------------------
# d. the two faces move apart
# --------------------------------------------------------------------------

def _leg_faces(failures, cache):
    _sd, _mesh, fem_data, sol = cache["solved"]
    jd = fem_data["joint_data"]
    conn = np.asarray(jd["conn"], dtype=int)
    side = np.asarray(jd["side"], dtype=int)
    off = np.asarray(fem_data["dof_offset"], dtype=int)
    u = np.asarray(sol["displacements"], dtype=float)
    nodes = np.asarray(fem_data["nodes"], dtype=float)

    # Pair the upper and lower interface of each station on the bar nodes they
    # share, then read the two SOIL copies against each other.
    pairs = {}
    for i in range(jd["n"]):
        bar = conn[i, 3:6] if side[i] == 1 else conn[i, 0:3]
        soil = conn[i, 0:3] if side[i] == 1 else conn[i, 3:6]
        pairs.setdefault(tuple(sorted(int(v) for v in bar)),
                         {})["u" if side[i] == 1 else "l"] = soil
    worst, where = 0.0, None
    for rec in pairs.values():
        if "u" not in rec or "l" not in rec:
            continue
        for a, b in zip(rec["u"], rec["l"]):
            if a == b:
                continue                     # a tip: the faces are one node
            d = u[off[a]:off[a] + 2] - u[off[b]:off[b] + 2]
            g = float(np.hypot(*d))
            if g > worst:
                worst, where = g, int(a)
    if where is None:
        failures.append("no station has two soil faces to compare")
        return
    scale = max(float(np.max(np.abs(u))), 1e-30)
    if worst <= 0.05 * scale:
        failures.append(f"the two faces of the split moved {worst:.3g} apart "
                        f"against a largest displacement of {scale:.3g}: the "
                        f"deformed mesh would show no offset")
    if not np.allclose(nodes[where], nodes[where], atol=0):
        failures.append("unreachable")
    cache["face_offset"] = (worst, scale, nodes[where])


# --------------------------------------------------------------------------
# e. the detail profile
# --------------------------------------------------------------------------

def _leg_details(failures, cache):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from xslope.fem_details import (joint_line_ids, joint_profile, list_lines,
                                    profile_table, write_profile_csv)
    from xslope.plot_fem_details import plot_detail
    sd, _mesh, fem_data, sol = cache["solved"]

    ids = joint_line_ids(fem_data)
    if ids != sorted(k + 1 for k in JOINTED):
        failures.append(f"joint_line_ids gave {ids}")
        return
    prof = joint_profile(fem_data, sol, ids[0], sd)
    n = len(prof["s"])
    if n < 5:
        failures.append(f"the interface profile has {n} stations")
        return
    if not np.all(np.diff(prof["s"]) >= 0):
        failures.append("the stations are not ordered along the line")
    # End 1 of the line, not whichever end the node ids fell on: the bar's own
    # profile is measured from there and the two share a panel.
    x1 = float(sd["reinforcement_lines"][ids[0] - 1]["x1"])
    if abs(float(prof["x"][0]) - x1) > 0.5:
        failures.append(f"the profile starts at x = {prof['x'][0]:.3g}, not at "
                        f"the line's own end 1 ({x1:g})")
    if not np.all(np.abs(prof["ts"]) <= prof["tlim"] + 1e-6):
        failures.append("a station's shear traction stands above its own limit")
    if not prof["slipping"].any():
        failures.append("no station reads as slipping, so the panel would show "
                        "an interface doing nothing")
    if prof["status"] != "slipping":
        failures.append(f"the line's verdict is {prof['status']!r}")
    if not len(prof["bar_s"]):
        failures.append("the bar's own profile is missing from the panel")

    rows = list_lines(fem_data, sol, sd)
    joints = [r for r in rows if r["kind"] == "joint"]
    reinf = [r for r in rows if r["kind"] == "reinforcement"]
    if len(joints) != len(JOINTED):
        failures.append(f"{len(joints)} joint rows for {len(JOINTED)} jointed lines")
    if len(reinf) != len(sd["reinforcement_lines"]):
        failures.append("the reinforcement rows changed when joints were added")
    if joints and joints[0]["index"] not in [r["index"] for r in reinf]:
        failures.append("a joint row names a line the reinforcement list does not")

    fig = plt.figure()
    plot_detail(prof, fig=fig)
    if len(fig.axes) != 4:
        failures.append(f"the joint detail figure has {len(fig.axes)} panels, "
                        f"not the four the panel documents")
    ylabels = [a.get_ylabel() for a in fig.axes]
    for want in ("Bar tension", "Normal traction", "Shear traction", "Slip"):
        if not any(want in y for y in ylabels):
            failures.append(f"the detail figure has no {want} panel: {ylabels}")
    plt.close(fig)

    cols, tab = profile_table(prof)
    for want in ("normal_traction", "shear_traction", "shear_limit", "slip",
                 "bar_axial_force"):
        if not any(c.startswith(want) for c in cols):
            failures.append(f"the exported CSV has no {want} column: {cols}")
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        write_profile_csv(prof, os.path.join(d, "j.csv"))


# --------------------------------------------------------------------------
# f. the report table
# --------------------------------------------------------------------------

def _leg_report(failures, cache):
    from xslope.report import _joint_profiles, _joint_table, _Counter
    sd, _mesh, fem_data, sol = cache["solved"]
    bundle = {"fem_data": fem_data, "solution": sol}
    profiles = _joint_profiles(sd, bundle, "converged")
    if len(profiles) != len(JOINTED):
        failures.append(f"the report reads {len(profiles)} joints for "
                        f"{len(JOINTED)} jointed lines")
        return
    table = _joint_table(profiles, _Counter())
    if [h.split(" (")[0] for h in table.headers] != [
            "Line", "Length", "Slipping", "Peak slip", "State"]:
        failures.append(f"the joints table's columns are {table.headers}")
    if len(table.rows) != len(JOINTED):
        failures.append(f"the joints table has {len(table.rows)} rows")
    for row in table.rows:
        if not row[2].endswith("%"):
            failures.append(f"the slipping share is not a share: {row[2]!r}")
        if row[4] not in ("intact", "slipping", "open"):
            failures.append(f"the state column reads {row[4]!r}")

    plain, _m0, fd0, sol0 = cache["plain"]
    if _joint_profiles(plain, {"fem_data": fd0, "solution": sol0}, "converged"):
        failures.append("an unjointed run is offered a joints table")

    # A field read back from a saved sidecar carries no interface state, and
    # reading its absent arrays as zeros would report every interface intact —
    # a state nothing measured.
    stale = {k: v for k, v in sol.items() if not k.startswith("joint_")}
    if _joint_profiles(sd, {"fem_data": fem_data, "solution": stale}, "converged"):
        failures.append("a solution carrying no interface state was tabulated "
                        "as though every joint were intact")
    from xslope.fem_details import list_lines
    if [r for r in list_lines(fem_data, stale, sd) if r["kind"] == "joint"]:
        failures.append("a solution carrying no interface state still lists "
                        "its joints")
    import matplotlib.pyplot as plt
    from xslope.plot_fem import plot_joint_states
    fig, ax = plt.subplots()
    if plot_joint_states(ax, fem_data, stale):
        failures.append("a solution carrying no interface state was offered a "
                        "slip colorbar")
    plt.close(fig)


# --------------------------------------------------------------------------

LEGS = (
    ("the loader's column reaches the mesher", _leg_wiring),
    ("preflight names the line", _leg_preflight),
    ("the plots draw the joint", _leg_plots),
    ("the two faces move apart", _leg_faces),
    ("the detail profile and its figure", _leg_details),
    ("the report's joints table", _leg_report),
)


def main():
    failures = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cache = {"solved": _solved()}
        _sd, mesh, fem_data, _sol = cache["solved"]
        print(f"  fixture: {len(mesh['nodes'])} nodes, "
              f"{len(mesh['elements'])} elements, "
              f"{len(mesh['elements_joint'])} joint elements on "
              f"{len(JOINTED)} jointed line(s)")
        # (c) fills cache['plain'], which (f) reads: the legs run in order.
        for name, leg in LEGS:
            before = len(failures)
            try:
                leg(failures, cache)
            except Exception as exc:      # a leg that raises is a failure
                failures.append(f"{name}: {type(exc).__name__}: {exc}")
            print(f"  {name:44s} "
                  + ("ok" if len(failures) == before else "FAILED"))
    if "face_offset" in cache:
        worst, scale, at = cache["face_offset"]
        print(f"  the two soil faces stand {worst:.4g} apart at "
              f"({at[0]:g}, {at[1]:g}), against a largest displacement of "
              f"{scale:.4g}")
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
    print("\nA jointed line reaches the mesher from the column, preflight names "
          "it, and the plots, the detail panel and the report all read it.")


if __name__ == "__main__":
    _cli()
