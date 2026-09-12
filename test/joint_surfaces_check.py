"""What a jointed line reaches: the mesher, preflight, the plots, the detail
panel and the report.

R1 built the mesh split and R2 the element. This is the path a MODEL takes to
them — the `Joint` column on the reinforce sheet — and what a solved jointed
model then shows.

Six legs, all on one solve of the shipped reinforcement sample with two of its
six lines made joints in memory:

  a. the wiring. ``mesh.extract_joint_options`` reads the column off a loaded
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
     the jointed lines over the mesh; the results panel draws every joint as a
     hairline colored by slip, with the slip colorbar as its only legend and no
     colorbar where nothing slipped, and draws nothing at all on a model with
     no joint. The displacement panel is the scaled deformed mesh, with the
     joint faces over its grid, on a jointed model and the arrow field on every
     other one — asserted in both layouts a results figure is drawn in, the
     stacked multi-panel one and the single panel Studio and the report render,
     which take different paths through plot_fem_results.
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
  h. the blocks. Cutting the element adjacency along the joint faces — which
     the mesh split has already done, since the two sides carry their own nodes
     — leaves the bodies the joints cut the section into. Read on two grids
     whose answer is known by inspection: a rectangle cut through is two blocks,
     and one whose joint stops inside is one, because the material wraps around
     the tip. That is what the deformed panel tints and outlines.
  g. the saved field keeps its joints. Slip and the opened record are the
     solve's own history and cannot be recovered from the displacements, so a
     field exported without them reloads as a model whose interfaces went
     quiet. The round trip is exact, and a file that is not this model's is
     refused whole.

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

#: The heaviest a joint span may be drawn (points), and the heaviest its white
#: under-stroke may be. The overlay was three-point bars over a four-and-a-half
#: point black under-stroke, which a model with hundreds of joints turns into a
#: solid mat over the field. It is thin lines now: a slipping span carries the
#: reading and takes the weight, an intact one stays out of the way, and only
#: the slipping one is backed by white so a dark green line still reads over a
#: dark blue or red patch of field. These are the bounds that keep it that way.
_SPAN_MAX = 1.6
_INTACT_MAX = 1.0
_HALO_MAX = _SPAN_MAX + 1.3


def _same_color(collection, color):
    """Whether every line of a collection is drawn in this color."""
    import matplotlib.colors as mcolors
    want = mcolors.to_rgba(color)
    got = collection.get_colors()
    return len(got) > 0 and all(tuple(c) == want for c in got)


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
                             extract_joint_options, extract_size_regions,
                             get_material_polygons)
    lines, _nr, _np = extract_constraint_line_geometry(sd)
    polys = get_material_polygons(sd, reinf_lines=lines)
    return _quiet(build_mesh_from_polygons, polys, target_size=TARGET_SIZE,
                  element_type="tri6", lines=lines, element_size_1d=SIZE_1D,
                  size_regions=extract_size_regions(sd),
                  joint_lines=extract_joint_options(sd))


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
    from xslope.mesh import extract_joint_options, line_is_jointed
    sd, mesh, fem_data, _sol = cache["solved"]

    got = extract_joint_options(sd)
    if got is None or sorted(got) != sorted(JOINTED):
        failures.append(f"extract_joint_options returned {got}, not the two "
                        f"flagged lines {JOINTED}")
    else:
        for k in JOINTED:
            # A reinforcement line flagged Joint carries its end anchorages and
            # says it HAS a bar; the joints sheet's own lines say the opposite.
            if set(got[k]) != {"tend1", "tend2", "bar"}:
                failures.append(f"line {k}'s options are {sorted(got[k])}, not "
                                f"the end anchorages and the bar flag the mesher "
                                f"reads")
            elif got[k].get("bar") is not True:
                failures.append(f"line {k} is a reinforcement line flagged Joint "
                                f"and must carry a bar; its option says "
                                f"{got[k].get('bar')!r}")

    plain = _model(jointed=())
    if extract_joint_options(plain) is not None:
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
    from xslope import plot_fem as _PF
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

    # The inputs overlay thins out as the count rises: a handful of joint lines
    # keep the fault symbol (dashes and ticks), a network gets plain traces you
    # can follow and a legend that counts them.
    from xslope.plot import JOINT_TICK_MAX_LINES, plot_joint_lines

    def _joint_only(n):
        return {"joint_lines": [{"x1": 0.0, "y1": float(k), "x2": 50.0,
                                 "y2": float(k) + 5.0} for k in range(n)]}

    for n, symbol in ((JOINT_TICK_MAX_LINES, True),
                      (JOINT_TICK_MAX_LINES + 1, False)):
        fig, ax = plt.subplots()
        plot_joint_lines(ax, _joint_only(n))
        labels = [str(t.get_label()) for t in ax.get_lines()
                  if not str(t.get_label()).startswith("_")]
        if f"Joint ({n} lines)" not in labels:
            failures.append(f"the inputs overlay does not count its lines at "
                            f"{n}: {labels}")
        ticked = [t for t in ax.get_lines()
                  if str(t.get_marker()) not in ("None", "", "none")]
        if symbol and not ticked:
            failures.append(f"{n} joint lines lost the fault symbol, which is "
                            f"readable at that count")
        if not symbol and ticked:
            failures.append(f"{n} joint lines still carry ticks, which past "
                            f"{JOINT_TICK_MAX_LINES} is all a reader sees")
        plt.close(fig)

    fig = plt.figure()
    _quiet(plot_mesh, mesh, materials=sd.get("materials"), fig=fig)
    ax = fig.axes[0]
    leg = ax.get_legend()
    texts = [t.get_text() for t in leg.get_texts()] if leg else []
    if "Reinforcement (joint)" not in texts:
        failures.append(f"the mesh plot does not draw the jointed line: {texts}")
    plt.close(fig)

    # A BAR-LESS joint over the mesh is drawn in the joint color and backed by
    # white. The element edges are black and the zone fills are the palette's,
    # so a near-black trace could not be told from an edge; the color that tells
    # them apart is the family the slip ramp ends on, and the white under-stroke
    # is what holds the trace off a dense grid. The fixture's own joints carry
    # bars, so the same mesh is read with its joint records made bar-less — the
    # station triple is [upper, bar, lower] and a bar-less line records -1 in
    # the middle of it.
    from xslope.plot import JOINT_COLOR, JOINT_LINEWIDTH
    import matplotlib.colors as _mc
    mesh_bl = dict(mesh)
    mesh_bl["joints"] = [
        dict(rec, bar=False,
             stations=[[int(st[0]), -1, int(st[-1])]
                       for st in (rec.get("stations") or [])])
        for rec in (mesh.get("joints") or [])]
    fig = plt.figure()
    _quiet(plot_mesh, mesh_bl, materials=sd.get("materials"), fig=fig)
    ax = fig.axes[0]
    want = _mc.to_rgba(JOINT_COLOR)
    # The trace runs through the line's stations; the ticks across it are the
    # two-point marks in the same color, and they keep their own style.
    traces = [ln for ln in ax.get_lines()
              if _mc.to_rgba(ln.get_color()) == want and len(ln.get_xdata()) > 2]
    if not traces:
        failures.append("the mesh panel draws no bar-less joint trace in the "
                        "joint color, so a joint is not told from a mesh edge")
    for ln in traces:
        if abs(ln.get_linewidth() - JOINT_LINEWIDTH) > 1e-9:
            failures.append(f"a mesh-panel joint trace is not the joint weight: "
                            f"{ln.get_linewidth()}")
        if not ln.get_path_effects():
            failures.append("a mesh-panel joint trace carries no white "
                            "under-stroke over the element edges")
    if want[1] <= max(want[0], want[2]):
        failures.append(f"the joint color {JOINT_COLOR} is not green-dominant, "
                        f"so it is not the family the slip ramp ends on and the "
                        f"panels no longer agree")
    plt.close(fig)

    # The results overlay: hairlines, no state legend, a colorbar only where
    # something slipped. A network of hundreds of traces is the case this is
    # drawn for, so weight and legend entries are both part of the contract.
    from matplotlib.collections import LineCollection
    fig, ax = plt.subplots()
    specs = plot_joint_states(ax, fem_data, sol)
    cols = [c for c in ax.collections if isinstance(c, LineCollection)]
    if not cols:
        failures.append("the results overlay draws no joint")
    widths = [w for c in cols for w in c.get_linewidths()]
    if widths and max(widths) > _HALO_MAX:
        failures.append(f"a joint is drawn heavier than a thin line and its "
                        f"under-stroke: {widths}")
    # The slipping spans carry the reading, so they take the weight and the
    # white backing; the intact ones stay lighter and carry none. The three
    # collections are told apart by color: white is the under-stroke, the
    # neutral gray is the intact hairline, whatever is left is the slip.
    import matplotlib.colors as mcolors
    _white = mcolors.to_rgba(_PF._JOINT_HALO_COLOR)
    _gray = mcolors.to_rgba(_PF._JOINT_INTACT_COLOR)
    halo = [c for c in cols if tuple(c.get_colors()[0]) == _white]
    faint = [c for c in cols if tuple(c.get_colors()[0]) == _gray]
    slid = [c for c in cols if c not in halo and c not in faint]
    if not halo:
        failures.append("the slipping spans carry no white under-stroke, so a "
                        "dark one over a dark patch of field has nothing to "
                        "stand off")
    if not faint:
        failures.append("the fixture's intact spans are not drawn in the "
                        "neutral gray")
    for c in faint:
        if max(c.get_linewidths()) > _INTACT_MAX:
            failures.append(f"an intact span is drawn as heavily as a slipping "
                            f"one: {c.get_linewidths()}")
    for c in slid:
        if max(c.get_linewidths()) > _SPAN_MAX:
            failures.append(f"a slipping span is heavier than the contract: "
                            f"{c.get_linewidths()}")
    for c in halo:
        if max(c.get_linewidths()) <= _SPAN_MAX:
            failures.append("the under-stroke is no wider than the line it is "
                            "meant to back")
        if slid and c.get_zorder() >= min(x.get_zorder() for x in slid):
            failures.append("the white under-stroke is drawn over its line")

    # The ramp shares no color with the field it is drawn over. The field is
    # coolwarm — blue through white to red — and the slip ramp is green: every
    # color on it has green as its strongest channel, which no coolwarm value
    # has. Measured on the ramps themselves rather than asserted of a name, so
    # swapping either for another cannot quietly lose the separation.
    ramp = _PF._joint_slip_cmap()
    if ramp.name != _PF._JOINT_SLIP_CMAP:
        failures.append(f"the slip ramp is not the one the module names: "
                        f"{ramp.name!r} vs {_PF._JOINT_SLIP_CMAP!r}")
    if specs and getattr(specs[0][0], "cmap", None) is not None:
        if specs[0][0].cmap.name != ramp.name:
            failures.append("the colorbar is drawn from a different ramp than "
                            "the spans, so the reading would not be honest")
    _ramp = ramp(np.linspace(0.0, 1.0, 64))[:, :3]
    _sep = float((_ramp[:, 1] - np.maximum(_ramp[:, 0], _ramp[:, 2])).min())
    if _sep <= 0.1:
        failures.append(f"the slip ramp is not green-dominant throughout "
                        f"(worst margin {_sep:.3f}), so it can be read as part "
                        f"of the coolwarm field under it")
    _field = plt.get_cmap('coolwarm')(np.linspace(0.0, 1.0, 256))[:, :3]
    if float((_field[:, 1] - np.maximum(_field[:, 0], _field[:, 2])).max()) > 0.0:
        failures.append("the field ramp is green-dominant somewhere, so green "
                        "no longer separates the joints from it")
    labels = [str(t.get_label()) for t in ax.get_lines()]
    if any(s.startswith("Joint (") for s in labels):
        failures.append(f"the overlay still carries state legend entries: "
                        f"{[s for s in labels if s.startswith('Joint (')]}")
    if not specs or "Joint slip" not in specs[0][1]:
        failures.append(f"no slip colorbar was offered: {specs}")
    plt.close(fig)

    # A jointed model where nothing slipped is drawn, but offered no colorbar:
    # an empty ramp would say a quantity was measured that was not.
    stuck = dict(sol)
    stuck["joint_slipping"] = np.zeros_like(np.asarray(sol["joint_slipping"]))
    stuck["joint_slip"] = np.zeros_like(np.asarray(sol["joint_slip"]))
    fig, ax = plt.subplots()
    if plot_joint_states(ax, fem_data, stuck):
        failures.append("a model with no slipping joint was offered a colorbar")
    if not [c for c in ax.collections if isinstance(c, LineCollection)]:
        failures.append("a model with no slipping joint drew no joint at all")
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
    if [c for c in ax.collections if isinstance(c, LineCollection)]:
        failures.append("an unjointed model had joint states drawn on it")
    plt.close(fig)
    cache["plain"] = (plain, mesh0, fd0, sol0)

    # The displacement panel. A jointed model's mechanism is blocks moving as
    # bodies on their joints, which an arrow field sampled at nodes does not
    # show, so that panel becomes the scaled deformed mesh with the joint faces
    # drawn. An unjointed model keeps the arrows.
    #
    # The rule belongs to the PANEL, not to a figure, so it is asserted on every
    # layout a results figure is drawn in: the stacked multi-panel figure the
    # driver scripts and the docs render, AND the one-panel-at-a-time figure
    # Studio's results view and the report render. Those take different layout
    # branches through plot_fem_results (deferred stacked colorbars vs. the
    # single-panel make_axes_locatable path), so one passing does not prove the
    # other.
    from matplotlib.quiver import Quiver
    from xslope.plot_fem import plot_fem_results
    layouts = (["shear_strain", "displace_vector"],   # stacked figure
               ["displace_vector"])                   # Studio / report: one panel
    for fd, sol, jointed in ((fem_data, sol, True), (fd0, sol0, False)):
        for panels in layouts:
            where = f"{len(panels)}-panel"
            fig, axes = _quiet(plot_fem_results, fd, sol, plot_type=list(panels),
                               figsize=(9, 7))
            # One panel returns the Axes itself, not a list of them.
            ax_d = axes if len(panels) == 1 else axes[-1]
            title = ax_d.get_title()
            arrows = [a for a in ax_d.collections if isinstance(a, Quiver)]
            faces = [c for c in ax_d.collections if isinstance(c, LineCollection)]
            if jointed and arrows:
                failures.append(f"a jointed model's {where} displacement panel "
                                f"still draws an arrow field")
            if jointed and "Deformation" not in title:
                failures.append(f"a jointed model's {where} displacement panel "
                                f"is not the deformed mesh: {title!r}")
            if jointed and "Scale" not in title:
                failures.append(f"the {where} deformed-mesh panel does not print "
                                f"its exaggeration: {title!r}")
            # The joint faces are the point of the substitution: a deformed grid
            # with no faces on it says nothing the arrows didn't.
            if jointed and len(faces) < 2:
                failures.append(f"the {where} deformed-mesh panel drew no joint "
                                f"faces over its grid: {len(faces)} collections")
            if not jointed and not arrows:
                failures.append(f"an unjointed model lost its {where} "
                                f"displacement vectors")
            if not jointed and "Displacement Vectors" not in title:
                failures.append(f"an unjointed model's {where} displacement "
                                f"panel changed: {title!r}")
            # The blocks: a faint tint per body, and the outside of the deformed
            # mesh as a line of its own. Without the tint two blocks that touch
            # are one gray field; without the boundary the only thing carrying
            # the deformed ground surface is the light element grid.
            from matplotlib.collections import PolyCollection
            # A Quiver is itself a PolyCollection, so the arrow field is not a
            # block tint however much it looks like one to isinstance.
            tints = [c for c in ax_d.collections
                     if isinstance(c, PolyCollection)
                     and not isinstance(c, Quiver)]
            edges = [c for c in ax_d.collections if isinstance(c, LineCollection)
                     and _same_color(c, _PF._DEFORMED_BOUNDARY_COLOR)]
            if jointed and not tints:
                failures.append(f"the {where} deformed panel tints no blocks, "
                                f"so two bodies that touch read as one")
            if jointed and not edges:
                failures.append(f"the {where} deformed panel draws no exterior "
                                f"boundary, so the moved ground surface is "
                                f"carried only by the element grid")
            if not jointed and tints:
                failures.append(f"an unjointed model's {where} panel was "
                                f"tinted by block")
            plt.close(fig)


    # A NETWORK's deformed panel draws no element edges at all. Past the count
    # at which the section drawings drop the joint ticks, the grid is all a
    # reader can see and the block outlines — the reading — are buried in it.
    # Measured on the fixture with its joints renumbered one per station, which
    # is what a generated network looks like to the drawing.
    from xslope.plot import JOINT_TICK_MAX_LINES as _JMAX
    from xslope.plot_fem import plot_deformed_mesh
    # The layout loop above rebound sol to the unjointed model's field; the
    # jointed pair is the fixture's own.
    _sd, _mesh, fem_data, sol = cache["solved"]
    n_j = int(fem_data["joint_data"]["n"])
    net = dict(fem_data)
    net["joint_data"] = dict(fem_data["joint_data"],
                             line_id=np.arange(n_j, dtype=int))
    for fd, dense in ((fem_data, False), (net, True)):
        n_lines = len(np.unique(np.asarray(fd["joint_data"]["line_id"])))
        if (n_lines > _JMAX) != dense:
            failures.append(f"the fixture for the {'dense' if dense else 'few'} "
                            f"case does not straddle the threshold: {n_lines}")
        fig, ax = plt.subplots()
        _quiet(plot_deformed_mesh, ax, fd, sol, 1000.0, joint_faces=True)
        grid = [c for c in ax.collections if isinstance(c, LineCollection)
                and _same_color(c, _PF._DEFORMED_GRID_UNDER_JOINTS)]
        if dense and grid:
            failures.append(f"a {n_lines}-line network's deformed panel still "
                            f"draws its element edges over the blocks")
        if not dense and not grid:
            failures.append(f"a {n_lines}-line model's deformed panel lost the "
                            f"light element edges, which it is still readable "
                            f"with")
        plt.close(fig)

    # The exaggeration is bounded, and a field with nothing in it says so rather
    # than being magnified into a shape. Both read off the scale the panel is
    # actually drawn at.
    from xslope.plot_fem import (deformation_below_resolution, deformation_scale,
                                 displacement_magnitude)
    # The cap binds where a section is no wider than it is tall — the height
    # rule alone would then swing the crest much further than the ceiling — so
    # it is read on the fixture AND on a copy squeezed to that shape, where the
    # two bounds disagree and the smaller has to win.
    narrow = dict(fem_data)
    narrow["nodes"] = (np.asarray(fem_data["nodes"], dtype=float)
                       * np.array([0.1, 1.0]))
    biggest = float(np.max(displacement_magnitude(fem_data, sol)))
    for tag, fd in (("the fixture", fem_data), ("a narrow section", narrow)):
        extent = _PF._section_extent(fd)
        height = float(np.ptp(np.asarray(fd["nodes"], dtype=float)[:, 1]))
        moved = deformation_scale(fd, sol) * biggest
        if moved > _PF._DEFORM_MAX_FRACTION * extent * 1.001:
            failures.append(f"on {tag} the deformed mesh is drawn with a point "
                            f"moved {moved / extent:.3f} of the section, past "
                            f"the {_PF._DEFORM_MAX_FRACTION} cap")
        if tag == "a narrow section" and \
                height * 0.15 <= _PF._DEFORM_MAX_FRACTION * extent:
            failures.append("the narrow copy does not make the two bounds "
                            "disagree, so the cap is never exercised")
    quiet_sol = dict(sol)
    tiny = _PF._DEFORM_NEGLIGIBLE_FRACTION * extent * 1e-2
    factor = tiny / float(np.max(displacement_magnitude(fem_data, sol)))
    for key in ("displacements", "displacements_elastic"):
        if quiet_sol.get(key) is not None:
            quiet_sol[key] = np.asarray(quiet_sol[key], dtype=float) * factor
    if not deformation_below_resolution(fem_data, quiet_sol):
        failures.append("a field whose largest displacement is a hundredth of "
                        "the negligible threshold is still called drawable")
    if deformation_scale(fem_data, quiet_sol) != 1.0:
        failures.append("a below-resolution field is still exaggerated")
    fig, ax = plt.subplots()
    _quiet(plot_deformed_mesh, ax, fem_data, quiet_sol,
           deformation_scale(fem_data, quiet_sol), joint_faces=True)
    if "below drawing resolution" not in ax.get_title():
        failures.append(f"an undrawable deformation does not say so: "
                        f"{ax.get_title()!r}")
    if "Scale" in ax.get_title():
        failures.append(f"an undrawable deformation still prints an "
                        f"exaggeration: {ax.get_title()!r}")
    plt.close(fig)
    if deformation_below_resolution(fem_data, sol):
        failures.append("the fixture's own mechanism reads as below drawing "
                        "resolution, so the threshold is set over real motion")


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
    if prof["status"] not in ("slipping", "open"):
        failures.append(f"the line's verdict is {prof['status']!r}")
    # A station may only be OPEN where the line reaches a free face. The sample's
    # sheets start ON the slope face, so the split gives that end two soil faces
    # and the crack can open there; the other end is buried, its two faces are
    # one node, and a tip pair never opens (see xslope/joint.py, "Tips").
    from shapely.geometry import Point as _Pt
    ring = sd["domain_polygon"].exterior
    for k in np.flatnonzero(prof["open"]):
        at_end = k in (0, n - 1)
        on_face = ring.distance(_Pt(float(prof["x"][k]),
                                    float(prof["y"][k]))) <= 1e-6
        if not (at_end and on_face):
            failures.append(
                f"station {k} at ({prof['x'][k]:.3g}, {prof['y'][k]:.3g}) reads "
                f"as open and is not an end of the line standing on the model's "
                f"external boundary")
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
# g. the saved field keeps its joints
# --------------------------------------------------------------------------

def _leg_sidecar(failures, cache):
    """A solved jointed field, exported and read back, still knows what its
    joints did — otherwise every figure of a reloaded run is drawn on a model
    whose interfaces have silently gone quiet, and a corpus figure can only be
    re-rendered by solving again.
    """
    import tempfile
    from xslope.fem import export_fem_solution, import_fem_solution
    from xslope.plot_fem import solution_has_joint_state
    _sd, _mesh, fem_data, sol = cache["solved"]

    with tempfile.TemporaryDirectory() as d:
        stem = os.path.join(d, "rt")
        _quiet(export_fem_solution, fem_data, sol, stem, meta={"FS": 1.0})
        joints_csv = os.path.join(d, "rt_fem_joints.csv")
        if not os.path.exists(joints_csv):
            failures.append("a solved jointed field writes no joint sidecar, "
                            "so its state cannot be read back")
            return
        back = _quiet(import_fem_solution, fem_data, stem)
        if not solution_has_joint_state(fem_data, back):
            failures.append("a reloaded jointed field carries no joint state")
        for key in ("joint_slip", "joint_open", "joint_slipping",
                    "joint_tn", "joint_ts", "joint_tlim"):
            a = np.asarray(sol[key]).astype(float)
            b = np.asarray(back.get(key, [])).astype(float)
            if a.shape != b.shape or not np.allclose(a, b):
                failures.append(f"{key} does not survive the round trip")

        # And a file that is not this model's is refused whole rather than
        # grafted: half a set of restored pairs reads as a solved interface.
        with open(joints_csv) as f:
            lines = f.readlines()
        with open(joints_csv, "w") as f:
            f.writelines(lines[:-1])
        hurt = _quiet(import_fem_solution, fem_data, stem)
        if not hurt.get("sidecar_notes"):
            failures.append("a joint sidecar that is not this model's is "
                            "grafted on without a word")
        if solution_has_joint_state(fem_data, hurt):
            failures.append("a refused joint sidecar still left joint state "
                            "on the solution")


# --------------------------------------------------------------------------
# h. the blocks a joint cuts
# --------------------------------------------------------------------------

def _grid_fixture(split_through):
    """A 2x2 grid of quads, cut along its middle row by one joint.

    ``split_through`` runs the joint the full width — the upper row gets its own
    copies of all three nodes on the line, and nothing connects the halves. Where
    it does not, the joint stops at the middle node: only the left node is
    copied, the middle one is the tip and is shared, and the material wraps
    around it. Nine nodes, then the copies; no mesher and no solve, so the
    component rule is read on a shape whose answer is known by inspection.
    """
    xy = [(x, y) for y in (0.0, 1.0, 2.0) for x in (0.0, 1.0, 2.0)]
    if split_through:
        xy += [(0.0, 1.0), (1.0, 1.0), (2.0, 1.0)]          # 9, 10, 11
        elements = [[0, 1, 4, 3], [1, 2, 5, 4],
                    [9, 10, 7, 6], [10, 11, 8, 7]]
        conn = [[3, 4, -1, 9, 10, -1], [4, 5, -1, 10, 11, -1]]
    else:
        xy += [(0.0, 1.0)]                                   # 9 only
        elements = [[0, 1, 4, 3], [1, 2, 5, 4],
                    [9, 4, 7, 6], [4, 5, 8, 7]]
        conn = [[3, 4, -1, 9, 4, -1]]
    return {"nodes": np.array(xy, dtype=float),
            "elements": np.array(elements, dtype=int),
            "element_types": np.array([4, 4, 4, 4], dtype=int),
            "joint_data": {"n": len(conn), "conn": np.array(conn, dtype=int),
                           "line_id": np.zeros(len(conn), dtype=int)}}


def _leg_blocks(failures, cache):
    from xslope.mesh import block_boundary_edges, block_components

    through = _grid_fixture(True)
    comp = block_components(through)
    if len(set(comp.tolist())) != 2:
        failures.append(f"a rectangle cut in two by a through-going joint came "
                        f"back as {len(set(comp.tolist()))} block(s)")
    elif set(comp[:2]) == set(comp[2:]):
        failures.append("the two halves of a through-cut rectangle are in the "
                        "same block")
    blocks, exterior = block_boundary_edges(through, comp)
    if sorted(len(v) for v in blocks.values()) != [2, 2]:
        failures.append(f"the two blocks do not carry two joint faces each: "
                        f"{ {k: len(v) for k, v in blocks.items()} }")
    if any(sorted(e[:2]) in ([3, 4], [4, 5]) and len(e) > 2 for e in exterior):
        failures.append("a joint face was counted as the outside of the mesh")
    if len(exterior) != 8:
        failures.append(f"a 2x2 grid has eight edges on its outside, not "
                        f"{len(exterior)}")

    inside = _grid_fixture(False)
    comp2 = block_components(inside)
    if len(set(comp2.tolist())) != 1:
        failures.append(f"a joint that stops inside the mass cut the section "
                        f"into {len(set(comp2.tolist()))} blocks; the material "
                        f"wraps around its tip and nothing can leave")
    blocks2, _ext2 = block_boundary_edges(inside, comp2)
    if sum(len(v) for v in blocks2.values()) != 2:
        failures.append(f"the joint that stops inside lost its faces: "
                        f"{ {k: len(v) for k, v in blocks2.items()} }")


# --------------------------------------------------------------------------

LEGS = (
    ("the loader's column reaches the mesher", _leg_wiring),
    ("preflight names the line", _leg_preflight),
    ("the plots draw the joint", _leg_plots),
    ("the two faces move apart", _leg_faces),
    ("the saved field keeps its joints", _leg_sidecar),
    ("the detail profile and its figure", _leg_details),
    ("the report's joints table", _leg_report),
    ("the blocks a joint cuts", _leg_blocks),
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
