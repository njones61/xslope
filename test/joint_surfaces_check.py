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
     thin line colored by slip with no white backing, with the slip colorbar
     and a key naming the three states placed off the section, no colorbar
     where nothing slipped, and nothing at all on a model with no joint. With
     Show joints on the strain panel of a jointed model draws EVERY joint,
     a jointed reinforcement sheet as its two faces either side of the bar,
     whether its soil can yield (over the strain field, strain title kept) or
     not (the Joint slip panel). The deformation panel draws the joints the
     same way except that a closed joint keeps the joint-face weight, and it
     returns the slip colorbar too. The displacement panel is the scaled deformed mesh, with the
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
     the tip. That is what the deformed panel tints and outlines. With Color
     by block on (``color_blocks=True``) the shipped toppling model's panel,
     drawn from its stored solution, fills more than one tint, neighboring
     blocks differing; off, it fills exactly one, its single material's.
  g. the saved field keeps its joints. Slip and the opened record are the
     solve's own history and cannot be recovered from the displacements, so a
     field exported without them reloads as a model whose interfaces went
     quiet. The round trip is exact, and a file that is not this model's is
     refused whole.
  i. the Run FEM dialog. A model with a joint in it opens the dialog on the
     Hybrid failure criterion and every other model opens it on
     Non-convergence, because a near-critical trial on a jointed model settles
     into slip the other criterion can only read as failure. A criterion the
     user has already chosen this session still wins.
  j. a bar-less joint trace is not a bar. A joints-sheet line has 1D elements
     of its own with no member between its faces, and both readers of
     ``elements_1d`` have to know it: the mesh plot, which drew the trace red
     and counted it as reinforcement, and the reinforcement sidecar, which
     wrote a per-bar force row for it with every force and every capacity zero.

Run directly:  PYTHONPATH=. python3 test/joint_surfaces_check.py
"""

import contextlib
import io
import math
import os
import sys
import warnings

import numpy as np

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

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

#: The heaviest a joint span may be drawn over a strain field (points). The
#: overlay was three-point bars over a four-and-a-half point black under-stroke,
#: which a model with hundreds of joints turns into a solid mat over the field.
#: It is thin lines now, with no white backing: a slipping span carries the
#: reading and takes the weight, a closed one stays out of the way. These are
#: the bounds that keep it that way on the strain panel; the block drawing's
#: closed joints are the joint-face weight instead (asserted below).
_SPAN_MAX = 1.6
_INTACT_MAX = 1.0


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
    got = fired.get("joint.bond_inputs_ignored")
    if got is None:
        failures.append("the sample's jointed lines carry Lp1/Lp2, and nothing "
                        "said they are not read")
    elif "fills Lp1 and Lp2," not in got.message:
        # Tres IS read on a jointed line (the bar softens to it), so the INFO
        # names the development lengths and nothing else as unread.
        failures.append(f"the INFO does not name exactly Lp1 and Lp2 as the "
                        f"unread inputs: {got.message[:120]}")

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

    # The results overlay: thin lines, no per-state legend entries on the
    # axes, a colorbar only where something slipped. A network of hundreds of
    # traces is the case this is drawn for, so weight and legend entries are
    # both part of the contract. Read on the fixture with nothing opened, so
    # every white collection would be a backing: white is kept for the gap
    # down an OPENED stretch, and nothing else carries it.
    from matplotlib.collections import LineCollection
    import matplotlib.colors as mcolors
    _white = mcolors.to_rgba("white")
    _gray = mcolors.to_rgba(_PF._JOINT_INTACT_COLOR)
    no_open = dict(sol)
    no_open["joint_open"] = np.zeros_like(np.asarray(sol["joint_open"]))

    def _by_color(ax):
        cols_ = [c for c in ax.collections if isinstance(c, LineCollection)
                 and 6.44 < c.get_zorder() < 6.61]
        white_ = [c for c in cols_ if len(c.get_colors())
                  and tuple(c.get_colors()[0]) == _white]
        gray_ = [c for c in cols_ if len(c.get_colors())
                 and tuple(c.get_colors()[0]) == _gray]
        return cols_, white_, gray_, [c for c in cols_
                                      if c not in white_ and c not in gray_]

    fig, ax = plt.subplots()
    specs = plot_joint_states(ax, fem_data, no_open)
    cols, white, faint, slid = _by_color(ax)
    if not cols:
        failures.append("the results overlay draws no joint")
    if white:
        failures.append(f"the strain overlay draws {len(white)} white "
                        f"collection(s) with nothing opened: a white backing "
                        f"under the slipping spans")
    if not slid:
        failures.append("the fixture's slipping spans are not drawn")
    if not faint:
        failures.append("the fixture's closed spans are not drawn in the "
                        "neutral gray")
    for c in faint:
        if max(c.get_linewidths()) > _INTACT_MAX:
            failures.append(f"a closed span on the strain panel is not the "
                            f"hairline: {c.get_linewidths()}")
    for c in slid:
        if max(c.get_linewidths()) > _SPAN_MAX:
            failures.append(f"a slipping span is heavier than the contract: "
                            f"{c.get_linewidths()}")
    plt.close(fig)

    # The block drawing: the same styling (no white backing under a slipping
    # span), except that a closed joint keeps the joint-face weight, because
    # the thick gray closed joints are what make the blocks read as blocks;
    # and it hands back the slip colorbar's spec for its caller to place.
    from xslope.plot_fem import JOINT_FACE_LINEWIDTH, plot_deformed_mesh as _pdm0
    fig, ax = plt.subplots()
    d_specs = _quiet(_pdm0, ax, fem_data, no_open, 1000.0, joint_faces=True,
                     single_panel=True) or []
    _cols, d_white, d_gray, d_slid = _by_color(ax)
    if not any("Joint slip" in lab for _sm, lab in d_specs):
        failures.append(f"the deformation panel returns no slip colorbar "
                        f"spec: {d_specs!r}")
    if d_white:
        failures.append(f"the block drawing draws {len(d_white)} white "
                        f"collection(s) with nothing opened: a white backing "
                        f"under the slipping spans")
    if not d_gray:
        failures.append("the block drawing draws no closed joint in gray")
    for c in d_gray:
        if abs(min(c.get_linewidths()) - JOINT_FACE_LINEWIDTH) > 1e-9:
            failures.append(f"a closed joint on the block drawing is not the "
                            f"joint-face weight {JOINT_FACE_LINEWIDTH}: "
                            f"{c.get_linewidths()}")
    for c in faint:
        if max(c.get_linewidths()) >= JOINT_FACE_LINEWIDTH:
            failures.append(f"a closed span on the strain panel is not "
                            f"thinner than on the block drawing: "
                            f"{c.get_linewidths()}")
    plt.close(fig)

    # The key is placed off the section: on each panel of the results figure
    # its box in data coordinates does not touch the ground the panel draws
    # (the shape the placer tests; a corner clear of the ground, or below the
    # axes when no corner is).
    from shapely.geometry import box as _sbox
    from xslope.plot_fem import plot_fem_results as _pfr
    fig, axes_k = _quiet(_pfr, fem_data, no_open,
                         plot_type=["deformation", "shear_strain"],
                         figsize=(9, 7))
    fig.canvas.draw()
    rr = fig.canvas.get_renderer()
    for ax_k in axes_k:
        leg = ax_k.get_legend()
        texts = [] if leg is None else [t.get_text() for t in leg.get_texts()]
        if "closed, no slip" not in texts:
            failures.append(f"the {ax_k.get_title()[:30]!r} panel draws no "
                            f"joint key: {texts!r}")
            continue
        shape = _PF._section_shape(ax_k)
        bb = leg.get_window_extent(rr)
        inv = ax_k.transData.inverted()
        (x0, y0), (x1, y1) = inv.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
        if shape is None or shape.intersects(
                _sbox(min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))):
            failures.append(f"the joint key on the {ax_k.get_title()[:30]!r} "
                            f"panel lies over the section")
    plt.close(fig)
    fig, ax = plt.subplots()
    specs = plot_joint_states(ax, fem_data, sol)

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

    # The panel the joints get when nothing else can happen. A jointed model
    # whose every material is linear elastic accumulates no viscoplastic strain
    # anywhere, so the strain panel would be a flat field with a scale invented
    # for it; that panel is the slip's instead, and it says so in its title. The
    # rule reads the MATERIALS, so it is asserted by flipping them on the same
    # solved model rather than by planting a small field.
    from xslope.plot_fem import (JOINT_SLIP_PANEL_LABEL, SHEAR_STRAIN_LABEL,
                                 plot_shear_strain_contours)

    fig, ax = plt.subplots()
    plot_shear_strain_contours(ax, fem_data, sol, single_panel=True)
    if SHEAR_STRAIN_LABEL not in ax.get_title():
        failures.append(f"a model that can yield lost the strain panel: "
                        f"{ax.get_title()!r}")
    plt.close(fig)

    fd_el = dict(fem_data)
    fd_el["elastic_materials"] = list(fem_data.get("material_names") or [])
    if not _PF._joint_slip_panel(fd_el):
        failures.append("an all-elastic jointed model is not read as a slip panel")
    if _PF._joint_slip_panel(fem_data):
        failures.append("a model that can yield is read as a slip panel")
    fig, ax = plt.subplots()
    mappable, _ = plot_shear_strain_contours(ax, fd_el, sol, single_panel=True)
    if JOINT_SLIP_PANEL_LABEL not in ax.get_title():
        failures.append(f"an all-elastic jointed model still names the strain "
                        f"field: {ax.get_title()!r}")
    if mappable is not None:
        failures.append("the slip panel contoured a field and offered its scale, "
                        "where the field it would scale is zero by construction")
    plt.close(fig)

    # Every joint draws on the strain panel with Show joints on: a jointed
    # reinforcement SHEET as its two faces either side of the bar, a bar-less
    # joint (wall contact, rock joint) as its line, each with the slip
    # colorbar and the key. The fixture's joints are all sheets, so a MIXED copy is made by
    # declaring one jointed line bar-less (``barless_1d_mask`` over its 1D
    # elements, which is what a joints-sheet line produces); the line chosen
    # is one with slipping spans, so its slip colorbar has something to scale.
    _states = ("closed, no slip", "slipping (color = slip)", "opened")
    n_1d = len(fem_data["elements_1d"])
    line_of = np.asarray(fem_data["element_materials_1d"], dtype=int)
    spans_all = _PF._joint_spans(fem_data, sol)
    slipping_lines = sorted({r["line"] for r in spans_all if r["slipping"]})
    joint_lines = sorted({r["line"] for r in spans_all})
    fd_mix = None
    if len(joint_lines) < 2 or not slipping_lines:
        failures.append(f"the fixture cannot make a mixed model: joint lines "
                        f"{joint_lines}, slipping {slipping_lines}")
    else:
        fd_mix = dict(fem_data)
        fd_mix["barless_1d_mask"] = line_of == int(slipping_lines[0])
        if not any(_PF._spans_on_bars(fd_mix, spans_all)):
            failures.append("the mixed fixture has no sheet left on a bar")
    fd_mix_el = None if fd_mix is None else dict(
        fd_mix, elastic_materials=list(fem_data.get("material_names") or []))

    def _overlay(ax):
        """The joint overlay's own collections: plot_joint_states draws them
        between z 6.45 and 6.6, and the strain panel draws nothing else there
        when its bars are left off."""
        return [c for c in ax.collections if isinstance(c, LineCollection)
                and 6.44 < c.get_zorder() < 6.61]

    def _strain(fd, **kw):
        fig, ax = plt.subplots()
        _m, specs_ = plot_shear_strain_contours(ax, fd, sol, single_panel=True,
                                                **kw)
        leg = ax.get_legend()
        key = [] if leg is None else [t.get_text() for t in leg.get_texts()]
        out = {"title": ax.get_title(), "mappable": _m,
               "key": [k for k in key if k in _states],
               "slip_bar": [lab for _sm, lab in specs_ if "Joint slip" in lab],
               "overlay": len(_overlay(ax)),
               "artists": len(ax.collections) + len(ax.lines)}
        plt.close(fig)
        return out

    # Sheets only: the bar AND the sheet faces, with the key and the slip
    # colorbar; with Show joints off, the bar alone.
    for fd, tag in ((fem_data, "yielding"), (fd_el, "all-elastic")):
        with_bar = _strain(fd, show_joints=True)
        no_bar = _strain(fd, show_joints=True, show_reinforcement=False)
        off = _strain(fd, show_joints=False, show_reinforcement=False)
        if not no_bar["overlay"]:
            failures.append(f"the {tag} sheet-only strain panel draws no sheet "
                            f"faces")
        if off["overlay"]:
            failures.append(f"with Show joints off the {tag} sheet-only strain "
                            f"panel still draws the sheet faces")
        if with_bar["artists"] <= no_bar["artists"]:
            failures.append(f"the {tag} sheet-only strain panel does not draw "
                            f"the sheets' bars")
        if not with_bar["key"]:
            failures.append(f"the {tag} sheet-only strain panel draws no joint "
                            f"key")
        if not with_bar["slip_bar"]:
            failures.append(f"the {tag} sheet-only strain panel offers no slip "
                            f"colorbar")
    yld = _strain(fem_data, show_joints=True)
    if SHEAR_STRAIN_LABEL not in yld["title"] or yld["mappable"] is None:
        failures.append(f"a yielding sheet-only model lost its strain title or "
                        f"colorbar: {yld['title']!r}")

    # Both kinds: the bar-less joint and the sheet both draw, with the key and
    # the slip colorbar. With Show joints off, none of it.
    if fd_mix is not None:
        from xslope.plot_fem import plot_joint_states as _pjs
        for show in (True, False):
            on = "on" if show else "off"
            m = _strain(fd_mix, show_joints=show)
            nb = _strain(fd_mix, show_joints=show, show_reinforcement=False)
            if SHEAR_STRAIN_LABEL not in m["title"] or m["mappable"] is None:
                failures.append(f"a yielding mixed model with Show joints {on} "
                                f"lost the strain title or colorbar")
            if show and not (m["key"] and m["slip_bar"] and nb["overlay"]):
                failures.append(f"the mixed strain panel does not draw its "
                                f"bar-less joint with the key and the slip "
                                f"colorbar: key {m['key']!r}, "
                                f"bars {m['slip_bar']!r}, "
                                f"overlay {nb['overlay']}")
            if not show and (m["key"] or m["slip_bar"] or nb["overlay"]):
                failures.append(f"with Show joints off the mixed strain panel "
                                f"still draws joints: {m!r}")
        # The sheet's faces are offset pairs either side of its bar, so the
        # strain overlay holds more segments than the bar-less line has spans:
        # every joint is drawn, the sheet's included.
        fig, ax = plt.subplots()
        _pjs(ax, fd_mix, sol)
        strain_segs = sum(len(c.get_segments()) for c in _overlay(ax))
        plt.close(fig)
        barless_spans = sum(1 for r, ob in zip(
            spans_all, _PF._spans_on_bars(fd_mix, spans_all)) if not ob)
        sheet_spans = len(spans_all) - barless_spans
        if strain_segs < barless_spans + 2 * sheet_spans:
            failures.append(f"the strain overlay leaves joints out: "
                            f"{strain_segs} segments for {barless_spans} "
                            f"bar-less spans and {sheet_spans} sheet spans "
                            f"(two faces each)")

    # The deformation panel draws the sheet's faces too.
    from xslope.plot_fem import plot_deformed_mesh as _pdm
    counts = []
    for faces_on in (True, False):
        fig, ax = plt.subplots()
        _quiet(_pdm, ax, fem_data, sol, 1000.0, joint_faces=faces_on)
        counts.append(len([c for c in ax.collections
                           if isinstance(c, LineCollection)]))
        plt.close(fig)
    if counts[0] <= counts[1]:
        failures.append(f"the deformation panel of a sheet-only model no longer "
                        f"draws the sheet faces: {counts[0]} collections with "
                        f"faces, {counts[1]} without")

    # An opened stretch is the one joint state with no colorbar to explain it.
    # It is drawn as its two faces apart, and the joint overlay's own key
    # (closed / slipping / opened, drawn by plot_joint_states off the section)
    # names it — on a figure that draws an opened stretch and on no other. The
    # overlay is drawn on the strain panel for every joint of a jointed model
    # with Show joints on, so the key is asserted on the
    # MIXED model's all-elastic slip panel and yielding strain panel, as the
    # results figure draws them (bars included), and on the overlay by itself.
    opened = dict(sol)
    opened["joint_open"] = np.ones_like(np.asarray(sol["joint_open"]))
    shut = dict(sol)
    shut["joint_open"] = np.zeros_like(np.asarray(sol["joint_open"]))

    def _key(ax):
        leg = ax.get_legend()
        return [] if leg is None else [t.get_text() for t in leg.get_texts()]

    cases_open = [("the joint overlay",
                   lambda ax, s_: plot_joint_states(ax, fem_data, s_))]
    if fd_mix is not None:
        cases_open += [
            ("the slip panel",
             lambda ax, s_: plot_shear_strain_contours(ax, fd_mix_el, s_,
                                                       single_panel=True)),
            ("the strain panel",
             lambda ax, s_: plot_shear_strain_contours(ax, fd_mix, s_,
                                                       single_panel=True))]
    for where, draw in cases_open:
        fig, ax = plt.subplots()
        draw(ax, opened)
        if "opened" not in _key(ax):
            failures.append(f"{where}: the opened stretch is drawn with nothing "
                            f"saying what it is: key {_key(ax)!r}")
        plt.close(fig)
        fig, ax = plt.subplots()
        draw(ax, shut)
        if "opened" in _key(ax):
            failures.append(f"{where}: nothing opened, yet the key names the "
                            f"opened state: key {_key(ax)!r}")
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

    # The block picture. A jointed model's mechanism is blocks moving as bodies
    # on their joints, so with Show joints on (the default) its deformation panel
    # draws the scaled deformed mesh as blocks with the joint faces over them. The
    # displacement-vector panel is the arrow field on every model, jointed or not.
    #
    # Each rule belongs to a PANEL, not to a figure, so it is asserted on every
    # layout a results figure is drawn in: the stacked multi-panel figure the
    # driver scripts and the docs render, AND the one-panel-at-a-time figure
    # Studio's results view and the report render. Those take different layout
    # branches through plot_fem_results (deferred stacked colorbars vs. the
    # single-panel make_axes_locatable path), so one passing does not prove the
    # other.
    from matplotlib.collections import PolyCollection
    from matplotlib.quiver import Quiver
    from xslope.plot_fem import plot_fem_results
    cases = ((fem_data, sol, True, "deformation"),
             (fem_data, sol, True, "displace_vector"),
             (fd0, sol0, False, "displace_vector"))
    for fd, sol, jointed, pt in cases:
        for panels in (["shear_strain", pt],   # stacked figure
                       [pt]):                  # Studio / report: one panel
            where = f"{len(panels)}-panel {pt}"
            fig, axes = _quiet(plot_fem_results, fd, sol, plot_type=list(panels),
                               figsize=(9, 7))
            # One panel returns the Axes itself, not a list of them.
            ax_d = axes if len(panels) == 1 else axes[-1]
            title = ax_d.get_title()
            arrows = [a for a in ax_d.collections if isinstance(a, Quiver)]
            faces = [c for c in ax_d.collections if isinstance(c, LineCollection)]
            # A Quiver is itself a PolyCollection, so the arrow field is not a
            # block tint however much it looks like one to isinstance.
            tints = [c for c in ax_d.collections
                     if isinstance(c, PolyCollection)
                     and not isinstance(c, Quiver)]
            edges = [c for c in ax_d.collections if isinstance(c, LineCollection)
                     and _same_color(c, _PF._DEFORMED_BOUNDARY_COLOR)]
            if pt == "displace_vector":
                # The arrow field on every model: a jointed model's third panel
                # is not the block picture a second time.
                if not arrows:
                    failures.append(f"{'a jointed' if jointed else 'an unjointed'}"
                                    f" model's {where} panel draws no "
                                    f"displacement vectors")
                if "Displacement Vectors" not in title:
                    failures.append(f"{'a jointed' if jointed else 'an unjointed'}"
                                    f" model's {where} panel is not the "
                                    f"vectors panel: {title!r}")
                if tints:
                    failures.append(f"{'a jointed' if jointed else 'an unjointed'}"
                                    f" model's {where} panel was tinted by "
                                    f"block")
            else:
                if arrows:
                    failures.append(f"a jointed model's {where} panel draws an "
                                    f"arrow field")
                if "Deformation" not in title:
                    failures.append(f"a jointed model's {where} panel is not "
                                    f"the deformed mesh: {title!r}")
                if "Scale" not in title:
                    failures.append(f"the {where} panel does not print its "
                                    f"exaggeration: {title!r}")
                # The joint faces are the point of the block picture: a deformed
                # grid with no faces on it says nothing about the joints.
                if len(faces) < 2:
                    failures.append(f"the {where} panel drew no joint faces "
                                    f"over its grid: {len(faces)} collections")
                # The blocks: a faint tint per body, and the outside of the
                # deformed mesh as a line of its own. Without the tint two blocks
                # that touch are one gray field; without the boundary the only
                # thing carrying the deformed ground surface is the light
                # element grid.
                if not tints:
                    failures.append(f"the {where} panel tints no blocks, so two "
                                    f"bodies that touch read as one")
                if not edges:
                    failures.append(f"the {where} panel draws no exterior "
                                    f"boundary, so the moved ground surface is "
                                    f"carried only by the element grid")
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
    for want in ("Bar tension", "Normal stress", "Shear stress", "Slip"):
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


    # Color by block, on the shipped toppling model's stored solution (no
    # solve): on, a tint per block with neighbors differing, so more than one;
    # off, one tint per material, and the model has one material.
    import matplotlib.figure as mplfig
    import matplotlib.colors as mcolors
    from matplotlib.collections import PolyCollection
    from xslope.fileio import load_slope_data
    from xslope.fem import build_fem_data, import_fem_solution
    from xslope.mesh import import_mesh_from_json
    from xslope.plot_fem import plot_fem_results
    stem = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "docs", "tutorials", "files", "xslope_rock_toppling")
    sd = _quiet(load_slope_data, stem + ".xlsx")
    fd = _quiet(build_fem_data, sd, _quiet(import_mesh_from_json,
                                          stem + "_mesh.json"))
    sol = _quiet(import_fem_solution, fd, stem)
    n_mats = len(np.unique(np.asarray(fd["element_materials"], dtype=int)))
    if n_mats != 1:
        failures.append(f"the toppling model carries {n_mats} materials; the "
                        f"Color by block assertion expects one")
    for on, want in ((True, "more than one"), (False, "exactly one")):
        fig = mplfig.Figure(figsize=(6.0, 4.0))
        _quiet(plot_fem_results, fd, sol, plot_type=["deformation"], fig=fig,
               show_joints=True, color_blocks=on)
        tints = {tuple(np.round(c[:3], 6))
                 for ax in fig.axes for coll in ax.collections
                 if isinstance(coll, PolyCollection) and coll.get_zorder() == 0.5
                 for c in coll.get_facecolors()}
        ok = len(tints) > 1 if on else len(tints) == 1
        if not ok:
            failures.append(f"with Color by block {'on' if on else 'off'} the "
                            f"toppling model's deformation panel fills "
                            f"{len(tints)} tint(s), not {want}")


# --------------------------------------------------------------------------
# i. the Run FEM dialog opens on the criterion the model needs
# --------------------------------------------------------------------------

def _leg_run_dialog(failures, cache):
    """A jointed model opens the Run FEM dialog on Hybrid.

    Under Non-convergence a near-critical trial on a jointed model counts as a
    failure whatever its joints are doing, so the bisection walks its lower
    bound to the floor and the run reports no factor of safety. The dialog
    reads the model the way the mesher does — either sheet's lines — and opens
    on the criterion that can rule on it.
    """
    from PySide6.QtWidgets import QApplication
    from studio.dialogs import RunFemDialog

    app = QApplication.instance() or QApplication([])   # noqa: F841
    sd, _mesh, _fem_data, _sol = cache["solved"]

    dlg = _quiet(RunFemDialog, defaults={}, slope_data=sd,
                 material_names=[m.get("name") for m in sd.get("materials") or []])
    if dlg.failure_criterion.currentData() != "hybrid":
        failures.append(f"a jointed model opens the Run FEM dialog on "
                        f"{dlg.failure_criterion.currentData()!r}, not 'hybrid'")
    dlg.deleteLater()

    # The same model with its Joint columns cleared is an ordinary bonded one.
    plain = _model(jointed=())
    dlg2 = _quiet(RunFemDialog, defaults={}, slope_data=plain,
                  material_names=[m.get("name") for m in plain.get("materials") or []])
    if dlg2.failure_criterion.currentData() != "non_convergence":
        failures.append(f"a model with no joint opens the Run FEM dialog on "
                        f"{dlg2.failure_criterion.currentData()!r}, not "
                        f"'non_convergence'")
    dlg2.deleteLater()

    # A criterion already chosen this session is not overwritten by the model.
    dlg3 = _quiet(RunFemDialog, defaults={"failure_criterion": "non_convergence"},
                  slope_data=sd,
                  material_names=[m.get("name") for m in sd.get("materials") or []])
    if dlg3.failure_criterion.currentData() != "non_convergence":
        failures.append("the model's own default overrode a criterion the user "
                        "had already chosen this session")
    dlg3.deleteLater()

    # The joints sheet is the other source, and is read the same way.
    sheet = _model(jointed=())
    sheet["joint_lines"] = [{"label": "bedding", "x1": 0.0, "y1": 0.0,
                             "x2": 10.0, "y2": 5.0, "phi": 35.0}]
    dlg4 = _quiet(RunFemDialog, defaults={}, slope_data=sheet,
                  material_names=[m.get("name") for m in sheet.get("materials") or []])
    if dlg4.failure_criterion.currentData() != "hybrid":
        failures.append("a model whose joints are on the joints sheet opens "
                        "the Run FEM dialog on "
                        f"{dlg4.failure_criterion.currentData()!r}, not 'hybrid'")
    dlg4.deleteLater()


# --------------------------------------------------------------------------
# j. a joints-sheet line's 1D elements are not a member
# --------------------------------------------------------------------------

def _leg_barless(failures, cache):
    """A joints-sheet line has 1D elements of its own and no bar between its
    faces, so nothing that reports on members may report on it.

    Both places that split ``elements_1d`` read the same ``barless_1d_mask``:
    the mesh plot, which drew the trace red and counted it as reinforcement, and
    the reinforcement sidecar, which wrote a per-bar force row for it with every
    force and every capacity zero. The fixture's own lines all carry bars, so
    the mask is flipped on a copy of its fem_data — which is exactly what a
    joints-sheet line produces — and the two readers are asked again.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from xslope.fem import _fem_reinforcement_dataframe as _reinf_df
    from xslope.plot_fem import plot_fem_data
    _sd, _mesh, fem_data, sol = cache["solved"]

    n_1d = len(fem_data["elements_1d"])
    if not n_1d:
        failures.append("the fixture carries no 1D elements to read")
        return

    # As built: every line carries a bar, so nothing is bar-less.
    df = _quiet(_reinf_df, fem_data, sol)
    if df is None or len(df) != n_1d:
        failures.append(f"the bonded fixture's sidecar lost rows: "
                        f"{0 if df is None else len(df)} of {n_1d}")
    _quiet(plot_fem_data, fem_data)
    title = plt.gcf().axes[0].get_title()
    if "reinforcement" not in title:
        failures.append(f"a bar on the mesh plot stopped being reinforcement: "
                        f"{title!r}")
    plt.close("all")

    # The same elements declared bar-less: a joints-sheet line.
    barless = dict(fem_data)
    barless["barless_1d_mask"] = np.ones(n_1d, dtype=bool)
    if _quiet(_reinf_df, barless, sol) is not None:
        failures.append("a bar-less joint line was written into the "
                        "reinforcement sidecar as a member with no force")
    _quiet(plot_fem_data, barless)
    title = plt.gcf().axes[0].get_title()
    if "reinforcement" in title:
        failures.append(f"the mesh plot counts a joint trace as reinforcement: "
                        f"{title!r}")
    if f"{n_1d} joint" not in title:
        failures.append(f"the mesh plot does not name the joint trace: "
                        f"{title!r}")
    labels = [h.get_label() for h in plt.gcf().legends[0].legend_handles] \
        if plt.gcf().legends else []
    if labels and not any(str(t).startswith("Joint (") for t in labels):
        failures.append(f"the mesh plot has no joint legend entry: {labels}")
    plt.close("all")


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
    ("the Run FEM dialog opens on Hybrid", _leg_run_dialog),
    ("a bar-less joint trace is not a bar", _leg_barless),
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
