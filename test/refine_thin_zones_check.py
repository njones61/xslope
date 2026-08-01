"""Standing checks on *Refine thin zones* — the Build-mesh dialog's toggle, the
mechanism it drives on each element family, and the resolution it promises.

A material zone too thin for the mesh to resolve does not fail. It solves, and it
returns a factor of safety that is too high, because a shear band cannot form
across one element. That is the failure this toggle exists to prevent, and it is
why the toggle is ON by default. What is guarded here:

  A. THE CONTROL — the dialog carries the checkbox, it is on when nothing says
     otherwise, ``options()`` reports what the box shows, and a stored preference
     reopens on the choice the user last made rather than reverting to the default.
  B. THE DERIVATION — ``thin_zone_size_regions`` (the quadrilateral mechanism)
     finds a thin band and asks for a size of one quarter of its width; it steps
     aside for a zone whose polygon already declares a Size, because that Size is
     the user saying how finely to mesh it; and it returns nothing at all on a
     section with no thin zone, which is what makes the option free there.
  C. THE WIRE — the toggle reaches ``build_mesh_from_polygons``, and it reaches it
     as a DIFFERENT argument per element family: ``thin_zones`` in
     ``refine_features`` for a triangular mesh, an extra ``size_regions`` entry for
     a quadrilateral one. Off sends neither. The builder is replaced by a recorder,
     so this leg is about the arguments and meshes nothing.
  D. THE RESOLUTION — the claim itself, measured on real meshes: with the toggle
     on, a thin band carries about four element rows across its width on BOTH
     families; with it off, it carries fewer than two. The split in C is not a
     style preference — a local size alone does not resolve a band on a triangular
     mesh — so the mechanism has to be measured, not asserted.
  E. THE THICKNESS IT MEASURES — a layer the section geometry splits into pieces
     is measured whole. A layer thick enough to resolve at the global size is left
     alone even when its individual pieces are not, and the refinement a zone gets
     does not move with the dialog's refinement factor, which is about a different
     feature class entirely.
  F. THE REPORT — a refinement the user did not ask for on this run says so: one
     Log line per zone naming it, the size it got and the rows that buys, and one
     line saying which control turns it off. A build that refines nothing is silent.

Rows are counted as the zone's width over the median element EDGE length inside it.
Not over the root of the element AREA: that is the same thing for a quad but
overstates a triangle by about 1.5x, and it was that proxy which read 4.7 rows on a
triangular mesh that was in fact delivering 3.

Skips cleanly (exit 0) when PySide6 is not installed.
"""
import contextlib
import io
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

#: The Griffiths thin-band section: a soft seam through an otherwise uniform slope,
#: and the case whose historical answer a quarter-of-the-thickness local size
#: reproduced. Its seam polygon carries a Size, which several legs strip so the
#: toggle is the only thing acting.
THIN = os.path.join(_REPO, "docs/fem/files/xslope_griffiths3_r0p2_thin.xlsx")
#: A section with no thin zone at all — the negative control for B.
PLAIN = os.path.join(_REPO, "docs/inputs/slope/xslope_simple1.xlsx")
#: A section whose top layer is ONE material stored as two polygons, because the
#: ground surface steps through it. Measured piecewise each polygon looks too thin
#: to resolve; measured as the layer it is, it is not thin at all. The case for E.
SPLIT = os.path.join(_REPO, "docs/verification/files/rocscience/rs2_51.xlsx")
#: Target for SPLIT — the size at which its 10-unit top layer reads as two 5-unit
#: bands. Below the 3-rows-across threshold for each piece, above it for the layer.
SPLIT_TARGET = 2.8

#: Target element size for the measured leg. Coarse enough that the seam is
#: unresolved without the toggle, which is the condition the check is about.
TARGET = 3.0


def _quiet(fn, *args, **kwargs):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        return fn(*args, **kwargs)


def _thin_model():
    """The thin-band model with its own zone Size stripped."""
    from xslope.fileio import load_slope_data
    sd = _quiet(load_slope_data, THIN)
    for p in sd.get("polygons") or []:
        p["size"] = None
    return sd


def _seam(sd):
    """(shapely polygon, local width) of the model's thin seam."""
    from shapely.geometry import Polygon
    from xslope.mesh import detect_thin_zones, get_material_polygons
    polys = get_material_polygons(sd)
    coords = [p["coords"] for p in polys]
    zones = [z for z in detect_thin_zones(coords, TARGET) if z["kind"] == "whole"]
    z = zones[0]
    return Polygon(coords[z["poly_index"]]), z["width"]


def _rows_across(mesh, seam, width):
    """Element rows the mesh puts across ``seam`` — its width over the median EDGE
    length of the elements whose centroid falls inside it.

    Edge length, not root-area: for a quad the two agree, but an equilateral triangle
    of edge h has area 0.433h², so its root-area is 0.66h and a triangular mesh reads
    1.5x more rows than it has. Both families are compared against one promise here,
    so the metric has to mean the same thing on both."""
    import numpy as np
    from shapely.geometry import Point
    nodes = np.asarray(mesh["nodes"])
    lengths = []
    for e, t in zip(np.asarray(mesh["elements"]), np.asarray(mesh["element_types"])):
        n = 3 if t in (3, 6) else 4
        if not seam.contains(Point(nodes[e[:t]].mean(axis=0))):
            continue
        c = nodes[e[:n]]
        for i in range(n):
            lengths.append(float(np.hypot(*(c[(i + 1) % n] - c[i]))))
    if not lengths:
        return 0.0
    return float(width) / float(np.median(lengths))


def _dialog(**defaults):
    from studio.dialogs import BuildMeshDialog
    return BuildMeshDialog(None, defaults=defaults)


# ------------------------------------------------------------------ A. control
def test_control():
    fails = []
    from studio.dialogs import REFINE_THIN_ZONES_TIP
    from studio.runners import REFINE_THIN_ZONES_DEFAULT

    dlg = _dialog()
    if dlg.refine_thin_zones.isChecked() is not REFINE_THIN_ZONES_DEFAULT:
        fails.append(f"a fresh dialog opened with the thin-zone toggle "
                     f"{dlg.refine_thin_zones.isChecked()}, want the declared "
                     f"default {REFINE_THIN_ZONES_DEFAULT}")
    if not REFINE_THIN_ZONES_DEFAULT:
        fails.append("the thin-zone toggle defaults OFF: an unresolved thin zone "
                     "returns a factor of safety that is too high and says nothing, "
                     "so the safe state is the default one")
    if dlg.options().get("refine_thin_zones") != dlg.refine_thin_zones.isChecked():
        fails.append("options() disagrees with the checkbox on open")
    if dlg.refine_thin_zones.toolTip() != REFINE_THIN_ZONES_TIP:
        fails.append("the thin-zone toggle carries no tooltip — the option is "
                     "on by default and unexplained")

    # The box must be readable in both directions, or the user cannot turn it off.
    for want in (False, True):
        dlg.refine_thin_zones.setChecked(want)
        if dlg.options()["refine_thin_zones"] is not want:
            fails.append(f"setting the toggle to {want} was not reported by options()")

    # Reopening on the session's last choice, including the non-default one.
    for want in (False, True):
        d = _dialog(refine_thin_zones=want)
        if d.options()["refine_thin_zones"] is not want:
            fails.append(f"a dialog defaulted to {want} opened on "
                         f"{d.options()['refine_thin_zones']}")

    # Not gated on the element type: both families have the failure.
    for et in ("tri3", "tri6", "quad4", "quad8", "quad9"):
        d = _dialog(element_type=et)
        idx = d.element_type.findData(et)
        d.element_type.setCurrentIndex(idx)
        if not d.refine_thin_zones.isEnabled():
            fails.append(f"the thin-zone toggle is dimmed for {et}, which has the "
                         f"same unresolved-band failure as every other type")
    return fails


# -------------------------------------------------------------- B. derivation
def test_derivation():
    fails = []
    from studio.runners import THIN_ZONE_ELEMS, thin_zone_size_regions
    from xslope.mesh import get_material_polygons

    sd = _thin_model()
    seam, width = _seam(sd)
    polys = get_material_polygons(sd)

    regions = thin_zone_size_regions(polys, TARGET)
    if not regions:
        fails.append("the seam section derived no size region at all")
    else:
        want = width / THIN_ZONE_ELEMS
        got = min(r["size"] for r in regions)
        if abs(got - want) > 1e-9:
            fails.append(f"the seam was sized at {got:.4g}, want its width over "
                         f"{THIN_ZONE_ELEMS} = {want:.4g}")
        if got >= TARGET:
            fails.append(f"the derived size {got:.4g} is not finer than the target "
                         f"{TARGET:g}, so it would be inert")
        for r in regions:
            if len(r.get("polygon") or []) < 3:
                fails.append("a derived region carries no usable ring")

    # A zone that states its own Size is the user's call, and stands.
    sd2 = _thin_model()
    for p in sd2["polygons"]:
        p["size"] = 1.25
    kept = thin_zone_size_regions(get_material_polygons(sd2), TARGET)
    if any(r["size"] != 1.25 for r in kept):
        fails.append("a zone that declares its own Size was given a derived one — "
                     "an automatic option overruled an explicit input")

    # Nothing to do on a section with no thin zone: this is what makes the option
    # free to leave on.
    from xslope.fileio import load_slope_data
    plain = _quiet(load_slope_data, PLAIN)
    if thin_zone_size_regions(get_material_polygons(plain), 1.0):
        fails.append("a section with no thin zone still derived a size region")
    return fails


# --------------------------------------------------------------------- C. wire
def test_wire():
    fails = []
    import xslope.mesh as mesh_mod
    from studio.runners import MeshWorker

    sd = _thin_model()
    seen = []
    real = mesh_mod.build_mesh_from_polygons

    def recorder(polygons, **kwargs):
        """Stands in for the builder: records the call and returns a mesh-shaped
        stub, so this leg is about the arguments and never meshes anything."""
        seen.append(kwargs)
        return {"nodes": [], "elements": [], "element_types": [],
                "element_materials": []}

    def call(element_type, thin):
        dlg = _dialog(element_type=element_type, auto_size=False,
                      target_size=TARGET, refine_thin_zones=thin)
        idx = dlg.element_type.findData(element_type)
        dlg.element_type.setCurrentIndex(idx)
        worker = MeshWorker()
        errors = []
        worker.failed.connect(errors.append)
        _quiet(worker.build, sd, dlg.options())
        return (seen[-1] if seen else None), errors

    n_plain = len(mesh_mod.extract_size_regions(sd))
    mesh_mod.build_mesh_from_polygons = recorder
    try:
        # Triangular: the refine-features path.
        on, err = call("tri6", True)
        if err:
            fails.append(f"the tri6 build failed: {err[0]}")
        elif "thin_zones" not in (on.get("refine_features") or ()):
            fails.append(f"tri6 with the toggle on did not ask for thin_zones "
                         f"(refine_features={on.get('refine_features')!r})")
        elif not on.get("refine_factor"):
            fails.append("tri6 asked for thin_zones with no refinement factor, "
                         "which turns the whole refine path off")

        off, err = call("tri6", False)
        if err:
            fails.append(f"the tri6 build with the toggle off failed: {err[0]}")
        elif "thin_zones" in (off.get("refine_features") or ()):
            fails.append("tri6 with the toggle OFF still asked for thin_zones — "
                         "the toggle is not wired, it is hardcoded on")
        elif off.get("refine_factor") is not None:
            fails.append("tri6 with every refinement off still set a refine factor")

        # Quadrilateral: a derived local size, not the refine-features path.
        qon, err = call("quad4", True)
        if err:
            fails.append(f"the quad4 build failed: {err[0]}")
        else:
            n = len(qon.get("size_regions") or [])
            if n <= n_plain:
                fails.append(f"quad4 with the toggle on passed {n} size region(s), "
                             f"the same {n_plain} the model declares — the derived "
                             f"size never reached the builder")
        qoff, err = call("quad4", False)
        if err:
            fails.append(f"the quad4 build with the toggle off failed: {err[0]}")
        elif len(qoff.get("size_regions") or []) != n_plain:
            fails.append("quad4 with the toggle OFF still carried a derived size "
                         "region — the toggle is not wired")

        # The two families must not be sent the same thing: that is the whole
        # finding behind the split, and a refactor that collapses them would pass
        # every other assertion here.
        if (qon.get("refine_features") or ()) == (on.get("refine_features") or ()) \
                and len(qon.get("size_regions") or []) == len(on.get("size_regions") or []):
            fails.append("the toggle sent quads and triangles identical arguments; "
                         "each family needs its own mechanism (see "
                         "studio.runners.thin_zone_size_regions)")
    finally:
        mesh_mod.build_mesh_from_polygons = real
    return fails


# --------------------------------------------------------------- D. resolution
def test_resolution():
    fails = []
    from studio.runners import MeshWorker, THIN_ZONE_ELEMS

    sd = _thin_model()
    seam, width = _seam(sd)
    seen = {}

    for et in ("tri6", "quad4"):
        rows = {}
        for thin in (False, True):
            dlg = _dialog(element_type=et, auto_size=False, target_size=TARGET,
                          refine_thin_zones=thin)
            idx = dlg.element_type.findData(et)
            dlg.element_type.setCurrentIndex(idx)
            worker = MeshWorker()
            built, errors = [], []
            worker.succeeded.connect(built.append)
            worker.failed.connect(errors.append)
            _quiet(worker.build, sd, dlg.options())
            if errors or not built:
                fails.append(f"{et} thin={thin}: the mesh build failed "
                             f"({errors[0] if errors else 'no mesh emitted'})")
                rows[thin] = None
                continue
            rows[thin] = _rows_across(built[-1], seam, width)

        if rows.get(True) is None or rows.get(False) is None:
            continue
        seen[et] = rows[True]
        # The promise, which is FOUR rows and the same four on either family. A
        # lower bound of three would pass a triangular mesh that quietly delivered
        # three while the quad path delivered four — which is what it did, hidden
        # behind a root-area row count that read the difference as 4.7 against 4.1.
        # Banded, because the ceiling is the half that catches over-refinement.
        if not (3.5 <= rows[True] <= 5.0):
            fails.append(f"{et} with the toggle on puts {rows[True]:.2f} element "
                         f"rows across the seam, want about "
                         f"{THIN_ZONE_ELEMS} (3.5 to 5.0)")
        # ...and it has to be the toggle doing it, not the target size.
        if rows[False] >= 2.0:
            fails.append(f"{et} already resolves the seam with the toggle off "
                         f"({rows[False]:.2f} rows) — this section no longer "
                         f"demonstrates the failure the toggle prevents")

    # The two families must land on the SAME resolution: one option, one promise.
    if all(isinstance(v, float) and v for v in seen.values()):
        spread = max(seen.values()) - min(seen.values())
        if spread > 0.5:
            fails.append("the families disagree on how finely a thin band is meshed "
                         + ", ".join(f"{k}={v:.2f} rows" for k, v in seen.items())
                         + " — the toggle promises one resolution, not one per "
                           "element type")
    return fails


# ------------------------------------------------------- E. the width it measures
def test_measured_thickness():
    """A layer is measured whole, and a layer thick enough to resolve is left alone.

    rs2_51's top layer is 10 units thick and stored as two 5-unit polygons, because
    the ground surface steps through it. At a 2.8 target each piece is under the
    three-rows-across threshold on its own while the layer is over it, so the two
    readings disagree and one of them is the material's. Measured piecewise the
    section came back with the layer meshed at roughly twice the rows its thickness
    asks for; measured whole it needs no refinement at all.
    """
    fails = []
    import numpy as np
    from xslope.fileio import load_slope_data
    from xslope.mesh import (build_mesh_from_polygons, detect_thin_zones,
                             extract_constraint_line_geometry,
                             get_material_polygons, thin_zone_plan)

    sd = _quiet(load_slope_data, SPLIT)
    lines, _, _ = extract_constraint_line_geometry(sd)
    polys = get_material_polygons(sd, reinf_lines=lines)
    coords = [p["coords"] for p in polys]
    mat_ids = [p["mat_id"] for p in polys]

    piecewise = [z for z in detect_thin_zones(coords, SPLIT_TARGET)
                 if z["kind"] == "whole"]
    if len(piecewise) < 2:
        fails.append("rs2_51 no longer reads as several thin pieces when measured "
                     "polygon by polygon — this section no longer demonstrates the "
                     "split-layer case")
    whole = detect_thin_zones(coords, SPLIT_TARGET, group_ids=mat_ids)
    if whole:
        fails.append(f"rs2_51's top layer is {10.0 / SPLIT_TARGET:.1f} element rows "
                     f"thick at a {SPLIT_TARGET} target and was still flagged thin "
                     f"({[z['kind'] for z in whole]}) — a layer split into polygons "
                     f"is being measured piece by piece")
    if thin_zone_plan(coords, SPLIT_TARGET, group_ids=mat_ids):
        fails.append("rs2_51 planned a thin-zone refinement it does not need")

    # ...and that has to show in the mesh, not only in the plan: refining nothing
    # must mean meshing exactly what the toggle-off build meshes.
    off = _quiet(build_mesh_from_polygons, polys, target_size=SPLIT_TARGET,
                 element_type="tri6", lines=lines or None)
    on = _quiet(build_mesh_from_polygons, polys, target_size=SPLIT_TARGET,
                element_type="tri6", lines=lines or None,
                refine_factor=3.0, refine_features=["thin_zones"])
    if len(on["elements"]) != len(off["elements"]):
        fails.append(f"rs2_51 meshes to {len(on['elements'])} elements with thin-zone "
                     f"refinement against {len(off['elements'])} without — a section "
                     f"with no thin zone must pay nothing for the option")

    # The refinement a zone gets is its own thickness over four, so the dialog's
    # refinement factor — which sizes the distance bands around lines and tips —
    # must not move it. It used to be read as the thing that set the size.
    sd2 = _thin_model()
    tpolys = get_material_polygons(sd2)
    counts = set()
    for factor in (2.0, 3.0, 8.0):
        m = _quiet(build_mesh_from_polygons, tpolys, target_size=TARGET,
                   element_type="tri6", refine_factor=factor,
                   refine_features=["thin_zones"])
        counts.add(len(m["elements"]))
    if len(counts) != 1:
        fails.append(f"the thin band meshed to {sorted(counts)} elements at "
                     f"refinement factors 2/3/8 — a thin zone is sized by its own "
                     f"thickness and must not move with the factor")
    return fails


# ------------------------------------------------------------------- F. the report
def test_report():
    """The Log names what was refined, or says nothing at all."""
    fails = []
    import re
    from studio.runners import MeshWorker, THIN_ZONE_ELEMS

    def log_for(sd, element_type, thin=True, target=TARGET):
        dlg = _dialog(element_type=element_type, auto_size=False, target_size=target,
                      refine_thin_zones=thin)
        dlg.element_type.setCurrentIndex(dlg.element_type.findData(element_type))
        buf = io.StringIO()
        worker = MeshWorker()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(io.StringIO()):
            worker.build(sd, dlg.options())
        return buf.getvalue()

    CLOSER = "('Refine thin zones' in the Build mesh dialog disables this)"
    #: name -> local size (~rows across width)
    LINE = re.compile(r"thin zone refined: '(?P<name>[^']+)' -> local size "
                      r"(?P<size>[\d.]+) \(~(?P<rows>\d+) rows across (?P<w>[\d.]+)")

    for et in ("tri6", "quad4"):
        text = log_for(_thin_model(), et)
        hits = [LINE.search(l) for l in text.splitlines()]
        hits = [h for h in hits if h]
        if not hits:
            fails.append(f"{et} refined a thin zone and the Log never said so:\n{text}")
            continue
        if len(hits) != 1:
            fails.append(f"{et} reported {len(hits)} thin zones on a section with one")
        m = hits[0]
        if not m.group("name").strip():
            fails.append(f"{et} reported a thin zone with no name")
        if m.group("rows") != str(THIN_ZONE_ELEMS):
            fails.append(f"{et} reported {m.group('rows')} rows, but the option "
                         f"asks for {THIN_ZONE_ELEMS} — the line does not describe "
                         f"what the mesher did")
        # The printed size and width have to be each other's quotient at the
        # promised row count, or the line is decoration rather than a report.
        size, w = float(m.group("size")), float(m.group("w"))
        if size <= 0 or abs(w / size - THIN_ZONE_ELEMS) > 0.05 * THIN_ZONE_ELEMS:
            fails.append(f"{et} printed size {size} across width {w}: that is "
                         f"{w / size:.2f} rows, not the {THIN_ZONE_ELEMS} the "
                         f"option promises")
        if CLOSER not in text:
            fails.append(f"{et} refined without saying which control turns it off")

    # Nothing refined -> nothing said. An option that narrates itself on every
    # build is an option users stop reading. (PLAIN at a 1-unit target, the size at
    # which it carries no thin zone — at 3 the erosion leaves a residue in its toe
    # wedge, which is a real detection and does get reported.)
    from xslope.fileio import load_slope_data
    quiet_text = log_for(_quiet(load_slope_data, PLAIN), "tri6", target=1.0)
    if "thin zone refined" in quiet_text or CLOSER in quiet_text:
        fails.append("a section with no thin zone still reported a refinement:\n"
                     + quiet_text)

    # ...and the toggle OFF is silent even where there IS a zone.
    off_text = log_for(_thin_model(), "tri6", thin=False)
    if "thin zone refined" in off_text or CLOSER in off_text:
        fails.append("the toggle was off and the Log still reported refinement:\n"
                     + off_text)
    return fails


CHECKS = [("the toggle + its default", test_control),
          ("the derived per-zone size (quads)", test_derivation),
          ("the wire, one mechanism per family", test_wire),
          ("element rows across a real thin band", test_resolution),
          ("the thickness it measures (split layers)", test_measured_thickness),
          ("what the Log reports", test_report)]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    try:
        import PySide6                                    # noqa: F401
    except Exception:
        print("refine thin zones: PySide6 not installed — checks skipped.")
        return []
    failures = []
    for name, fn in CHECKS:
        try:
            fs = fn()
        except Exception as exc:
            import traceback
            traceback.print_exc()
            fs = [f"{name} raised: {exc!r}"]
        print(f"  {name:44s} {'ok' if not fs else f'FAIL ({len(fs)})'}")
        failures += fs
    return failures


def main():
    print("thin-zone mesh refinement (Studio):")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll thin-zone refinement checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication
        _app = QApplication.instance() or QApplication([])
    except Exception:
        pass
    main()
