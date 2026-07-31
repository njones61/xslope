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
     on, a thin band carries at least three element rows across its width on both
     families; with it off, it carries fewer than two. The split in C is not a
     style preference — a local size alone does not resolve a band on a triangular
     mesh, and the refine-feature alone does not resolve one on a quad mesh — so
     the mechanism has to be measured, not asserted.

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
    """Element rows the mesh puts across ``seam`` — its width over the median size
    of the elements whose centroid falls inside it."""
    import numpy as np
    from shapely.geometry import Point
    nodes = np.asarray(mesh["nodes"])
    sizes = []
    for e, t in zip(np.asarray(mesh["elements"]), np.asarray(mesh["element_types"])):
        n = 3 if t in (3, 6) else 4
        if not seam.contains(Point(nodes[e[:t]].mean(axis=0))):
            continue
        c = nodes[e[:n]]
        area = 0.5 * abs(sum(c[i][0] * c[(i + 1) % n][1] - c[(i + 1) % n][0] * c[i][1]
                             for i in range(n)))
        if area > 0:
            sizes.append(area ** 0.5)
    if not sizes:
        return 0.0
    return float(width) / float(np.median(sizes))


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
    from studio.runners import MeshWorker

    sd = _thin_model()
    seam, width = _seam(sd)

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
        # The promise. Three rows is the mesher's own definition of resolved; the
        # toggle asks for four, so a floor of three leaves room for the mesher's
        # own smoothing without letting a regression through.
        if rows[True] < 3.0:
            fails.append(f"{et} with the toggle on puts {rows[True]:.2f} element "
                         f"rows across the seam, want at least 3")
        # ...and it has to be the toggle doing it, not the target size.
        if rows[False] >= 2.0:
            fails.append(f"{et} already resolves the seam with the toggle off "
                         f"({rows[False]:.2f} rows) — this section no longer "
                         f"demonstrates the failure the toggle prevents")
    return fails


CHECKS = [("the toggle + its default", test_control),
          ("the derived per-zone size (quads)", test_derivation),
          ("the wire, one mechanism per family", test_wire),
          ("element rows across a real thin band", test_resolution)]


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
