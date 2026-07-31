"""Standing checks on the quadrilateral meshing style as the user sets it — the
Studio *Build mesh* dialog's style control, and the wire from that control to the
mesh builder.

The two styles are a real discretization choice: 'free' recombines a triangulation
into quadrilaterals everywhere, 'structured' additionally sweeps regular grids
through the zones that can carry one. They are offered per run and stored nowhere,
so nothing in the input file would catch it if the control stopped reaching the
builder. What is guarded here:

  A. THE CONTROL — the dialog carries a two-item radio group, the labels are the
     two styles, exactly one is selected, and the selection is 'free' on a dialog
     opened with no preference. A dialog that quietly opened on 'structured' would
     change every mesh built from it.
  B. THE GATE — the group is live for the quadrilateral element types and dimmed,
     with the reason as its tooltip, for the triangular ones, which the builder
     does not apply a quadrilateral style to at all. The gate has to follow the
     element type live, because the type is chosen in the same dialog.
  C. THE WIRE — the selected style reaches ``build_mesh_from_polygons`` as
     ``quad_style``, and a different selection sends a different value. This is the
     check that fails if the option is added to the dialog and never passed on: the
     builder is replaced by a recorder, so the mesh itself is not built and the
     assertion is about the argument, not about the mesh.
  D. THE REOPEN — the dialog opens on a style handed to it as a default, which is
     how the session's last choice greets the next build rather than the choice
     silently reverting.

The mesh builder's own behaviour under each style — element quality, delivered
size, which zones sweep — is checked in test/quad_mesh_check.py; nothing here
builds a mesh.

Skips cleanly (exit 0) when PySide6 is not installed.
"""
import io
import contextlib
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

#: A section with quadrilateral-friendly zones — the runner needs a real model to
#: reach the builder call, but the builder is mocked, so nothing is meshed.
LEVEE = os.path.join(_REPO, "docs/seep/files/xslope_levee_poly.xlsx")

TRI_TYPES = ("tri3", "tri6")
QUAD_TYPES = ("quad4", "quad8", "quad9")


def _quiet(fn, *args, **kwargs):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        return fn(*args, **kwargs)


def _set_type(dlg, element_type):
    idx = dlg.element_type.findData(element_type)
    dlg.element_type.setCurrentIndex(idx)
    return idx


# ------------------------------------------------------------------ A. control
def test_control():
    fails = []
    from studio.dialogs import QUAD_STYLES, BuildMeshDialog

    dlg = BuildMeshDialog(None)
    radios = dlg.quad_style_buttons.buttons()
    if len(radios) != 2:
        fails.append(f"the quad style group has {len(radios)} buttons, want 2")
    if not dlg.quad_style_buttons.exclusive():
        fails.append("the quad style group is not exclusive — both styles could be set")

    labels = [b.text() for b in radios]
    for _, want, _ in QUAD_STYLES:
        if want not in labels:
            fails.append(f"no radio labelled {want!r} (have {labels})")

    checked = [b.text() for b in radios if b.isChecked()]
    if len(checked) != 1:
        fails.append(f"{len(checked)} styles selected on open, want exactly 1")

    _set_type(dlg, "quad4")
    if dlg.quad_style() != "free":
        fails.append(f"a fresh dialog opened on {dlg.quad_style()!r}, want 'free'")
    if dlg.options().get("quad_style") != "free":
        fails.append("options() disagrees with the dialog's own selection on open")

    # The tooltips are the entire explanation the user gets, so their absence is a
    # failure of the control, not a cosmetic one.
    for key, _, tip in QUAD_STYLES:
        if dlg._quad_style_radios[key].toolTip() != tip:
            fails.append(f"the {key!r} radio carries no tooltip on a quad mesh")

    # Selecting the other style must be what the dialog then reports.
    dlg._quad_style_radios["structured"].setChecked(True)
    if dlg.quad_style() != "structured" or dlg.options()["quad_style"] != "structured":
        fails.append("selecting Structured did not change the reported style")
    dlg._quad_style_radios["free"].setChecked(True)
    if dlg.quad_style() != "free":
        fails.append("selecting Free back did not change the reported style")
    return fails


# --------------------------------------------------------------------- B. gate
def test_gate():
    fails = []
    from studio.dialogs import QUAD_STYLE_TRI_TIP, BuildMeshDialog

    dlg = BuildMeshDialog(None)
    for et in QUAD_TYPES:
        _set_type(dlg, et)
        if not dlg.quad_style_group.isEnabled():
            fails.append(f"the quad style group is dimmed for {et}")
        if dlg.quad_style_group.toolTip():
            fails.append(f"{et} still carries the not-applicable tooltip")
    for et in TRI_TYPES:
        _set_type(dlg, et)
        if dlg.quad_style_group.isEnabled():
            fails.append(f"the quad style group is live for {et}")
        if dlg.quad_style_group.toolTip() != QUAD_STYLE_TRI_TIP:
            fails.append(f"{et} dims the styles without saying why "
                         f"(tooltip {dlg.quad_style_group.toolTip()!r})")

    # The gate has to track the type live, in both directions, within one dialog.
    _set_type(dlg, "quad8")
    if not dlg.quad_style_group.isEnabled():
        fails.append("switching tri -> quad did not re-enable the style group")
    _set_type(dlg, "tri3")
    if dlg.quad_style_group.isEnabled():
        fails.append("switching quad -> tri did not dim the style group")

    # Dimmed is not the same as forgotten: the dialog still reports a style, so a
    # user who sets Structured, glances at tri, and comes back has not lost it.
    _set_type(dlg, "quad4")
    dlg._quad_style_radios["structured"].setChecked(True)
    _set_type(dlg, "tri6")
    _set_type(dlg, "quad4")
    if dlg.quad_style() != "structured":
        fails.append("a trip through a triangular element type reset the style")
    return fails


# --------------------------------------------------------------------- C. wire
def test_wire():
    fails = []
    import xslope.mesh as mesh_mod
    from studio.dialogs import BuildMeshDialog
    from studio.runners import MeshWorker
    from xslope.fileio import load_slope_data

    slope_data = _quiet(load_slope_data, LEVEE)

    seen = []
    real = mesh_mod.build_mesh_from_polygons

    def recorder(polygons, **kwargs):
        """Stands in for the builder: records the call and returns a mesh-shaped
        stub, so this check is about the argument and never meshes anything."""
        seen.append(kwargs)
        return {"nodes": [], "elements": [], "element_types": [],
                "element_materials": []}

    mesh_mod.build_mesh_from_polygons = recorder
    try:
        for want in ("free", "structured"):
            dlg = BuildMeshDialog(None)
            _set_type(dlg, "quad4")
            dlg._quad_style_radios[want].setChecked(True)
            worker = MeshWorker()
            errors = []
            worker.failed.connect(errors.append)
            _quiet(worker.build, slope_data, dlg.options())
            if errors:
                fails.append(f"the {want} build failed: {errors[0]}")
                continue
            if not seen:
                fails.append(f"the {want} build never reached the mesh builder")
                continue
            got = seen[-1].get("quad_style")
            if got != want:
                fails.append(f"the dialog selected {want!r} and the builder was "
                             f"called with quad_style={got!r}")
    finally:
        mesh_mod.build_mesh_from_polygons = real

    if len(seen) == 2 and seen[0].get("quad_style") == seen[1].get("quad_style"):
        fails.append("both styles reached the builder as the same value — "
                     "re-meshing after toggling would rebuild the same mesh")

    # A triangular mesh still passes a legal value: the builder rejects anything
    # that is not one of the two styles, whatever the element type.
    mesh_mod.build_mesh_from_polygons = recorder
    try:
        dlg = BuildMeshDialog(None)
        _set_type(dlg, "tri6")
        worker = MeshWorker()
        _quiet(worker.build, slope_data, dlg.options())
    finally:
        mesh_mod.build_mesh_from_polygons = real
    if seen[-1].get("quad_style") not in ("free", "structured"):
        fails.append(f"a triangular build passed quad_style="
                     f"{seen[-1].get('quad_style')!r}, which the builder refuses")
    return fails


# ------------------------------------------------------------------- D. reopen
def test_reopen():
    fails = []
    from studio.dialogs import BuildMeshDialog

    for want in ("free", "structured"):
        dlg = BuildMeshDialog(None, defaults={"element_type": "quad4",
                                              "quad_style": want})
        if dlg.quad_style() != want:
            fails.append(f"a dialog defaulted to {want!r} opened on "
                         f"{dlg.quad_style()!r}")
    # An unreadable stored value falls back to the default rather than leaving the
    # group with nothing selected.
    dlg = BuildMeshDialog(None, defaults={"quad_style": "paved"})
    if dlg.quad_style() != "free":
        fails.append(f"an unknown stored style opened on {dlg.quad_style()!r}, "
                     "want the default 'free'")
    if not any(b.isChecked() for b in dlg.quad_style_buttons.buttons()):
        fails.append("an unknown stored style left no style selected at all")
    return fails


CHECKS = [("the style control + its default", test_control),
          ("the quad-only gate + its reason", test_gate),
          ("the wire to build_mesh_from_polygons", test_wire),
          ("reopening on a remembered style", test_reopen)]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    try:
        import PySide6                                    # noqa: F401
    except Exception:
        print("quad style dialog: PySide6 not installed — checks skipped.")
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
    print("quadrilateral mesh style (Studio):")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll quadrilateral mesh style checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication
        _app = QApplication.instance() or QApplication([])
    except Exception:
        pass
    main()
