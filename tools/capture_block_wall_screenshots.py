"""Studio dialog captures for Tutorial FEM-3, headlessly (offscreen Qt).

Same practice as ``tools/capture_tutorial_screenshots.py``, whose helpers this
module reuses rather than copies: build the real Studio widget, ``.show()`` it
under the ``offscreen`` QPA platform, let the layout settle, then
``QWidget.grab()`` it to a PNG.  Nothing here touches a live display, and nothing
here is restyled — every shot is the widget as the reader meets it.

Four shots, one for each thing FEM-3 asks the reader to do that a worksheet
cannot show: type the seven contacts of a block wall into the joints editor, set
the Joint column on three geogrid layers, run the strength reduction at the
settings a jointed model needs, and read the model checks that tell a bonded
sheet it should have been a joint.

The plot figures are a separate producer: ``tools/make_block_wall_figures.py``.

Run:  python3 tools/capture_block_wall_screenshots.py --all      # every shot
      python3 tools/capture_block_wall_screenshots.py joints     # by name
"""

from __future__ import annotations

import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import matplotlib
matplotlib.use("Agg")

try:
    from PySide6.QtWidgets import QApplication                     # noqa: F401
except Exception:                       # engine-only install — no studio layer
    print("capture_block_wall_screenshots: PySide6 not installed — skipped.")
    sys.exit(0)

# The offscreen QApplication, the modal-dialog muzzle, the settle pump, the grab,
# the table/list sizing idioms and the app-defaults guard all belong to the
# tutorial capture pipeline and are used unchanged, so a change to how a tutorial
# dialog is photographed reaches FEM-3's shots too.
from tools.capture_tutorial_screenshots import (OUT_DIR,             # noqa: E402
                                                _app_defaults, _fem_only,
                                                _grab, _line_table, _load,
                                                _settle)

FILES = os.path.join(REPO_ROOT, "docs", "tutorials", "files")
FEM03_WALL = os.path.join(FILES, "xslope_block_wall.xlsx")
FEM03_GRID = os.path.join(FILES, "xslope_block_wall_grid.xlsx")
FEM03_LINER_BONDED = os.path.join(FILES, "xslope_liner_bonded.xlsx")

#: The mesh the page builds for part 1: quadratic triangles at a 0.8 m global
#: target, with each block polygon's own 0.3 m Size read from the file.
FEM03_ELEMENT_TYPE = "tri6"
FEM03_TARGET_SIZE = 0.8
#: Run FEM, as the page has the reader set it.  The failure criterion is NOT in
#: this list, and that is the point of the shot: the dialog opens a jointed model
#: on Hybrid by itself, so what is photographed is the choice the reader does not
#: have to make.  The bracket and the tolerance are the dialog's own; the sweep
#: budget is raised to 100,000, which is the model checks' floor for a jointed
#: model and the box's own ceiling.
FEM03_BRACKET = (1.0, 2.0)
FEM03_TOLERANCE = 0.01
FEM03_MAX_ITERATIONS = 100000
#: How wide the model-checks panel is photographed.  A finding's line is elided to
#: the panel's width, and this is the width at which all three of the bonded
#: liner's warnings still carry the line they are about before the ellipsis — a
#: panel narrow enough to cut that would photograph three findings about nothing
#: in particular.
FEM03_CHECKS_WIDTH = 1180


def _meshed(path):
    """The model with a SPLIT mesh attached — the state Build Mesh leaves behind,
    and the only state Studio's Run FEM action is reachable in.

    ``extract_joint_options`` is what makes it a jointed mesh rather than a mesh
    with lines drawn on it: it tells the mesher which constraint lines to split
    along.  Run FEM's own dialog reads the same function to decide which failure
    criterion to open on, so a shot taken on a mesh built without it would
    photograph the wrong default.
    """
    import contextlib
    import io

    from xslope.mesh import (build_mesh_from_polygons,
                             extract_constraint_line_geometry,
                             extract_joint_options,
                             extract_point_constraints,
                             extract_size_regions, get_material_polygons)

    data = _load(path)
    lines, _n_reinf, _n_pile = extract_constraint_line_geometry(data)
    with contextlib.redirect_stdout(io.StringIO()):
        data["mesh"] = build_mesh_from_polygons(
            get_material_polygons(data, reinf_lines=lines),
            FEM03_TARGET_SIZE, FEM03_ELEMENT_TYPE, lines=lines or None,
            element_size_1d=data.get("element_size_1d"),
            point_constraints=extract_point_constraints(data),
            size_regions=extract_size_regions(data),
            joint_lines=extract_joint_options(data))
    return data


def fem03_joints_editor():
    """The joints editor, table view, with the seven contacts of part 1 in it.

    The table is taken out to Jred, its last column, rather than to the last
    column the seven rows fill.  Everything past t_cut is BLANK on every row and
    the page is telling the reader to leave it blank — the residual strengths, the
    dilation, and the normal and shear stiffnesses derived from the material
    around the contact — so a shot cut off at t_cut would photograph a table whose
    empty cells are the instruction.

    No usage band is set: every field of a joint line is FEM-only, so the editor
    shows one band and there is nothing to untick.
    """
    from studio.editors import JointsEditor

    dlg = JointsEditor().build(_load(FEM03_WALL), None)
    return _grab(_line_table(dlg, through="jred"),
                 "fem03_studio_joints_editor.png")


def fem03_reinforce_editor():
    """The reinforcement editor, table view, with part 2's three geogrid layers.

    The FEM band alone, because this section runs nothing but the finite element
    engine: the band carries the Joint column, which is the one cell the page is
    about, together with the Adhesion and Delta the interface takes its strength
    from and the E and Area the bar takes its stiffness from.  Out to Jred, the
    band's last column, for the reason the joints shot gives — the stiffnesses are
    left blank on purpose.
    """
    from studio.editors import ReinforcementEditor

    dlg = _fem_only(ReinforcementEditor().build(_load(FEM03_GRID), None))
    return _grab(_line_table(dlg, through="jred"),
                 "fem03_studio_reinforce_editor.png")


def fem03_run_fem():
    """Run FEM on the meshed wall, at the page's settings.

    Only the bracket, the tolerance and the sweep budget are passed.  The failure
    criterion is left to the dialog, which opens a jointed model on Hybrid and
    every other model on Non-convergence, so the shot shows the reader the box
    already reading what their model needs.
    """
    from studio.dialogs import RunFemDialog

    data = _meshed(FEM03_WALL)
    dlg = RunFemDialog(defaults={"analysis": "ssrm",
                                 "F_min": FEM03_BRACKET[0],
                                 "F_max": FEM03_BRACKET[1],
                                 "tolerance": FEM03_TOLERANCE,
                                 "max_iterations": FEM03_MAX_ITERATIONS},
                       material_names=[m.get("name") for m in data["materials"]],
                       slope_data=data)
    dlg.resize(dlg.sizeHint())
    return _grab(dlg, "fem03_studio_run_fem.png")


def fem03_preflight_liner():
    """The model checks on the bonded liner — three warnings, one sheet.

    The panel is photographed on its own rather than inside a Run dialog because
    the three lines ARE the subject: a sheet on a material boundary, a sheet flat
    and long enough to slide on, and an interface friction angle well below the
    soil's own.  Each is a different reading of the same row, and together they
    are the argument the page's part 3 then runs twice to settle.

    The panel opens at a height meant for a model with a screenful of findings,
    so on a three-warning model most of it is empty.  The width is set first —
    wide enough that every line still names its own subject, ``Reinforcement line
    1 ('liner')``, before the elision — and the height is then taken from the
    selected entry's own wrapped text at that width, the ``heightForWidth`` idiom,
    so the figure is the three lines and their detail and nothing else.  The first
    line is selected, which is the one the page quotes.
    """
    from PySide6.QtWidgets import QLabel

    from studio.preflight_panel import PreflightPanel

    data = _load(FEM03_LINER_BONDED)
    panel = PreflightPanel("fem", slope_data=data)
    panel.resize(FEM03_CHECKS_WIDTH, panel.height())
    panel.show()
    _settle()
    panel._list.setCurrentRow(0)
    _settle()
    # Measured rather than computed once: the detail is a scrolled stack whose
    # page rewraps as the panel is resized, so the fit is iterated until the
    # viewport is the height the wrapped text asks for.
    for _ in range(4):
        page = panel._detail.currentWidget()
        width = panel._detail_scroll.viewport().width()
        want = max([lbl.heightForWidth(width)
                    for lbl in page.findChildren(QLabel)] + [0])
        margins = page.layout().contentsMargins() if page.layout() else None
        if margins is not None:
            want += margins.top() + margins.bottom()
        have = panel._detail_scroll.viewport().height()
        if abs(want - have) <= 2:
            break
        panel.resize(panel.width(), panel.height() + (want - have))
        _settle()
    report = panel.report
    print("   checks      %d error(s), %d warning(s) on %s"
          % (len(report.errors), len(report.warnings),
             os.path.basename(FEM03_LINER_BONDED)))
    for f in report.warnings:
        print("   warning     %s" % f.rule_id)
    return _grab(panel, "fem03_preflight_liner.png")


SHOTS = {
    "fem03_joints_editor": fem03_joints_editor,
    "fem03_reinforce_editor": fem03_reinforce_editor,
    "fem03_run_fem": fem03_run_fem,
    "fem03_preflight_liner": fem03_preflight_liner,
}


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    sweep = "--all" in argv
    argv = [a for a in argv if not a.startswith("-")]
    # A bare run writes nothing, for the reason the tutorial capture script gives:
    # every shot is regenerated from live Qt, so a sweep rewrites every committed
    # PNG whether or not it changed.
    if not argv and not sweep:
        print(__doc__.split("Run:")[0].strip())
        print("\nName a shot (substring match), or --all to sweep.\n")
        print("shots: %s" % ", ".join(sorted(SHOTS)))
        return 1
    os.makedirs(OUT_DIR, exist_ok=True)
    names = [n for n in SHOTS if not argv or any(a in n for a in argv)]
    if not names:
        print("no shot matching %s; known shots: %s"
              % (argv, ", ".join(sorted(SHOTS))))
        return 1
    with _app_defaults():
        for name in names:
            SHOTS[name]()
    print("\nwrote %d screenshot(s) to docs/tutorials/images/" % len(names))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
