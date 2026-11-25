"""Capture XSLOPE Studio's FEM reliability results for docs/reliability/fem.md.

    python tools/capture_reliability_fem.py [--work DIR]

The model is the FEM samples page's reliability example
(docs/fem/files/xslope_simple_mult_layers_fem.xlsx), copied with its mesh file
(xslope_simple_mult_layers_fem_mesh.json, 985 tri6 elements) into a scratch
folder (--work, else a new temporary one), so the repo's copy is never written.
Studio is driven the way a user runs it: open the workbook (Studio attaches the
mesh file), switch to FEM mode, and click Reliability... with the Taylor series,
F min 0.7, F max 1.6, and an SSRM tolerance of 0.001, the run the samples page's
REL-FEM-page tag locks. Thirteen strength-reduction solves: 8 to 16 minutes on
one process.

Studio's settings go to a scratch ini (tools/studio_settings.py) for the whole
run, so the machine's own Studio settings are left as they are.

Writes docs/reliability/images/, at Studio's native resolution:

- fem_reliability_window.png: the main window at the end of the run, the FEM
  Results view at the most-likely values over the Log pane's table (the pane
  a whole number of lines tall, so no line is cut off at its foot).
- fem_reliability_summary.png: the Reliability results box Studio shows.

and prints the summary's numbers and the Log pane's per-parameter table.
"""
import contextlib
import io
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
OUT = REPO / "docs" / "reliability" / "images"
MODEL = REPO / "docs" / "fem" / "files" / "xslope_simple_mult_layers_fem.xlsx"
MESH = MODEL.with_name(MODEL.stem + "_mesh.json")
REL = {"F_min": 0.7, "F_max": 1.6, "reliability_tol": 0.001}   # the REL-FEM-page tag's run
WINDOW = (1280, 1000)
LOG_HEIGHT = 400                     # the Log pane, tall enough for the whole table


def capture(work):
    from PySide6.QtWidgets import QApplication, QMessageBox
    app = QApplication.instance() or QApplication([])

    def settle(n=40):
        for _ in range(n):
            app.processEvents()
            time.sleep(0.02)

    def wait(done, label, timeout=3600):
        t0 = time.time()
        while not done():
            app.processEvents()
            time.sleep(0.05)
            if time.time() - t0 > timeout:
                raise TimeoutError(label)
        settle()
        return time.time() - t0

    import studio
    import xslope
    for mod in (xslope, studio):
        assert Path(mod.__file__).resolve().is_relative_to(REPO), f"{mod.__name__} from {mod.__file__}"
    print(f"xslope {xslope.__version__} from {Path(xslope.__file__).parent}")
    from studio.dialogs import ReliabilityDialog
    from studio.main_window import MainWindow, MODES

    shots = {}

    def information(parent, title, text, *a, **k):
        print(f"[{title}]\n{text}")
        if title == "Reliability results":
            box = QMessageBox(QMessageBox.Information, title, text, QMessageBox.Ok, parent)
            box.show()
            settle()
            shots["summary"] = box.grab()
            box.close()
        return QMessageBox.Ok

    def warning(parent, title, text, *a, **k):
        raise RuntimeError(f"{title}: {text}")

    QMessageBox.information = staticmethod(information)
    QMessageBox.warning = staticmethod(warning)
    QMessageBox.critical = staticmethod(warning)

    def rel_exec(dlg):
        dlg.engine.setCurrentIndex(dlg.engine.findData("taylor"))
        dlg.f_min.setValue(REL["F_min"])
        dlg.f_max.setValue(REL["F_max"])
        dlg.reliability_tol.setValue(REL["reliability_tol"])
        app.processEvents()
        shots["rel_opts"] = dlg.options()
        dlg.accept()
        return dlg.result()

    ReliabilityDialog.exec = rel_exec

    model = Path(work) / MODEL.name
    shutil.copy2(MODEL, model)
    shutil.copy2(MESH, Path(work) / MESH.name)
    mw = MainWindow()
    mw._add_recent = lambda path: None          # not a recent file
    mw.resize(*WINDOW)
    for name in ("assistant_dock", "chat_dock", "ai_dock"):
        dock = getattr(mw, name, None)
        if dock is not None:
            dock.hide()
    mw.show()
    settle()
    with contextlib.redirect_stdout(io.StringIO()):
        mw.open_path(str(model))
    settle()
    mw.set_mode_index([m for _, m in MODES].index("fem"))
    settle()
    mesh = mw.doc.slope_data.get("mesh")
    assert mesh is not None, "Studio did not attach the mesh file"
    print(f"mesh: {len(mesh['elements'])} elements, {len(mesh['nodes'])} nodes (the mesh file's)")

    n0 = len(mw.log.toPlainText())
    mw.run_reliability()
    dt = wait(lambda: mw._rel_runner is None, "reliability")
    log = mw.log.toPlainText()[n0:]
    assert "summary" in shots, "no Reliability results box"
    print(f"reliability: {dt:.0f} s")

    OUT.mkdir(parents=True, exist_ok=True)
    shots["summary"].save(str(OUT / "fem_reliability_summary.png"))
    view = mw.fem_results_canvas
    assert view is not None and mw.view_tabs.currentWidget() is view, "the FEM Results view is not showing"
    # the Log pane tall enough for the table, scrolled so the table starts at its top
    from PySide6.QtCore import Qt
    mw.resizeDocks([mw.log_dock], [LOG_HEIGHT], Qt.Vertical)
    settle()
    text = mw.log.toPlainText()
    i = text.rfind("=== RELIABILITY ANALYSIS RESULTS ===")
    assert i >= 0, "no table in the Log pane"
    line = mw.log.document().findBlock(i).blockNumber()

    def scroll():
        mw.log.verticalScrollBar().setValue(line)
        mw.log.horizontalScrollBar().setValue(0)
        settle()

    def cut_off():
        """How far the last line the Log pane shows runs past its bottom (px)."""
        h = mw.log.viewport().height()
        block = mw.log.firstVisibleBlock()
        y = mw.log.blockBoundingGeometry(block).translated(mw.log.contentOffset()).top()
        while block.isValid():
            bottom = y + mw.log.blockBoundingRect(block).height()
            if bottom > h:
                return h - y if y < h else 0.0      # the part of a line showing
            y, block = bottom, block.next()
        return 0.0
    scroll()
    # the Log pane a whole number of lines tall: no line cut off at its bottom
    for _ in range(4):
        shown = cut_off()
        if shown <= 0:
            break
        mw.resizeDocks([mw.log_dock], [mw.log_dock.height() - max(1, round(shown))],
                       Qt.Vertical)
        settle()
        scroll()
    assert cut_off() <= 0, f"the Log pane still cuts a line ({cut_off():.1f} px of it shows)"
    view._render_timer.stop()
    view._render_current()
    settle()
    mw.grab().save(str(OUT / "fem_reliability_window.png"))
    from PIL import Image
    for name in ("fem_reliability_window.png", "fem_reliability_summary.png"):
        print(f"  {name}: {Image.open(OUT / name).size}")

    r = mw.doc.results["reliability"]["reliability"]
    print("summary: " + ", ".join(f"{k} {r[k]:.4f}" for k in (
        "F_MLV", "sigma_F", "COV_F", "beta_normal", "beta_ln", "pf_normal", "prob_failure")))
    print(f"dialog: {shots['rel_opts']}")
    j = log.rfind("=== RELIABILITY ANALYSIS RESULTS ===")
    print(log[j:] if j >= 0 else log[-3000:])


def main():
    work = sys.argv[sys.argv.index("--work") + 1] if "--work" in sys.argv else tempfile.mkdtemp(prefix="rel_fem_")
    os.makedirs(work, exist_ok=True)
    from tools import studio_settings
    import matplotlib
    matplotlib.use("Agg")
    with studio_settings.redirected(prefix="docs_rel_fem") as ini:
        print(f"settings: {ini}")
        capture(work)
    sys.stdout.flush()
    os._exit(0)


if __name__ == "__main__":
    main()
