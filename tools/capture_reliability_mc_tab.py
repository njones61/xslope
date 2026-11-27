"""Studio's Reliability · MC tab, for docs/studio/analysis.md.

Runs Monte Carlo through Studio's own runner on the shipped reliability sample
(docs/inputs/slope/xslope_reliability.xlsx: Bishop on its circle, 2,000 Latin
hypercube realizations, a fixed seed) and grabs the result tabs offscreen:

    docs/studio/images/analysis_reliability_mc_histogram.png    (Plot sub-tab)
    docs/studio/images/analysis_reliability_mc_convergence.png  (Convergence)

Settings are sandboxed (tools.studio_settings) and the workbook is a copy.

Run from the repo root:  PYTHONPATH=. python tools/capture_reliability_mc_tab.py
"""

import contextlib
import io
import os
import shutil
import sys
import tempfile

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from tools import studio_settings                      # noqa: E402

MODEL = os.path.join(_REPO, "docs", "inputs", "slope", "xslope_reliability.xlsx")
OUT = os.path.join(_REPO, "docs", "studio", "images")
OPTS = {"app_mode": "lem", "engine": "mc", "method": "bishop", "search": False,
        "n_samples": 2000, "rng_seed": 20240117, "distribution": "normal",
        "sampling": "lhs", "converge_rel": None, "num_slices": 40}
WINDOW = (1400, 900)


def _flush():
    from PySide6.QtCore import QCoreApplication, QEvent
    from PySide6.QtWidgets import QApplication
    QCoreApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    QApplication.instance().processEvents()


def main():
    from PySide6.QtWidgets import QApplication, QMessageBox
    app = QApplication.instance() or QApplication(sys.argv)
    tmp = tempfile.mkdtemp(prefix="mc_tab_")
    try:
        with studio_settings.redirected(prefix="mc_tab_capture"):
            from studio.main_window import MainWindow
            from studio.runners import ReliabilityRunner
            path = os.path.join(tmp, os.path.basename(MODEL))
            shutil.copy(MODEL, path)
            mw = MainWindow()
            mw.resize(*WINDOW)
            mw.show()
            with contextlib.redirect_stdout(io.StringIO()):
                mw.open_path(path)
            _flush()
            got = []
            runner = ReliabilityRunner(mw.doc.slope_data, dict(OPTS))
            runner.succeeded.connect(got.append)
            runner.failed.connect(lambda m: got.append(RuntimeError(m)))
            with contextlib.redirect_stdout(io.StringIO()):
                runner.run()
            if not got or isinstance(got[0], Exception):
                raise SystemExit(f"the Monte Carlo run failed: {got}")
            real = QMessageBox.information
            QMessageBox.information = staticmethod(lambda *a, **k: None)
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    mw._on_rel_succeeded(got[0])
            finally:
                QMessageBox.information = real
            _flush()
            view = mw.reliability_hist_canvas
            mw.view_tabs.setCurrentWidget(view)
            for label, name in (("Plot", "analysis_reliability_mc_histogram.png"),
                                ("Convergence", "analysis_reliability_mc_convergence.png")):
                idx = [view.tabs.tabText(i) for i in range(view.tabs.count())].index(label)
                view.tabs.setCurrentIndex(idx)
                _flush()
                canvas = view.tabs.widget(idx)
                if hasattr(canvas, "render_now"):
                    canvas.render_now()
                _flush()
                out = os.path.join(OUT, name)
                mw.view_tabs.grab().save(out)
                print("wrote", os.path.relpath(out, _REPO))
            rel = got[0]["reliability"]
            print(f"Pf {rel['pf_empirical'] * 100:.2f}%")
            mw.doc._dirty = False
            mw.close()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    del app


if __name__ == "__main__":
    main()
