"""Standing check on dropping a file onto the Studio window: an input workbook or a
project package dropped anywhere on the window opens as File → Open would, and
anything else is refused.

The canvases and the text panes take drags of their own, so the window answers
drags through an application-wide event filter. Each case here sends the drag to
one of those child widgets, which is where a real drop lands.

Every window is built against a scratch settings file (Studio's QSettings
constructor is replaced for the duration), so nothing of the user's is read or
written. Skips cleanly (exit 0) when PySide6 is not installed.
"""
import os
import sys
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("XSLOPE_NO_UPDATE_CHECK", "1")

SAMPLE = os.path.join(_REPO, "docs", "inputs", "slope", "xslope_simple1.xlsx")


def run():
    try:
        from PySide6.QtCore import QEvent, QMimeData, QPoint, QPointF, QSettings, Qt, QUrl
        from PySide6.QtGui import QDragEnterEvent, QDropEvent
        from PySide6.QtWidgets import QApplication, QDialog
    except ImportError:
        print("  PySide6 not installed; skipped")
        return []
    app = QApplication.instance() or QApplication([])
    from studio import main_window

    tmp = tempfile.mkdtemp(prefix="xslope-drop-")
    ini = os.path.join(tmp, "studio.ini")
    package = os.path.join(tmp, "project.xslz")
    text = os.path.join(tmp, "notes.txt")
    for path in (package, text):
        open(path, "w").close()

    def mime(*paths):
        data = QMimeData()
        data.setUrls([QUrl.fromLocalFile(p) for p in paths])
        return data

    def enter(widget, data):
        ev = QDragEnterEvent(QPoint(5, 5), Qt.CopyAction, data,
                             Qt.LeftButton, Qt.NoModifier)
        QApplication.sendEvent(widget, ev)
        return ev.isAccepted()

    def drop(widget, data):
        # Qt hands a drop only to the widget that took the drag as it entered,
        # so a drop is sent the way a real one arrives: the drag first.
        enter(widget, data)
        ev = QDropEvent(QPointF(5, 5), Qt.CopyAction, data, Qt.LeftButton,
                        Qt.NoModifier, QEvent.Drop)
        QApplication.sendEvent(widget, ev)
        app.processEvents()             # the open runs after the drop returns
        return ev.isAccepted()

    fails = []
    real = main_window.QSettings
    main_window.QSettings = lambda *a, **k: QSettings(ini, QSettings.IniFormat)
    mw = main_window.MainWindow()
    dialog = QDialog()
    dialog.setAcceptDrops(True)
    try:
        mw.show()
        app.processEvents()
        canvas = mw.canvas.view.viewport()
        targets = {"the Inputs canvas": canvas, "the Log pane": mw.log.viewport(),
                   "the window itself": mw}

        # 1. A workbook or a package is taken wherever on the window it lands.
        for where, widget in targets.items():
            for path in (SAMPLE, package):
                if not enter(widget, mime(path)):
                    fails.append(f"a {os.path.splitext(path)[1]} dragged over "
                                 f"{where} is refused")

        # 2. Other drops are left alone: another kind of file, several files, a
        #    file that does not exist, and a drop on a window other than this one.
        #    (The canvas's own scene takes any drag as it enters, as Qt's does, so
        #    what is checked is that nothing opens.)
        for what, widget, data in (
                ("a .txt file", canvas, mime(text)),
                ("two workbooks", canvas, mime(SAMPLE, SAMPLE)),
                ("a missing workbook", canvas, mime(os.path.join(tmp, "gone.xlsx"))),
                ("a workbook on another window", dialog, mime(SAMPLE))):
            drop(widget, data)
            if mw.doc.path is not None:
                fails.append(f"{what} dropped opened {mw.doc.path}")
                mw.doc.path = None

        # 3. A workbook dropped on the canvas opens.
        if not drop(canvas, mime(SAMPLE)):
            fails.append("the drop on the canvas is not accepted")
        elif os.path.normcase(str(mw.doc.path)) != os.path.normcase(SAMPLE):
            fails.append(f"the dropped workbook did not open (document: {mw.doc.path})")
    finally:
        main_window.QSettings = real
        mw.doc._dirty = False
        mw.close()
        dialog.close()
        app.processEvents()
    return fails


def main():
    print("Dropping a file on the Studio window:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("  a workbook or package is taken anywhere on the window and opens; other "
          "drags are left alone")
    print("\nAll drop-to-open checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    main()
