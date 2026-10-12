"""Standing check on the Log pane's font size: its starting size per platform, and
that the size the user picks is kept.

The size used to live only in the window. On Windows a file opened from Explorer,
or an Open in Studio link, starts a new copy of Studio, so every such open came
back at the default, which there is also drawn a third larger than on a Mac
(96 pixels per inch against 72). The size is now saved in Studio's settings, and
the default is 9 pt on Windows and 12 pt elsewhere.

Every window here is built against a scratch settings file (Studio's QSettings
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
        from PySide6.QtCore import QSettings
        from PySide6.QtWidgets import QApplication, QSpinBox
    except ImportError:
        print("  PySide6 not installed; skipped")
        return []
    app = QApplication.instance() or QApplication([])
    from studio import main_window

    fails = []
    want = 9 if sys.platform == "win32" else 12
    if main_window.LOG_FONT_DEFAULT_PT != want:
        fails.append(f"the default is {main_window.LOG_FONT_DEFAULT_PT} pt on "
                     f"{sys.platform}, not {want}")

    tmp = tempfile.mkdtemp(prefix="xslope-logfont-")
    ini = os.path.join(tmp, "studio.ini")
    real = main_window.QSettings
    main_window.QSettings = lambda *a, **k: QSettings(ini, QSettings.IniFormat)
    windows = []

    def window():
        mw = main_window.MainWindow()
        windows.append(mw)
        spin = mw.log_dock.titleBarWidget().findChild(QSpinBox)
        return mw, spin

    try:
        # 1. A first window starts at the platform default.
        mw, spin = window()
        if mw.log.font().pointSize() != want or spin.value() != want:
            fails.append(f"a first window starts at {mw.log.font().pointSize()} pt "
                         f"(spinner {spin.value()}), not {want}")

        # 2. A size picked in the spinner applies, is saved, and survives an open.
        spin.setValue(10)
        app.processEvents()
        if mw.log.font().pointSize() != 10:
            fails.append("the spinner does not resize the Log font")
        if int(QSettings(ini, QSettings.IniFormat).value("log_font_pt", 0)) != 10:
            fails.append("the size picked is not saved")
        mw.open_path(SAMPLE)
        app.processEvents()
        if mw.log.font().pointSize() != 10 or spin.value() != 10:
            fails.append(f"opening a file changed the Log font to "
                         f"{mw.log.font().pointSize()} pt")

        # 3. A new window (a new copy of Studio) starts at the saved size.
        mw2, spin2 = window()
        if mw2.log.font().pointSize() != 10 or spin2.value() != 10:
            fails.append(f"a new window starts at {mw2.log.font().pointSize()} pt, "
                         f"not the saved 10")

        # 4. A saved value that is not a usable size falls back or is held in range.
        store = QSettings(ini, QSettings.IniFormat)
        for saved, expect in (("garbage", want), (99, 28), (2, 6)):
            store.setValue("log_font_pt", saved)
            store.sync()
            mw2.settings.sync()
            got = mw2._saved_log_font_pt()
            if got != expect:
                fails.append(f"a saved {saved!r} gives {got} pt, not {expect}")
    finally:
        main_window.QSettings = real
        for mw in windows:
            mw.doc._dirty = False
            mw.close()
        app.processEvents()
    return fails


def main():
    print("Log pane font size:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("  starts at the platform default, keeps the size picked, and a new "
          "window starts at it")
    print("\nAll log-font checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    main()
