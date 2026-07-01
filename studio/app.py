"""Application entry point for XSlope Studio (the ``xslope-studio`` command)."""

from __future__ import annotations

import os

# Bind Matplotlib's Qt backend to PySide6 before any backend import.
os.environ.setdefault("QT_API", "pyside6")

import sys

from PySide6.QtGui import QIcon
from PySide6.QtWidgets import QApplication

from .main_window import APP_NAME, ORG_NAME, MainWindow

ICON_PATH = os.path.join(os.path.dirname(__file__), "resources", "icon.png")


def _apply_icon(app):
    """Give the app a custom icon across platforms. ``setWindowIcon`` covers the
    window title bar (and the Windows/Linux taskbar); the Dock on macOS and the
    taskbar grouping on Windows need a couple of extra platform-specific nudges."""
    if not os.path.exists(ICON_PATH):
        return
    app.setWindowIcon(QIcon(ICON_PATH))
    if sys.platform == "darwin":
        # setWindowIcon doesn't drive the Dock tile for a non-bundled process —
        # set it natively. Best-effort; needs pyobjc (the gui extra pulls it on macOS).
        try:
            from AppKit import NSApplication, NSImage
            img = NSImage.alloc().initWithContentsOfFile_(ICON_PATH)
            if img is not None:
                NSApplication.sharedApplication().setApplicationIconImage_(img)
        except Exception:
            pass
    elif sys.platform.startswith("win"):
        # Tie the taskbar grouping to our own AppUserModelID so Windows uses the
        # window icon rather than the python launcher's.
        try:
            import ctypes
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("xslope.studio")
        except Exception:
            pass


def _silence_aspect_limit_noise():
    """Drop Matplotlib's "Ignoring fixed x/y limits to fulfill fixed data aspect
    with adjustable data limits" warning. Many engine plots intentionally use
    ``set_aspect('equal', adjustable='datalim')`` with explicit limits; Matplotlib
    logs that warning (via ``matplotlib.axes._base``) on *every* redraw, which
    floods Studio's Log pane. It's purely informational and not actionable by the
    user, so filter just that message (leaving all other Matplotlib warnings)."""
    import logging

    class _DropAspectLimitNoise(logging.Filter):
        def filter(self, record):
            return "to fulfill fixed data aspect" not in record.getMessage()

    logging.getLogger("matplotlib.axes._base").addFilter(_DropAspectLimitNoise())


def main(argv=None):
    _silence_aspect_limit_noise()
    argv = list(sys.argv if argv is None else argv)
    app = QApplication.instance() or QApplication(argv)
    app.setApplicationName(APP_NAME)
    app.setOrganizationName(ORG_NAME)
    _apply_icon(app)

    win = MainWindow()
    win.setWindowIcon(app.windowIcon())
    win.show()

    # Optionally open a file passed on the command line.
    rest = argv[1:]
    if rest and os.path.exists(rest[0]):
        win.open_path(rest[0])

    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
