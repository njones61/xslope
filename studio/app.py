"""Application entry point for XSlope Studio (the ``xslope-studio`` command)."""

from __future__ import annotations

import os

# Bind Matplotlib's Qt backend to PySide6 before any backend import.
os.environ.setdefault("QT_API", "pyside6")

import sys

from PySide6.QtWidgets import QApplication

from .main_window import APP_NAME, ORG_NAME, MainWindow


def main(argv=None):
    argv = list(sys.argv if argv is None else argv)
    app = QApplication.instance() or QApplication(argv)
    app.setApplicationName(APP_NAME)
    app.setOrganizationName(ORG_NAME)

    win = MainWindow()
    win.show()

    # Optionally open a file passed on the command line.
    rest = argv[1:]
    if rest and os.path.exists(rest[0]):
        win.open_path(rest[0])

    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
