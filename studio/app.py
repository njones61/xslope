"""Application entry point for XSLOPE Studio (the ``xslope-studio`` command).

Studio runs one process per launch: there is no single-instance server, so on
Windows and Linux a double-clicked ``.xslz`` or a clicked ``xslope://`` link starts
a new copy of Studio rather than reaching one that is already running. macOS is the
exception by nature — it routes both to the running application as a
``QEvent.FileOpen`` (see :class:`StudioApplication`).
"""

from __future__ import annotations

import os

# Bind Matplotlib's Qt backend to PySide6 before any backend import.
os.environ.setdefault("QT_API", "pyside6")

import sys

from PySide6.QtCore import QEvent
from PySide6.QtGui import QIcon
from PySide6.QtWidgets import QApplication

from . import urlscheme, welcome
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


class StudioApplication(QApplication):
    """The application object, with the one platform event Studio has to answer.

    macOS does not pass a double-clicked document or a clicked ``xslope://`` link on
    the command line: it posts a ``QEvent.FileOpen`` to the running application,
    which may be a copy that has been open for hours or one the click just launched.
    Windows and Linux use ``argv`` instead (:func:`main` reads it).

    An event that arrives before the window exists — the usual case when the click
    is what launched the app — is held and replayed by :meth:`set_window`, so the
    launch and the already-running case end in the same call.
    """

    def __init__(self, argv):
        super().__init__(argv)
        self._window = None
        self._pending = []

    def set_window(self, window):
        self._window = window
        pending, self._pending = self._pending, []
        for request in pending:
            self._deliver(request)

    def event(self, event):
        if event.type() == QEvent.FileOpen:
            url = event.url()
            request = url.toString() if url.scheme() == urlscheme.SCHEME else event.file()
            if request:
                self._deliver(request)
                return True
        return super().event(event)

    def _deliver(self, request):
        if self._window is None:
            self._pending.append(request)
            return
        open_request(self._window, request)


def open_request(window, request):
    """Open whatever the OS just handed us: an ``xslope://`` link, or a file path.

    One place, because the two arrival routes (``argv`` on Windows and Linux,
    ``QEvent.FileOpen`` on macOS) carry exactly the same two kinds of thing and must
    treat them the same way. A link goes through the scheme handler, which confirms
    and checks the host before any network request; a path goes straight to the
    normal open. A path that does not exist is ignored, as it was before.

    A request that is acted on marks the window ``opened_by_request``, which is how
    :func:`main` knows the launch arrived with work in hand — here rather than at
    each caller, because this function is already the one place both arrival routes
    meet, and a request that was ignored is not work.
    """
    if urlscheme.is_scheme_url(request):
        handled = window.open_scheme_url(request)
    elif os.path.exists(request):
        handled = window.open_path(request)
    else:
        return None
    window.opened_by_request = True
    return handled


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
    app = QApplication.instance() or StudioApplication(argv)
    app.setApplicationName(APP_NAME)
    app.setOrganizationName(ORG_NAME)
    _apply_icon(app)

    win = MainWindow()
    win.setWindowIcon(app.windowIcon())
    win.show()
    if isinstance(app, StudioApplication):
        # macOS delivers a double-clicked file or an xslope:// link as an event,
        # possibly before this point; the app has been holding it for the window.
        app.set_window(win)

    # A file or an xslope:// link passed on the command line — how Windows and
    # Linux deliver both a double-clicked document and a clicked scheme link.
    rest = argv[1:]
    if rest:
        open_request(win, rest[0])

    # Did this launch arrive with work in hand? Both arrival routes have run by now
    # — the command line just above, and the events macOS delivers instead, replayed
    # by set_window before it — and both mark the window through open_request.
    with_work = bool(getattr(win, "opened_by_request", False))

    # The welcome window, on a launch whose settings still ask for one — the first
    # one, and afterwards only while the reader leaves it switched on. Deferred so
    # the main window paints behind it, and it lives HERE rather than in MainWindow
    # for the same reason the update check does: constructing a window, which every
    # Studio test does, must never raise a modal dialog. A launch that was handed a
    # document or a link is skipped whatever the preference says: the user asked for
    # that project, and a greeting over it is a dialog in the way of an answer — the
    # next bare launch still shows it. Set XSLOPE_NO_WELCOME=1 to suppress it
    # entirely (a capture, a kiosk, a classroom image).
    if (not with_work and welcome.show_at_launch(win.settings)
            and not os.environ.get("XSLOPE_NO_WELCOME")):
        from PySide6.QtCore import QTimer
        QTimer.singleShot(0, win.show_welcome)

    # The silent startup update check (Help → Check for Updates at Startup,
    # default on, at most once a day). Deferred so the window paints first, and
    # it lives HERE rather than in MainWindow so that constructing a window —
    # which every Studio test does — never touches the network. Set
    # XSLOPE_NO_UPDATE_CHECK=1 to suppress it entirely (CI, offline classrooms).
    if not os.environ.get("XSLOPE_NO_UPDATE_CHECK"):
        from PySide6.QtCore import QTimer
        QTimer.singleShot(2000, win.check_updates_at_startup)

    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
