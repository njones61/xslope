"""Studio's machine-wide settings store, replaced by a scratch file for the
length of a check.

Studio keeps its preferences — the recent-files list, the editor toggles, the
update-check stamp, the assistant's provider and model — in one machine-wide
store, opened as ``QSettings("XSlope", "XSlope Studio")`` from a dozen call sites,
none of which take a path. A check that builds a real ``MainWindow`` and hands it
a document goes through the real open path, and the real open path appends to
that list: run the suite and the person's own recent files are reordered by it.

**A store that cannot be redirected has to be replaced.** On macOS the
two-argument form is served by ``cfprefsd``, which resolves the domain for the
logged-in user whatever ``HOME`` says, and comes back NativeFormat at the real
plist whatever ``setDefaultFormat`` / ``setPath`` say (measured — see
``test/welcome_window_check.py``). So :class:`StoreRedirector` stands in for the
*class*, hands out a scratch ini wherever the machine store is asked for, and
passes every other construction through untouched.

The scope is one check, in one process. :func:`redirected` installs the stand-in
on entry and puts the real class back on exit, so a check that runs beside
``test/assistant_capture_check.py`` in the same suite worker leaves that check's
own before/after comparison of the real store reading the real store — which is
the guard that found this. The ini is named for the pid, so parallel workers
never share one.

Used by every check that opens a ``MainWindow``; the recorded-assistant-session
harness (``tools/assistant_sessions.py``) builds its own pinned store on the same
mechanism.
"""
from __future__ import annotations

import contextlib
import functools
import os
import sys
import tempfile

#: The organization and application Studio opens its store under.
ORG_NAME = "XSlope"
SETTINGS_APP = "XSlope Studio"


class StoreRedirector:
    """Stands in for ``QSettings`` and hands out a scratch store instead.

    Every construction that asks for Studio's machine-wide store — ``QSettings()``
    or any call carrying ``"XSlope"`` — comes back as the scratch ini. Every other
    construction (a path with an explicit format) is passed through untouched, and
    so is every attribute (``IniFormat``, ``UserScope``, …), so the name behaves as
    the class everywhere else.
    """

    def __init__(self, real, path):
        self._real = real
        self._path = path

    def __call__(self, *args, **kwargs):
        if self._is_machine_store(args):
            return self._real(self._path, self._real.IniFormat)
        return self._real(*args, **kwargs)

    @staticmethod
    def _is_machine_store(args):
        """Whether this construction asks for Studio's own machine-wide store."""
        if not args:
            return True          # QSettings() — the application-wide default store
        return any(a == ORG_NAME for a in args if isinstance(a, str))

    def __getattr__(self, name):
        return getattr(self._real, name)


def patch_targets():
    """Every already-imported module holding its own reference to ``QSettings``.

    Studio's modules bind the class at import time (``from PySide6.QtCore import
    QSettings``), so patching ``PySide6.QtCore`` alone reaches only the modules
    imported afterwards. Both are done: this sweep catches the ones already in,
    and the QtCore patch catches every later import (``studio.editors`` and
    ``studio.welcome`` are loaded lazily, mid-check).
    """
    import PySide6.QtCore
    real = PySide6.QtCore.QSettings
    targets = [(PySide6.QtCore, real)]
    for name, mod in list(sys.modules.items()):
        if mod is None or not (name == "studio" or name.startswith("studio.")):
            continue
        if getattr(mod, "QSettings", None) is real:
            targets.append((mod, real))
    return targets


def scratch_settings_path(prefix="studio_check", env=None):
    """The ini this process writes instead of the user's store.

    One file per process, named for the pid, so two workers running Studio checks
    at the same time cannot see or overwrite each other's. ``env`` names an
    environment variable that overrides the location outright, for a caller that
    has to hand the path to a second process.
    """
    override = os.environ.get(env) if env else None
    if override:
        os.makedirs(os.path.dirname(os.path.abspath(override)) or ".",
                    exist_ok=True)
        return os.path.abspath(override)
    root = os.path.join(tempfile.gettempdir(), "xslope_studio_settings")
    os.makedirs(root, exist_ok=True)
    return os.path.join(root, "%s_%d.ini" % (prefix, os.getpid()))


@contextlib.contextmanager
def redirected(path=None, prefix="studio_check"):
    """Every Studio store opened inside the block lands in a scratch ini.

    Yields the path. Restores the real ``QSettings`` on the way out, including on
    modules imported *during* the block: those bound the stand-in as their own
    ``QSettings`` and are not in the list taken on entry, so the sweep is repeated
    on exit for exactly that set. Without it a lazily-imported panel would keep
    sending the machine store to this ini for the rest of the process — and the
    next check in the same worker would be measuring a scratch file.

    A no-op when PySide6 is absent (an engine-only install has no Studio layer to
    isolate), so a check can wrap its whole body without a second code path.
    """
    try:
        import PySide6.QtCore                                     # noqa: F401
    except Exception:
        yield None
        return

    path = path or scratch_settings_path(prefix)
    targets = patch_targets()
    real_class = targets[0][1]
    for module, real in targets:
        setattr(module, "QSettings", StoreRedirector(real, path))
    try:
        yield path
    finally:
        for module, real in targets:
            setattr(module, "QSettings", real)
        for name, mod in list(sys.modules.items()):
            if mod is None or not (name == "studio" or name.startswith("studio.")):
                continue
            if isinstance(getattr(mod, "QSettings", None), StoreRedirector):
                setattr(mod, "QSettings", real_class)


def isolated(fn):
    """Decorator form of :func:`redirected`, for a check's entry point.

    A check's ``run()`` is what the suite calls, and it is the scope the isolation
    belongs to: the store is replaced for every leg the check runs and put back
    before the next check in the same worker starts. Decorating the entry point
    rather than each ``MainWindow()`` call site is what makes it hard to add a leg
    that escapes — a new window opened anywhere under ``run()`` is already inside.
    """
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        with redirected():
            return fn(*args, **kwargs)
    return wrapper
