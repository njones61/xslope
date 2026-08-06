"""Standing checks on what Studio's background threads may touch, and on how they
end.

A Qt object's C++ half belongs to the GUI thread. Release the last Python
reference to one anywhere else and shiboken does not delete it there — it queues
the deletion onto the main thread (``mainThreadDeletionHandler``) to run at the
next Python bytecode the GUI thread executes. A window sitting in ``exec()`` runs
no Python at all while it is idle, so that "next bytecode" is whatever the user's
next mouse move reaches: an event filter, halfway through Qt's own enter/leave
dispatch. Under PySide 6.11 on macOS that deferred destructor has been seen to
fault there, taking the process down with a segmentation fault seconds after a
run the user had already watched finish.

Two things have to hold for that to be impossible, and neither is visible in any
result:

  A. NOTHING QT REACHES A RUNNER. Every runner is built on the GUI thread and then
     left to itself, so whatever it carries is dropped on the worker when the run
     ends. Its options, its slope_data, its closures and its defaults are walked
     here — transitively, through dicts, lists, closure cells and bound methods —
     and a single ``QObject`` among them is a failure. (Its own Qt internals are
     not: the runner IS a QObject, and its parent is the window.)
  B. NOTHING QT IS COLLECTED ON A WORKER. A reference does not have to be handed
     over to be released in the wrong place: Python's cyclic collector runs on
     whichever thread trips its allocation threshold, and during a run that is the
     worker — it imports the engine and then allocates for seconds while the GUI
     thread is idle. Before ``runners.suspend_cycle_gc`` this was 33 of 114
     collections on the worker in one LEM run, any of which could have finalized a
     dialog and its combo boxes. This drives a real LEM run through the real window
     under a real ``exec()`` loop and asserts the collector never runs off the GUI
     thread, that no Qt object is finalized off it, and that the suspension is
     symmetric — held for the whole run, returned at the end, with the collector
     enabled again afterwards.

  C. EVERY RUNNER IS ONE. The guard lives in ``RunnerThread``, so a new runner that
     derives from ``QThread`` directly would silently opt out of it.
  D. QUITTING STOPS THE THREADS. ``QThread``'s destructor calls ``qFatal`` on a
     thread that is still running. Closing the window has always stopped them;
     quitting is not always a close — on macOS Cmd-Q ends ``exec()`` without
     sending one — and the window is then destroyed at interpreter shutdown with
     the persistent mesh thread still in its event loop, which aborts the process
     on the way out. Both paths are exercised here.

Nothing here is a solver check: one small LEM run is solved, for its allocations.

Skips cleanly (exit 0) when PySide6 is not installed.
"""
import contextlib
import gc
import io
import os
import sys
import threading
import types
import weakref

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

MAIN = threading.main_thread().name

#: A small model with a circular surface — solved once, single surface, for the
#: allocations a real run makes rather than for its answer.
MODEL = os.path.join(_REPO, "docs/inputs/slope/xslope_reinf.xlsx")


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


#: Every window these checks open, so the run can stop their threads on the way
#: out. A window left with its mesh thread running aborts the process at
#: interpreter shutdown — which is one of the things being guarded here, and would
#: otherwise fail the suite from outside any check.
_WINDOWS = []


def _window():
    """A real MainWindow with MODEL open, and the message boxes stubbed out."""
    from PySide6.QtWidgets import QApplication, QMessageBox
    QApplication.instance() or QApplication([])
    for name in ("warning", "information", "critical"):
        setattr(QMessageBox, name, staticmethod(lambda *a, **k: QMessageBox.Ok))
    from studio.main_window import MainWindow
    mw = MainWindow()
    _WINDOWS.append(mw)
    _quiet(mw.open_path, MODEL)
    return mw


# --- A. nothing Qt reaches a runner -----------------------------------------

def _qt_reachable(root, max_nodes=200000):
    """Every ``QObject`` reachable from ``root`` through Python containers, each
    with the attribute path that found it. Qt objects are not walked INTO: a
    widget's own children are Qt's business, and one widget is already a failure."""
    from PySide6.QtCore import QObject
    seen, out, stack, n = set(), [], [(root, "")], 0
    while stack and n < max_nodes:
        obj, path = stack.pop()
        n += 1
        if id(obj) in seen:
            continue
        seen.add(id(obj))
        if isinstance(obj, QObject):
            out.append((type(obj).__name__, path))
            continue
        if isinstance(obj, (str, bytes, bytearray, int, float, bool, type(None))):
            continue
        if isinstance(obj, dict):
            for k, v in list(obj.items())[:5000]:
                stack.append((v, f"{path}[{k!r}]"))
        elif isinstance(obj, (list, tuple, set, frozenset)):
            for i, v in enumerate(list(obj)[:5000]):
                stack.append((v, f"{path}[{i}]"))
        elif isinstance(obj, types.FunctionType):
            for i, cell in enumerate(obj.__closure__ or ()):
                try:
                    stack.append((cell.cell_contents, f"{path}.<cell {i}>"))
                except ValueError:
                    pass
            stack.append((obj.__defaults__ or (), f"{path}.<defaults>"))
        elif isinstance(obj, types.MethodType):
            stack.append((obj.__self__, f"{path}.__self__"))
        elif hasattr(obj, "__dict__"):
            for k, v in list(vars(obj).items()):
                stack.append((v, f"{path}.{k}"))
    return out


def _every_runner(mw):
    """One of each runner, built the way the window builds it."""
    from studio import runners
    sd = mw.doc.slope_data
    return {
        "LemRunner": runners.LemRunner(sd, {"method": "bishop"}, parent=mw),
        "SeepRunner": runners.SeepRunner(sd, {"bc": 1}, parent=mw),
        "FemRunner": runners.FemRunner(sd, {"analysis": "ssrm"}, parent=mw),
        "SensitivityRunner": runners.SensitivityRunner(
            sd, {"mode": "sensitivity"}, parent=mw),
        "ReliabilityRunner": runners.ReliabilityRunner(
            sd, {"engine": "taylor"}, parent=mw),
        "ReportRunner": runners.ReportRunner(
            sd, mw.report_solutions(), {}, os.path.join(_HERE, "_unused.docx"),
            parent=mw),
        "UpdateCheckRunner": runners.UpdateCheckRunner(parent=mw),
        "UpdateDownloadRunner": runners.UpdateDownloadRunner(
            {"filename": "x"}, _HERE, parent=mw),
        "AssistantModelsRunner": runners.AssistantModelsRunner(
            "anthropic", parent=mw),
        "MeshWorker": runners.MeshWorker(),
    }


def test_no_qt_in_runners():
    """A. No runner carries a Qt object — nor anything that reaches one."""
    failures = []
    mw = _window()
    for name, runner in _every_runner(mw).items():
        hits = []
        for attr, value in vars(runner).items():
            hits += [(tn, f"{attr}{path}") for tn, path in _qt_reachable(value)]
        if hits:
            failures.append(f"{name} carries Qt: " + ", ".join(
                f"{tn} at .{p}" for tn, p in hits[:4]))
    return failures


def test_runner_base():
    """C. Every runner derives from RunnerThread, so none can skip the guard."""
    from PySide6.QtCore import QThread
    from studio import runners
    failures = []
    for name in dir(runners):
        obj = getattr(runners, name)
        if (isinstance(obj, type) and issubclass(obj, QThread)
                and obj not in (QThread, runners.RunnerThread)):
            if not issubclass(obj, runners.RunnerThread):
                failures.append(f"{name} derives from QThread, not RunnerThread — "
                                "it would run without the collector guard")
    if not any(isinstance(getattr(runners, n), type)
               and getattr(runners, n) is not runners.RunnerThread
               and issubclass(getattr(runners, n), runners.RunnerThread)
               for n in dir(runners)):
        failures.append("no runner derives from RunnerThread at all")
    return failures


# --- B. nothing Qt is collected on a worker ---------------------------------

def test_gc_stays_on_the_gui_thread():
    """B. A real LEM run under a real exec() loop: the collector never runs off
    the GUI thread, nothing Qt is finalized off it, and the suspension is
    symmetric."""
    from PySide6.QtCore import QObject, QTimer
    from PySide6.QtWidgets import QApplication
    from studio import runners

    app = QApplication.instance() or QApplication([])
    mw = _window()
    failures = []

    passes = []            # (thread, generation)
    finalized = []         # (type name, thread)
    tracked = set()

    def _pass(phase, info):
        if phase == "start":
            passes.append((threading.current_thread().name, info.get("generation")))

    def _note(name):
        finalized.append((name, threading.current_thread().name))

    def track():
        for obj in gc.get_objects():
            if isinstance(obj, QObject) and id(obj) not in tracked:
                tracked.add(id(obj))
                try:
                    weakref.finalize(obj, _note, type(obj).__name__)
                except TypeError:
                    pass

    track()
    gc.callbacks.append(_pass)
    held = []              # was the collector suspended while the run was live?
    enabled_before = gc.isenabled()
    state = {}

    def start():
        # The real run path, minus the modal dialog: _start_lem builds and starts
        # the LemRunner exactly as run_lem does once the dialog is accepted.
        with contextlib.redirect_stdout(io.StringIO()):
            mw._start_lem({"method": "bishop", "analysis": "single_surface",
                           "surface": "circular", "num_slices": 20})
        held.append(runners.cycle_gc_suspended())
        state["started"] = mw._runner is not None
        QTimer.singleShot(50, poll)

    def poll():
        if mw._runner is not None:
            held.append(runners.cycle_gc_suspended())
            QTimer.singleShot(50, poll)
            return
        QTimer.singleShot(200, done)

    def done():
        state["after"] = runners.cycle_gc_suspended()
        state["enabled_after"] = gc.isenabled()
        app.quit()

    QTimer.singleShot(50, start)
    QTimer.singleShot(180000, app.quit)          # never hang the suite
    app.exec()
    gc.callbacks.remove(_pass)

    if not state.get("started"):
        return ["the LEM run never started — nothing was measured"]
    off = [p for p in passes if p[0] != MAIN]
    if off:
        failures.append(f"{len(off)} of {len(passes)} cyclic collections ran off "
                        f"the GUI thread ({sorted({p[0] for p in off})}) — a Qt "
                        "object finalized there is deleted from a mouse move")
    bad = [f for f in finalized if f[1] != MAIN]
    if bad:
        failures.append("Qt objects finalized off the GUI thread: "
                        + ", ".join(f"{n} on {t}" for n, t in bad[:4]))
    if not held or not all(held):
        failures.append("the cyclic collector was not suspended for the whole run "
                        f"({held.count(False)} of {len(held)} samples unheld)")
    if state.get("after"):
        failures.append("the collector was still suspended after the run — "
                        "suspend/resume is not symmetric")
    if enabled_before and not state.get("enabled_after"):
        failures.append("the collector was left disabled after the run")
    if not passes:
        failures.append("no cyclic collection happened at all — the measurement "
                        "cannot distinguish a fix from a dead instrument")
    return failures


# --- D. quitting stops the threads ------------------------------------------

def test_quit_stops_threads():
    """D. Both endings stop every thread: the application's aboutToQuit and the
    window's own close. A QThread still running when its wrapper is destroyed
    aborts the process."""
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    failures = []

    mw = _window()
    if not mw._mesh_thread.isRunning():
        return ["the mesh thread was not running — nothing was measured"]
    app.aboutToQuit.emit()
    if mw._mesh_thread.isRunning():
        failures.append("aboutToQuit left the mesh thread running — quitting with "
                        "Cmd-Q would destroy it mid-run and abort the process")
    if not hasattr(mw, "stop_threads"):
        failures.append("MainWindow has no stop_threads")
    mw.stop_threads()                    # idempotent: a second call must be safe

    mw2 = _window()
    if mw2.doc.dirty:                    # a freshly opened file is not dirty; if
        return failures                  # it were, close() would prompt
    mw2.close()
    if mw2._mesh_thread.isRunning():
        failures.append("closeEvent left the mesh thread running")
    return failures


CHECKS = [("nothing Qt reaches a runner", test_no_qt_in_runners),
          ("every runner is a RunnerThread", test_runner_base),
          ("the collector stays on the GUI thread", test_gc_stays_on_the_gui_thread),
          ("quitting stops every thread", test_quit_stops_threads)]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    try:
        import PySide6                                    # noqa: F401
    except Exception:
        print("thread safety: PySide6 not installed — checks skipped.")
        return []
    failures = []
    try:
        for name, fn in CHECKS:
            try:
                fs = fn()
            except Exception as exc:
                import traceback
                traceback.print_exc()
                fs = [f"{name} raised: {exc!r}"]
            print(f"  {name:44s} {'ok' if not fs else f'FAIL ({len(fs)})'}")
            failures += fs
    finally:
        for mw in _WINDOWS:              # see _WINDOWS
            try:
                mw.stop_threads()
            except Exception:
                pass
    return failures


def main():
    print("Studio thread safety:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll thread-safety checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication
        _app = QApplication.instance() or QApplication([])
    except Exception:
        pass
    main()
