"""Recorded assistant sessions — the capture harness behind the AI-assistant tutorial.

Every other tutorial figure photographs a dialog the producer drives itself, so the
shot is reproducible by construction. A tutorial about the assistant cannot work
that way: its subject is what the model says and runs, which only a real turn
produces. So this module runs the real thing — the real offscreen ``MainWindow``,
the real ``Assistant``, the real provider call — and records both halves of it:

  * a PNG of the Assistant dock after each turn, grown so the whole exchange is
    inside the frame rather than scrolled out of it, and
  * a plain-text transcript carrying the user's words, the assistant's words,
    every ``run_python`` snippet with its printed result and MODEL CHECKS block,
    and the token usage the turn cost.

A session that changed the model also writes the workbook it left behind, so the
tutorial can hand the reader the same file to open.

The turns cost real money against the key in the OS keychain, so nothing here runs
by accident: sessions live in their own registry, are never swept by an unnamed
capture run, and ``--dry-run`` exercises the whole path — window, dock, grab,
transcript, workbook — with a stub reply and no provider call at all.

Run:  python3 tools/capture_tutorial_screenshots.py w1 --dry-run   # plumbing only
      python3 tools/capture_tutorial_screenshots.py w1             # real turns
"""

from __future__ import annotations

import contextlib
import io
import os
import sys
import threading
import time

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import matplotlib
matplotlib.use("Agg")

from PySide6.QtCore import QEventLoop, QSettings, Qt, QTimer
from PySide6.QtGui import QImage
from PySide6.QtWidgets import QApplication, QMessageBox

IMAGES_DIR = os.path.join(REPO_ROOT, "docs", "tutorials", "images")
FILES_DIR = os.path.join(REPO_ROOT, "docs", "tutorials", "files")

#: Prefix for everything this harness writes. The tutorial is W-1 ("working with
#: the assistant"), and the prefix keeps its recorded sessions together in a
#: directory shared with every other tutorial's figures.
PREFIX = "w1"

_app = QApplication.instance() or QApplication([])
# Modal dialogs must never block a headless run (the same guard the screenshot
# producer installs — this module is importable on its own, so it repeats it).
QMessageBox.warning = staticmethod(lambda *a, **k: QMessageBox.Ok)
QMessageBox.information = staticmethod(lambda *a, **k: QMessageBox.Ok)
QMessageBox.critical = staticmethod(lambda *a, **k: QMessageBox.Ok)

#: Registry of recorded sessions, kept SEPARATE from the screenshot producer's
#: SHOTS so that a full capture sweep (no arguments) cannot spend API credit.
#: ``capture_tutorial_screenshots.main`` runs these only when named.
SESSIONS = {}


def _settle(cycles=12):
    """Pump the event loop so deferred layout and the canvas's debounced render fire."""
    for _ in range(cycles):
        _app.processEvents()
        time.sleep(0.02)


# --------------------------------------------------------------------------- #
# Token accounting
# --------------------------------------------------------------------------- #
#: Usage fields worth recording, in the order they are printed. The first three
#: are the Assistant's own accumulator — the same split the dock shows the reader,
#: with the prompt-cached share named INSIDE the input count rather than beside it.
#: The rest are LiteLLM's raw field names, used when a build has no accumulator and
#: the harness falls back to its own probe.
_USAGE_FIELDS = ("input", "cached_input", "output",
                 "prompt_tokens", "completion_tokens", "total_tokens",
                 "cache_creation_input_tokens", "cache_read_input_tokens",
                 "reasoning_tokens", "calls")


def _assistant_usage(assistant, scope):
    """One scope of the Assistant's own token accumulator, or ``None``.

    The dock's readout is the number to record: it is what the reader sees beside
    the conversation, and it is already split the way the tutorial talks about cost
    — ``input`` (with the cached part named inside it) and ``output``. ``turn``
    counts the completions since the last send, so a turn's cost is read straight
    off it with no differencing. Read through ``getattr`` because a build without
    the accumulator must still record a session, on the probe below.
    """
    usage = getattr(assistant, "usage", None)
    if not isinstance(usage, dict) or scope not in usage:
        return None
    return {k: v for k, v in (usage[scope] or {}).items()
            if isinstance(v, (int, float)) and not isinstance(v, bool)}


def _numeric_fields(obj):
    """The numeric public attributes of a usage object / mapping, flattened."""
    if obj is None:
        return {}
    if isinstance(obj, dict):
        items = obj.items()
    else:
        items = [(k, getattr(obj, k, None))
                 for k in dir(obj) if not k.startswith("_")]
    out = {}
    for key, value in items:
        if isinstance(value, bool) or callable(value):
            continue
        if isinstance(value, (int, float)):
            out[key] = value
        elif value is not None and key.endswith("_details"):
            for sub, subval in _numeric_fields(value).items():
                out.setdefault(sub, subval)
    return out


def _usage_delta(before, after):
    """What one turn added to a running accumulator."""
    return {k: v - before.get(k, 0) for k, v in after.items()
            if v - before.get(k, 0)}


def _format_usage(usage):
    if not usage:
        return "not recorded"
    ordered = [k for k in _USAGE_FIELDS if k in usage]
    ordered += [k for k in sorted(usage) if k not in _USAGE_FIELDS]
    return " · ".join("%s %s" % (k, format(usage[k], ",")
                                 if isinstance(usage[k], int) else usage[k])
                      for k in ordered)


@contextlib.contextmanager
def _usage_probe():
    """Accumulate token usage across the session's provider calls.

    The Assistant grew its own ``usage`` accumulator; this probe is the fallback
    for a build that does not have one yet, and it also cross-checks the one that
    does. It wraps ``litellm.completion`` — the single call the agent worker makes
    — rather than reaching into the assistant, so nothing in ``studio/ai`` has to
    cooperate. The worker calls it from its own thread, hence the lock.
    """
    try:
        import litellm
    except Exception:
        yield None
        return
    totals = {}
    lock = threading.Lock()
    real = litellm.completion

    def wrapped(*args, **kwargs):
        resp = real(*args, **kwargs)
        try:
            fields = _numeric_fields(getattr(resp, "usage", None))
        except Exception:
            fields = {}
        with lock:
            totals["calls"] = totals.get("calls", 0) + 1
            for key, value in fields.items():
                totals[key] = totals.get(key, 0) + value
        return resp

    litellm.completion = wrapped
    try:
        yield totals
    finally:
        litellm.completion = real


# --------------------------------------------------------------------------- #
# Settings pinning
# --------------------------------------------------------------------------- #
#: The organization / application pair Studio's own store is opened under —
#: ``QSettings("XSlope", "XSlope Studio")``, constructed in a dozen places across
#: the GUI. A call carrying this pair is the machine-wide store, and is what the
#: redirector below sends to the session's own ini instead.
ORG_NAME = "XSlope"
SETTINGS_APP = "XSlope Studio"

#: Environment override for the session store's path. The suite sets it so every
#: session of one run shares a store; unset, each PROCESS gets its own file.
SETTINGS_ENV = "XSLOPE_SESSION_SETTINGS"

#: Settings copied from the machine store into the session's own, so a session
#: still runs against the endpoint the person configured (a custom Ollama base,
#: a Z.ai coding-plan URL, the model list the dialog last fetched). Prefix match.
#: API KEYS ARE NOT COPIED: they live in the OS keychain, which the session reads
#: directly and unchanged, and the QSettings fallback copy of one would be a
#: secret written into a temp file for no gain.
_CARRIED_PREFIXES = ("ai/",)
_NEVER_CARRIED = ("ai/key_fallback/",)


def _session_settings_path():
    """The ini this process's sessions read and write — never the user's store.

    One file per process, named for the pid, so two harness processes recording
    at the same time cannot see or overwrite each other's pins. The old failure
    was exactly that: the pins went into the machine-wide store and were restored
    at the end of each session, so a second session starting mid-first one read
    (or restored) the first one's values — and a run whose ``ai/confirm`` had been
    put back to True stopped dead on the "Run code?" modal with nobody to click it.
    """
    override = os.environ.get(SETTINGS_ENV)
    if override:
        os.makedirs(os.path.dirname(os.path.abspath(override)) or ".",
                    exist_ok=True)
        return os.path.abspath(override)
    import tempfile
    root = os.path.join(tempfile.gettempdir(), "xslope_assistant_sessions")
    os.makedirs(root, exist_ok=True)
    return os.path.join(root, "settings_%d.ini" % os.getpid())


class _StoreRedirector:
    """Stands in for ``QSettings`` and hands out the session's store instead.

    Studio opens its store as ``QSettings("XSlope", "XSlope Studio")`` from a
    dozen call sites, none of which take a path, and on macOS that store CANNOT be
    redirected by ``setDefaultFormat``/``setPath`` (measured — see
    ``test/welcome_window_check.py``): the two-argument form comes back
    NativeFormat at the real plist whatever those say. A store that cannot be
    redirected has to be REPLACED, so the *class* is what this replaces, for the
    length of the session and inside this process only.

    Every other construction — a scratch ini, a path with an explicit format —
    is passed through untouched, and so is every attribute (``IniFormat``,
    ``UserScope``, …), so the name behaves as the class everywhere else.
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


def _patch_targets():
    """Every already-imported module holding its own reference to ``QSettings``.

    Studio's modules bind the class at import time (``from PySide6.QtCore import
    QSettings``), so patching ``PySide6.QtCore`` alone reaches only the modules
    imported afterwards. Both are done: this sweep catches the ones already in,
    and the QtCore patch catches every later import (``studio.editors`` and
    ``studio.welcome`` are loaded lazily, mid-session).
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


@contextlib.contextmanager
def _pinned_settings(values):
    """Run the session against a store of its own, carrying ``values``.

    The dock's caption and the completion call both read the stored provider /
    model selection, and the confirm-before-running gate decides whether an
    unattended turn stops on a modal — so a session has to state all three. It
    used to state them by writing the machine-wide store and putting it back
    afterwards, which made two concurrent sessions incompatible (see
    :func:`_session_settings_path`) and put the person's own preferences one
    crashed run away from being left on a recording's values.

    Now nothing of theirs is written or restored: the session gets a private ini,
    the class Studio constructs its store from is replaced for the length of the
    run so every part of the app lands in that ini, and the machine store is only
    ever READ — once, to carry the non-secret ``ai/`` settings across so the
    session still talks to the endpoint the person configured. API keys are not
    copied; they come from the OS keychain, which is untouched by any of this.
    """
    path = _session_settings_path()
    store = QSettings(path, QSettings.IniFormat)
    machine = QSettings(ORG_NAME, SETTINGS_APP)
    for key in machine.allKeys():
        if (key.startswith(_CARRIED_PREFIXES)
                and not key.startswith(_NEVER_CARRIED)
                and key not in values):
            store.setValue(key, machine.value(key))
    for key, value in values.items():
        store.setValue(key, value)
    store.sync()

    targets = _patch_targets()
    real_class = targets[0][1]
    for module, real in targets:
        setattr(module, "QSettings", _StoreRedirector(real, path))
    try:
        yield store
    finally:
        for module, real in targets:
            setattr(module, "QSettings", real)
        # A module imported DURING the session bound the redirector as its own
        # ``QSettings`` and is not in the list above; left alone it would keep
        # sending the machine store to this session's ini for the rest of the
        # process. The sweep is repeated on the way out for exactly that set.
        for name, mod in list(sys.modules.items()):
            if mod is None or not (name == "studio" or name.startswith("studio.")):
                continue
            if isinstance(getattr(mod, "QSettings", None), _StoreRedirector):
                setattr(mod, "QSettings", real_class)
        store.sync()


# --------------------------------------------------------------------------- #
# The watchdog
# --------------------------------------------------------------------------- #
#: Seconds a turn is allowed to overrun its own timeout before the watchdog
#: stops the process. The per-turn timeout is a QTimer, and a QTimer is a message
#: on the event loop — so it cannot fire while the event loop is the thing that is
#: stuck, which is exactly what a runaway snippet does: ``run_python`` executes on
#: the GUI thread, and a snippet that never returns takes the timer with it. The
#: kernel's own snippet limit (``ai/run_timeout``, pinned below) is what normally
#: ends such a run; this is the backstop for the case that limit cannot reach — a
#: C-level call between two bytecodes, where no Python-level interrupt lands.
WATCHDOG_GRACE_S = 120.0


@contextlib.contextmanager
def _watchdog(seconds, message, note=None):
    """Kill the process if the block has not finished ``seconds`` from now.

    Runs on its own thread, which keeps running while the GUI thread is blocked
    (a busy Python loop yields the GIL; a C call releases it), so it is the one
    part of the harness that still works when the event loop is gone. It writes
    the reason where the run's record is — the transcript, through ``note`` — and
    then ends the process hard, because there is nothing left to unwind through:
    the stack it would have to return along is the one that is stuck.
    """
    if not seconds or seconds <= 0:
        yield
        return
    done = threading.Event()

    def watch():
        if done.wait(seconds):
            return
        try:
            if note is not None:
                note(message)
        except Exception:
            pass
        sys.stderr.write("  ! %s\n" % message)
        sys.stderr.flush()
        os._exit(3)

    thread = threading.Thread(target=watch, daemon=True,
                              name="assistant-session-watchdog")
    thread.start()
    try:
        yield
    finally:
        done.set()


# --------------------------------------------------------------------------- #
# The dock: sizing and grabbing
# --------------------------------------------------------------------------- #
def _fit_transcript(win, max_height):
    """Size the window so the dock's viewport is exactly the transcript's height.

    The dock is a column of the main window, so the transcript's height is the
    window's height minus everything else; the way to show a long exchange whole is
    to make the window taller, not to resize a widget the layout will overrule. It
    works the other way too — a short first turn in a window-height dock is a figure
    that is mostly empty box — so the fit runs in both directions.

    Every number is measured off the laid-out document rather than guessed, and the
    measurement is repeated because each resize changes the wrapping (a viewport
    that no longer needs a scrollbar is wider, so the same text takes fewer lines).
    The bounds are the window's own minimum and the caller's cap; an exchange too
    long for the cap keeps the cap, and the caller is told rather than being handed
    a figure no page can carry.
    """
    view = win.chat_dock.widget().transcript
    floor = win.minimumSizeHint().height()
    clamped = False
    for _ in range(5):
        _settle(4)
        doc = view.document()
        doc.setTextWidth(view.viewport().width())
        wanted = int(doc.size().height() + doc.documentMargin() * 2) + 4
        deficit = wanted - view.viewport().height()
        if abs(deficit) <= 2:
            break
        target = win.height() + deficit
        if target > max_height:
            target = max_height
            clamped = True
        target = max(target, floor)
        if target == win.height():
            break
        win.resize(win.width(), target)
    _settle()
    bar = view.verticalScrollBar()
    bar.setValue(bar.maximum())      # a clamped transcript ends on its last turn
    _settle(4)
    return clamped


def _grab_dock(win, path, max_height):
    clamped = _fit_transcript(win, max_height)
    pix = win.chat_dock.grab()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    pix.save(path)
    print("   -> %s  (%dx%d%s)" % (os.path.basename(path), pix.width(),
                                   pix.height(),
                                   ", clamped" if clamped else ""))
    return path


# --------------------------------------------------------------------------- #
# One turn
# --------------------------------------------------------------------------- #
class _Recorder:
    """Collects one turn's assistant output off the Assistant's signals."""

    def __init__(self, assistant):
        self.assistant = assistant
        self.lines = []
        self.error = None
        self.done = False
        self._conns = [
            (assistant.assistant_text, self._on_text),
            (assistant.tool_ran, self._on_tool_ran),
            (assistant.tool_declined, self._on_declined),
            (assistant.failed, self._on_failed),
            (assistant.finished, self._on_finished),
        ]
        for signal, slot in self._conns:
            signal.connect(slot)

    def disconnect(self):
        for signal, slot in self._conns:
            with contextlib.suppress(RuntimeError, TypeError):
                signal.disconnect(slot)

    # -- signal slots --
    def _on_text(self, text):
        self.lines.append("Assistant: " + (text or "").strip())

    def _on_tool_ran(self, code, output, outputs):
        self.lines.append("Ran code:\n" + _indent(code))
        self.lines.append("Output:\n" + _indent(output or "(no output)"))
        for path in outputs or []:
            self.lines.append("[file: %s]" % os.path.basename(path))

    def _on_declined(self, code):
        self.lines.append("Declined to run the code.")

    def _on_failed(self, message):
        self.error = message
        self.lines.append("Error: " + message)
        self.done = True

    def _on_finished(self):
        self.done = True


def _indent(text, pad="    "):
    return "\n".join(pad + line for line in (text or "").rstrip().splitlines())


def _wait(recorder, timeout_s):
    """Block on the event loop until the turn finishes, fails, or times out."""
    if recorder.done:
        return False
    loop = QEventLoop()
    timed_out = {"v": False}

    def quit_loop():
        if loop.isRunning():
            loop.quit()

    def on_timeout():
        timed_out["v"] = True
        quit_loop()

    recorder.assistant.finished.connect(quit_loop)
    recorder.assistant.failed.connect(quit_loop)
    timer = QTimer()
    timer.setSingleShot(True)
    timer.timeout.connect(on_timeout)
    timer.start(int(timeout_s * 1000))
    try:
        loop.exec()
    finally:
        timer.stop()
        with contextlib.suppress(RuntimeError, TypeError):
            recorder.assistant.finished.disconnect(quit_loop)
        with contextlib.suppress(RuntimeError, TypeError):
            recorder.assistant.failed.disconnect(quit_loop)
    return timed_out["v"]


def _install_stub(assistant, reply="…"):
    """Replace the provider call with an immediate canned reply (``--dry-run``).

    Bound on the instance, so ``studio/ai`` is untouched and the ChatDock's own
    send path — input box, user block, busy state, finished signal — is the one
    exercised. The reply is emitted from the event loop rather than inline so the
    caller's wait behaves exactly as it does against a real turn.
    """
    def send(user_text, images=None):
        assistant._messages.append({"role": "user", "content": user_text})

        def reply_now():
            assistant.assistant_text.emit(reply)
            assistant.finished.emit()
        QTimer.singleShot(0, reply_now)

    assistant.send = send


# --------------------------------------------------------------------------- #
# The session
# --------------------------------------------------------------------------- #
def run_assistant_session(name, model_path, turns, *, provider="anthropic",
                          model="claude-opus-5", images=None,
                          out_dir=IMAGES_DIR, files_dir=FILES_DIR,
                          timeout_s=600, dry_run=False, dock_width=560,
                          window_size=(1500, 950), max_height=6400,
                          save_after=None, save_each_turn=False, settings=None,
                          prefix=None):
    """Record one assistant conversation against a real model, offscreen.

    Environment overrides, for playing the same sessions against another model
    without touching the tutorial's own artifacts: ``XSLOPE_SESSIONS_PROVIDER``,
    ``XSLOPE_SESSIONS_MODEL`` replace ``provider`` / ``model``;
    ``XSLOPE_SESSIONS_OUT`` is a directory that receives BOTH the images and the
    saved workbooks / transcripts in place of ``out_dir`` and ``files_dir``.

    Opens ``model_path`` in an offscreen ``MainWindow``, pins ``provider`` /
    ``model`` and switches the confirm-before-running gate off for the length of
    the run (there is nobody to click Run), then plays ``turns`` through the chat
    dock's own send path, one at a time, waiting for each to finish before the
    next is typed.

    Parameters
    ----------
    name : str
        Session name. Every artifact is named from it and from nothing else, so a
        re-run overwrites its own files: ``w1_{name}_{i}.png``,
        ``w1_{name}_transcript.md``, ``w1_{name}_after.xlsx``.
    model_path : str
        The .xlsx the session opens. ``None`` starts from an empty project.
    turns : list
        Each item is either the prompt text, or ``(text, image_path)`` to attach a
        picture to that turn the way a reader would paste one into the dock.
    images : list, optional
        Images for the FIRST turn, for the common single-turn case where writing a
        one-item tuple list would be noise. A turn's own images win.
    timeout_s : float
        Per turn. A turn that overruns is recorded as an error and the session
        stops there rather than hanging a producer run. It is also the limit the
        kernel is pinned to for ONE SNIPPET (``ai/run_timeout``), because the
        per-turn timeout is a QTimer and a snippet runs on the GUI thread: a
        runaway snippet stops the timer that was supposed to end it. Past that
        limit and a grace period, the watchdog thread notes the turn in the
        transcript and stops the process (:data:`WATCHDOG_GRACE_S`).
    dry_run : bool
        Do everything except call the provider: the reply is a stub "…", so the
        whole path — window, dock sizing, grab, transcript, workbook — is
        exercised without spending anything.
    max_height : int
        Ceiling on the window height the fit is allowed to reach. A long
        conversation whose dock does not fit inside it is grabbed showing its END,
        and the run says ``clamped`` so the figure is not mistaken for the whole
        exchange.
    save_after : bool, optional
        Whether to write the workbook the session left behind. Default (``None``)
        writes it exactly when the session actually changed the model.
    settings : dict, optional
        Extra QSettings keys pinned for the length of the run, beside the provider
        pins — anything the session needs the app to be set to (``report/finalize``
        off, say, so an unattended capture never drives Word).
    prefix : str, optional
        Overrides :data:`PREFIX` in every filename this session writes. The
        scenario suite (``tools/assistant_suite.py``) records dozens of sessions
        into a scratch directory of its own and needs them named for the suite,
        not for the W-1 tutorial; nothing else passes it.
    save_each_turn : bool
        Also write the workbook AFTER EACH TURN, as ``w1_{name}_after_{i}.xlsx``.
        A conversation that edits the model over several turns is only checkable
        turn by turn — the end state cannot show what the second of three edits
        did — so a multi-edit session records one file per turn.

    Returns
    -------
    dict
        ``{name, images, transcript, workbook, turns, usage, seconds, error}`` —
        ``turns`` is one record per turn (prompt, text, usage, seconds, error).
    """
    _env_provider = os.environ.get("XSLOPE_SESSIONS_PROVIDER")
    _env_model = os.environ.get("XSLOPE_SESSIONS_MODEL")
    _env_out = os.environ.get("XSLOPE_SESSIONS_OUT")
    if _env_provider:
        provider = _env_provider
    if _env_model:
        model = _env_model
    if _env_out:
        os.makedirs(_env_out, exist_ok=True)
        out_dir = files_dir = _env_out

    from studio.main_window import MainWindow

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(files_dir, exist_ok=True)
    stem = "%s_%s" % (prefix or PREFIX, name)
    transcript_path = os.path.join(files_dir, "%s_transcript.md" % stem)
    workbook_path = os.path.join(files_dir, "%s_after.xlsx" % stem)

    pins = {"ai/provider": provider, "ai/model/%s" % provider: model,
            "ai/confirm": False,
            # One snippet may not outlast the turn it belongs to. The kernel reads
            # this and stops a snippet that overruns it, which is what makes the
            # per-turn timeout reachable at all: the timeout is a QTimer, and a
            # snippet running on the GUI thread is what stops timers firing.
            "ai/run_timeout": float(timeout_s)}
    pins.update(settings or {})
    result = {"name": name, "images": [], "transcript": transcript_path,
              "workbook": None, "turns": [], "usage": {}, "seconds": 0.0,
              "error": None}
    started = time.time()
    win = None
    totals = {}          # the session's own accumulator, read before the window goes

    print("assistant session %r: %d turn(s), %s/%s%s"
          % (name, len(turns), provider, model, "  [dry run]" if dry_run else ""))

    with _pinned_settings(pins), _usage_probe() as probe:
        try:
            win = MainWindow()
            win.resize(*window_size)
            if model_path:
                win.open_path(model_path)
            else:
                # File -> New: the reader's own starting point for a build-from
                # scratch session. Without it the document is closed, the dock's
                # MODEL SUMMARY has nothing to describe, and the first snippet
                # would be the thing that creates the project.
                win.new_project()
            win.show()
            _settle()
            win.resizeDocks([win.chat_dock], [dock_width], Qt.Horizontal)
            if win.doc.is_open:
                win.canvas.render_inputs(win.doc.slope_data)
            _settle()

            assistant = win.assistant
            chat = win.chat_dock.widget()
            # Belt and braces: the pinned setting is what the config reads, and
            # this is the same value stated through the config's own setter, so a
            # build that caches the flag still runs unattended.
            assistant.config.set_confirm_before_run(False)
            if dry_run:
                _install_stub(assistant)

            _write_header(transcript_path, stem, model_path, provider, model,
                          dry_run)

            for i, turn in enumerate(turns, start=1):
                text, turn_images = _turn_parts(turn, images if i == 1 else None)
                print("  turn %d: %s" % (i, _one_line(text)))
                # The watchdog's own transcript note, written from its thread
                # while this one is stuck, so a killed run still reads as a
                # record of a turn rather than a file that stops mid-session.
                def stopped(message, _i=i, _text=text, _imgs=turn_images):
                    _append_turn(transcript_path, _i,
                                 {"prompt": _text, "images": list(_imgs),
                                  "body": "Error: " + message, "usage": {},
                                  "seconds": timeout_s + WATCHDOG_GRACE_S,
                                  "error": message})
                    _append_footer(transcript_path,
                                   dict(result, error=message,
                                        seconds=time.time() - started))

                with _watchdog(timeout_s + WATCHDOG_GRACE_S,
                               "turn %d did not return within %gs of its %gs "
                               "timeout — the run is stuck on the GUI thread, so "
                               "the harness stopped the process."
                               % (i, WATCHDOG_GRACE_S, timeout_s), stopped):
                    record = _play_turn(win, chat, assistant, text, turn_images,
                                        timeout_s, probe)
                png = _grab_dock(win, os.path.join(out_dir, "%s_%d.png"
                                                   % (stem, i)), max_height)
                record["image"] = png
                if save_each_turn and win.doc.is_open:
                    step = os.path.join(files_dir, "%s_after_%d.xlsx" % (stem, i))
                    with contextlib.redirect_stdout(io.StringIO()):
                        win.doc.save(step)
                    record["workbook"] = step
                    result.setdefault("step_workbooks", []).append(step)
                    print("   -> %s" % os.path.basename(step))
                result["images"].append(png)
                result["turns"].append(record)
                _append_turn(transcript_path, i, record)
                if record["error"]:
                    result["error"] = record["error"]
                    print("  ! %s" % record["error"])
                    break

            totals = _assistant_usage(assistant, "session") or {}
            changed = bool(win.doc.is_open and win.doc.dirty)
            if save_after if save_after is not None else changed:
                with contextlib.redirect_stdout(io.StringIO()):
                    win.doc.save(workbook_path)
                result["workbook"] = workbook_path
                print("   -> %s" % os.path.basename(workbook_path))
        finally:
            if win is not None:
                # The window asks about unsaved changes on close, with a modal
                # box that has nobody to click it. The session's own copy is
                # already written above, so the document is marked clean first
                # and the box never appears.
                with contextlib.suppress(Exception):
                    win.doc._set_dirty(False)
                win.close()
            _settle(4)

    # The Assistant's own totals where it keeps them, the probe's where it does
    # not; the probe's call count is added either way, since the number of
    # completions a turn took is not a token count and only it knows.
    result["usage"] = dict(totals or probe or {})
    if probe is not None:
        result["usage"]["calls"] = probe.get("calls", 0)
    result["seconds"] = time.time() - started
    _append_footer(transcript_path, result)
    print("   session %r: %.1fs, %s" % (name, result["seconds"],
                                        _format_usage(result["usage"])))
    return result


def _turn_parts(turn, default_images):
    """Split a turn item into ``(text, [image paths])``."""
    if isinstance(turn, (tuple, list)):
        text = turn[0]
        extra = turn[1] if len(turn) > 1 else None
        paths = [extra] if isinstance(extra, str) else list(extra or [])
    else:
        text, paths = turn, []
    if not paths and default_images:
        paths = [default_images] if isinstance(default_images, str) \
            else list(default_images)
    return text, paths


def _play_turn(win, chat, assistant, text, image_paths, timeout_s, probe):
    """Type one prompt into the dock, send it, and wait for the turn to finish."""
    before_usage = _numeric_fields(getattr(assistant, "usage", None))
    before_probe = dict(probe or {})
    calls_before = before_probe.get("calls", 0)
    started = time.time()

    for path in image_paths:
        img = QImage(path)
        if img.isNull():
            raise ValueError("attachment is not a readable image: %s" % path)
        chat._add_attachment(img)
    chat.input.setPlainText(text)

    recorder = _Recorder(assistant)
    try:
        chat._send()
        _settle(2)
        timed_out = _wait(recorder, timeout_s)
    finally:
        recorder.disconnect()

    error = recorder.error
    if timed_out and not error:
        error = "timed out after %gs" % timeout_s
        with contextlib.suppress(Exception):
            assistant.cancel()
        _settle(8)

    usage = _assistant_usage(assistant, "turn")
    if usage is None:
        # No accumulator on this build (or a flat one): fall back to what the
        # session's own probe saw pass through litellm between the two reads.
        after_usage = _numeric_fields(getattr(assistant, "usage", None))
        usage = (_usage_delta(before_usage, after_usage)
                 or _usage_delta(before_probe, dict(probe or {})))
    if probe is not None:
        usage = dict(usage, calls=probe.get("calls", 0) - calls_before)
    return {"prompt": text, "images": list(image_paths),
            "body": "\n\n".join(recorder.lines).strip(),
            "usage": usage, "seconds": time.time() - started, "error": error}


def _one_line(text, cap=70):
    line = " ".join((text or "").split())
    return line if len(line) <= cap else line[:cap - 1] + "…"


# --------------------------------------------------------------------------- #
# The transcript file
# --------------------------------------------------------------------------- #
def _write_header(path, stem, model_path, provider, model, dry_run):
    rel = os.path.relpath(model_path, REPO_ROOT) if model_path else "(empty project)"
    lines = [
        "# %s — recorded assistant session" % stem,
        "",
        "Captured by `tools/assistant_sessions.py` against the live Studio "
        "assistant, offscreen.",
        "",
        "- Project: `%s`" % rel,
        "- Provider / model: `%s` / `%s`%s" % (provider, model,
                                               "  (dry run — stub reply)"
                                               if dry_run else ""),
        "- Confirm-before-running: off (unattended capture)",
        "",
        "",
    ]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


def _append_turn(path, i, record):
    """One turn, verbatim, as a plain fenced block.

    Fenced as ``text`` rather than rendered: the point of the file is that the
    tutorial can quote exactly what the model said and ran, and the snippets it
    carries are already Python inside a transcript that also carries prose,
    printed output and a MODEL CHECKS block.
    """
    body = ["You: " + (record["prompt"] or "").strip()]
    for img in record["images"]:
        body.append("[attached: %s]" % os.path.basename(img))
    if record["body"]:
        body.append(record["body"])
    block = "\n\n".join(body)
    fence = "```"
    while fence in block:
        fence += "`"
    parts = [
        "## Turn %d — %s" % (i, _one_line(record["prompt"], 60)),
        "",
        fence + "text",
        block,
        fence,
        "",
        "Tokens: %s · %.1fs" % (_format_usage(record["usage"]),
                                record["seconds"]),
        "",
        "",
    ]
    with open(path, "a", encoding="utf-8") as fh:
        fh.write("\n".join(parts))


def _append_footer(path, result):
    if not os.path.exists(path):
        return
    lines = ["## Session total", "",
             "- Tokens: %s" % _format_usage(result["usage"]),
             "- Wall time: %.1fs" % result["seconds"]]
    if result["workbook"]:
        lines.append("- Workbook written: `%s`"
                     % os.path.relpath(result["workbook"], REPO_ROOT))
    if result["error"]:
        lines.append("- Ended on an error: %s" % result["error"])
    lines.append("")
    with open(path, "a", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


# --------------------------------------------------------------------------- #
# W-1 — Working with the assistant
# --------------------------------------------------------------------------- #
W1_SMOKE_MODEL = os.path.join(REPO_ROOT,
                              "docs/tutorials/files/xslope_reinforced_slope_start.xlsx")


def w1_smoke(dry_run=False):
    """One cheap question against a model that is already built.

    The end-to-end proof for everything above: a real provider call, a real
    ``run_python`` solve inside the live document, a dock grab and a transcript.
    Its prompt asks for a number the reader can check against the tutorial's own
    figures, and it changes nothing, so the session writes no workbook.
    """
    return run_assistant_session("smoke", W1_SMOKE_MODEL,
                                 ["What is the factor of safety of this model "
                                  "with Spencer?"],
                                 dry_run=dry_run)


SESSIONS["w1_smoke"] = w1_smoke


# --- the recorded conversations the tutorial is built from ------------------ #
#: The reinforced slope of LEM-8, finished — the project every session but the
#: first and the seventh opens.
W1_MODEL = os.path.join(REPO_ROOT,
                        "docs/tutorials/files/xslope_reinforced_slope.xlsx")
#: The same slope with three transcription errors written into it (material 2's
#: friction angle 3 instead of 37, the crest surcharge 2400 instead of 240, and
#: the bottom of the model at -100 instead of -10). Built by hand for the
#: diagnosis session; Spencer returns 0.07 on it.
W1_BROKEN = os.path.join(REPO_ROOT, "docs/tutorials/files/w1_diagnose_start.xlsx")
#: The drawing of the problem, the same figure LEM-8 opens on. The build session
#: pastes it into the dock and gives the assistant nothing else.
W1_SKETCH = os.path.join(IMAGES_DIR, "lem08_problem_sketch.png")


def w1_build_from_image(dry_run=False):
    """From an empty project and a drawing, to a solved model."""
    return run_assistant_session(
        "build_from_image", None,
        [("Build this model. Use the dimensions and properties on the drawing. "
          "Unit system: US customary (ft, psf, pcf). Add a starting circle and "
          "run Spencer with a search.", W1_SKETCH)],
        timeout_s=900, max_height=20000, dry_run=dry_run)


def w1_modify(dry_run=False):
    """Three edits to a finished model, in one conversation, each rerun."""
    return run_assistant_session(
        "modify", W1_MODEL,
        ["Change the slope face to 2:1 and rerun the search.",
         "Add a distributed load of 500 psf on the crest from x = 60 to x = 90 "
         "and rerun.",
         "Extend all the reinforcement lines 5 ft to the right and rerun."],
        timeout_s=900, max_height=20000, save_each_turn=True, dry_run=dry_run)


def w1_sweep_builtin(dry_run=False):
    """A sweep the kernel already has a mode for."""
    return run_assistant_session(
        "sweep_builtin", W1_MODEL,
        ["Sweep the geogrid Tmax from 500 to 3000 lb/ft in 6 steps with a search "
         "at each step and plot FS against Tmax."],
        timeout_s=900, max_height=20000, dry_run=dry_run)


def w1_sweep_adhoc(dry_run=False):
    """A study with no mode behind it — a loop the assistant has to write."""
    return run_assistant_session(
        "sweep_adhoc", W1_MODEL,
        ["Run the analysis with 2, 3, 4, 5 and 6 geogrid layers (removing the top "
         "layers first), searching each time, and tabulate FS against the number "
         "of layers."],
        timeout_s=900, max_height=20000, dry_run=dry_run)


def w1_elastic_fem(dry_run=False):
    """Asking for stiffnesses, then running the finite element analysis on them."""
    return run_assistant_session(
        "elastic_fem", W1_MODEL,
        ["Suggest values of Young's modulus and Poisson's ratio for these "
         "materials so I can run a finite element analysis, and explain your "
         "choice.",
         "Enter them, build a quadratic mesh at 2 ft, and run the strength "
         "reduction analysis."],
        timeout_s=1200, max_height=20000, dry_run=dry_run)


def w1_conceptual(dry_run=False):
    """Two questions that change nothing — what the assistant knows, not runs."""
    return run_assistant_session(
        "conceptual", W1_MODEL,
        ["How does a reliability analysis work in XSLOPE?",
         "How do I decide standard deviations for a reliability analysis if I "
         "only have a few tests?"],
        timeout_s=600, max_height=20000, save_after=False, dry_run=dry_run)


def w1_diagnose(dry_run=False):
    """A broken model, and no hint about where the breakage is."""
    return run_assistant_session(
        "diagnose", W1_BROKEN,
        ["This model gives a factor of safety below 1. Can you find what is "
         "wrong?"],
        timeout_s=900, max_height=20000, dry_run=dry_run)


def w1_report(dry_run=False):
    """The analysis report, asked for in a sentence.

    Two turns because the report documents what the session solved: the first
    turn is the search whose result the second turn writes up. ``report/finalize``
    is pinned off for the length of the run — the finish drives Word, and an
    unattended capture must never take over the machine's copy.
    """
    return run_assistant_session(
        "report", W1_MODEL,
        ["Run Spencer with a search.",
         "Generate the analysis report for this model."],
        timeout_s=900, max_height=20000, settings={"report/finalize": False},
        dry_run=dry_run)


for _fn in (w1_build_from_image, w1_modify, w1_sweep_builtin, w1_sweep_adhoc,
            w1_elastic_fem, w1_conceptual, w1_diagnose, w1_report):
    SESSIONS["w1_" + _fn.__name__[3:]] = _fn
del _fn


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    dry_run = "--dry-run" in argv
    argv = [a for a in argv if not a.startswith("-")]
    names = [n for n in SESSIONS if not argv or any(a in n for a in argv)]
    if not names:
        print("no session matching %s; known: %s"
              % (argv, ", ".join(sorted(SESSIONS))))
        return 1
    for name in names:
        SESSIONS[name](dry_run=dry_run)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
