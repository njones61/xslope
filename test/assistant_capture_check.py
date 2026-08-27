"""Standing checks on the recorded-assistant-session harness (tools/assistant_sessions.py).

The harness that produces the AI-assistant tutorial's figures is the one producer
in the pipeline that cannot be re-run casually: each session makes billed provider
calls, so a defect in it is discovered at the cost of the run that finds it. That
is what ``--dry-run`` is for — the whole path with a canned reply and no provider
call — and this file is what proves the dry run really is the whole path.

What is asserted, all offscreen, all against the real ``MainWindow``, the real
``Assistant`` and the real ``ChatDock``:

  A. THE PATH RUNS END TO END — a dry session opens the window on a real project,
     plays its turns through the dock's own send path, grabs the dock to a PNG per
     turn and writes the transcript. The PNG is a real image with the dock's own
     size, not an empty file; the transcript carries the prompt, the reply and a
     fence around them.
  B. NOTHING IS CALLED — ``litellm.completion`` is not entered by a dry run. This
     is the leg that keeps ``--dry-run`` free: it is asserted by counting calls
     through the same wrapper the harness itself uses for token accounting.
  C. THE SETTINGS COME BACK — provider, model and the confirm-before-running gate
     are pinned for the length of the session and restored afterwards, including
     the case of a key that was not set at all (restored by removal, not by an
     empty string). The store cannot be redirected on macOS (see
     ``test/welcome_window_check.py``), so this is asserted against the real one,
     which is precisely the store a leak would damage.
  D. A SWEEP CANNOT SPEND — the sessions registry is disjoint from the screenshot
     producer's ``SHOTS``, and an argument-less ``main`` selects no session. A
     session that leaked into ``SHOTS`` would be run, and billed, by every routine
     figure regeneration.
  E. THE WORKBOOK — a session asked to save writes a loadable .xlsx under the
     session's own deterministic name.
  F. A TURN THAT NEVER ANSWERS — the wait is bounded: a stub that never replies is
     recorded as a timeout error on that turn, the session stops there, and the
     harness returns rather than hanging the producer.
  G. THE NAMES ARE DETERMINISTIC — a second run of the same session overwrites its
     own files and adds none, so a re-record never leaves an orphan behind for the
     docs build to pick up.

Skips cleanly (exit 0) when PySide6 is not installed (engine-only install — no
Studio layer to capture).
"""
import contextlib
import io
import os
import shutil
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

#: The project the dry sessions open — a real tutorial model, so the window goes
#: through the real load/render path rather than an empty document.
MODEL = os.path.join(REPO_ROOT,
                     "docs/tutorials/files/xslope_reinforced_slope_start.xlsx")

PROMPT = "What is the factor of safety of this model with Spencer?"


def _fail(failures, name, message):
    failures.append("%s: %s" % (name, message))
    print("  FAIL  %s — %s" % (name, message))


def _ok(name, detail=""):
    print("  ok    %s%s" % (name, ("  (%s)" % detail) if detail else ""))


@contextlib.contextmanager
def _recent_files_preserved():
    """Keep the user's recent-files list out of this run.

    The harness opens its project through ``MainWindow.open_path``, the real path,
    which records the file in ``recent_files``. That store is the user's own (it
    cannot be redirected on macOS), so the list is stashed and put back — the check
    exercises the real open path without reordering anybody's File menu.
    """
    from PySide6.QtCore import QSettings

    settings = QSettings("XSlope", "XSlope Studio")
    had = settings.contains("recent_files")
    stashed = settings.value("recent_files") if had else None
    try:
        yield
    finally:
        if had:
            settings.setValue("recent_files", stashed)
        else:
            settings.remove("recent_files")
        settings.sync()


@contextlib.contextmanager
def _call_counter():
    """Count entries into ``litellm.completion`` for the length of the block."""
    calls = []
    try:
        import litellm
    except Exception:
        yield calls               # not installed: a dry run cannot call it anyway
        return
    real = litellm.completion

    def counted(*a, **k):
        calls.append(k.get("model"))
        return real(*a, **k)

    litellm.completion = counted
    try:
        yield calls
    finally:
        litellm.completion = real


def _run(sessions, tmp, name, turns, **kw):
    """One dry session into a scratch output directory."""
    kw.setdefault("dry_run", True)
    images = os.path.join(tmp, "images")
    files = os.path.join(tmp, "files")
    with contextlib.redirect_stdout(io.StringIO()):
        return sessions.run_assistant_session(name, MODEL, turns,
                                              out_dir=images, files_dir=files, **kw)


# --------------------------------------------------------------------------- #
# A + B + G — the path, the silence, the names
# --------------------------------------------------------------------------- #
def check_dry_session(sessions, failures):
    from PySide6.QtGui import QImage

    tmp = tempfile.mkdtemp(prefix="xslope_assistant_capture_")
    try:
        with _call_counter() as calls:
            result = _run(sessions, tmp, "probe", [PROMPT, "And with Bishop?"])

        name = "A. dry session runs end to end"
        if result["error"]:
            _fail(failures, name, "session reported an error: %s" % result["error"])
        elif len(result["images"]) != 2 or len(result["turns"]) != 2:
            _fail(failures, name, "expected 2 turns and 2 images, got %d and %d"
                  % (len(result["turns"]), len(result["images"])))
        else:
            bad = [p for p in result["images"] if not os.path.exists(p)]
            if bad:
                _fail(failures, name, "missing PNG: %s" % bad)
            else:
                img = QImage(result["images"][0])
                if img.isNull() or img.width() < 200 or img.height() < 200:
                    _fail(failures, name, "the grab is not a usable image (%s)"
                          % ("null" if img.isNull()
                             else "%dx%d" % (img.width(), img.height())))
                else:
                    _ok(name, "%dx%d dock grab, 2 turns"
                        % (img.width(), img.height()))

        name = "A. the transcript carries both halves of every turn"
        path = result["transcript"]
        text = open(path, encoding="utf-8").read() if os.path.exists(path) else ""
        missing = [w for w in (PROMPT, "And with Bishop?", "## Turn 1", "## Turn 2",
                               "```text", "You:", "Assistant:", "anthropic",
                               "claude-opus-5") if w not in text]
        if not text:
            _fail(failures, name, "no transcript written at %s" % path)
        elif missing:
            _fail(failures, name, "transcript is missing %s" % missing)
        elif text.count("```") % 2:
            _fail(failures, name, "unbalanced code fences in the transcript")
        else:
            _ok(name, "%d chars, %d turn heading(s)"
                % (len(text), text.count("\n## Turn ")))

        name = "B. a dry run contacts no provider"
        if calls:
            _fail(failures, name, "litellm.completion was called %d time(s): %s"
                  % (len(calls), calls))
        else:
            _ok(name)

        name = "G. a re-record overwrites its own files and adds none"
        before = sorted(os.listdir(os.path.join(tmp, "images")))
        with _call_counter():
            again = _run(sessions, tmp, "probe", [PROMPT, "And with Bishop?"])
        after = sorted(os.listdir(os.path.join(tmp, "images")))
        if before != after:
            _fail(failures, name, "the second run changed the file set: %s -> %s"
                  % (before, after))
        elif again["images"] != result["images"]:
            _fail(failures, name, "image paths are not deterministic")
        elif again["transcript"] != result["transcript"]:
            _fail(failures, name, "transcript path is not deterministic")
        else:
            _ok(name, ", ".join(before))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- #
# C — the pinned settings come back
# --------------------------------------------------------------------------- #
def check_settings_restored(sessions, failures):
    from PySide6.QtCore import QSettings

    settings = QSettings("XSlope", "XSlope Studio")
    keys = ["ai/provider", "ai/model/anthropic", "ai/confirm",
            "ai/model/openai"]                     # the last one is the unset case
    settings.remove("ai/model/openai")
    settings.sync()
    before = {k: (settings.value(k) if settings.contains(k) else None) for k in keys}

    tmp = tempfile.mkdtemp(prefix="xslope_assistant_capture_")
    try:
        _run(sessions, tmp, "settings", [PROMPT], provider="anthropic",
             model="claude-haiku-4-5")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    fresh = QSettings("XSlope", "XSlope Studio")
    after = {k: (fresh.value(k) if fresh.contains(k) else None) for k in keys}
    name = "C. the pinned settings are restored"
    drift = {k: (before[k], after[k]) for k in keys if before[k] != after[k]}
    if drift:
        _fail(failures, name, "left behind %s" % drift)
    elif after["ai/model/openai"] is not None:
        _fail(failures, name, "an unset key came back set")
    else:
        _ok(name, "provider/model/confirm unchanged, unset key still unset")


# --------------------------------------------------------------------------- #
# D — an unnamed sweep runs no session
# --------------------------------------------------------------------------- #
def check_sweep_cannot_spend(sessions, failures):
    import tools.capture_tutorial_screenshots as capture

    name = "D. sessions are not screenshots"
    overlap = sorted(set(sessions.SESSIONS) & set(capture.SHOTS))
    if overlap:
        _fail(failures, name, "a session is registered as a shot: %s" % overlap)
    else:
        _ok(name, "%d session(s), %d shot(s), disjoint"
            % (len(sessions.SESSIONS), len(capture.SHOTS)))

    name = "D. an argument-less capture run selects no session"
    ran = []
    stashed = dict(sessions.SESSIONS)
    original = {k: capture.SESSIONS.get(k) for k in list(capture.SESSIONS)}
    try:
        for key in list(sessions.SESSIONS):
            sessions.SESSIONS[key] = lambda dry_run=False, _k=key: ran.append(_k)
        for key in list(capture.SESSIONS):
            capture.SESSIONS[key] = sessions.SESSIONS[key]
        # The selection is what is under test, not the shots — which are what the
        # sweep would otherwise spend minutes rendering — so the shot table is
        # emptied for the call and put back.
        shots = dict(capture.SHOTS)
        capture.SHOTS.clear()
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                capture.main([])
        finally:
            capture.SHOTS.update(shots)
    finally:
        sessions.SESSIONS.update(stashed)
        capture.SESSIONS.update({k: v for k, v in original.items() if v})
    if ran:
        _fail(failures, name, "a sweep ran %s" % ran)
    else:
        _ok(name)

    name = "D. naming a session runs it"
    ran = []
    stashed = dict(sessions.SESSIONS)
    original = dict(capture.SESSIONS)
    try:
        for key in list(capture.SESSIONS):
            capture.SESSIONS[key] = lambda dry_run=False, _k=key: ran.append(
                (_k, dry_run))
        shots = dict(capture.SHOTS)
        capture.SHOTS.clear()
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                capture.main(["w1", "--dry-run"])
        finally:
            capture.SHOTS.update(shots)
    finally:
        sessions.SESSIONS.update(stashed)
        capture.SESSIONS.update(original)
    if not ran:
        _fail(failures, name, "'w1' selected nothing")
    elif not all(dry for _, dry in ran):
        _fail(failures, name, "--dry-run did not reach the session: %s" % ran)
    else:
        _ok(name, ", ".join(k for k, _ in ran))


# --------------------------------------------------------------------------- #
# E — the workbook a session leaves behind
# --------------------------------------------------------------------------- #
def check_workbook_written(sessions, failures):
    from xslope.fileio import load_slope_data

    tmp = tempfile.mkdtemp(prefix="xslope_assistant_capture_")
    name = "E. a session that saves writes a loadable workbook"
    try:
        result = _run(sessions, tmp, "saved", [PROMPT], save_after=True)
        path = result["workbook"]
        if not path or not os.path.exists(path):
            _fail(failures, name, "no workbook written")
            return
        if os.path.basename(path) != "w1_saved_after.xlsx":
            _fail(failures, name, "unexpected name %s" % os.path.basename(path))
            return
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                sd = load_slope_data(path)
        except Exception as exc:
            _fail(failures, name, "the workbook does not load back: %s" % exc)
            return
        if not sd.get("materials"):
            _fail(failures, name, "the workbook came back empty")
        else:
            _ok(name, "%s, %d material(s)"
                % (os.path.basename(path), len(sd["materials"])))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def check_no_workbook_when_unchanged(sessions, failures):
    tmp = tempfile.mkdtemp(prefix="xslope_assistant_capture_")
    name = "E. a session that changed nothing writes no workbook"
    try:
        result = _run(sessions, tmp, "readonly", [PROMPT])
        if result["workbook"]:
            _fail(failures, name, "wrote %s from a read-only session"
                  % result["workbook"])
        else:
            _ok(name)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- #
# F — the wait is bounded
# --------------------------------------------------------------------------- #
def check_timeout(sessions, failures):
    """A turn whose reply never arrives must end the session, not hang it."""
    tmp = tempfile.mkdtemp(prefix="xslope_assistant_capture_")
    name = "F. a turn that never answers times out"
    real = sessions._install_stub
    try:
        sessions._install_stub = lambda assistant, reply="…": setattr(
            assistant, "send", lambda user_text, images=None: None)
        result = _run(sessions, tmp, "stuck", [PROMPT, "second turn"],
                      timeout_s=1.0)
        if not result["error"] or "timed out" not in result["error"]:
            _fail(failures, name, "no timeout recorded (error=%r)" % result["error"])
        elif len(result["turns"]) != 1:
            _fail(failures, name, "the session continued past the stuck turn (%d turns)"
                  % len(result["turns"]))
        elif not os.path.exists(result["transcript"]):
            _fail(failures, name, "the stuck turn was not written to the transcript")
        else:
            text = open(result["transcript"], encoding="utf-8").read()
            if "timed out" not in text:
                _fail(failures, name, "the transcript does not record the timeout")
            else:
                _ok(name, result["error"])
    finally:
        sessions._install_stub = real
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- #
def run():
    failures = []
    import tools.assistant_sessions as sessions

    if not os.path.exists(MODEL):
        return ["the check's project is missing: %s" % MODEL]

    with _recent_files_preserved():
        check_dry_session(sessions, failures)
        check_settings_restored(sessions, failures)
        check_sweep_cannot_spend(sessions, failures)
        check_workbook_written(sessions, failures)
        check_no_workbook_when_unchanged(sessions, failures)
        check_timeout(sessions, failures)
    return failures


def main():
    print("assistant capture harness (tools/assistant_sessions.py):")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print("  - %s" % f)
        raise SystemExit(1)
    print("\nAll assistant capture checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication      # noqa: F401
    except Exception:
        print("assistant_capture_check: PySide6 not installed — skipped.")
        raise SystemExit(0)
    main()
