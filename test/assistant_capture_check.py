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
  C. THE MACHINE'S SETTINGS ARE NEVER TOUCHED — provider, model and the
     confirm-before-running gate are stated for the length of the session in a
     store of the session's OWN, and the user's is only read. It used to be
     written and restored, which cost a measured corpus run: two sessions running
     at once pinned the same machine-wide store, the second restored the first's
     ``ai/confirm``, and the first stopped dead on the "Run code?" modal with
     nobody to click it. So what is asserted is the whole store, key by key,
     unchanged across a session — and the pins present in the session's own file.
  H. TWO SESSIONS AT ONCE DO NOT COLLIDE — two dry sessions in separate processes,
     held inside their pinned region at the same time by a file handshake, each
     read back the provider and model THEY asked for, out of two different stores,
     and the user's store is unchanged by both. Under the restore-the-machine-store
     design this is the leg that fails.
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


def _machine_settings():
    """Every key and value in Studio's own machine-wide store, as a dict.

    The comparison this check is built on. A session must leave it identical —
    including ``recent_files``, which the real open path would otherwise reorder,
    and which is now written to the session's own store like everything else.
    """
    from PySide6.QtCore import QSettings

    settings = QSettings("XSlope", "XSlope Studio")
    settings.sync()
    return {key: settings.value(key) for key in settings.allKeys()}


def _machine_settings_unchanged(before, failures, name, what):
    after = _machine_settings()
    drift = {k: (before.get(k), after.get(k))
             for k in set(before) | set(after) if before.get(k) != after.get(k)}
    if drift:
        _fail(failures, name, "%s changed the user's own settings: %s"
              % (what, drift))
        return False
    return True


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
# C — the session states its settings in a store of its own
# --------------------------------------------------------------------------- #
def check_settings_are_the_sessions_own(sessions, failures):
    from PySide6.QtCore import QSettings

    before = _machine_settings()
    tmp = tempfile.mkdtemp(prefix="xslope_assistant_capture_")
    seen = {}
    real_stub = sessions._install_stub

    def stub(assistant, reply="…"):
        # Read what the running window actually resolved, from inside it.
        seen.update(provider=assistant.config.provider(),
                    model=assistant.config.model(),
                    confirm=assistant.config.confirm_before_run(),
                    store=assistant.config._s.fileName())
        real_stub(assistant, reply)

    sessions._install_stub = stub
    try:
        _run(sessions, tmp, "settings", [PROMPT], provider="anthropic",
             model="claude-haiku-4-5")
    finally:
        sessions._install_stub = real_stub
        shutil.rmtree(tmp, ignore_errors=True)

    name = "C. the user's own settings are untouched"
    if _machine_settings_unchanged(before, failures, name, "a dry session"):
        _ok(name, "%d key(s) identical" % len(before))

    name = "C. the session reads its own store"
    machine = QSettings("XSlope", "XSlope Studio").fileName()
    if (seen.get("model") != "claude-haiku-4-5"
            or seen.get("provider") != "anthropic"):
        _fail(failures, name, "the window resolved %s/%s, not the session's pins"
              % (seen.get("provider"), seen.get("model")))
    elif seen.get("confirm") is not False:
        _fail(failures, name, "confirm-before-running is %r, so an unattended "
              "turn would stop on the modal" % seen.get("confirm"))
    elif (not seen.get("store")
            or os.path.abspath(seen["store"]) == os.path.abspath(machine)):
        _fail(failures, name, "the window is reading the machine store (%s)"
              % seen.get("store"))
    elif str(os.getpid()) not in os.path.basename(seen["store"]):
        _fail(failures, name, "the session's store is not per-process: %s"
              % seen["store"])
    else:
        _ok(name, os.path.basename(seen["store"]))


# --------------------------------------------------------------------------- #
# H — two sessions at once
# --------------------------------------------------------------------------- #
#: Driver for one session in a process of its own. It holds the session inside its
#: pinned region until the other process is in there too (the handshake through
#: ``sync``), and only then reads back what the running window resolved — so a
#: store shared between the two is read AFTER the other one has written it, which
#: is the collision this leg exists to catch.
_PARALLEL_DRIVER = '''
import json, os, sys, time

sys.path.insert(0, %(repo)r)
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import matplotlib
matplotlib.use("Agg")
from PySide6.QtWidgets import QApplication
_app = QApplication.instance() or QApplication([])

import tools.assistant_sessions as sessions

model, tag, other, sync, out_dir = sys.argv[1:6]
seen = {}
real = sessions._install_stub


def stub(assistant, reply="."):
    open(os.path.join(sync, tag), "w").close()
    deadline = time.time() + 120
    while not os.path.exists(os.path.join(sync, other)) and time.time() < deadline:
        time.sleep(0.05)
    seen.update(provider=assistant.config.provider(),
                model=assistant.config.model(),
                confirm=assistant.config.confirm_before_run(),
                store=assistant.config._s.fileName(),
                both_in=os.path.exists(os.path.join(sync, other)))
    real(assistant, reply)


sessions._install_stub = stub
result = sessions.run_assistant_session(
    tag, %(model_path)r, ["ping"], dry_run=True, model=model,
    out_dir=os.path.join(out_dir, "images"),
    files_dir=os.path.join(out_dir, "files"))
seen["error"] = result["error"]
sys.stderr.write("RESULT " + json.dumps(seen) + "\\n")
'''


def check_parallel_sessions_do_not_collide(sessions, failures):
    import json
    import subprocess
    import threading

    name = "H. two sessions at once do not collide"
    before = _machine_settings()
    tmp = tempfile.mkdtemp(prefix="xslope_assistant_parallel_")
    sync = os.path.join(tmp, "sync")
    os.makedirs(sync)
    driver = os.path.join(tmp, "driver.py")
    with open(driver, "w", encoding="utf-8") as fh:
        fh.write(_PARALLEL_DRIVER % {"repo": REPO_ROOT, "model_path": MODEL})
    runs = {}

    def go(tag, model, other):
        proc = subprocess.run(
            [sys.executable, driver, model, tag, other, sync,
             os.path.join(tmp, tag)],
            capture_output=True, text=True, timeout=900)
        line = [ln for ln in proc.stderr.splitlines() if ln.startswith("RESULT ")]
        runs[tag] = (json.loads(line[-1][len("RESULT "):]) if line
                     else {"error": "no result line; stderr tail: %s"
                                    % proc.stderr[-400:]})

    threads = [threading.Thread(target=go, args=a) for a in
               (("alpha", "claude-haiku-4-5", "beta"),
                ("beta", "claude-sonnet-5", "alpha"))]
    try:
        for t in threads:
            t.start()
        for t in threads:
            t.join()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    wanted = {"alpha": "claude-haiku-4-5", "beta": "claude-sonnet-5"}
    problems = []
    for tag, want in wanted.items():
        got = runs.get(tag) or {}
        if got.get("error"):
            problems.append("%s failed: %s" % (tag, got["error"]))
            continue
        if not got.get("both_in"):
            problems.append("%s never overlapped the other session, so nothing "
                            "about concurrency was measured" % tag)
        if got.get("model") != want:
            problems.append("%s ran as %r, not %r — the two sessions share a store"
                            % (tag, got.get("model"), want))
        if got.get("confirm") is not False:
            problems.append("%s had confirm-before-running back on, so it would "
                            "hang on the modal" % tag)
    stores = {tag: (runs.get(tag) or {}).get("store") for tag in wanted}
    if len(set(stores.values())) != len(stores):
        problems.append("both sessions used the same store: %s" % stores)
    if problems:
        for problem in problems:
            _fail(failures, name, problem)
    else:
        _ok(name, ", ".join("%s=%s" % (t, os.path.basename(s or ""))
                            for t, s in stores.items()))
    _machine_settings_unchanged(before, failures, name, "two parallel sessions")


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

    check_dry_session(sessions, failures)
    check_settings_are_the_sessions_own(sessions, failures)
    check_parallel_sessions_do_not_collide(sessions, failures)
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
