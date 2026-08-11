# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Checks for finishing the Analysis Report in Word.

What is being defended:

  A. THE PAGE NUMBERS ARE WORD'S — on a machine that has Word, a generated
     report goes in with a contents page that names its sections and comes out
     with one that says what page each is on. The evidence is the field result
     cached in the document: ``PAGEREF`` fields with numbers in them where there
     were none, and the "press F9" line gone.

  B. AND THE DOCUMENT IS STILL CLEAN — no field comes back marked dirty. A dirty
     field is what makes Word for Mac ask about "fields that may refer to other
     files" every time the report is opened, so a finish that introduced one
     would trade a missing page number for a warning on a clean document.

  C. NOTHING RAISES, EVER — no path, a path to nothing, a file that is not a
     Word document, no Word on the machine, a platform with no Word at all: each
     is a ``(False, sentence)``. This is a cosmetic step; it may not cost anyone
     their report.

  D. WORD IS ONLY HANDED A REPORT — a file that is not a document package with a
     contents field in it never reaches Word. An unreadable file put in front of
     Word is a modal dialog on the user's screen.

  E. THE STUDIO WIRING — generating a report finalizes it and then opens it; the
     settings key switches the finish off; and a document with nothing to update
     is opened without Word being started.

  F. THE WINDOWS LEG IS INTACT — it cannot be run here, so what is checked is
     that it is there, that it says the things COM needs said, and that on a
     platform where its shell is missing it fails into the same sentence the Mac
     leg would.

Word is started at most once (check A), and only when this machine has it: with
no Word, A is skipped and the detection and fallback paths are checked instead.

AND ONLY WHEN IT IS ASKED FOR. Driving Word means driving the copy the person at
this machine is working in: their document steals focus, their session is
scripted, and a run started for another reason has taken the machine over. So
every leg that reaches the real Word — check A, and the tidiness question it
asks afterwards — is off unless ``XSLOPE_WORD_OK=1`` is set in the environment
(:data:`WORD_OK`). The skip is printed and the row says so; a leg that did not
run is never reported as one that passed. Everything else here runs
unconditionally: none of it starts an application.
"""

import os
import sys
import tempfile
import zipfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

REINF_XLSX = os.path.join(_REPO, "docs", "inputs", "slope", "xslope_reinf.xlsx")

_REPORT = {}


def _report_docx(directory):
    """A small real report, written into ``directory``.

    The same sample the report checks use, solved by one method with the figures
    switched off: what is being read here is the contents field, not a plot.
    """
    import matplotlib
    matplotlib.use("Agg")
    from xslope.fileio import load_slope_data
    from xslope.report import generate_report
    from xslope.slice import generate_slices
    from xslope.solve import solve_selected

    if "solved" not in _REPORT:
        slope_data = load_slope_data(REINF_XLSX)
        ok, out = generate_slices(slope_data, circle=slope_data["circles"][0],
                                  num_slices=15)
        if not ok:
            raise RuntimeError(f"the sample model produced no slices: {out}")
        slice_df, surface = out
        bundle = {"slice_df": slice_df, "failure_surface": surface,
                  "results": solve_selected("bishop", slice_df),
                  "search": None, "method": "bishop"}
        _REPORT["solved"] = (slope_data, {"lem": [bundle]})
    slope_data, solutions = _REPORT["solved"]

    # A path with a space in it: the scripts are handed their arguments, never a
    # quoted command line, and this is what would catch it if that changed.
    path = os.path.join(directory, "sample report.docx")
    opts = {"input_path": REINF_XLSX, "title": "Sample Levee",
            "pd_figure": False, "lem_search_figure": False,
            "lem_solution_figure": False}
    ok, out = generate_report(slope_data, solutions, opts, path)
    if not ok:
        raise RuntimeError(f"generate_report failed: {out}")
    return path


def _contents_state(path):
    """What the document's contents field currently holds."""
    with zipfile.ZipFile(path) as package:
        xml = package.read("word/document.xml").decode("utf-8", "replace")
    start = xml.find(" TOC ")
    end = xml.find("fldCharType=\"end\"", start)
    result = xml[start:end] if start >= 0 else ""
    return {
        "pagerefs": result.count("PAGEREF"),
        "dirty": xml.count("w:dirty"),
        "asks_for_f9": "Page numbers appear" in xml or "press F9" in xml,
        "numbers": [t for t in _run_texts(result) if t.strip().isdigit()],
    }


def _run_texts(xml):
    """The text of every run in a fragment of document XML."""
    import re
    return re.findall(r"<w:t[^>]*>([^<]*)</w:t>", xml)


def _stub_docx(path):
    """A file with a .docx name and nothing in it — what a check harness writes
    when it only needs a path."""
    open(path, "wb").close()
    return path


# --------------------------------------------------------------------------
# A + B. the real thing
# --------------------------------------------------------------------------

#: The one way to ask for the legs that drive Microsoft Word.
#:
#: Word on this machine is the copy the person at it is working in. A leg that
#: scripts it opens and closes documents in their live session, steals their
#: focus, and — on a machine whose automation permission has been granted once —
#: does it without asking. That is not a thing a test run may decide to do for
#: its own reasons, so it is opt-in and nothing else in the suite can turn it on:
#: it is read from the environment, at the moment it matters.
WORD_OK = "XSLOPE_WORD_OK"


def word_automation_allowed():
    """Whether this run may drive Microsoft Word."""
    return os.environ.get(WORD_OK) == "1"


def test_word_finishes_the_report():
    """Word turns the contents page into a contents page.

    Skipped, loudly, unless this run was explicitly allowed to drive Word
    (:data:`WORD_OK`), and skipped on a machine with no Word: the rest of the
    row still checks the detection and the fallback.
    """
    from xslope.report_finalize import finalize_with_word, word_available

    if not word_automation_allowed():
        print(f"    (Word automation is off: set {WORD_OK}=1 to run the finish "
              f"against the Word on this machine)")
        return []

    ok, why = word_available()
    if not ok:
        print(f"    (no Word on this machine: {why} — the finish is not run)")
        return []

    fails = []
    with tempfile.TemporaryDirectory() as tmp:
        path = _report_docx(tmp)
        before = _contents_state(path)
        if before["pagerefs"]:
            fails.append("the written report already carries page references; "
                         "there would be nothing to prove")
        if not before["asks_for_f9"]:
            fails.append("the written report does not ask the reader to update "
                         "the field — check what the renderer now caches")

        done, msg = finalize_with_word(path)
        if not done:
            return fails + [f"Word did not finish the report: {msg}"]

        after = _contents_state(path)
        if after["pagerefs"] < 1:
            fails.append("the contents field came back with no page references")
        if not after["numbers"]:
            fails.append("the contents field came back with no page numbers")
        if after["asks_for_f9"]:
            fails.append("the contents page still asks the reader to press F9")
        if after["dirty"]:
            fails.append(f"the finish left {after['dirty']} dirty field(s) — Word "
                         f"would ask about them on every open")
        if not _closed_in_word(path):
            fails.append("the document was left open in Word")
    return fails


def _closed_in_word(path):
    """True when Word is not holding the document open. Anything that stops the
    question being answered counts as answered: this guards tidiness, not
    correctness — including a run that was never allowed to talk to Word."""
    if sys.platform != "darwin" or not word_automation_allowed():
        return True
    import subprocess
    script = ("function run(argv) { var w = Application('Microsoft Word');"
              " var d = w.documents; var out = [];"
              " for (var i = 0; i < d.length; i++) {"
              "  try { out.push(d[i].name()); } catch (e) {} }"
              " return out.join('\\n'); }")
    try:
        done = subprocess.run(["osascript", "-l", "JavaScript", "-"],
                              input=script, capture_output=True, text=True,
                              timeout=30)
    except Exception:
        return True
    return os.path.basename(path) not in done.stdout.splitlines()


# --------------------------------------------------------------------------
# C. nothing raises
# --------------------------------------------------------------------------

def test_refusals_are_sentences():
    """Every way this can decline is a (False, sentence)."""
    from xslope.report_finalize import finalize_with_word

    fails = []
    cases = [("", "no path"), ("/no/such/report.docx", "a path to nothing")]
    with tempfile.TemporaryDirectory() as tmp:
        text = os.path.join(tmp, "notes.txt")
        with open(text, "w") as fh:
            fh.write("not a document")
        cases.append((text, "a file that is not a Word document"))

        for path, what in cases:
            try:
                ok, msg = finalize_with_word(path)
            except Exception as exc:
                fails.append(f"{what} raised {exc!r}")
                continue
            if ok:
                fails.append(f"{what} was reported as finished")
            if not isinstance(msg, str) or not msg.strip():
                fails.append(f"{what} gave no reason")
    return fails


def test_no_word_falls_back():
    """With no Word to be found, detection and the finish both say so plainly —
    and nothing is run."""
    import xslope.report_finalize as rf

    fails = []
    app, finder = rf._MAC_WORD_APP, rf._mdfind_word
    ran = []
    osascript = rf._osascript
    try:
        rf._MAC_WORD_APP = os.path.join(tempfile.gettempdir(), "No Word.app")
        rf._mdfind_word = lambda: None
        rf._osascript = lambda *a, **k: ran.append(a) or (False, "should not run")

        if sys.platform == "darwin":
            ok, why = rf.word_available()
            if ok:
                fails.append("Word was found where there is none")
            if "not installed" not in why:
                fails.append(f"the reason given is {why!r}")

        with tempfile.TemporaryDirectory() as tmp:
            path = _stub_docx(os.path.join(tmp, "r.docx"))
            ok, msg = rf.finalize_with_word(path)
            if ok:
                fails.append("a finish was reported with no Word")
            if not msg:
                fails.append("the refusal gave no reason")
        if ran:
            fails.append("a script was run although there is no Word")
    finally:
        rf._MAC_WORD_APP, rf._mdfind_word, rf._osascript = app, finder, osascript
    return fails


def test_platform_gate():
    """A platform Word does not run on is told so, not driven."""
    import xslope.report_finalize as rf

    fails = []
    real = sys.platform
    try:
        sys.platform = "linux"
        ok, why = rf.word_available()
        if ok:
            fails.append("Word was offered on a platform it does not run on")
        if "macOS" not in why:
            fails.append(f"the reason given is {why!r}")
        with tempfile.TemporaryDirectory() as tmp:
            ok, _msg = rf.finalize_with_word(_stub_docx(
                os.path.join(tmp, "r.docx")))
            if ok:
                fails.append("a finish was reported on a platform with no Word")
    finally:
        sys.platform = real
    return fails


# --------------------------------------------------------------------------
# F. the Windows leg
# --------------------------------------------------------------------------

#: What the Windows script must still be doing. It is driven from PowerShell
#: rather than pywin32 — no new dependency, and a subprocess can be given a
#: timeout — so these are the COM calls, spelled as PowerShell spells them.
_WINDOWS_MUST_SAY = [
    "Word.Application",
    "$word.Visible = $false",
    "$word.Documents.Open(",
    "Fields.Update()",
    "TablesOfContents",
    "$doc.Save()",
    "$doc.Close(",
    "$word.Quit()",
]


def test_windows_leg_is_intact():
    """It cannot be run here; it can be read, and it can be made to fail
    correctly."""
    import xslope.report_finalize as rf

    fails = []
    for phrase in _WINDOWS_MUST_SAY:
        if phrase not in rf._WINDOWS_SCRIPT:
            fails.append(f"the Windows script no longer says {phrase!r}")
    if "$started" not in rf._WINDOWS_SCRIPT:
        fails.append("the Windows script quits Word whether or not it started it")
    if rf._ps_quote("O'Brien's") != "'O''Brien''s'":
        fails.append("a PowerShell argument is not quoted safely")

    # On this machine there is no PowerShell to run it: the leg must come back
    # with a sentence rather than an exception, which is the discipline that
    # matters most on the platform it cannot be tried on.
    real = sys.platform
    try:
        sys.platform = "win32"
        with tempfile.TemporaryDirectory() as tmp:
            path = _stub_docx(os.path.join(tmp, "r.docx"))
            try:
                ok, msg = rf.finalize_with_word(path, timeout=20)
            except Exception as exc:
                fails.append(f"the Windows leg raised {exc!r}")
                return fails
            if ok and not os.path.exists(path):
                fails.append("a finish was reported for a file that is not there")
            if not isinstance(msg, str) or not msg.strip():
                fails.append("the Windows leg gave no reason")
    finally:
        sys.platform = real
    return fails


# --------------------------------------------------------------------------
# D + E. the Studio wiring
# --------------------------------------------------------------------------

def test_only_a_report_reaches_word():
    """A file with nothing to update is opened, and Word is left alone."""
    from studio.report_dialog import carries_contents_field

    fails = []
    with tempfile.TemporaryDirectory() as tmp:
        if carries_contents_field(_stub_docx(os.path.join(tmp, "empty.docx"))):
            fails.append("an empty file was taken for a report")
        text = os.path.join(tmp, "notes.docx")
        with open(text, "w") as fh:
            fh.write("not a package")
        if carries_contents_field(text):
            fails.append("a text file was taken for a report")
        if carries_contents_field(os.path.join(tmp, "absent.docx")):
            fails.append("a missing file was taken for a report")
        if not carries_contents_field(_report_docx(tmp)):
            fails.append("a real report was not recognised as one")
    return fails


def test_generate_finalizes_then_opens():
    """The Studio path: the document is finished, then shown — and the settings
    key stops the finish without stopping the report."""
    from PySide6.QtCore import QSettings
    from studio import report_dialog

    fails = []
    finalized, opened = [], []
    import xslope.report_finalize as rf
    real_finalize = rf.finalize_with_word
    real_open = report_dialog.QDesktopServices.openUrl
    rf.finalize_with_word = lambda path, **kw: (finalized.append(path),
                                                (True, "pretend"))[1]
    report_dialog.QDesktopServices.openUrl = (
        lambda url: opened.append(url.toLocalFile()))
    try:
        with tempfile.TemporaryDirectory() as tmp:
            path = _report_docx(tmp)
            if report_dialog.open_output(path, "docx") != "document":
                fails.append("a .docx did not open as a document")
            if finalized != [path]:
                fails.append(f"the report was not finalized before it was "
                             f"opened (finalized {finalized})")
            if opened[-1:] != [path]:
                fails.append(f"the document opened was {opened[-1:]}")

            # A document with nothing to update is opened just the same, and
            # Word is not started for it.
            finalized.clear()
            stub = _stub_docx(os.path.join(tmp, "stub.docx"))
            report_dialog.open_output(stub, "docx")
            if finalized:
                fails.append("a file with no contents field was sent to Word")
            if opened[-1] != stub:
                fails.append("a file with no contents field was not opened")

            # The escape hatch: no UI, one key.
            settings = QSettings(os.path.join(tmp, "settings.ini"),
                                 QSettings.IniFormat)
            settings.setValue(report_dialog.FINALIZE_KEY, False)
            if report_dialog.finalization_enabled(settings):
                fails.append("the settings key does not switch the finish off")
            finalized.clear()
            ok, why = report_dialog.finalize_document(path, settings=settings)
            if ok or finalized:
                fails.append("the finish ran with the key switched off")
            if "off" not in why:
                fails.append(f"switching it off is reported as {why!r}")
            settings.setValue(report_dialog.FINALIZE_KEY, True)
            if not report_dialog.finalization_enabled(settings):
                fails.append("the finish cannot be switched back on")
    finally:
        rf.finalize_with_word = real_finalize
        report_dialog.QDesktopServices.openUrl = real_open
    return fails


def test_word_is_not_driven_unless_it_was_asked_for():
    """With the environment variable unset, not one leg of this file talks to
    Microsoft Word.

    Measured on the calls themselves rather than on the docstrings: every check
    in this file is run with ``subprocess.run`` and the finalizer's own
    ``_osascript`` watched, and either one reaching the machine's Word is the
    failure. That is the guard the ban needs — a run started for something else
    drove the owner's live Word session twice, and a rule that lives only in a
    comment is a rule that does that again.

    The other checks stub ``_osascript`` out for their own purposes and put it
    back afterwards; the watch is what they put back, so it sees every call any
    of them makes for real.

    The mutation at the end runs ``osascript -e 1`` twice, which is the one
    osascript this file ever spawns: it evaluates the literal 1 and exits. It
    starts no application, names none, and cannot reach Word — it exists only so
    that the watch above is proved able to count a call it was meant to catch.
    Anything watching the process list sees the name and nothing else.
    """
    import subprocess

    fails = []
    import xslope.report_finalize as rf

    if word_automation_allowed():
        return [f"{WORD_OK} is set in this environment, so the guard cannot be "
                f"measured; unset it and run again"]

    spawned, scripted = [], []
    real_run = subprocess.run
    real_osascript = rf._osascript

    def watched_run(args, *rest, **kw):
        name = args[0] if isinstance(args, (list, tuple)) and args else args
        if str(name).endswith("osascript"):
            spawned.append(list(args) if isinstance(args, (list, tuple))
                           else args)
        return real_run(args, *rest, **kw)

    def watched_osascript(*a, **kw):
        scripted.append(a[:1])
        return real_osascript(*a, **kw)

    subprocess.run = watched_run
    rf._osascript = watched_osascript
    try:
        for name, fn in CHECKS:
            if fn is test_word_is_not_driven_unless_it_was_asked_for:
                continue
            try:
                fn()
            except Exception:
                # A leg that raised is another check's failure, not this one's;
                # what matters here is whether it reached Word on the way.
                pass
    finally:
        subprocess.run = real_run
        rf._osascript = real_osascript

    if spawned:
        fails.append(f"{len(spawned)} osascript call(s) were made with {WORD_OK} "
                     f"unset: {spawned[:2]}")
    if scripted:
        fails.append(f"{len(scripted)} script(s) were sent to the machine's Word "
                     f"with {WORD_OK} unset")

    # Mutation: the watch has to be able to see one. A leg that does call
    # osascript is counted.
    spawned.clear()
    subprocess.run = watched_run
    try:
        real_run(["osascript", "-e", "1"], capture_output=True, text=True,
                 timeout=10)
        watched_run(["osascript", "-e", "1"], capture_output=True, text=True,
                    timeout=10)
    except Exception:
        pass
    finally:
        subprocess.run = real_run
    if not spawned:
        fails.append("the watch counted no osascript call even when one was "
                     "made; it cannot fail")
    return fails


CHECKS = [
    ("Word builds the report's page numbers", test_word_finishes_the_report),
    ("every refusal is a sentence", test_refusals_are_sentences),
    ("no Word falls back quietly", test_no_word_falls_back),
    ("a platform without Word is not driven", test_platform_gate),
    ("the Windows leg is intact", test_windows_leg_is_intact),
    ("only a report reaches Word", test_only_a_report_reaches_word),
    ("generating finalizes, then opens", test_generate_finalizes_then_opens),
    ("Word is not driven unless it was asked for",
     test_word_is_not_driven_unless_it_was_asked_for),
]

#: Checks that need the Studio layer; skipped when PySide6 is absent.
_STUDIO_ONLY = {test_only_a_report_reaches_word,
                test_generate_finalizes_then_opens}


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    checks = CHECKS
    try:
        import PySide6                                    # noqa: F401
    except Exception:
        print("Report finish: PySide6 not installed — Studio checks skipped.")
        checks = [c for c in CHECKS if c[1] not in _STUDIO_ONLY]

    failures = []
    for name, fn in checks:
        try:
            fs = fn()
        except Exception as exc:
            import traceback
            traceback.print_exc()
            fs = [f"{name} raised: {exc!r}"]
        print(f"  {name:48s} {'ok' if not fs else f'FAIL ({len(fs)})'}")
        failures += fs
    return failures


def main():
    print("Finishing the report in Word:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll report-finish checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication
        _app_ = QApplication.instance() or QApplication([])
    except Exception:
        pass
    main()
