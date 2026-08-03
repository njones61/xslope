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

"""Finish a written report through Word itself, so its page numbers are real.

The Analysis Report carries a real table-of-contents field. The field is
deliberately not marked dirty — a dirty field is what makes Word for Mac ask
"this document contains fields that may refer to other files" every time the
report is opened — and the consequence is that its page numbers only appear
once someone presses F9. Page numbers cannot be guessed: only a page layout
engine knows what page a heading lands on.

So the report is handed to the one program that already knows: Word opens it,
updates every field and every table of contents, saves, and closes. What is left
on disk is a document whose contents page is already right, opened by nobody.

What Word is handed is a COPY, which takes the report's place when it comes
back. A finish that fails, hangs or is killed halfway therefore leaves the
report exactly as it was written, and Word's lock file is never written beside
the user's document.

**Nothing here is required.** :func:`finalize_with_word` returns the package's
``(bool, message)`` and never raises: no Word, no automation permission, a Word
that does not answer inside :data:`TIMEOUT` — all of them are a ``False`` and a
sentence. The report is complete without this step; it is merely nicer with it.

**Nothing here is automatic.** :func:`xslope.report.generate_report` does not
call it: a headless script that writes fifty reports in a loop must not open
fifty documents in Word. Studio calls it after a generate, and a script that
wants the same finish calls it too::

    from xslope.report_finalize import finalize_with_word
    ok, msg = finalize_with_word("north_levee_report.docx")

**macOS** drives Word over Apple events (``osascript``). Word 16 for Mac does
not answer the ``close`` message in AppleScript but does in JavaScript for
Automation, so the document is opened and updated by one script and closed by
the other — see :data:`_MAC_UPDATE` and :data:`_MAC_CLOSE`.

**Windows** drives ``Word.Application`` over COM from PowerShell rather than
from ``win32com``: pywin32 is not a dependency of this package, and a COM call
made from a subprocess can be given a timeout, which an in-process one cannot.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile

__all__ = ["TIMEOUT", "word_available", "finalize_with_word"]

#: Seconds any one call into Word may take before it is abandoned. Generous:
#: Word's first launch of the day is slow, and the alternative to waiting is a
#: report with no page numbers.
TIMEOUT = 60

#: Where a Mac keeps Word. Checked before ``mdfind``, which is slower and can be
#: switched off.
_MAC_WORD_APP = "/Applications/Microsoft Word.app"
_MAC_BUNDLE_ID = "com.microsoft.Word"

# ---------------------------------------------------------------------------
# The scripts
# ---------------------------------------------------------------------------

#: Open the document, update every field in it — the body's, and the headers'
#: and footers' one section at a time, which is where the page count lives —
#: then rebuild every table of contents with its page numbers, and save.
#:
#: The document is offered to Word twice, and both offers matter. Word for Mac
#: is sandboxed: it may only read a file it has been granted, and an Apple event
#: carrying an ALIAS is what carries that grant — without the alias line, asking
#: Word to open a file it has never seen blocks until the call is abandoned.
#: Whether that first event also opens the document depends on how long Word has
#: been running, so the document is looked for among the open ones before being
#: opened by name. ``file name`` is given an HFS path (Word 16 hangs on a POSIX
#: one), converted by AppleScript itself from the argument, so a path with
#: spaces, quotes or non-ASCII characters is never pasted into a script.
#:
#: There is no ``revert``: asking Word to re-read a file it already has open
#: raises a confirmation Word never shows and never returns from. Freshness is
#: guaranteed instead by :func:`_finalize_macos`, which hands Word a copy at a
#: path it cannot already have open.
_MAC_UPDATE = r'''
on updateFieldsOf(theRange)
    -- Every field in one story, each update on its own: a field Word declines
    -- to rebuild must not cost the document the rest of them.
    set updated to 0
    tell application "Microsoft Word"
        try
            repeat with aField in (get fields of theRange)
                try
                    update field aField
                    set updated to updated + 1
                end try
            end repeat
        end try
    end tell
    return updated
end updateFieldsOf

on run argv
    set thePath to item 1 of argv
    set theAlias to (POSIX file thePath) as alias
    set hfsPath to (POSIX file thePath) as text
    tell application "Microsoft Word"
        try
            open theAlias
        end try
        set theDoc to missing value
        try
            repeat with aDoc in (get documents)
                try
                    if (full name of aDoc) is hfsPath then set theDoc to aDoc
                end try
            end repeat
        end try
        if theDoc is missing value then
            set theDoc to open file name hfsPath add to recent files false
        end if
        set fieldCount to my updateFieldsOf(theDoc)
        -- A section need not have a header or a footer of every kind, and asking
        -- for one it has not got is an error, not an empty answer.
        try
            repeat with aSection in (get sections of theDoc)
                repeat with anIndex in {header footer primary, header footer first page, header footer even pages}
                    try
                        set aPart to (get header aSection index anIndex)
                        set fieldCount to fieldCount + my updateFieldsOf(text object of aPart)
                    end try
                    try
                        set aPart to (get footer aSection index anIndex)
                        set fieldCount to fieldCount + my updateFieldsOf(text object of aPart)
                    end try
                end repeat
            end repeat
        end try
        set tocCount to 0
        try
            repeat with aToc in (get tables of contents of theDoc)
                try
                    update aToc
                    update page numbers aToc
                    set tocCount to tocCount + 1
                end try
            end repeat
        end try
        -- Word occasionally answers a save with "the object does not exist" and
        -- saves the same document happily a moment later, so the one step that
        -- puts the work on disk is the one step that is tried twice.
        set fullPath to ""
        try
            set fullPath to (full name of theDoc) as text
        end try
        try
            save theDoc
        on error
            delay 1
            save theDoc
        end try
        set saved of theDoc to true
        return fullPath & "|" & fieldCount & "|" & tocCount
    end tell
end run
'''

#: Close the window of the document that was just finalized, and nothing else.
#: The document is matched on its full name so that a report is never mistaken
#: for something of the user's that happens to be open beside it.
#:
#: This is JavaScript because Word 16 does not answer the AppleScript ``close``
#: message at all, and it tries three routes because it does not always answer
#: this one either — whether it does appears to depend on how long Word has been
#: running and how much has been asked of it. Failing to close is not a failure
#: of the finish: the document is saved, and it is the document the user is
#: about to be shown anyway.
_MAC_CLOSE = r'''
function run(argv) {
    var full = argv[0];
    var word = Application('Microsoft Word');
    var docs = word.documents;
    for (var i = 0; i < docs.length; i++) {
        var name;
        try { name = docs[i].fullName(); } catch (e) { continue; }
        if (name !== full) { continue; }
        try { docs[i].saved = true; } catch (e) {}
        var title = null;
        try { title = docs[i].windows[0].name(); } catch (e) {}
        try { docs[i].windows[0].close(); return 'closed'; } catch (e) {}
        try { docs[i].close(); return 'closed'; } catch (e) {}
        if (title !== null) {
            var wins = word.windows;
            for (var j = 0; j < wins.length; j++) {
                try {
                    if (wins[j].name() !== title) { continue; }
                    wins[j].close();
                    return 'closed';
                } catch (e) {}
            }
        }
        return 'left open';
    }
    return 'gone';
}
'''

#: The Windows leg. ``$Path`` arrives as an argument, never spliced into the
#: script.
#:
#: ``Word.Application`` is a single-instance COM server, so asking for one when
#: the user already has Word open hands back THEIRS. Two things follow: the
#: window is only forced hidden when this script started Word — hiding the
#: user's own Word would be startling — and Word is only quit on the way out
#: when this script started it. Each story range is walked through its
#: ``NextStoryRange`` chain, which is where the headers and footers of every
#: section past the first live.
_WINDOWS_SCRIPT = r'''
param([Parameter(Mandatory=$true)][string]$Path)
$ErrorActionPreference = 'Stop'
if ([Type]::GetTypeFromProgID('Word.Application') -eq $null) {
    Write-Output 'no-word'
    exit 3
}
$started = -not [bool](Get-Process -Name WINWORD -ErrorAction SilentlyContinue)
$word = New-Object -ComObject Word.Application
if ($started) { $word.Visible = $false }
$word.DisplayAlerts = 0
$doc = $null
try {
    # Open(FileName, ConfirmConversions, ReadOnly, AddToRecentFiles)
    $doc = $word.Documents.Open($Path, $false, $false, $false)
    $fields = 0
    foreach ($story in $doc.StoryRanges) {
        $range = $story
        while ($range -ne $null) {
            try { $range.Fields.Update() | Out-Null; $fields += $range.Fields.Count } catch {}
            $range = $range.NextStoryRange
        }
    }
    $tocs = 0
    foreach ($toc in $doc.TablesOfContents) {
        try { $toc.Update(); $toc.UpdatePageNumbers(); $tocs += 1 } catch {}
    }
    $doc.Save()
    $doc.Saved = $true
    Write-Output ("{0}|{1}|{2}" -f $doc.FullName, $fields, $tocs)
} finally {
    if ($doc -ne $null) { try { $doc.Close(0) } catch {} }
    if ($started) { try { $word.Quit() } catch {} }
    try { [void][Runtime.InteropServices.Marshal]::ReleaseComObject($word) } catch {}
}
'''


# ---------------------------------------------------------------------------
# Is there a Word to finalize with?
# ---------------------------------------------------------------------------

def word_available():
    """``(True, path_or_name)`` when this machine has a Word to drive, else
    ``(False, reason)``.

    The reason is written to be shown to a person: it says what is missing, not
    which call failed.
    """
    if sys.platform == "darwin":
        if os.path.isdir(_MAC_WORD_APP):
            return True, _MAC_WORD_APP
        found = _mdfind_word()
        if found:
            return True, found
        return False, "Microsoft Word is not installed."
    if sys.platform.startswith("win"):
        ok, out = _powershell(
            "if ([Type]::GetTypeFromProgID('Word.Application')) "
            "{ 'yes' } else { 'no' }", timeout=20)
        if not ok:
            return False, out
        if out.strip() == "yes":
            return True, "Word.Application"
        return False, "Microsoft Word is not installed."
    return False, ("Reports are finished in Word, which runs on macOS and "
                   "Windows only.")


def _mdfind_word():
    """Word somewhere other than /Applications, or None. Spotlight can be off,
    which is a None and not an error."""
    try:
        out = subprocess.run(
            ["mdfind", f"kMDItemCFBundleIdentifier == '{_MAC_BUNDLE_ID}'"],
            capture_output=True, text=True, timeout=15)
    except (OSError, subprocess.SubprocessError):
        return None
    for line in out.stdout.splitlines():
        if line.strip().endswith(".app"):
            return line.strip()
    return None


# ---------------------------------------------------------------------------
# The finish itself
# ---------------------------------------------------------------------------

def finalize_with_word(path, timeout=None):
    """Update the fields of the document at ``path`` in Word and save it.

    Parameters
    ----------
    path : str
        A ``.docx`` that has already been written.
    timeout : float, optional
        Seconds to wait for Word, per call. Defaults to :data:`TIMEOUT`.

    Returns
    -------
    (bool, str)
        ``(True, "…")`` when Word saved the document, else ``(False, reason)``.
        The reason is a sentence for a person; the caller decides whether it is
        worth showing at all — for Studio it is a log line, because a report
        without page numbers is still a report.
    """
    timeout = TIMEOUT if timeout is None else timeout
    if not path:
        return False, "No document was named."
    path = os.path.abspath(path)
    if not os.path.isfile(path):
        return False, f"There is no document at {path}."
    if not path.lower().endswith(".docx"):
        return False, "Only a Word document can be finished in Word."

    ok, why = word_available()
    if not ok:
        return False, why

    if sys.platform == "darwin":
        return _finalize_macos(path, timeout)
    if sys.platform.startswith("win"):
        return _finalize_windows(path, timeout)
    return False, why                       # unreachable: word_available said no


def _finalize_macos(path, timeout):
    """Word for Mac, over Apple events — and never on the user's own file.

    Word is given a copy at a path of our making, and the finished copy takes
    the report's place when it comes back. Three things follow from that, all of
    them the reason for it:

    * Word cannot already have that document open, so it cannot hand back an
      older copy of a report that has just been regenerated — and it is never
      asked to revert, which is a confirmation Word never shows and never
      returns from.
    * The lock file Word writes beside anything it opens is written beside the
      copy, and goes when the copy's directory goes. Nothing appears next to the
      user's report.
    * A Word that fails, hangs or is killed halfway leaves the report exactly as
      it was written, because it never had it.

    The copy is made in the report's own directory so that putting it back is a
    rename rather than a second write.
    """
    folder = os.path.dirname(path) or "."
    try:
        work_dir = tempfile.mkdtemp(prefix=".xslope_finalize_", dir=folder)
    except OSError as exc:
        return False, f"The report's folder is not writable: {exc}."
    work = os.path.join(work_dir, os.path.basename(path))
    try:
        shutil.copy2(path, work)
        ok, out = _osascript(_MAC_UPDATE, [work], timeout=timeout)
        if ok:
            parts = out.strip().split("|")
            full = parts[0] if parts else ""
            fields = parts[1] if len(parts) > 1 else "?"
            tocs = parts[2] if len(parts) > 2 else "?"
            if full:
                # Closing is a courtesy, not the result: the copy is saved.
                _osascript(_MAC_CLOSE, [full], timeout=timeout,
                           language="JavaScript")
            os.replace(work, path)
            return True, (f"Word updated {fields} fields and {tocs} table(s) of "
                          f"contents.")
        return False, out
    except OSError as exc:
        return False, f"The finished report could not be put back: {exc}."
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def _finalize_windows(path, timeout):
    """Word for Windows, over COM from PowerShell.

    Untested on macOS by definition; written so that every step that can fail
    fails into the same ``(False, sentence)`` the Mac leg returns.
    """
    ok, out = _powershell(_WINDOWS_SCRIPT, args=[path], timeout=timeout)
    if not ok:
        return False, out
    text = out.strip().splitlines()[-1] if out.strip() else ""
    if text == "no-word":
        return False, "Microsoft Word is not installed."
    parts = text.split("|")
    fields = parts[1] if len(parts) > 1 else "?"
    tocs = parts[2] if len(parts) > 2 else "?"
    return True, (f"Word updated {fields} fields and {tocs} table(s) of "
                  f"contents.")


# ---------------------------------------------------------------------------
# Talking to the two shells
# ---------------------------------------------------------------------------

def _osascript(script, args=(), timeout=TIMEOUT, language=None):
    """Run ``script`` through ``osascript``, with ``args`` as its arguments.

    The script arrives on standard input and the arguments as argv, so nothing
    is ever quoted into a command line. Returns ``(True, stdout)`` or
    ``(False, sentence)``: a refused automation permission and a Word that never
    answers are both ordinary outcomes here, not exceptions.
    """
    argv = ["osascript"]
    if language:
        argv += ["-l", language]
    argv += ["-"] + [str(a) for a in args]
    try:
        done = subprocess.run(argv, input=script, capture_output=True,
                              text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return False, (f"Word did not answer within {timeout:g} seconds — the "
                       f"report was left as written.")
    except OSError as exc:
        return False, f"Word could not be reached: {exc}."
    if done.returncode != 0:
        return False, _apple_event_reason(done.stderr)
    return True, done.stdout


def _apple_event_reason(stderr):
    """Turn osascript's error text into something a person can act on."""
    text = (stderr or "").strip()
    if "-1743" in text or "Not authorized to send Apple events" in text:
        return ("Word could not be controlled: allow XSLOPE Studio to control "
                "Microsoft Word under System Settings → Privacy & Security → "
                "Automation.")
    if "-600" in text or "isn't running" in text:
        return "Microsoft Word did not start."
    return text.splitlines()[-1] if text else "Word reported an unknown error."


def _powershell(script, args=(), timeout=TIMEOUT):
    """Run ``script`` through PowerShell, with ``args`` after it.

    ``-File -`` is not a thing, so the script is passed with ``-Command`` and
    reads its argument from ``$args``; the argument is a separate argv entry and
    is never spliced into the text.
    """
    argv = ["powershell", "-NoProfile", "-NonInteractive",
            "-ExecutionPolicy", "Bypass", "-Command", "-"]
    body = script
    if args:
        # The script declares param(...); calling it as a script block keeps the
        # argument out of the script text.
        quoted = " ".join(_ps_quote(str(a)) for a in args)
        body = f"& {{{script}}} {quoted}"
    try:
        done = subprocess.run(argv, input=body, capture_output=True, text=True,
                              timeout=timeout)
    except subprocess.TimeoutExpired:
        return False, (f"Word did not answer within {timeout:g} seconds — the "
                       f"report was left as written.")
    except OSError as exc:
        return False, f"Word could not be reached: {exc}."
    if done.returncode not in (0, 3):
        text = (done.stderr or done.stdout or "").strip()
        return False, (text.splitlines()[-1] if text
                       else "Word reported an unknown error.")
    return True, done.stdout


def _ps_quote(value):
    """A PowerShell single-quoted literal: the one escape is a doubled quote."""
    return "'" + value.replace("'", "''") + "'"
