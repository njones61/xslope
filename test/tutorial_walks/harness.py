"""The shared machinery of the standing tutorial walk.

A page walk reads its inputs off the page and checks its outputs against the
page, so the page is what is under test:

* **Inputs come from the page.** A table the reader pastes is taken from the
  page's own markdown (:meth:`Page.table`), a count the reader is told ("Add row
  three times") is read from the sentence (:meth:`Page.count`), and a value the
  reader types is read from where the page prints it. A page that tells the reader
  to do something that no longer works fails the walk at that step.
* **Outputs are checked at the precision the page prints.** :meth:`Walk.expect`
  finds the number in the page's text and passes only if the run rounds to it as
  written (``1.282`` is three decimals: 1.2816 passes, 1.2826 does not), the rule
  the docs tags run under.
* **A step that cannot be done as written stops the page.** :meth:`Walk.step` is
  a block: a missing control or label, a refused paste, a dialog that does not
  open or does not behave as the page says, raises inside it, and the walk
  records "cannot be done as written" and goes no further on that page, because
  every later step stands on it.

Studio runs offscreen in a settings sandbox (``tools/studio_settings.py``), on
copies of the page's files in a temporary folder (Studio and the solvers write
companion files beside a model). Dialogs are taken at ``exec()``: the walk queues
what it does in each dialog before the click that opens it, and a dialog nobody
queued an action for is itself a failure.
"""
from __future__ import annotations

import contextlib
import glob
import io
import os
import re
import shutil
import sys
import time
import traceback
from decimal import Decimal, InvalidOperation, ROUND_HALF_EVEN

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DOCS = os.path.join(REPO, "docs")
TUTORIALS = os.path.join(DOCS, "tutorials")
#: Where the page text is read from. The files a page links to are always taken
#: from the repository; the text may be read from a copy, which is how a check of
#: this check plants a stale number or a dead step without touching docs/.
PAGE_DIR = os.environ.get("XSLOPE_WALK_PAGE_DIR") or TUTORIALS


def repo_path(rel):
    """A path given relative to the repository root."""
    return os.path.join(REPO, rel)


class StepFailed(Exception):
    """A page step that cannot be done as written."""


# --------------------------------------------------------------------------- #
# The page
# --------------------------------------------------------------------------- #
_WORDS = {"once": 1, "one": 1, "twice": 2, "two": 2, "three": 3, "four": 4,
          "eleven": 11, "twelve": 12, "thirteen": 13, "fourteen": 14, "fifteen": 15,
          "sixteen": 16, "seventeen": 17, "eighteen": 18, "nineteen": 19, "twenty": 20,
          "five": 5, "six": 6, "seven": 7, "eight": 8, "nine": 9, "ten": 10}

_TABLE_RULE = re.compile(r":?-{3,}:?")


def number_text(s):
    """A number as a page prints it, as plain decimal text.

    Handles thousands separators, the typographic minus, a trailing ``%``, and
    ``a × 10<sup>b</sup>``. Returns None for text that is not a number.
    """
    t = (s or "").strip().replace("−", "-").replace("**", "").replace("`", "")
    t = t.replace("{,}", "").replace(",", "").rstrip("%").strip()
    m = re.fullmatch(r"(-?\d+(?:\.\d+)?)\s*(?:×|x)\s*10\s*<sup>\s*(-?\d+)\s*</sup>", t)
    if m:
        return "%se%s" % (m.group(1), m.group(2))
    try:
        Decimal(t)
    except (InvalidOperation, ValueError):
        return None
    return t


def rounds_to(value, written):
    """Does ``value`` round to ``written`` at the precision it is written?

    Rounding as Python's own formatting rounds (the exact binary value, half to
    even), which is what Studio and the Log print: the same rule the docs tags
    are checked under (run_tests._as_written_miss). Significant figures for
    ``1.246e-06``, decimal places otherwise."""
    want = Decimal(written)
    if "e" in written.lower():
        mant = written.lower().split("e")[0]
        digits = len(mant.replace("-", "").replace(".", "").lstrip("0"))
        got = Decimal("%.*e" % (max(digits - 1, 0), float(value)))
        return got == want, "%.*e" % (max(digits - 1, 0), float(value))
    got = Decimal(float(value)).quantize(want, rounding=ROUND_HALF_EVEN)
    return got == want, str(got)


class Page:
    """One tutorial page's markdown, searched as flowing text.

    Searches run on the text with its line breaks turned into spaces, so a phrase
    the page wraps over two lines is still one phrase; offsets are kept, so every
    match still reports its line.
    """

    def __init__(self, name):
        self.name = name
        self.path = os.path.join(PAGE_DIR, name)
        with open(self.path, encoding="utf-8") as fh:
            self.text = fh.read()
        # Every run of whitespace (a line break and the next line's indent
        # included) becomes one space, so a pattern is written with single spaces
        # whatever the page's wrapping; _orig maps each offset back to the text.
        flat, orig, prev_ws = [], [], False
        for i, ch in enumerate(self.text):
            if ch.isspace():
                if prev_ws:
                    continue
                flat.append(" ")
                prev_ws = True
            else:
                flat.append(ch)
                prev_ws = False
            orig.append(i)
        self.flat = "".join(flat)
        self._orig = orig + [len(self.text)]

    def line_of(self, offset):
        """The page line of a ``flat`` offset."""
        return self.text.count("\n", 0, self._orig[min(offset, len(self._orig) - 1)]) + 1

    @staticmethod
    def lit(text):
        """A regex for literal ``text``, any run of whitespace matching any other
        (the page may wrap it anywhere)."""
        return r"\s+".join(re.escape(p) for p in text.split())

    def finditer(self, pattern, start=0):
        return list(re.compile(pattern).finditer(self.flat, start))

    def find(self, pattern, nth=None, start=0, what=None):
        """The one match of ``pattern`` (or the ``nth``), else StepFailed."""
        hits = self.finditer(pattern, start)
        if nth is None:
            if len(hits) != 1:
                raise StepFailed("%s: expected the page to say %s once; it does %d "
                                 "time(s)" % (self.name, what or repr(pattern), len(hits)))
            return hits[0]
        if nth >= len(hits):
            raise StepFailed("%s: the page does not say %s (looked for match %d, "
                             "found %d)" % (self.name, what or repr(pattern), nth + 1,
                                            len(hits)))
        return hits[nth]

    def says(self, text):
        """Is ``text`` on the page (whitespace-insensitive)?"""
        return bool(self.finditer(self.lit(text)))

    def number(self, pattern, nth=None, group=1, start=0):
        """``(written, line)`` of the number captured by ``pattern``."""
        m = self.find(pattern, nth=nth, start=start)
        raw = m.group(group)
        written = number_text(raw)
        if written is None:
            raise StepFailed("%s:%d: %r is not a number" % (self.name, self.line_of(m.start()), raw))
        return written, self.line_of(m.start(group))

    def count(self, pattern, nth=None, start=0):
        """A count the page gives in words or digits ("three times", "twice")."""
        m = self.find(pattern, nth=nth, start=start)
        word = m.group(1).lower()
        if word.isdigit():
            return int(word)
        if word not in _WORDS:
            raise StepFailed("%s:%d: cannot read %r as a count"
                             % (self.name, self.line_of(m.start()), word))
        return _WORDS[word]

    def table(self, after, nth=0, start=0):
        """The first markdown table after the ``nth`` match of ``after`` (a regex
        on the flowing text), as rows of cells, header row first, with the cell
        markup (bold, code) taken off."""
        m = self.find(after, nth=nth, start=start)
        line0 = self.line_of(m.end())
        lines = self.text.split("\n")
        rows, started = [], False
        for line in lines[line0 - 1:]:
            s = line.strip()
            if s.startswith("|") and s.endswith("|"):
                started = True
                cells = [c.strip() for c in s.strip("|").split("|")]
                if all(_TABLE_RULE.fullmatch(c) for c in cells):
                    continue
                rows.append([c.replace("**", "").strip("`").strip() for c in cells])
            elif started:
                break
        if not rows:
            raise StepFailed("%s:%d: no table after %r" % (self.name, line0, after))
        return rows

    def code_block(self, after, nth=0, start=0):
        """The lines of the first fenced code block after the ``nth`` match of
        ``after`` (a Log excerpt the page quotes)."""
        m = self.find(after, nth=nth, start=start)
        line0 = self.line_of(m.end())
        lines = self.text.split("\n")
        out, inside = None, False
        for line in lines[line0 - 1:]:
            if line.strip().startswith("```"):
                if inside:
                    return out
                inside, out = True, []
                continue
            if inside:
                out.append(line)
        raise StepFailed("%s:%d: no code block after %r" % (self.name, line0, after))

    def tables(self):
        """Every markdown table on the page, as (line, rows)."""
        out, cur, first = [], None, None
        for i, line in enumerate(self.text.split("\n"), 1):
            s = line.strip()
            if s.startswith("|") and s.endswith("|"):
                cells = [c.strip() for c in s.strip("|").split("|")]
                if cur is None:
                    cur, first = [], i
                if not all(_TABLE_RULE.fullmatch(c) for c in cells):
                    cur.append([c.replace("**", "").strip("`").strip() for c in cells])
            elif cur is not None:
                out.append((first, cur))
                cur = None
        if cur is not None:
            out.append((first, cur))
        return out


# --------------------------------------------------------------------------- #
# The walk
# --------------------------------------------------------------------------- #
class Walk:
    """One page walked once: its failures, its notes, its files, its Studio."""

    def __init__(self, module, workdir, verbose=False):
        self.module = module
        self.page = Page(module.PAGE)
        self.workdir = workdir
        self.verbose = verbose
        self.failures = []
        self.notes = []
        self.checked = 0
        self.current = None
        self._studio = None

    # -- reporting ---------------------------------------------------------- #
    def say(self, *a):
        if self.verbose:
            print("   ", *a, file=sys.__stdout__, flush=True)

    def fail(self, message):
        msg = "%s: %s" % (self.page.name, message) if not message.startswith(self.page.name) else message
        if self.current:
            msg = "%s [%s]" % (msg, self.current)
        self.failures.append(msg)
        self.say("FAIL", msg)

    def note(self, message):
        self.notes.append(message)
        self.say("note:", message)

    @contextlib.contextmanager
    def step(self, name):
        """A page step. Anything raised inside fails the page at this step and
        stops the walk (later steps stand on this one)."""
        prev, self.current = self.current, name
        self.say("--", name)
        try:
            yield
        except StepFailed as exc:
            self.failures.append("%s: cannot be done as written at '%s': %s"
                                 % (self.page.name, name, exc))
            raise _Abort()
        except _Abort:
            raise
        except Exception as exc:
            tb = traceback.format_exc(limit=6)
            self.failures.append("%s: cannot be done as written at '%s': %s: %s"
                                 % (self.page.name, name, type(exc).__name__, exc))
            self.say(tb)
            raise _Abort()
        finally:
            self.current = prev

    def check(self, cond, message):
        """A behavior the page describes. Records a failure and goes on."""
        self.checked += 1
        if not cond:
            self.fail(message)
        return bool(cond)

    def require(self, cond, message):
        """A behavior the next steps need: fails the step."""
        self.checked += 1
        if not cond:
            raise StepFailed(message)

    def expect(self, label, value, pattern, nth=None, group=1, start=0, scale=1.0):
        """The page prints ``label`` as the number ``pattern`` captures; the run
        gives ``value``. Passes only if ``value`` rounds to it as written."""
        self.checked += 1
        try:
            written, line = self.page.number(pattern, nth=nth, group=group, start=start)
        except StepFailed as exc:
            self.fail("%s: %s" % (label, exc))
            return False
        if value is None:
            self.fail("%s:%d: the page prints %s = %s; the run gives nothing"
                      % (self.page.name, line, label, written))
            return False
        ok, got = rounds_to(float(value) * scale, written)
        if not ok:
            self.fail("%s:%d: the page prints %s = %s; the run gives %.6g, which "
                      "rounds to %s" % (self.page.name, line, label, written,
                                        float(value) * scale, got))
        else:
            self.say("ok %s = %s (run %.6g)" % (label, written, float(value) * scale))
        return ok

    def expect_written(self, label, value, written):
        """``value`` against a number already read off the page (one of several a
        single sentence prints), at the precision it is written."""
        self.checked += 1
        txt = number_text(written)
        if txt is None or value is None:
            self.fail("%s: the page prints %r; the run gives %r" % (label, written, value))
            return False
        ok, got = rounds_to(float(value), txt)
        if not ok:
            self.fail("the page prints %s = %s; the run gives %.6g, which rounds to %s"
                      % (label, txt, float(value), got))
        else:
            self.say("ok %s = %s (run %.6g)" % (label, txt, float(value)))
        return ok

    def expect_count(self, label, value, pattern, nth=None, start=0):
        """An integer the page prints, exactly (node counts, frame counts)."""
        return self.expect(label, value, pattern, nth=nth, start=start)

    def expect_about(self, label, value, pattern, within, nth=None, start=0):
        """A count the page gives as "about N" (iteration and step counts, which
        the pages say differ by a few between computers): passes within
        ``within`` of N."""
        self.checked += 1
        try:
            written, line = self.page.number(pattern, nth=nth, start=start)
        except StepFailed as exc:
            self.fail("%s: %s" % (label, exc))
            return False
        if value is None or abs(float(value) - float(written)) > within:
            self.fail("%s:%d: the page gives %s as about %s; the run gives %s (more "
                      "than %s away)" % (self.page.name, line, label, written, value, within))
            return False
        self.say("ok %s about %s (run %s)" % (label, written, value))
        return True

    def check_log(self, block, log, label):
        """Every line of a Log excerpt the page quotes is in the Log the run wrote.

        Lines are compared with runs of spaces collapsed; an elision line (``…``
        or ``...``) is skipped, and so is the time a search took (``elapsed
        time=4.17 seconds``), which no two runs share. Everything else, every
        printed number included, must be there as the page prints it.
        """
        have = {" ".join(l.split()) for l in log.splitlines()}
        flat = " ".join(log.split())
        for line in block:
            s = " ".join(line.split())
            if not s or s in ("…", "...", "⋮"):
                continue
            self.checked += 1
            if "elapsed time=" in s:
                s = s.split("elapsed time=")[0].strip()
                if s in flat:
                    continue
                self.fail("%s: the Log does not print the page's line %r" % (label, line.strip()))
                continue
            if s not in have and s not in flat:
                self.fail("%s: the Log does not print the page's line %r" % (label, line.strip()))

    # -- files ---------------------------------------------------------------- #
    def copy_model(self, rel, sub=None):
        """Copy a repository workbook and every companion sharing its stem
        (``_mesh.json``, ``_seep.csv``, ``_tseep*`` ...) into the work folder;
        return the copy's path. The repository file is never opened in place."""
        src = repo_path(rel)
        if not os.path.exists(src):
            raise StepFailed("%s does not exist" % rel)
        dst_dir = os.path.join(self.workdir, sub or "models")
        os.makedirs(dst_dir, exist_ok=True)
        stem = os.path.splitext(src)[0]
        for f in glob.glob(glob.escape(stem) + "*"):
            if os.path.isfile(f):
                rest = os.path.basename(f)[len(os.path.basename(stem)):]
                # Only companions: the stem followed by '_' or '.', not a longer name.
                if rest.startswith((".", "_")):
                    shutil.copy2(f, dst_dir)
        return os.path.join(dst_dir, os.path.basename(src))

    def page_link(self, name):
        """The repository path of a file the page links to (``files/x.xlsx``,
        ``../lem/files/x.xlsx``, ``../inputs/x.xlsx``); fails if the page does not
        link it."""
        hits = re.findall(r"\]\(((?:\.\./)*[^)\s:#]*?(?<![\w.-])%s)\)" % re.escape(name),
                          self.page.text)
        if not hits:
            raise StepFailed("the page does not link %s" % name)
        return os.path.relpath(os.path.normpath(os.path.join(TUTORIALS, hits[0])), REPO)

    # -- engine --------------------------------------------------------------- #
    @staticmethod
    def load(path):
        from xslope.fileio import load_slope_data
        with contextlib.redirect_stdout(io.StringIO()):
            return load_slope_data(path)

    @staticmethod
    def quiet(fn, *a, **k):
        with contextlib.redirect_stdout(io.StringIO()):
            return fn(*a, **k)

    def lem(self, sd, method, analysis="auto_search", num_slices=40, **kw):
        """``run_lem_analysis``, the call Studio's Run LEM makes."""
        import copy
        from xslope.search import run_lem_analysis
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            b = run_lem_analysis(copy.deepcopy(sd), method, analysis=analysis,
                                 num_slices=num_slices, **kw)
        b = dict(b)
        b["_log"] = buf.getvalue()
        return b

    # -- Studio --------------------------------------------------------------- #
    def studio(self):
        """The page's Studio window (one per page), built on first use."""
        if self._studio is None:
            self._studio = Studio(self)
        return self._studio

    def close(self):
        if self._studio is not None:
            self._studio.close()
            self._studio = None


class _Abort(Exception):
    """Unwinds a page walk after a failed step."""


# --------------------------------------------------------------------------- #
# Studio, offscreen
# --------------------------------------------------------------------------- #
def _norm(s):
    return " ".join((s or "").replace("&", "").replace("…", "").replace("...", "")
                    .rstrip(":").split()).lower()


class DialogQueue:
    """Dialogs taken at ``exec()``: each pops the next queued action.

    An action returns None to accept the dialog (as clicking its OK/Run), or the
    string ``"reject"`` to cancel it. A dialog with nothing queued for it is an
    unexpected dialog, which fails the step that opened it. Message boxes shown
    through the static helpers are recorded (``boxes``) and answered OK / Yes.
    """

    def __init__(self):
        self.queue = []
        self.seen = []
        self.errors = []
        self.boxes = []
        self._saved = {}

    def install(self):
        from PySide6.QtWidgets import QDialog, QMessageBox
        self._saved["exec"] = QDialog.exec
        self._saved["exec_"] = getattr(QDialog, "exec_", None)
        q = self

        def _exec(dlg, *a, **k):
            q.seen.append((type(dlg).__name__, dlg.windowTitle()))
            if not q.queue:
                q.errors.append("an unexpected dialog opened: %s '%s'"
                                % (type(dlg).__name__, dlg.windowTitle()))
                dlg.reject()
                return 0
            action = q.queue.pop(0)
            try:
                out = action(dlg)
            except Exception as exc:
                q.errors.append("in %s '%s': %s: %s" % (type(dlg).__name__,
                                                         dlg.windowTitle(),
                                                         type(exc).__name__, exc))
                q.errors.append(traceback.format_exc(limit=5))
                dlg.reject()
                return 0
            if out == "reject":
                dlg.reject()
                return 0
            if out == "keep":       # the action already closed the dialog itself
                return dlg.result()
            dlg.accept()
            return dlg.result()

        QDialog.exec = _exec
        QDialog.exec_ = _exec
        # Nothing a walk does may open a document or a web page on the desktop:
        # a generated report, a docs link. Each request is recorded instead.
        from PySide6.QtGui import QDesktopServices
        self._saved["openUrl"] = QDesktopServices.openUrl
        self.urls = []
        QDesktopServices.openUrl = staticmethod(lambda url, *a: q.urls.append(url) or True)
        for name in ("information", "warning", "critical", "question"):
            self._saved[name] = getattr(QMessageBox, name)

            def _box(*a, _n=name, **k):
                title = a[1] if len(a) > 1 else ""
                text = a[2] if len(a) > 2 else ""
                q.boxes.append((_n, str(title), str(text)))
                if _n != "question":
                    return QMessageBox.Ok
                buttons = a[3] if len(a) > 3 else k.get("buttons")
                try:
                    if buttons is not None and (buttons & QMessageBox.Discard):
                        return QMessageBox.Discard     # "Unsaved changes": discard
                except TypeError:
                    pass
                return QMessageBox.Yes
            setattr(QMessageBox, name, staticmethod(_box))

    def uninstall(self):
        from PySide6.QtWidgets import QDialog, QMessageBox
        if "exec" in self._saved:
            QDialog.exec = self._saved["exec"]
            if self._saved.get("exec_") is not None:
                QDialog.exec_ = self._saved["exec_"]
        for name in ("information", "warning", "critical", "question"):
            if name in self._saved:
                setattr(QMessageBox, name, self._saved[name])
        if "openUrl" in self._saved:
            from PySide6.QtGui import QDesktopServices
            QDesktopServices.openUrl = self._saved["openUrl"]
        self._saved.clear()


class Studio:
    """The real Studio main window, driven the way the page tells the reader."""

    def __init__(self, walk):
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        os.environ.setdefault("XSLOPE_NO_UPDATE_CHECK", "1")
        # A generated report's contents page is finished headlessly, never by
        # driving Word on the desktop the walk runs on.
        os.environ.setdefault("XSLOPE_FINALIZE", "libreoffice")
        from PySide6.QtWidgets import QApplication
        self.app = QApplication.instance() or QApplication([])
        self.walk = walk
        self.dialogs = DialogQueue()
        self.dialogs.install()
        from studio.main_window import MainWindow
        self.mw = MainWindow()
        self.mw.resize(1600, 1000)
        self.mw.show()
        self.settle()

    def close(self):
        try:
            self.mw.doc.mark_clean() if hasattr(self.mw.doc, "mark_clean") else None
        except Exception:
            pass
        try:
            self.mw.close()
        except Exception:
            pass
        self.settle(3)
        self.dialogs.uninstall()

    # -- event loop ----------------------------------------------------------- #
    def settle(self, n=10):
        for _ in range(n):
            self.app.processEvents()
            time.sleep(0.01)

    def wait(self, cond, timeout, label):
        t0 = time.time()
        while not cond():
            self.app.processEvents()
            time.sleep(0.02)
            if time.time() - t0 > timeout:
                raise StepFailed("%s did not finish in %d s" % (label, timeout))
        self.settle(5)
        return time.time() - t0

    def _drain_errors(self):
        errs = [e for e in self.dialogs.errors if not e.startswith("Traceback")]
        tbs = [e for e in self.dialogs.errors if e.startswith("Traceback")]
        self.dialogs.errors.clear()
        if tbs:
            self.walk.say("\n".join(tbs))
        return errs

    def run_queued(self, actions, trigger, label):
        """Queue one action per dialog the click will open, click, and require
        that every queued dialog opened and its action ran cleanly."""
        self.dialogs.queue[:] = list(actions)
        trigger()
        self.settle()
        left = len(self.dialogs.queue)
        self.dialogs.queue[:] = []
        errs = self._drain_errors()
        if errs:
            raise StepFailed("%s: %s" % (label, "; ".join(errs)))
        if left:
            raise StepFailed("%s: %d dialog(s) the page describes did not open"
                             % (label, left))

    # -- the window ----------------------------------------------------------- #
    @property
    def doc(self):
        return self.mw.doc

    @property
    def data(self):
        return self.mw.doc.slope_data

    def log_text(self):
        return self.mw.log.toPlainText()

    def mark(self):
        return len(self.log_text())

    def log_since(self, mark):
        return self.log_text()[mark:]

    def status(self):
        return self.mw.statusBar().currentMessage()

    def open(self, path):
        self.run_queued([], lambda: self.mw.open_path(path), "open %s" % os.path.basename(path))
        self.settle(20)
        if not self.doc.is_open:
            raise StepFailed("%s did not open" % os.path.basename(path))

    def menu_action(self, menu, text):
        from PySide6.QtWidgets import QMenu
        want = _norm(text)
        for m in self.mw.menuBar().findChildren(QMenu):
            if _norm(m.title()) != _norm(menu):
                continue
            for a in m.actions():
                if _norm(a.text()) == want:
                    return a
        raise StepFailed("there is no %s → %s menu item" % (menu, text))

    def trigger(self, menu, text, *actions):
        a = self.menu_action(menu, text)
        if not a.isEnabled():
            raise StepFailed("%s → %s is grayed out" % (menu, text))
        self.run_queued(actions, a.trigger, "%s → %s" % (menu, text))

    def set_mode(self, name):
        from studio.main_window import MODES
        idx = [m[0] for m in MODES].index(name)
        self.mw.set_mode_index(idx)
        self.settle()

    @property
    def mode(self):
        return self.mw._mode

    def tree_items(self):
        from PySide6.QtWidgets import QTreeWidgetItemIterator
        it = QTreeWidgetItemIterator(self.mw.inputs_tree)
        out = {}
        while it.value():
            item = it.value()
            out[item.text(0)] = item
            it += 1
        return out

    def tree_count(self, label):
        items = self.tree_items()
        if label not in items:
            raise StepFailed("the Inputs tree has no '%s' (it has %s)" % (label, sorted(items)))
        return items[label].text(1)

    def edit(self, label, *actions):
        """Click ``label`` in the Inputs tree; ``actions`` run in the dialog(s)
        it opens."""
        items = self.tree_items()
        if label not in items:
            raise StepFailed("the Inputs tree has no '%s' (it has %s)" % (label, sorted(items)))
        item = items[label]
        self.run_queued(actions, lambda: self.mw.inputs_tree.itemClicked.emit(item, 0),
                        "Inputs → %s" % label)
        self.settle()

    # -- runs ----------------------------------------------------------------- #
    def _run(self, menu_text, runner_attr, action, timeout, label):
        mark = self.mark()
        t0 = time.time()
        self.trigger("Run", menu_text, action)
        self.wait(lambda: getattr(self.mw, runner_attr) is None, timeout, label)
        return time.time() - t0, self.log_since(mark)

    def run_lem(self, action, timeout=1800):
        return self._run("Run LEM", "_runner", action, timeout, "the LEM run")

    def run_seep(self, action, timeout=3600):
        return self._run("Run Seep", "_seep_runner", action, timeout, "the seepage run")

    def run_fem(self, action, timeout=4 * 3600):
        return self._run("Run FEM", "_fem_runner", action, timeout, "the FEM run")

    def run_parametric(self, action, timeout=4 * 3600):
        return self._run("Parametric", "_sens_runner", action, timeout, "the parametric run")

    def run_reliability(self, action, timeout=4 * 3600):
        return self._run("Reliability", "_rel_runner", action, timeout, "the reliability run")

    def build_mesh(self, action, timeout=1800):
        mark = self.mark()
        t0 = time.time()
        self.trigger("Run", "Build Mesh", action)
        self.wait(lambda: not self.mw._mesh_busy, timeout, "the mesh build")
        return time.time() - t0, self.log_since(mark)

    @property
    def results(self):
        return self.doc.results


# --------------------------------------------------------------------------- #
# Dialog helpers: controls found by the words the page uses for them
# --------------------------------------------------------------------------- #
def _qt():
    from PySide6 import QtWidgets
    return QtWidgets


def labels(dlg, visible=True):
    W = _qt()
    return [l for l in dlg.findChildren(W.QLabel)
            if l.text() and (l.isVisibleTo(dlg) or not visible)]


def field(dlg, label, nth=0):
    """The input widget the visible label ``label`` names: a form row's field, a
    label's buddy, or the widget beside it in its layout."""
    W = _qt()
    key = _norm(label)
    hits = [l for l in labels(dlg) if _norm(l.text()) == key]
    if not hits:
        hits = [l for l in labels(dlg) if _norm(l.text()).startswith(key)]
    if not hits:
        raise StepFailed("%s has no field labeled '%s' (labels: %s)"
                         % (_title(dlg), label, sorted({_norm(l.text()) for l in labels(dlg)})[:60]))
    lab = hits[nth]
    if lab.buddy() is not None:
        return lab.buddy()
    parent = lab.parentWidget()
    lay = parent.layout() if parent is not None else None

    def find(l):
        for i in range(l.count()):
            it = l.itemAt(i)
            if it.widget() is lab:
                return l, i
            if it.layout() is not None:
                r = find(it.layout())
                if r:
                    return r
        return None
    found = find(lay) if lay is not None else None
    if found is None:
        raise StepFailed("cannot find the control beside '%s'" % label)
    l, i = found
    if isinstance(l, W.QFormLayout):
        r, _role = l.getItemPosition(i)
        item = l.itemAt(r, W.QFormLayout.FieldRole)
        w = item.widget() if item is not None else None
        if w is None and item is not None and item.layout() is not None:
            w = item.layout().itemAt(0).widget()
        if w is not None:
            return w
    if isinstance(l, W.QGridLayout):
        r, c, _rs, _cs = l.getItemPosition(i)
        item = l.itemAtPosition(r, c + 1)
        if item is not None and item.widget() is not None:
            return item.widget()
    for j in range(i + 1, l.count()):
        w = l.itemAt(j).widget()
        if w is not None:
            return w
    raise StepFailed("nothing beside '%s'" % label)


def _title(dlg):
    return "%s '%s'" % (type(dlg).__name__, dlg.windowTitle())


def button(dlg, text, kind=None, visible=True):
    W = _qt()
    kind = kind or W.QAbstractButton
    key = _norm(text)
    cands = [b for b in dlg.findChildren(kind) if (b.isVisibleTo(dlg) or not visible)]
    hits = [b for b in cands if _norm(b.text()) == key]
    if not hits:
        hits = [b for b in cands if _norm(b.text()).startswith(key)]
    if len(hits) != 1:
        raise StepFailed("%s: %d button(s) read '%s' (buttons: %s)"
                         % (_title(dlg), len(hits), text,
                            [b.text() for b in cands if b.text()][:40]))
    return hits[0]


def press(dlg, text, settle=True):
    b = button(dlg, text)
    if not b.isEnabled():
        raise StepFailed("%s: '%s' is grayed out" % (_title(dlg), text))
    b.click()
    if settle:
        _process()
    return b


def tick(dlg, text, on=True):
    W = _qt()
    b = button(dlg, text, W.QAbstractButton)
    if not b.isEnabled():
        raise StepFailed("%s: '%s' is grayed out" % (_title(dlg), text))
    if b.isChecked() != bool(on):
        b.click()
        _process()
    return b


def set_value(w, value):
    W = _qt()
    if isinstance(w, W.QComboBox):
        i = w.findText(str(value))
        if i < 0:
            i = next((j for j in range(w.count())
                      if _norm(w.itemText(j)) == _norm(str(value))), -1)
        if i < 0:
            i = next((j for j in range(w.count())
                      if _norm(w.itemText(j)).startswith(_norm(str(value)))), -1)
        if i < 0:
            raise StepFailed("the list has no '%s' (it has %s)"
                             % (value, [w.itemText(j) for j in range(w.count())]))
        w.setCurrentIndex(i)
    elif isinstance(w, W.QDoubleSpinBox):
        w.setValue(float(value))
    elif isinstance(w, W.QSpinBox):
        w.setValue(int(float(value)))
    elif isinstance(w, W.QLineEdit):
        w.setFocus()
        w.selectAll()
        w.clear()
        w.insert(str(value))
        w.editingFinished.emit()
    elif isinstance(w, W.QAbstractButton):
        w.setChecked(bool(value))
    else:
        raise StepFailed("cannot type into a %s" % type(w).__name__)
    _process()


def value(w):
    W = _qt()
    if isinstance(w, W.QComboBox):
        return w.currentText()
    if isinstance(w, (W.QSpinBox, W.QDoubleSpinBox)):
        return w.value()
    if isinstance(w, W.QLineEdit):
        return w.text()
    if isinstance(w, W.QAbstractButton):
        return w.isChecked()
    if isinstance(w, W.QLabel):
        return w.text()
    return None


def set_field(dlg, label, val, nth=0):
    w = field(dlg, label, nth=nth)
    if not w.isEnabled():
        raise StepFailed("%s: '%s' is grayed out" % (_title(dlg), label))
    set_value(w, val)
    return w


def _process(n=3):
    from PySide6.QtWidgets import QApplication
    for _ in range(n):
        QApplication.processEvents()


def show_parameters_for(dlg, wanted):
    """Set an editor's **Show parameters for:** toggles to exactly ``wanted``
    (the words on the checkboxes: LEM, Seepage, FEM, Reliability)."""
    toggles = getattr(dlg, "_toggles", None)
    if not toggles:
        raise StepFailed("%s has no Show parameters for: toggles" % _title(dlg))
    names = {cb.text(): cb for cb in toggles.values()}
    for w in wanted:
        if w not in names:
            raise StepFailed("%s has no '%s' toggle (it has %s)" % (_title(dlg), w, sorted(names)))
    for text, cb in names.items():
        if cb.isChecked() != (text in wanted):
            cb.click()
            _process()


def list_rows(lw):
    """The texts of a list widget's rows."""
    return [lw.item(i).text() for i in range(lw.count())]


def select_row(lw, prefix):
    """Select the list row that starts with ``prefix`` (as a reader clicks it)."""
    rows_ = list_rows(lw)
    hits = [i for i, t in enumerate(rows_) if _norm(t).startswith(_norm(prefix))]
    if len(hits) != 1:
        raise StepFailed("the list has %d row(s) starting '%s' (rows: %s)"
                         % (len(hits), prefix, rows_))
    lw.setCurrentRow(hits[0])
    _process()
    return hits[0]


def seep_bc_set(dlg, index=0):
    """The Seep BC editor's set ``index`` (0 = Set 1), its tab made current."""
    W = _qt()
    tabs = dlg.findChildren(W.QTabWidget)
    if not tabs:
        raise StepFailed("%s has no set tabs" % _title(dlg))
    tabs[0].setCurrentIndex(index)
    _process()
    return tabs[0].widget(index)


def editable(dlg, nth=0):
    """The visible editable table (Studio's ``_EditableTable``) in a dialog."""
    from studio.editors import _EditableTable
    ts = [t for t in dlg.findChildren(_EditableTable) if t.isVisibleTo(dlg)]
    if len(ts) <= nth:
        raise StepFailed("%s shows no table" % _title(dlg))
    return ts[nth]


def headers(table):
    t = getattr(table, "table", table)
    hh = t.horizontalHeader()
    cols = sorted([c for c in range(t.columnCount()) if not t.isColumnHidden(c)],
                  key=hh.visualIndex)
    return [(t.horizontalHeaderItem(c).text() if t.horizontalHeaderItem(c) else "")
            for c in cols], cols


def column(table, name):
    """The (logical) column whose header reads ``name``, or starts with
    ``name (`` (a unit)."""
    hdr, cols = headers(table)
    for h, c in zip(hdr, cols):
        if _norm(h) == _norm(name):
            return c
    for h, c in zip(hdr, cols):
        if _norm(h).startswith(_norm(name) + " ("):
            return c
    raise StepFailed("the table has no column '%s' (columns: %s)" % (name, hdr))


def rows(table):
    """Every row as text, visible columns in visual order (combo cells read)."""
    W = _qt()
    t = getattr(table, "table", table)
    _hdr, cols = headers(t)
    out = []
    for r in range(t.rowCount()):
        row = []
        for c in cols:
            w = t.cellWidget(r, c)
            if w is not None:
                cb = w if isinstance(w, W.QComboBox) else w.findChild(W.QComboBox)
                if cb is not None:
                    row.append(cb.currentText())
                    continue
                le = w.findChild(W.QLineEdit)
                row.append(le.text() if le is not None else "")
                continue
            it = t.item(r, c)
            row.append(it.text() if it is not None else "")
        out.append(row)
    return out


def paste(table, block, row=0, col=0):
    """Paste ``block`` (rows of cells) with a real Cmd/Ctrl+V at (row, col), as a
    reader pastes a table copied off the page."""
    from PySide6.QtCore import QEvent, Qt
    from PySide6.QtGui import QKeyEvent
    from PySide6.QtWidgets import QApplication
    t = getattr(table, "table", table)
    QApplication.clipboard().setText("\n".join("\t".join(str(c) for c in r) for r in block))
    t.setCurrentCell(row, col)
    QApplication.sendEvent(t, QKeyEvent(QEvent.KeyPress, Qt.Key_V, Qt.ControlModifier, "v"))
    _process()


def type_cells(table, block, row=0, col=0):
    """Type ``block`` into existing cells, as a reader types row by row. A row the
    table does not have fails the step."""
    t = getattr(table, "table", table)
    _hdr, cols = headers(t)
    for i, r in enumerate(block):
        if row + i >= t.rowCount():
            raise StepFailed("the table has %d row(s); the page has the reader type "
                             "%d" % (t.rowCount(), row + len(block)))
        for j, v in enumerate(r):
            c = cols[cols.index(col) + j] if col in cols else col + j
            it = t.item(row + i, c)
            if it is None:
                raise StepFailed("cell (%d, %d) cannot be typed into" % (row + i, c))
            it.setText(str(v))
    _process()


def page_rows(table_rows, drop_header=True):
    """A page table as the block a reader copies: header dropped, cells as text."""
    return [list(r) for r in (table_rows[1:] if drop_header else table_rows)]


def preflight_items(dlg):
    """What a Run dialog's Model checks panel lists, as text."""
    W = _qt()
    pf = getattr(dlg, "preflight", None)
    if pf is None:
        return None
    out = []
    lst = getattr(pf, "_list", None)
    if lst is not None:
        for i in range(lst.count()):
            it = lst.item(i)
            t = it.text()
            wid = lst.itemWidget(it)
            if wid is not None:
                labs = [l.text() for l in wid.findChildren(W.QLabel) if l.text()]
                if isinstance(wid, W.QLabel) and wid.text():
                    labs.insert(0, wid.text())
                t = (t + " " + " | ".join(labs)).strip()
            out.append(t)
    try:
        out.insert(0, pf.title())
    except Exception:
        pass
    # Each entry's full text, which the panel shows beneath the list when the
    # entry is selected (the list line may be a shortened form of it).
    rep = getattr(pf, "report", None)
    for f in (getattr(rep, "findings", None) or ()):
        msg = getattr(f, "message", None)
        if msg:
            out.append("%s: %s" % (getattr(f, "severity", ""), msg))
    return out


def run_button_enabled(dlg):
    W = _qt()
    for bb in dlg.findChildren(W.QDialogButtonBox):
        b = bb.button(W.QDialogButtonBox.Ok)
        if b is not None:
            return b.isEnabled()
    return True
