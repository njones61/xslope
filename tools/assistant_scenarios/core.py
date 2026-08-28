"""The pieces a scored assistant scenario is built from.

A scenario is a starting workbook, an optional set of planted faults, one or more
prompts, and a list of CRITERIA. Running it records a real conversation through
``tools/assistant_sessions.run_assistant_session``; scoring it reads the workbook
the conversation left behind and the transcript it produced, and asks each
criterion one yes/no question with a one-line reason.

Nothing here calls a provider. The runner (``tools/assistant_suite.py``) does that;
this module is what the runner scores WITH, which is why replaying a saved run
re-scores it exactly — every criterion reads files, and re-solves the model itself
rather than believing anything the assistant said about it.

The ground truth is always an independent re-run: the workbook is reloaded from
disk and put through :func:`xslope.search.run_lem_analysis` (or the seepage / FEM
entry points ``run_tests`` uses) in this process, with no assistant in the loop.
Where the file also carries a published lock — a ``<!-- test: ... -->`` tag on a
docs page — that lock is checked against the independent run as a second opinion,
so a model the assistant quietly corrupted cannot pass by agreeing with itself.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
import os
import re
import sys
import warnings

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)


def repo(*parts):
    """An absolute path inside the repository."""
    return os.path.join(REPO_ROOT, *parts)


# --------------------------------------------------------------------------- #
# The scenario record
# --------------------------------------------------------------------------- #
class Scenario:
    """One scored conversation.

    Parameters
    ----------
    name : str
        Unique; every artifact the run writes is named from it.
    family : str
        The group it is reported under (``build``, ``edit``, ``diagnose``, ...).
    model : str or None
        The workbook the session opens, absolute. ``None`` starts from File → New.
    turns : list
        Prompts, in order; an item may be ``(text, image_path)``.
    criteria : list of Criterion
    faults : list of Fault, optional
        Written into the copy the session opens. Every scenario opens a copy of
        ``model`` under the run's own directory — faults or none — so a session
        can neither corrupt a repository model nor write a file beside one.
    what : str
        One line saying what the scenario is for; printed in the scorecard.
    """

    def __init__(self, name, family, model, turns, criteria, faults=(), what="",
                 timeout_s=900, save_each_turn=False, settings=None,
                 images=None, save_after=None):
        self.name = name
        self.family = family
        self.model = model
        self.turns = list(turns)
        self.criteria = list(criteria)
        self.faults = list(faults)
        self.what = what
        self.timeout_s = timeout_s
        self.save_each_turn = save_each_turn
        self.settings = dict(settings or {})
        self.images = images
        self.save_after = save_after

    def __repr__(self):
        return "<Scenario %s (%s), %d turn(s), %d criteria>" % (
            self.name, self.family, len(self.turns), len(self.criteria))


class Criterion:
    """One yes/no question about a finished session.

    ``check(ctx)`` returns ``(ok, reason)``. The reason is written whether the
    criterion passed or failed — a passing criterion that says WHAT it measured is
    what makes the scorecard readable ("FS 1.5867 vs re-run 1.5867").
    """

    def __init__(self, name, check, kind="behavior"):
        self.name = name
        self.check = check
        self.kind = kind

    def __call__(self, ctx):
        try:
            ok, reason = self.check(ctx)
        except Exception as exc:                      # a scorer must never abort a run
            return False, "scorer raised %s: %s" % (type(exc).__name__, exc)
        return bool(ok), str(reason)


class Fault:
    """One deliberate error written into a copy of a model.

    ``apply(slope_data)`` mutates the model in place. ``names`` is the pattern that
    recognizes the assistant NAMING this fault in its prose — the fault is found
    only when the answer says so, not when a snippet happens to print the value.
    """

    def __init__(self, ident, describe, apply, names):
        self.id = ident
        self.describe = describe
        self.apply = apply
        self.names = names if hasattr(names, "search") else re.compile(names, re.I)


# --------------------------------------------------------------------------- #
# Reading a finished session
# --------------------------------------------------------------------------- #
#: Markers the recorder writes at column 0 in a turn's body.
_MARKERS = ("Assistant: ", "Ran code:", "Output:", "[file: ",
            "Declined to run the code.", "Error: ")


def parse_body(body):
    """Split one turn's recorded body into ``(texts, codes, outputs, files)``.

    The recorder (``tools/assistant_sessions._Recorder``) writes the assistant's
    prose, every ``run_python`` snippet, and every snippet's printed result — the
    MODEL CHECKS block included — into one block of text. Scoring needs the three
    apart: what was CLAIMED (prose), what was RUN (code), and what was MEASURED
    (output).
    """
    texts, codes, outputs, files = [], [], [], []
    current, sink = None, None
    for line in (body or "").splitlines():
        marker = next((m for m in _MARKERS if line.startswith(m)), None)
        if marker is not None:
            if current is not None and sink is not None:
                sink.append("\n".join(current).strip())
            current, sink = None, None
            if marker == "Assistant: ":
                current, sink = [line[len(marker):]], texts
            elif marker == "Ran code:":
                current, sink = [], codes
            elif marker == "Output:":
                current, sink = [], outputs
            elif marker == "[file: ":
                files.append(line[len("[file: "):].rstrip("]"))
            continue
        if current is not None:
            current.append(line[4:] if line.startswith("    ") else line)
    if current is not None and sink is not None:
        sink.append("\n".join(current).strip())
    return texts, codes, outputs, files


class Turn:
    """One recorded turn, parsed."""

    def __init__(self, record):
        self.prompt = record.get("prompt") or ""
        self.body = record.get("body") or ""
        self.usage = dict(record.get("usage") or {})
        self.seconds = record.get("seconds") or 0.0
        self.error = record.get("error")
        self.workbook = record.get("workbook")
        self.image = record.get("image")
        self.texts, self.codes, self.outputs, self.files = parse_body(self.body)

    @property
    def prose(self):
        return "\n\n".join(self.texts)

    @property
    def code(self):
        return "\n\n".join(self.codes)

    @property
    def output(self):
        return "\n\n".join(self.outputs)

    @property
    def completions(self):
        return int(self.usage.get("calls") or 0)


class ScoreCtx:
    """Everything a criterion may look at, and the engine calls it scores with."""

    def __init__(self, scenario, session, outdir):
        self.scenario = scenario
        self.session = session
        self.outdir = outdir
        self.turns = [Turn(r) for r in session.get("turns") or []]
        self.usage = dict(session.get("usage") or {})
        self.seconds = session.get("seconds") or 0.0
        self.error = session.get("error")
        self.start_model = session.get("start_model")
        # Every session opens a COPY of its model, so the path a published lock
        # is indexed under is not the path that was opened.
        self.lock_model = session.get("lock_model") or self.start_model
        self.step_workbooks = list(session.get("step_workbooks") or [])
        # None = nothing changed — EXCEPT on a session that saved per turn, where
        # the last turn's file already holds the end state and the end-of-session
        # save is skipped because the document is no longer dirty.
        self.workbook = session.get("workbook") or (self.step_workbooks[-1]
                                                    if self.step_workbooks
                                                    else None)
        self.transcript = session.get("transcript")

    # -- the words and the code -------------------------------------------
    @property
    def prose(self):
        return "\n\n".join(t.prose for t in self.turns)

    @property
    def code(self):
        return "\n\n".join(t.code for t in self.turns)

    @property
    def output(self):
        return "\n\n".join(t.output for t in self.turns)

    @property
    def final_prose(self):
        for turn in reversed(self.turns):
            if turn.prose.strip():
                return turn.prose
        return ""

    # -- the model, before and after --------------------------------------
    def saved(self, turn=None):
        """The workbook this session left behind (or after ``turn``), or None."""
        if turn is None:
            return self.workbook
        if 1 <= turn <= len(self.turns):
            return self.turns[turn - 1].workbook
        return None

    def before(self):
        return load_model(self.start_model) if self.start_model else None

    def after(self, turn=None):
        path = self.saved(turn)
        return load_model(path) if path else None

    def model_changed(self):
        return bool(self.workbook)


# --------------------------------------------------------------------------- #
# The engine, driven independently of the assistant
# --------------------------------------------------------------------------- #
_SOLVE_CACHE_PATH = None
_solve_cache = {}


def use_solve_cache(path):
    """Point the independent-run cache at ``path`` and load what it already holds.

    Re-solving is the expensive half of scoring — a search is seconds to minutes —
    and a replay re-scores the same files. Keying on the file's own SHA-256 means a
    cache entry can never be served for a workbook that changed.
    """
    global _SOLVE_CACHE_PATH, _solve_cache
    _SOLVE_CACHE_PATH = path
    if path and os.path.exists(path):
        try:
            with open(path, encoding="utf-8") as fh:
                _solve_cache = json.load(fh)
        except Exception:
            _solve_cache = {}


def _cache_put(key, value):
    _solve_cache[key] = value
    if _SOLVE_CACHE_PATH:
        with contextlib.suppress(Exception):
            os.makedirs(os.path.dirname(_SOLVE_CACHE_PATH), exist_ok=True)
            with open(_SOLVE_CACHE_PATH, "w", encoding="utf-8") as fh:
                json.dump(_solve_cache, fh, indent=1, sort_keys=True)


def digest(path):
    """SHA-256 of a file, or ``''``."""
    try:
        with open(path, "rb") as fh:
            return hashlib.sha256(fh.read()).hexdigest()
    except Exception:
        return ""


@contextlib.contextmanager
def _quiet():
    """Engine chatter off — the solvers announce themselves on stdout."""
    buf = io.StringIO()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with contextlib.redirect_stdout(buf):
            yield buf


def load_model(path):
    """``load_slope_data`` on ``path``, quietly; ``None`` if it will not load."""
    from xslope.fileio import load_slope_data
    try:
        with _quiet():
            return load_slope_data(path)
    except Exception:
        return None


def declared_method(slope_data):
    """The method the MODEL declares, or ``spencer`` where it declares none —
    the same rule ``run_lem`` follows, so a scorer asks the same question the
    assistant was answering."""
    method = (slope_data or {}).get("lem_method")
    method = str(method or "").strip().lower()
    return method or "spencer"


def solve(path, method=None, search=True, num_slices=40, rapid=False):
    """Run one LEM analysis on ``path`` independently and return a small dict.

    ``{'FS', 'Xo', 'Yo', 'R', 'method', 'error'}``. This is the ground truth every
    numeric criterion compares against: it reloads the file from disk and drives
    :func:`xslope.search.run_lem_analysis` — the same entry point Studio's Run LEM
    dialog and the assistant's ``run_lem`` both use — with nothing carried over
    from the session.
    """
    return solve_variant(path, None, "", method=method, search=search,
                         num_slices=num_slices, rapid=rapid)


def solve_variant(path, apply=None, key="", method=None, search=True,
                  num_slices=40, rapid=False, want_slices=False):
    """:func:`solve`, with ``apply(slope_data)`` made to the model first.

    This is how an EDIT or a SWEEP is checked. The assistant is asked to change one
    input and re-run; the independent answer to that is the same file, reloaded
    from disk, changed the same way here, and re-solved — with nothing of the
    session's carried over. ``key`` names the change and rides in the cache key
    beside the file's own digest, so two variants of one file never collide and a
    cached answer can never be served for a workbook that moved.

    ``want_slices`` returns the solved slice table as well; it is not cached,
    because a DataFrame is not JSON.
    """
    from xslope.search import run_lem_analysis

    sd = load_model(path)
    if sd is None:
        return {"FS": None, "error": "the workbook does not load"}
    method = (method or declared_method(sd)).lower()
    cache_key = "|".join(["lem", digest(path), key, method, str(bool(search)),
                          str(num_slices), str(bool(rapid))])
    if cache_key in _solve_cache and not want_slices:
        return dict(_solve_cache[cache_key])
    out = {"FS": None, "Xo": None, "Yo": None, "R": None,
           "method": method, "error": None}
    try:
        if apply is not None:
            apply(sd)
        with _quiet():
            bundle = run_lem_analysis(
                sd, method,
                analysis="auto_search" if search else "single_surface",
                surface="noncircular" if (not sd.get("circular")
                                          and sd.get("non_circ")) else "circular",
                num_slices=num_slices, rapid=rapid, announce=False)
        results = bundle.get("results") or {}
        out["FS"] = results.get("FS")
        cache = ((bundle.get("search") or {}).get("fs_cache") or [None])[0]
        if isinstance(cache, dict):
            out["Xo"], out["Yo"] = cache.get("Xo"), cache.get("Yo")
            if cache.get("Yo") is not None and cache.get("Depth") is not None:
                out["R"] = cache["Yo"] - cache["Depth"]
        if out["FS"] is None:
            out["error"] = str(results.get("failure") or "no solution")
        if want_slices:
            out["slice_df"] = bundle.get("slice_df")
    except Exception as exc:
        out["error"] = "%s: %s" % (type(exc).__name__, exc)
    if want_slices:
        return out
    _cache_put(cache_key, out)
    return dict(out)


def seep_flow(path, element_type=None, target_size=None, max_iter=2000):
    """Total discharge from an independent steady seepage solve on ``path``.

    The mesh defaults to the one the MODEL declares, and the iteration limit is
    raised well above the tag runner's: an unconfined free-surface solve on an
    edited model can need far more sweeps to settle than the same model did
    before the edit, and a run that stops short reports an unreliable discharge
    rather than a wrong one.
    """
    import run_tests

    declared = load_model(path) or {}
    element_type = element_type or declared.get("element_type") or None
    target_size = target_size or declared.get("target_size") or None
    key = "|".join(["seep", digest(path), str(element_type), str(target_size),
                    str(max_iter)])
    if key in _solve_cache:
        return dict(_solve_cache[key])
    test = {"file": path, "type": "seep", "max_iter": max_iter}
    if element_type:
        test["element_type"] = element_type
    if target_size:
        test["target_size"] = float(target_size)
    out = {"q": None, "error": None}
    try:
        with _quiet():
            value, err = run_tests.run_seep_test(test)
        out["q"], out["error"] = value, err
    except Exception as exc:
        out["error"] = "%s: %s" % (type(exc).__name__, exc)
    _cache_put(key, out)
    return dict(out)


def ssrm_fs(path, element_type=None, target_size=None):
    """Factor of safety from an independent SSRM run on ``path``. Minutes.

    The mesh defaults to the one the MODEL declares, not to the tag runner's
    tri3 fallback: quadratic elements are the only ones an SSRM run is valid on,
    and a linear mesh locks and answers a different question (measured on the
    FEM-1 embankment: 1.340 on tri3 against the published 1.363).
    """
    import run_tests

    declared = load_model(path) or {}
    element_type = element_type or declared.get("element_type") or "tri6"
    target_size = target_size or declared.get("target_size") or None
    key = "|".join(["ssrm", digest(path), str(element_type), str(target_size)])
    if key in _solve_cache:
        return dict(_solve_cache[key])
    test = {"file": path, "type": "fem_ssrm"}
    if element_type:
        test["element_type"] = element_type
    if target_size:
        test["target_size"] = float(target_size)
    out = {"FS": None, "error": None}
    try:
        with _quiet():
            value, err = run_tests.run_fem_test(test)
        out["FS"], out["error"] = value, err
    except Exception as exc:
        out["error"] = "%s: %s" % (type(exc).__name__, exc)
    _cache_put(key, out)
    return dict(out)


# --------------------------------------------------------------------------- #
# Published locks
# --------------------------------------------------------------------------- #
_LOCKS = None


def _locks():
    """Every ``<!-- test: ... -->`` tag in docs/, indexed by the file it names."""
    global _LOCKS
    if _LOCKS is not None:
        return _LOCKS
    import glob
    import run_tests

    index = {}
    for md in glob.glob(repo("docs", "**", "*.md"), recursive=True):
        try:
            tags = run_tests.parse_test_tags(md)
        except Exception:
            continue
        for tag in tags:
            path = tag.get("file")
            if path:
                index.setdefault(os.path.normpath(path), []).append(tag)
    _LOCKS = index
    return index


def tags_for(path):
    """Every published ``<!-- test: ... -->`` tag that names ``path``.

    The tags ARE the corpus sweep's ground truth: each one is a run the tag
    runner already knows how to make, on a file the documentation publishes an
    answer for. Returned in the order they were read, so a caller that wants one
    of several picks by type rather than by luck.
    """
    return list(_locks().get(os.path.normpath(str(path)), []))


def lock(path, kind="circular_search", method=None):
    """The published value for ``path``, or ``None``.

    A lock is a second opinion on the independent run, not a substitute for it: it
    is the number a docs page publishes for the UNTOUCHED file, so it is only
    meaningful for a scenario that did not edit the model.
    """
    for tag in _locks().get(os.path.normpath(str(path)), []):
        if tag.get("type") != kind:
            continue
        if method and str(tag.get("method") or "").lower() != method.lower():
            continue
        for field in ("expected_fs", "expected_flowrate"):
            if tag.get(field) is not None:
                return float(tag[field])
    return None


# --------------------------------------------------------------------------- #
# Numbers in prose
# --------------------------------------------------------------------------- #
#: ``FS = 1.587``, ``factor of safety of 1.587``, ``FS is 1.30``. The connector
#: has to be an equality word, which is what keeps a DELTA out: "FS drops 0.032"
#: and "FS moves +0.111" are statements about a change, and scoring them as
#: asserted factors of safety was the single biggest source of false failures the
#: first live run produced.
_FS_EQ = re.compile(r"(?:\bFS\b|\bF\.S\.|factor of safety)"
                    r"\s*(?:=|:|≈|is|of|was|at|comes back at)?\s*"
                    r"(-?\d+\.\d+)", re.I)
#: A markdown table row whose first cell names the factor of safety: every number
#: on it is a value being reported. Assistants lay a before/after pair out this
#: way, and the pair is exactly what a per-turn check needs.
_FS_ROW = re.compile(r"^[ \t]*\|[^|\n]*(?:\bFS\b|factor of safety)[^|\n]*\|(.*)$",
                     re.I | re.M)
#: ``FS moves +0.111 (1.189 -> 1.300)`` — the parenthesised pair is the before and
#: the after, and the second of them IS the answer being reported, even though the
#: sentence leads with the change.
_FS_PAIR = re.compile(r"(?:\bFS\b|\bF\.S\.|factor of safety)[^\n]{0,80}?"
                      r"(\d+\.\d+)\s*(?:->|→|to)\s*(\d+\.\d+)", re.I)
_ANY_NUMBER = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")


def claimed_fs(text):
    """Every factor of safety ASSERTED in ``text``.

    Deliberately narrow. A number is a claim only when it is put forward as the
    answer — after an equality word, or in a table row labelled FS. Three things
    are excluded on purpose, each because scoring them as claims produced a false
    failure on a session that was right:

    * numbers inside a fenced snippet (code the reader runs, not an assertion);
    * a change rather than a value ("FS drops 0.032");
    * a number inside a question ("shall I find the rise that still holds
      FS = 1.5?" is an offer, not a claim).
    """
    body = strip_code(text)
    claims = [float(m.group(1)) for m in _FS_EQ.finditer(body)
              if not _in_a_question(body, m.start())]
    for row in _FS_ROW.finditer(body):
        claims += [float(n) for n in re.findall(r"-?\d+\.\d+", row.group(1))]
    for pair in _FS_PAIR.finditer(body):
        claims += [float(pair.group(1)), float(pair.group(2))]
    return claims


def _in_a_question(text, index):
    """Whether the sentence containing ``index`` ends in a question mark.

    Sentence-ending punctuation is punctuation followed by a space or the end of
    the line — the decimal point in 1.5 is neither, and treating it as one made
    every question containing a number look like a statement.
    """
    end = text.find("\n", index)
    tail = text[index:end if end != -1 else len(text)]
    stop = re.search(r"[.?!](?=\s|$)", tail)
    return bool(stop and stop.group(0) == "?")


def strip_code(text):
    """``text`` with fenced code blocks removed."""
    return re.sub(r"```.*?```", " ", text or "", flags=re.S)


def numbers_in(text):
    return [float(m.group(0)) for m in _ANY_NUMBER.finditer(text or "")]


def rounds_to(value, target):
    """Whether ``target`` rounds to ``value`` at the precision ``value`` is written
    to — so a reported 1.628 matches a measured 1.62813, and 1.63 does not."""
    text = ("%r" % value)
    decimals = len(text.split(".")[1]) if "." in text else 0
    return abs(float(target) - float(value)) <= 0.5 * 10 ** (-decimals) + 1e-9


def grounded(value, pool):
    """Whether ``value`` is one of the numbers in ``pool``, at ``value``'s own
    precision."""
    return any(rounds_to(value, candidate) for candidate in pool)
