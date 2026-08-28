"""Criteria — the questions a finished session is scored against.

Each factory returns a :class:`~tools.assistant_scenarios.core.Criterion`. The rule
every one of them follows: **measure, never believe**. A criterion about a number
re-runs the engine; a criterion about an edit reloads the workbook from disk; a
criterion about a claim compares the claim with what the session's own snippets
printed. Nothing is scored from the assistant's description of what it did.
"""

from __future__ import annotations

import os
import re

from .core import (Criterion, claimed_fs, declared_method, grounded, lock,
                   load_model, numbers_in, repo, rounds_to, seep_flow, solve,
                   ssrm_fs, strip_code)


# --------------------------------------------------------------------------- #
# Round-trip discipline
# --------------------------------------------------------------------------- #
_LATEX = re.compile(r"\$\$|\\frac|\\tan\b|\\sqrt|\\begin\{|\\text\{|\\phi\b|"
                    r"\\alpha\b|\\beta\b|\$[^\s$][^$\n]{0,80}\$")


def no_latex():
    """The dock renders GitHub markdown and no math. LaTeX reaches the reader as
    literal backslashes."""
    def check(ctx):
        hits = sorted({m.group(0)[:24] for m in _LATEX.finditer(ctx.prose)})
        if hits:
            return False, "LaTeX in the answer: %s" % ", ".join(hits[:4])
        return True, "no LaTeX"
    return Criterion("no LaTeX", check, kind="discipline")


_EXPLORE = re.compile(r"\bhelp\s*\(|\bdir\s*\(|\binspect\b|\bpkgutil\b|"
                      r"getsource|signature\s*\(|__doc__")


def no_exploration():
    """The helpers are documented in the brief. A turn that opens by reading the
    engine's source has spent a completion on something it was already told."""
    def check(ctx):
        hits = sorted({m.group(0).strip() for m in _EXPLORE.finditer(ctx.code)})
        if hits:
            return False, "exploration calls in snippets: %s" % ", ".join(hits[:4])
        return True, "no help()/dir()/inspect"
    return Criterion("no exploration calls", check, kind="discipline")


def completions_at_most(n):
    """Completions per turn — the unit the user is billed in."""
    def check(ctx):
        worst = max([t.completions for t in ctx.turns] or [0])
        counts = ", ".join(str(t.completions) for t in ctx.turns)
        if worst > n:
            return False, "completions per turn %s (cap %d)" % (counts, n)
        return True, "completions per turn %s (cap %d)" % (counts, n)
    return Criterion("completions <= %d/turn" % n, check, kind="discipline")


def tokens_at_most(n):
    def check(ctx):
        total = int(ctx.usage.get("input") or 0) + int(ctx.usage.get("output") or 0)
        ok = total <= n
        return ok, "%s tokens (budget %s)" % (format(total, ","), format(n, ","))
    return Criterion("tokens <= %s" % format(n, ","), check, kind="discipline")


def finished():
    """No turn ended on an error or a timeout."""
    def check(ctx):
        for i, turn in enumerate(ctx.turns, start=1):
            if turn.error:
                return False, "turn %d: %s" % (i, turn.error)
        if not ctx.turns:
            return False, "no turns recorded"
        return True, "%d turn(s), no errors" % len(ctx.turns)
    return Criterion("session completed", check, kind="discipline")


# --------------------------------------------------------------------------- #
# Numbers
# --------------------------------------------------------------------------- #
def fs_matches(tol=0.005, method=None, search=True, on="saved", num_slices=40,
               rapid=False):
    """The factor of safety the answer ends on matches an independent re-run.

    ``on='saved'`` re-solves the workbook the session left behind (the right
    ground truth after an edit); ``on='start'`` re-solves the file it opened.
    """
    def check(ctx):
        path = ctx.saved() if on == "saved" else ctx.start_model
        if on == "saved" and not path:
            path = ctx.start_model
        if not path:
            return False, "no workbook to re-solve"
        claims = claimed_fs(ctx.final_prose) or claimed_fs(ctx.prose)
        if not claims:
            return False, "no factor of safety stated in the answer"
        run = solve(path, method=method, search=search, num_slices=num_slices,
                    rapid=rapid)
        if run.get("FS") is None:
            return False, "independent re-run failed: %s" % run.get("error")
        best = min(claims, key=lambda v: abs(v - run["FS"]))
        ok = abs(best - run["FS"]) <= tol
        return ok, ("stated FS %s vs independent %s (%s)"
                    % (best, round(run["FS"], 4), run["method"]))
    return Criterion("FS matches an independent re-run", check, kind="truth")


def claims_grounded(allow=0):
    """Every factor of safety the prose asserts was printed by a snippet.

    This is the check the W-1 audit's failures all sit upstream of: the arithmetic
    was never wrong, but a number can still be asserted that no run produced.
    """
    def check(ctx):
        pool = numbers_in(ctx.output)
        loose = [v for v in claimed_fs(ctx.prose) if not grounded(v, pool)]
        if len(loose) > allow:
            return False, ("%d stated FS value(s) never printed by a snippet: %s"
                           % (len(loose), ", ".join(str(v) for v in loose[:5])))
        return True, "every stated FS appears in the printed output"
    return Criterion("stated numbers were computed", check, kind="truth")


def lock_holds(kind="circular_search", method=None, tol=0.01):
    """The starting file still reproduces its published value.

    A fixture guard reported as a criterion: if this fails, the scenario's other
    numeric criteria are measuring against a model that moved.
    """
    def check(ctx):
        path = ctx.start_model
        published = lock(ctx.lock_model, kind=kind, method=method)
        if published is None:
            return True, "no published lock for this file"
        if kind == "seep":
            run = seep_flow(path)
            got = run.get("q")
        elif kind == "fem_ssrm":
            run = ssrm_fs(path)
            got = run.get("FS")
        else:
            run = solve(path, method=method, search=True)
            got = run.get("FS")
        if got is None:
            return False, "independent run failed: %s" % run.get("error")
        ok = abs(got - published) <= tol * max(1.0, abs(published))
        return ok, "published %s vs independent %s" % (published, round(got, 4))
    return Criterion("published lock reproduces", check, kind="truth")


def flow_matches(tol=0.02):
    """The discharge the answer reports matches an independent seepage solve."""
    def check(ctx):
        path = ctx.saved() or ctx.start_model
        run = seep_flow(path)
        if run.get("q") is None:
            return False, "independent seepage run failed: %s" % run.get("error")
        pool = [abs(v) for v in numbers_in(strip_code(ctx.prose))]
        want = abs(run["q"])
        near = [v for v in pool if want and abs(v - want) <= tol * want]
        if not near:
            return False, ("independent q = %.6g; the answer states no number "
                           "within %d%%" % (run["q"], int(tol * 100)))
        return True, "stated %.6g vs independent %.6g" % (near[0], run["q"])
    return Criterion("discharge matches an independent solve", check, kind="truth")


def ssrm_matches(tol=0.03):
    """The SSRM factor of safety the answer reports matches the published one.

    Where the model carries a published ``fem_ssrm`` lock, that IS the independent
    run — made by the suite on every full pass, at the mesh the docs page names —
    and it is used, because an SSRM solve costs minutes and re-running one per
    scored scenario would put the cost of scoring above the cost of the turns.
    A model with no lock is re-run here instead. That the assistant really solved
    rather than recalled is a separate criterion (``claims_grounded``): the number
    it states has to appear in a snippet's own printed output.
    """
    def check(ctx):
        claims = claimed_fs(ctx.prose)
        if not claims:
            return False, "no factor of safety stated in the answer"
        published = lock(ctx.lock_model, kind="fem_ssrm")
        if published is not None:
            best = min(claims, key=lambda v: abs(v - published))
            ok = abs(best - published) <= tol * max(1.0, abs(published))
            return ok, "stated %s vs published SSRM %s" % (best, published)
        run = ssrm_fs(ctx.saved() or ctx.start_model)
        if run.get("FS") is None:
            return False, "independent SSRM run failed: %s" % run.get("error")
        best = min(claims, key=lambda v: abs(v - run["FS"]))
        ok = abs(best - run["FS"]) <= tol * max(1.0, abs(run["FS"]))
        return ok, "stated %s vs independent SSRM %s" % (best, round(run["FS"], 4))
    return Criterion("SSRM FS matches the published run", check, kind="truth")


def method_is(expected=None):
    """The method that ran is the model's own (or the one the user named), and the
    answer says which method it was."""
    def check(ctx):
        sd = load_model(ctx.start_model) if ctx.start_model else None
        want = (expected or (declared_method(sd) if sd else "spencer")).lower()
        named = re.search(r"\b(oms|bishop|janbu|corps|lowe|spencer|"
                          r"morgenstern[- ]price|mprice)\b", ctx.prose, re.I)
        if not named:
            return False, "the answer never names the method it ran"
        said = named.group(1).lower().replace("morgenstern-price", "mprice") \
                                     .replace("morgenstern price", "mprice")
        forced = re.search(r"method\s*=\s*['\"](\w+)['\"]", ctx.code)
        ran = (forced.group(1).lower() if forced else want)
        ok = (said == want) and (ran == want)
        return ok, "model declares %s; answer says %s; snippet ran %s" % (
            want, said, ran)
    return Criterion("ran the model's method", check, kind="truth")


# --------------------------------------------------------------------------- #
# The model on disk
# --------------------------------------------------------------------------- #
def edited(label, predicate, turn=None):
    """An edit, verified on the RELOADED workbook.

    ``predicate(before, after)`` returns ``(ok, reason)``; both are ``slope_data``
    dicts read back from disk, so an edit that the canvas silently rebuilt away is
    seen here as not made.

    In a cumulative conversation ``before`` is the PREVIOUS turn's workbook, not
    the file the session opened: turn 3 of "change the face, add a load, extend
    the reinforcement" has to be measured against turn 2, or the face change of
    turn 1 — which moved every reinforcement line with it — is charged to turn 3.
    """
    def check(ctx):
        path = ctx.saved(turn)
        if not path:
            return False, "the session saved no workbook (nothing changed)"
        after = load_model(path)
        if after is None:
            return False, "the saved workbook does not load"
        previous = ctx.saved(turn - 1) if turn and turn > 1 else None
        before = load_model(previous) if previous else ctx.before()
        return predicate(before, after)
    return Criterion(label, check, kind="edit")


def model_unchanged():
    """The session left the model alone — the right answer to a question, and to a
    request it should have asked about instead of guessing."""
    def check(ctx):
        if ctx.model_changed():
            return False, "the model was edited and saved"
        return True, "the model was not changed"
    return Criterion("model left unchanged", check, kind="edit")


def geometry_source_kept():
    """A profile-line model is still profile-line native, a polygon model still
    polygon native — the edit went into the source, not the derived copy."""
    def check(ctx):
        before, after = ctx.before(), ctx.after()
        if after is None:
            return False, "the session saved no workbook"
        was = bool((before or {}).get("profile_lines"))
        now = bool(after.get("profile_lines"))
        if was and not now:
            return False, "the model was profile-line native and now has none"
        if not was and now:
            return False, "polygons were the source; profile lines were added"
        return True, "geometry source unchanged (%s native)" % (
            "profile-line" if was else "polygon")
    return Criterion("edited the geometry source", check, kind="edit")


def solves(method=None):
    """The model the session leaves behind can be solved at all."""
    def check(ctx):
        path = ctx.saved() or ctx.start_model
        run = solve(path, method=method, search=True)
        if run.get("FS") is None:
            return False, "the saved model does not solve: %s" % run.get("error")
        return True, "solves at FS = %s (%s)" % (round(run["FS"], 4), run["method"])
    return Criterion("the model solves", check, kind="edit")


def model_matches(label, predicate):
    """A property of the saved model (or the starting one, when nothing changed)."""
    def check(ctx):
        sd = ctx.after() or load_model(ctx.start_model)
        if sd is None:
            return False, "no model to read"
        return predicate(sd)
    return Criterion(label, check, kind="edit")


# --------------------------------------------------------------------------- #
# Which tools it reached for
# --------------------------------------------------------------------------- #
def helper_used(*names):
    """Every named helper was actually called in a snippet."""
    def check(ctx):
        missing = [n for n in names if not re.search(r"\b%s\s*\(" % re.escape(n),
                                                     ctx.code)]
        if missing:
            return False, "never called: %s" % ", ".join(missing)
        return True, "called %s" % ", ".join(names)
    return Criterion("uses %s" % "/".join(names), check, kind="behavior")


def helper_not_used(*names):
    def check(ctx):
        used = [n for n in names if re.search(r"\b%s\s*\(" % re.escape(n), ctx.code)]
        if used:
            return False, "called %s" % ", ".join(used)
        return True, "did not call %s" % ", ".join(names)
    return Criterion("avoids %s" % "/".join(names), check, kind="behavior")


def snippets_at_most(n):
    def check(ctx):
        total = sum(len(t.codes) for t in ctx.turns)
        return total <= n, "%d snippet(s) (cap %d)" % (total, n)
    return Criterion("snippets <= %d" % n, check, kind="discipline")


def asks_rather_than_guesses():
    """An ambiguous request is answered with a question and no edit."""
    def check(ctx):
        if ctx.model_changed():
            return False, "it edited the model instead of asking"
        if "?" not in ctx.final_prose:
            return False, "no question asked"
        return True, "asked instead of guessing"
    return Criterion("asks when the request is ambiguous", check, kind="behavior")


def warns_about(pattern, label=None):
    """The answer raises a named concern (a units mismatch, an overlap, a risk)."""
    rx = pattern if hasattr(pattern, "search") else re.compile(pattern, re.I)
    def check(ctx):
        found = rx.search(ctx.prose)
        return bool(found), ("said %r" % found.group(0)[:48]) if found \
            else "never raised it"
    return Criterion(label or "raises the concern", check, kind="behavior")


def reports_warnings():
    """Every ``WARNING:`` a snippet came back with is passed on to the user."""
    def check(ctx):
        lines = [ln.strip() for ln in ctx.output.splitlines()
                 if ln.strip().startswith("WARNING:")]
        if not lines:
            return True, "no warnings were raised"
        prose = ctx.prose.lower()
        silent = []
        for line in lines:
            words = [w for w in re.findall(r"[a-z_]{5,}", line.lower())
                     if w not in ("warning", "instead", "because")]
            if not any(w in prose for w in words[:6]):
                silent.append(line[:60])
        if silent:
            return False, "warning not passed on: %s" % silent[0]
        return True, "%d warning(s), all reported" % len(lines)
    return Criterion("reports the kernel's warnings", check, kind="behavior")


# --------------------------------------------------------------------------- #
# Citations
# --------------------------------------------------------------------------- #
_DOC_URL = re.compile(r"xslope\.readthedocs\.io/en/latest/([A-Za-z0-9_/\-]*)")


def cites_docs(at_least=1):
    """The pages named exist in this repository's ``docs/`` tree."""
    def check(ctx):
        pages = [p.strip("/") for p in _DOC_URL.findall(ctx.prose)]
        if len(pages) < at_least:
            return False, "cited %d page(s), wanted %d" % (len(pages), at_least)
        bad = []
        for page in pages:
            candidates = [repo("docs", page + ".md"),
                          repo("docs", page, "index.md"),
                          repo("docs", page.rstrip("/") + "/index.md")]
            if page in ("", "/"):
                continue
            if not any(os.path.exists(c) for c in candidates):
                bad.append(page)
        if bad:
            return False, "page(s) that do not exist: %s" % ", ".join(sorted(set(bad)))
        return True, "%d real page(s): %s" % (len(pages), ", ".join(pages[:4]))
    return Criterion("cites real documentation pages", check, kind="behavior")


def cites_corpus():
    """A worked example was looked up, not remembered: ``corpus_index`` was called
    and something it returned is named in the answer."""
    def check(ctx):
        if not re.search(r"\bcorpus_index\s*\(", ctx.code):
            return False, "corpus_index was never called"
        rows = ctx.output
        titles = re.findall(r"'title':\s*'([^']{6,60})'", rows) \
            + re.findall(r"\b(VP\d+|RS2-\d+|GS2-\d+|LEM-\d+|SEEP-\d+|FEM-\d+)\b",
                         rows)
        named = [t for t in titles if t and t.lower() in ctx.prose.lower()]
        if not named:
            return False, "corpus_index ran but nothing it returned is cited"
        return True, "cites %s" % named[0]
    return Criterion("cites a corpus example it looked up", check, kind="behavior")


# --------------------------------------------------------------------------- #
# Diagnosis
# --------------------------------------------------------------------------- #
def faults_named():
    """Every planted fault is named in the answer."""
    def check(ctx):
        prose = ctx.prose
        missed = [f.id for f in ctx.scenario.faults if not f.names.search(prose)]
        if missed:
            return False, "missed: %s" % ", ".join(missed)
        return True, "named all %d planted fault(s)" % len(ctx.scenario.faults)
    return Criterion("names every planted fault", check, kind="truth")


def no_false_accusation(*patterns):
    """Nothing correct in the model is reported as the bug."""
    compiled = [(p if hasattr(p, "search") else re.compile(p, re.I))
                for p in patterns]
    def check(ctx):
        for rx in compiled:
            found = rx.search(ctx.prose)
            if found:
                return False, "accused a sound input: %r" % found.group(0)[:60]
        return True, "no sound input was accused"
    return Criterion("accuses nothing sound", check, kind="truth")


_MUTATES = re.compile(r"(slope_data\s*\[[^\]]+\]\s*(\[[^\]]+\]\s*)*=|"
                      r"\bmat\w*\s*\[[^\]]+\]\s*=|\['(phi|c|gamma|Normal|"
                      r"max_depth|Tmax)'\]\s*=)")
_SOLVES = re.compile(r"\brun_lem\s*\(|\bcircular_search\s*\(|\bsolve_\w+\s*\(")


def varied_inputs(at_least=2):
    """A diagnosis is measured, not narrated: at least N inputs were changed and
    re-solved on.

    Counted in VARIATIONS rather than snippets. One snippet that sets four
    different fields and re-solves after each is the efficient way to do this, and
    counting snippets scored it as one probe — which failed a session that had in
    fact tested more than the criterion asked for.
    """
    def check(ctx):
        probes = [code for t in ctx.turns for code in t.codes
                  if _MUTATES.search(code) and _SOLVES.search(code)]
        changes = sum(len(_MUTATES.findall(code)) for code in probes)
        ok = changes >= at_least
        return ok, "%d input change(s) re-solved on, in %d snippet(s) (wanted %d)" % (
            changes, len(probes), at_least)
    return Criterion("varies inputs and measures", check, kind="truth")


def states_the_measurement():
    """The answer reports a measured change, not only a reading of the inputs."""
    def check(ctx):
        prose = strip_code(ctx.prose)
        rx = re.compile(r"(\d+\.\d+)\s*(?:->|→|to|→)\s*(\d+\.\d+)", re.I)
        if rx.search(prose):
            return True, "reports a before/after pair"
        if len(claimed_fs(prose)) >= 2:
            return True, "reports more than one measured factor of safety"
        # A table with a factor-of-safety column and two or more rows is the same
        # measurement laid out the other way: what was changed, and what it gave.
        header = re.search(r"^\|.*(?:\bFS\b|factor of safety).*\|\s*$",
                           prose, re.I | re.M)
        if header:
            rows = [ln for ln in prose[header.end():].splitlines()
                    if ln.strip().startswith("|")
                    and re.search(r"\d+\.\d+", ln)]
            if len(rows) >= 2:
                return True, "a factor-of-safety table with %d measured rows" % len(rows)
        return False, "no before/after measurement in the answer"
    return Criterion("reports what the change was worth", check, kind="truth")


# --------------------------------------------------------------------------- #
# Deliverables
# --------------------------------------------------------------------------- #
def wrote_file(suffix, contains=()):
    """A file of the given kind was produced, and (for a .docx) carries the strings
    the answer says it carries."""
    def check(ctx):
        # A file offered in the chat, OR a path a snippet printed. Only files
        # written into the kernel's own output folder are offered in the chat, and
        # `generate_report` writes beside the project instead — so a report that
        # exists and is correct is invisible to the first list.
        made = [f for t in ctx.turns for f in t.files if f.lower().endswith(suffix)]
        made += re.findall(r"[\w./\\ -]+%s" % re.escape(suffix), ctx.output)
        made = [m.strip() for m in made if m.strip()]
        if not made:
            return False, "no %s was written" % suffix
        path = None
        for candidate in reversed(made):
            if os.path.isabs(candidate) and os.path.exists(candidate):
                path = candidate
                break
        if path is None:
            for root, _dirs, names in os.walk(ctx.outdir):
                for name in names:
                    if name == os.path.basename(made[-1]):
                        path = os.path.join(root, name)
        if path is None:
            return False, "%s is named but no such file exists" % made[-1]
        wanted = [c(ctx) if callable(c) else c for c in contains]
        if wanted:
            text = _docx_text(path) if suffix == ".docx" else _read(path)
            missing = [c for c in wanted if str(c).lower() not in text.lower()]
            if missing:
                return False, "%s does not contain %s" % (
                    os.path.basename(path), str(missing[0])[:40])
        return True, "wrote %s%s" % (os.path.basename(path),
                                     " carrying %d checked string(s)" % len(wanted)
                                     if wanted else "")
    return Criterion("produced the %s" % suffix, check, kind="behavior")


def _read(path):
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            return fh.read()
    except Exception:
        return ""


def _docx_text(path):
    """The document's text, straight out of ``word/document.xml`` — no Word, no
    converter, and no dependency on how the paragraph was styled."""
    import zipfile
    try:
        with zipfile.ZipFile(path) as zf:
            xml = zf.read("word/document.xml").decode("utf-8", "replace")
    except Exception:
        return ""
    return re.sub(r"<[^>]+>", " ", xml)


# --------------------------------------------------------------------------- #
# More of the same, for the scenarios that need them
# --------------------------------------------------------------------------- #
def helper_any(*names):
    """At least one of the named helpers was called — for a job two of them do."""
    def check(ctx):
        used = [n for n in names if re.search(r"\b%s\s*\(" % re.escape(n), ctx.code)]
        if not used:
            return False, "none of %s was called" % ", ".join(names)
        return True, "called %s" % ", ".join(used)
    return Criterion("uses one of %s" % "/".join(names), check, kind="behavior")


def per_turn_fs(tol=0.005, method=None):
    """Each turn's stated factor of safety matches a re-run of THAT TURN's workbook.

    A conversation that edits the model three times has three answers, and the end
    state cannot show whether the second one was right — so each turn is re-solved
    against the file it left behind.
    """
    def check(ctx):
        rows, bad = [], []
        for i, turn in enumerate(ctx.turns, start=1):
            path = turn.workbook or (ctx.workbook if i == len(ctx.turns) else None)
            claims = claimed_fs(turn.prose)
            if not path or not claims:
                continue
            run = solve(path, method=method, search=True)
            if run.get("FS") is None:
                bad.append("turn %d: re-run failed (%s)" % (i, run.get("error")))
                continue
            best = min(claims, key=lambda v: abs(v - run["FS"]))
            rows.append("t%d %s/%s" % (i, best, round(run["FS"], 4)))
            if abs(best - run["FS"]) > tol:
                bad.append("turn %d: stated %s, re-run %s"
                           % (i, best, round(run["FS"], 4)))
        if not rows and not bad:
            return False, "no turn stated a factor of safety on a saved workbook"
        if bad:
            return False, "; ".join(bad)
        return True, "stated/re-run " + ", ".join(rows)
    return Criterion("every turn's FS matches its own workbook", check, kind="truth")


def fs_close_to(published, rel=0.10, label=None):
    """The model the session built lands near a published answer for the same
    problem — a fidelity check on what was MODELLED, not on the arithmetic."""
    def check(ctx):
        path = ctx.saved() or ctx.start_model
        run = solve(path, search=True)
        if run.get("FS") is None:
            return False, "the model does not solve: %s" % run.get("error")
        off = abs(run["FS"] - published) / max(1e-9, abs(published))
        return off <= rel, "built model solves at %s; published %s (%.1f%% off)" % (
            round(run["FS"], 4), published, 100 * off)
    return Criterion(label or "within %d%% of the published answer" % int(rel * 100),
                     check, kind="truth")


_CAUSAL = re.compile(r"\b(because|the reason|is capped by|is governed by|"
                     r"is controlled by|is driven by|which is why|due to)\b", re.I)


def mechanism_claims_tested():
    """A cause offered for a number was tested by varying something.

    The one failure class the W-1 audit found that no interface change reaches: a
    confident account of a mechanism that no measurement supports. A snippet that
    changed an input and re-solved is what makes the account falsifiable.
    """
    def check(ctx):
        claims = _CAUSAL.findall(strip_code(ctx.prose))
        if not claims:
            return True, "no causal claim made"
        probes = [c for t in ctx.turns for c in t.codes
                  if _MUTATES.search(c) and _SOLVES.search(c)]
        if not probes:
            return False, ("%d causal claim(s) (%r) with no snippet that varied "
                           "an input" % (len(claims), claims[0]))
        return True, "%d causal claim(s), %d probe(s) behind them" % (
            len(claims), len(probes))
    return Criterion("mechanisms offered were tested", check, kind="truth")


def stated_values_grounded(pattern, label):
    """Every number the prose attaches to ``pattern`` was printed by a snippet.

    ``claims_grounded`` does this for the factor of safety; this is the same rule
    for a reliability index, a probability of failure, a discharge — anything the
    brief says must come from a helper rather than from arithmetic in the answer.
    """
    rx = re.compile(pattern + r"[^0-9\n]{0,24}(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)",
                    re.I)
    # "the standard error on P_f is about 0.84%" states an uncertainty, not a
    # probability of failure, and grading it as one failed a session that was
    # right about every number it did report.
    about_error = re.compile(r"(error|interval|uncertaint|noise|tolerance|"
                             r"spread|±|\+/-)", re.I)
    def check(ctx):
        pool = numbers_in(ctx.output)
        stated = [float(m.group(1)) for m in rx.finditer(strip_code(ctx.prose))
                  if not about_error.search(
                      strip_code(ctx.prose)[max(0, m.start() - 60):m.start()])]
        if not stated:
            return False, "no %s stated" % label
        loose = [v for v in stated if not grounded(v, pool)]
        if loose:
            return False, "%s not printed by any snippet: %s" % (
                label, ", ".join(str(v) for v in loose[:4]))
        return True, "%d %s value(s), all computed" % (len(stated), label)
    return Criterion("%s comes from the helper" % label, check, kind="truth")


def crossing_verified(pattern, apply, target_fs, tol=0.03, label=None):
    """The value the answer reports really does put the slope at the target.

    The answer names a cohesion, a friction angle, a load. This takes that number,
    writes it into a fresh copy of the starting model, solves, and checks the
    factor of safety lands on the target — so a sweep that was read off the wrong
    axis, or extrapolated past its bracket, fails here.
    """
    rx = re.compile(pattern, re.I)
    def check(ctx):
        import os
        from .faults import plant
        from .core import Fault

        found = rx.search(strip_code(ctx.prose))
        if not found:
            return False, "the answer names no value"
        # A thousands separator is how an engineer writes 1,050 psf; without this
        # the number parsed out of the answer would be 50.
        value = float(found.group(1).replace(",", ""))
        scratch = os.path.join(ctx.outdir, "verify")
        fault = Fault("verify", "apply the reported value",
                      lambda sd, v=value: apply(sd, v), r"(?!)")
        path = plant(ctx.start_model, [fault], scratch)
        run = solve(path, search=True)
        if run.get("FS") is None:
            return False, "the verifying run failed: %s" % run.get("error")
        off = abs(run["FS"] - target_fs)
        return off <= tol, "reported %g gives FS %s (target %s)" % (
            value, round(run["FS"], 4), target_fs)
    return Criterion(label or "the reported value hits the target", check,
                     kind="truth")
