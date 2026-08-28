"""The corpus sweep — every shipped workbook, one turn each.

The scored scenarios in :mod:`~tools.assistant_scenarios.definitions` are thirty
conversations chosen to cover the FAMILIES of thing a person asks for. They say
nothing about coverage of the INPUTS: a model with piles, a model with a tension
crack and a model with a Hoek-Brown material all reach the assistant through the
same helper, and a defect that only shows on one of them is invisible to a suite
that never opens one.

This sweep is the other axis. It enumerates every workbook the repository ships,
opens each in its own session, asks one cheap question — run it as the file
declares it, and say what came back — and scores the answer against the same tag
runner ``run_tests.py`` drives. Each file is labelled with the input columns it
actually exercises, read off the loaded model rather than off its name, so the
failures group by input class and a pattern ("every model with a tension crack")
is visible where thirty hand-written scenarios would show one anecdote.

Six things are asked of every session:

* the number it reports matches an independent run at the tag runner's tolerance
  — and where it reports NO number, no independent run produces one either;
* the run used the method the file declares, not one of the assistant's choosing;
* no snippet opened by reading the engine's source;
* the turn cost at most four completions;
* no snippet was broken code — a solver reporting that it has no solution is a
  finding about the model and is not counted against the session;
* and — the point of the sweep — any workbook the session SAVED still loads
  through ``load_slope_data`` and still reproduces every published value its
  source file carries.

The last one is why the sweep is worth its money. A session that answers a
question and quietly rewrites the model it was asked about is the failure mode a
transcript cannot show, and 390 workbooks is a wide enough net to catch a save
path that only misbehaves on one input column.
"""

from __future__ import annotations

import glob
import os
import re

from .core import (Criterion, Scenario, claimed_fs, declared_method, digest,
                   load_model, numbers_in, repo, strip_code, tags_for)

# --------------------------------------------------------------------------- #
# Which files
# --------------------------------------------------------------------------- #
#: Directories swept flat, in the order the scorecard lists them.
CORPUS_DIRS = ("docs/inputs/slope", "docs/tutorials/files", "docs/lem/files",
               "docs/seep/files", "docs/fem/files", "docs/parametric/files")
#: Directories swept recursively — the verification corpus is one tree with a
#: directory per vendor.
CORPUS_TREES = ("docs/verification/files",)

#: The one prompt. Every session gets exactly this, so the only thing that varies
#: between sessions is the model — which is what makes a per-class failure rate
#: mean something.
PROMPT = ("Run this model's analysis as the file declares it (its method, its "
          "surface, its seepage/FEM settings if it is a seepage or FEM model) "
          "and report the result.")

#: Name suffixes that mark a companion copy of another workbook. One is skipped
#: only when it is byte-identical to the sibling it is named after: a `_start`
#: file that really is the tutorial's starting point differs from the finished
#: model and is a model in its own right, and a `_KEY` that has drifted from its
#: problem file is exactly the kind of thing a sweep should be opening.
_COMPANION_SUFFIXES = ("_start", "_KEY")


def discover():
    """Every shipped workbook, deduplicated, in a stable order."""
    paths = []
    for folder in CORPUS_DIRS:
        paths += sorted(glob.glob(repo(folder, "*.xlsx")))
    for tree in CORPUS_TREES:
        paths += sorted(glob.glob(repo(tree, "**", "*.xlsx"), recursive=True))
    return [p for p in paths
            if not os.path.basename(p).startswith("~$")
            and not duplicate_of_sibling(p)]


def duplicate_of_sibling(path):
    """Whether ``path`` is a ``_start``/``_KEY`` copy identical to its sibling."""
    stem, ext = os.path.splitext(path)
    for suffix in _COMPANION_SUFFIXES:
        if stem.endswith(suffix):
            sibling = stem[:-len(suffix)] + ext
            if os.path.exists(sibling) and digest(sibling) == digest(path):
                return True
    return False


# --------------------------------------------------------------------------- #
# What a file exercises
# --------------------------------------------------------------------------- #
#: Input classes in the order a case's PRIMARY class is chosen — rarest first, so
#: a model that carries piles is filed under piles rather than under the polygon
#: geometry every model has. The sweep is ordered by round-robin over these
#: buckets, which is what makes a run stopped at the budget still say something
#: about every column rather than everything about the first directory.
CLASS_ORDER = (
    "transient seepage", "rapid drawdown", "piles", "reinforcement",
    "Hoek-Brown", "power-curve strength", "c/p ratio", "elastic (FEM)",
    "seismic", "tension crack", "line loads", "staged (set 2)",
    "seepage BCs", "seepage pore pressure", "ru pore pressure",
    "piezometric line", "non-circular", "distributed loads",
    "saturated unit weight", "metric units", "declared mesh", "max depth",
    "profile lines", "polygons",
)


def input_classes(sd):
    """The input columns ``sd`` actually exercises, read off the loaded model."""
    if sd is None:
        return ["does not load"]
    mats = sd.get("materials") or []
    u_kinds = {str(m.get("u") or "none").strip().lower() for m in mats}
    options = {str(m.get("option") or "").strip().lower() for m in mats}
    found = set()

    def mark(flag, label):
        if flag:
            found.add(label)

    mark(sd.get("tseep"), "transient seepage")
    mark(any(float(m.get("psi") or 0.0) for m in mats), "rapid drawdown")
    mark(sd.get("pile_lines"), "piles")
    mark(sd.get("reinforcement_lines"), "reinforcement")
    mark("hb" in options, "Hoek-Brown")
    mark("pow" in options, "power-curve strength")
    mark("cp" in options, "c/p ratio")
    mark("elastic" in options, "elastic (FEM)")
    mark(float(sd.get("k_seismic") or 0.0), "seismic")
    mark(float(sd.get("tcrack_depth") or 0.0), "tension crack")
    mark(sd.get("line_loads"), "line loads")
    mark(sd.get("dloads2") or sd.get("has_seepage_bc2")
         or sd.get("piezo_line2"), "staged (set 2)")
    mark((sd.get("seepage_bc") or {}).get("specified_heads")
         or (sd.get("seepage_bc") or {}).get("specified_fluxes"), "seepage BCs")
    mark("seep" in u_kinds, "seepage pore pressure")
    mark("ru" in u_kinds, "ru pore pressure")
    mark(sd.get("piezo_line"), "piezometric line")
    mark(sd.get("non_circ"), "non-circular")
    mark(sd.get("dloads"), "distributed loads")
    mark(any(m.get("gamma_sat") for m in mats), "saturated unit weight")
    mark(str(sd.get("unit_system") or "").lower() == "metric", "metric units")
    mark(sd.get("element_type") or sd.get("target_size"), "declared mesh")
    mark(float(sd.get("max_depth") or 0.0), "max depth")
    mark(sd.get("profile_lines"), "profile lines")
    mark(not sd.get("profile_lines") and sd.get("polygons"), "polygons")
    ordered = [c for c in CLASS_ORDER if c in found]
    return ordered or ["plain"]


def primary_class(classes):
    """The one bucket a case is ordered under."""
    for label in CLASS_ORDER:
        if label in classes:
            return label
    return classes[0] if classes else "plain"


# --------------------------------------------------------------------------- #
# One file's ground truth
# --------------------------------------------------------------------------- #
#: Tag types the sweep can re-run as an independent answer, and which engine each
#: one belongs to. Everything else a docs page tags a file with (a roundtrip, a
#: mesh count, an importer conformance test) says nothing about "run this model's
#: analysis" and is left out of the ground truth rather than scored against it.
_LEM_TAGS = ("circular_search", "single_circle", "noncircular_search",
             "single_noncirc")
_SEEP_TAGS = ("seep",)
_FEM_TAGS = ("fem_ssrm",)


class Case:
    """One workbook, with everything the scorers need to judge a session on it."""

    def __init__(self, path):
        self.path = path
        self.name = _case_name(path)
        self.model = load_model(path)
        self.classes = input_classes(self.model)
        self.primary = primary_class(self.classes)
        self.loads = self.model is not None
        self.method = declared_method(self.model) if self.model else None
        self.declares_method = bool((self.model or {}).get("lem_method"))
        self.tags = [t for t in tags_for(path)
                     if t.get("type") in _LEM_TAGS + _SEEP_TAGS + _FEM_TAGS]
        self.kind = _kind(self.model, self.tags)
        # A material reading pore pressure off a seepage solution needs that
        # solve made before any stability run — the flow first, then the slices.
        # Unless the field already travels with the workbook: a `_seep.csv`
        # sidecar IS the solved field, the loader attaches it as `seep_u`, and a
        # run that solved the flow again would answer a different question from
        # the one the file's published value was made on.
        self.stages_seep = (
            any(str(m.get("u") or "").strip().lower() == "seep"
                for m in (self.model or {}).get("materials") or [])
            and (self.model or {}).get("seep_u") is None)

    @property
    def timeout_s(self):
        return {"fem": 2400, "seep": 1200}.get(self.kind, 600)

    def __repr__(self):
        return "<Case %s (%s, %s)>" % (self.name, self.kind, self.primary)


def _case_name(path):
    """A scenario name from the path — unique across the whole corpus."""
    rel = os.path.relpath(path, repo("docs"))
    stem = os.path.splitext(rel)[0].replace(os.sep, "_")
    return re.sub(r"[^A-Za-z0-9_]+", "_", stem)


def _engine(tag_type):
    if tag_type in _LEM_TAGS:
        return "lem"
    if tag_type in _SEEP_TAGS:
        return "seep"
    return "fem"


def _kind(sd, tags):
    """Which engine "the file's analysis" means: ``lem``, ``seep`` or ``fem``.

    The MODEL answers this, not the tags. Much of the verification corpus carries
    an ``fem_ssrm`` value beside its stability ones — the same slope answered two
    ways — and reading the kind off whichever tag came first filed the LEM-8
    tutorial model as a FEM case, which would have scored a correct Spencer
    answer against a published SSRM one. A model that defines a failure surface
    is a model whose analysis is a stability analysis; one that defines specified
    heads and no surface is a seepage model; what is left is FEM.
    """
    types = {t.get("type") for t in tags}
    if sd is None:
        return _engine(next(iter(types), "circular_search"))
    if sd.get("circles") or sd.get("non_circ"):
        return "lem"
    bc = sd.get("seepage_bc") or {}
    if bc.get("specified_heads") or bc.get("specified_fluxes"):
        return "seep"
    if types & set(_FEM_TAGS) or sd.get("element_type") or sd.get("ssr_zones"):
        return "fem"
    return "lem"


def truth_candidates(case, ran=None):
    """Independent answers a correct session could legitimately report, in the
    order they are tried.

    Every entry is ``(label, callable)`` returning ``(value, error, tolerance)``
    — the tolerance travels with the answer because it belongs to the run that
    produced it, and a seepage discharge and a factor of safety are not compared
    to the same closeness. They are evaluated lazily and in order, so a session
    that agrees with the first one costs one solve; only a session that agrees
    with none pays for them all.

    Why more than one for LEM: the prompt asks for the analysis "as the file
    declares it", and the file declares a method and a surface family but not
    whether to search. The docs page that publishes an answer for the file does
    say — that tag is tried first — and where the file carries no tag, both
    readings are offered, because scoring one of two defensible runs as wrong
    would measure the ambiguity rather than the assistant. The same reasoning
    covers the method: where the declared method has no solution on the declared
    surface, another method on that surface is what an engineer reports, and the
    candidates follow.

    ``ran`` is the engine the session actually drove, read off its snippets. A
    model can define a failure surface and still be a page's finite-element
    example, and comparing an SSRM answer with a limit-equilibrium one measures
    nothing; the engine the session used is tried first, the one the file reads
    as next.
    """
    groups = {"lem": _lem_candidates(case),
              "seep": _seep_candidates(case),
              "fem": _fem_candidates(case)}
    order = []
    for key in (ran, case.kind, "lem", "seep", "fem"):
        if key in groups and key not in order:
            order.append(key)
    out = []
    for key in order:
        out += groups[key]
    return out


def _tagged(case, types):
    """The file's published tags of the given types, declared method first."""
    def key(tag):
        same = 0 if str(tag.get("method") or "").lower() == (case.method or "") \
            else 1
        return (same, types.index(tag.get("type")))
    return sorted((t for t in case.tags if t.get("type") in types), key=key)


def _lem_candidates(case):
    import run_tests
    from .core import solve

    path, out = case.path, []
    for tag in _tagged(case, _LEM_TAGS):
        out.append(("re-run %s tag (%s)" % (tag.get("type"), tag.get("method")),
                    lambda t=tag: _run_tag(run_tests, t, path)))
    circular = not (not (case.model or {}).get("circular")
                    and (case.model or {}).get("non_circ"))
    surface = "circular" if circular else "noncircular"
    methods = [case.method] + [str(t.get("method") or "").lower()
                               for t in _tagged(case, _LEM_TAGS)] + ["bishop"]
    seen, ordered = set(), []
    for method in methods:
        if method and method not in seen:
            seen.add(method)
            ordered.append(method)
    plain = {"type": "circular_search"}
    for method in ordered:
        for search in (True, False):
            how = "search" if search else "single surface"
            if case.stages_seep:
                # A model whose materials read pore pressure from a seepage
                # solution has no field until that solve is made. The tag runner
                # stages it (`seep=steady`); driving the LEM entry point straight
                # would be refused by preflight, which is a fact about the run
                # this scorer is making rather than about the answer.
                test = {"type": ("circular_search" if search else "single_circle")
                        if circular else
                        ("noncircular_search" if search else "single_noncirc"),
                        "method": method, "num_slices": 40, "seep": "steady"}
                out.append(("independent %s %s (%s, seepage staged)"
                            % (surface, how, method),
                            lambda t=test: _run_tag(run_tests, t, path)))
            else:
                out.append(("independent %s %s (%s)" % (surface, how, method),
                            lambda m=method, s=search:
                            _unpack(solve(path, method=m, search=s), "FS", plain)))
    return out


def _seep_candidates(case):
    import run_tests
    from .core import seep_flow

    path, out = case.path, []
    for tag in _tagged(case, _SEEP_TAGS):
        out.append(("re-run seep tag",
                    lambda t=tag: _run_tag(run_tests, t, path)))
    out.append(("independent steady seepage solve",
                lambda: _unpack(seep_flow(path), "q", {"type": "seep"})))
    return out


def _fem_candidates(case):
    from .core import ssrm_fs

    path, out = case.path, []
    for tag in _tagged(case, _FEM_TAGS):
        if tag.get("expected_fs") is not None:
            # An SSRM solve is minutes. Where the tag publishes the answer, that
            # published number IS an independent run — made by the regression
            # suite, at the mesh the page names — and re-making it here would put
            # the cost of scoring a sweep above the cost of the sweep.
            out.append(("published fem_ssrm tag",
                        lambda t=tag: (float(t["expected_fs"]), None,
                                       tolerance_for(t, float(t["expected_fs"])))))
    out.append(("independent SSRM run",
                lambda: _unpack(ssrm_fs(path), "FS", {"type": "fem_ssrm"})))
    return out


def _run_tag(run_tests, tag, path):
    """One published tag, re-run through the tag runner, against ``path``."""
    test = dict(tag)
    test["file"] = path
    try:
        value, error, _note = run_tests.run_test(test)
    except Exception as exc:
        return None, "%s: %s" % (type(exc).__name__, exc), None
    return value, error, tolerance_for(tag, value)


def _unpack(run, field, tag):
    value = run.get(field)
    return value, run.get("error"), tolerance_for(tag, value)


def tolerance_for(tag, value):
    """The tag runner's own tolerance for one answer.

    ``run_tests._expected_and_tol`` is the single place the regression suite
    decides how close is close enough, so a number this sweep accepts is a number
    that suite would accept. A seepage tolerance there is a fraction of the
    PUBLISHED discharge; where the run produced the number itself there is no
    published value to take a fraction of, and the same fraction is taken of the
    run. An SSRM answer is resolved no finer than the bracket it was closed to,
    which is the tag's own ``tolerance`` and 0.05 by default.
    """
    import run_tests

    if value is None:
        return None
    test = dict(tag or {})
    kind = test.get("type") or "circular_search"
    if kind in _SEEP_TAGS:
        return float(test.get("tolerance", 0.05)) * abs(float(value))
    if kind in _FEM_TAGS:
        return float(test.get("tolerance", 0.05))
    _expected, tol = run_tests._expected_and_tol(test, 0.01)
    return float(tol) if tol else 0.01


# --------------------------------------------------------------------------- #
# The criteria
# --------------------------------------------------------------------------- #
#: The exception types a SNIPPET is at fault for — a name that does not exist, an
#: argument that does not fit, a key that was never there. These are mistakes in
#: the code the assistant wrote.
_SNIPPET_BUG = re.compile(
    r"^\s*(NameError|AttributeError|TypeError|KeyError|IndexError|ImportError|"
    r"ModuleNotFoundError|SyntaxError|IndentationError|UnboundLocalError|"
    r"ZeroDivisionError|FileNotFoundError)\s*:", re.M)
#: What the ENGINE raises when it has no answer: a method with no admissible
#: solution, a surface that cannot be built, inputs preflight refuses. A run that
#: comes back this way is a finding about the model, and an answer that passes it
#: on is doing its job — so it is not counted against the session.
_ENGINE_REFUSAL = re.compile(
    r"^\s*(RuntimeError|AnalysisError|PreflightError|ValueError|"
    r"AnalysisCancelled)\s*:", re.M)


def no_snippet_errors():
    """No snippet was broken code.

    The distinction that matters is between a snippet that does not work and a
    model that has no answer. ``run_lem`` raising "only solutions with anomalous
    base tension found" is the engine reporting on the model, and a session that
    passes that on has done the right thing; ``AttributeError: module
    'xslope.solve' has no attribute 'generate_failure_surface'`` is a snippet
    calling something that does not exist. Only the second is counted.
    """
    def check(ctx):
        bugs, refusals = [], 0
        for i, turn in enumerate(ctx.turns, start=1):
            for out in turn.outputs:
                found = _SNIPPET_BUG.search(out)
                if found:
                    bugs.append("turn %d: %s" % (i, found.group(1)))
                elif _ENGINE_REFUSAL.search(out):
                    refusals += 1
            if turn.error:
                bugs.append("turn %d: %s" % (i, str(turn.error)[:48]))
        if bugs:
            return False, "; ".join(bugs[:3])
        if refusals:
            return True, "no broken snippet; %d run(s) the engine refused" % refusals
        return True, "no snippet raised"
    return Criterion("no exceptions in snippets", check, kind="discipline")


_RAN_ENGINE = ((r"\brun_fem\s*\(|\bsolve_ssrm\s*\(", "fem"),
               (r"\brun_seep\s*\(|run_seepage_analysis\s*\(|"
                r"run_transient_seepage\s*\(", "seep"),
               (r"\brun_lem\s*\(|\bcircular_search\s*\(|\bnoncircular_search\s*\(",
                "lem"))


def engine_ran(ctx):
    """Which engine the session actually drove, or ``None``."""
    for pattern, name in _RAN_ENGINE:
        if re.search(pattern, ctx.code):
            return name
    return None


def reported_value_matches(case):
    """The number the answer reports matches an independent run of the same file.

    ``run_tests`` makes that run: the published tag where the file carries one —
    literally the same call the regression suite makes — and the engine entry
    point otherwise. The comparison is at the tag runner's own tolerance, so a
    number the sweep accepts is a number the regression suite would accept.

    An answer that reports NO number passes only when no independent run
    produces one either. Several shipped workbooks carry a starting circle that
    daylights once, and no engine has an answer for them as they stand; a session
    that says so and asks how to proceed is right, and scoring it as a miss would
    measure the corpus rather than the assistant.
    """
    def check(ctx):
        ran = engine_ran(ctx)
        stated = _stated_values(case, ctx, ran)
        tried, answered = [], False
        for label, run in truth_candidates(case, ran):
            value, error, tol = run()
            if value is None:
                tried.append("%s failed (%s)" % (label, _short(error)))
                continue
            answered = True
            tol = tol or 0.01
            tried.append("%s = %s" % (label, _round(value)))
            if not stated:
                # Nothing to compare: the only question left is whether ANY run
                # answers this file, and one that does settles it.
                break
            best = min(stated, key=lambda v: abs(v - value))
            if abs(best - value) <= tol:
                return True, "stated %s vs %s = %s (tol %.4g)" % (
                    best, label, _round(value), tol)
        if not stated:
            if not answered:
                return True, ("no run of this model answers as the file declares "
                              "it, and the answer states no number: %s"
                              % (tried[0] if tried else "nothing to run"))
            return False, "the answer states no number; %s" % "; ".join(tried[:3])
        return False, "stated %s; %s" % (
            ", ".join(str(v) for v in stated[:4]), "; ".join(tried[:4])
            or "no independent run could be made")
    return Criterion("reported value matches an independent run", check,
                     kind="truth")


def _short(text, limit=90):
    """One line of an error message, for a scorecard cell."""
    line = " ".join(str(text or "").split())
    return line[:limit] + ("…" if len(line) > limit else "")


#: A markdown table with a factor-of-safety COLUMN. ``claimed_fs`` reads a table
#: whose ROW is labelled FS — the before/after shape an edit is reported in — and
#: misses the other layout entirely: a method-by-method comparison, where FS is a
#: column and each row is a method. Both are answers.
def _table_fs(prose):
    """Numbers in the body of any markdown table with an FS column."""
    lines = strip_code(prose).splitlines()
    out, in_table = [], False
    for line in lines:
        stripped = line.strip()
        if not stripped.startswith("|"):
            in_table = False
            continue
        cells = [c.strip().lower() for c in stripped.strip("|").split("|")]
        if any(re.fullmatch(r"\**\s*(fs|f\.s\.|factor of safety)[^|]*\**", c)
               for c in cells):
            in_table = True
            continue
        if in_table:
            out += numbers_in(stripped)
    return out


def _stated_values(case, ctx, ran=None):
    """The numbers the answer puts forward as its result."""
    prose = ctx.final_prose or ctx.prose
    if (ran or case.kind) == "seep":
        # A discharge has no "FS =" to anchor on; every number in the answer is a
        # candidate, and the tolerance is what discriminates.
        return [abs(v) for v in numbers_in(strip_code(prose))]
    return (claimed_fs(prose) or _table_fs(prose)
            or claimed_fs(ctx.prose) or _table_fs(ctx.prose))


def _round(value):
    try:
        return round(float(value), 6)
    except Exception:
        return value


_FORCED_METHOD = re.compile(r"method\s*=\s*['\"](\w+)['\"]")
_NAMES_METHOD = re.compile(r"\b(oms|bishop|janbu|corps|lowe|spencer|"
                           r"morgenstern[- ]price|mprice)\b", re.I)


def ran_declared_method(case):
    """The run was the file's own method, and a file that declares one has it
    named in the answer.

    Asked only of a session that made a stability solve. A model can define a
    failure surface and still be a page's finite-element example, and a session
    that answered it with SSRM has not ignored a declared method — it ran a
    different analysis, which the value criterion judges on its own terms.
    """
    def check(ctx):
        if case.kind != "lem":
            return True, "%s model — no LEM method to declare" % case.kind
        forced = {m.group(1).lower() for m in _FORCED_METHOD.finditer(ctx.code)}
        wrong = sorted(f for f in forced if f != case.method)
        if wrong:
            return False, "model declares %s; a snippet forced %s" % (
                case.method, ", ".join(wrong))
        if not re.search(r"\brun_lem\s*\(", ctx.code):
            return True, "no stability solve was made"
        if not case.declares_method:
            return True, "the file declares no method; nothing else was forced"
        said = _NAMES_METHOD.search(strip_code(ctx.prose))
        name = (said.group(1).lower().replace("morgenstern-price", "mprice")
                .replace("morgenstern price", "mprice")) if said else None
        if name != case.method:
            return False, "model declares %s; the answer says %s" % (
                case.method, name or "nothing")
        return True, "ran and named %s" % case.method
    return Criterion("ran the method the file declares", check, kind="truth")


def saved_model_still_holds(case):
    """A workbook the session saved reloads, and still reproduces every published
    value its source file carries.

    Most sessions save nothing — the question does not ask for an edit — and this
    passes by saying so. The ones that DO save are the sweep's whole reason for
    existing: an answer can be right about the number and still leave a model
    behind that no longer loads, or that loads and now solves at something else.
    """
    def check(ctx):
        saved = ctx.saved()
        if not saved:
            return True, "the session saved nothing"
        reloaded = load_model(saved)
        if reloaded is None:
            return False, "the saved workbook does not load"
        published = [t for t in case.tags
                     if t.get("expected_fs") is not None
                     or t.get("expected_flowrate") is not None]
        # An SSRM lock costs minutes to re-check and says the same thing about a
        # corrupted save that a stability or seepage lock on the same file says in
        # seconds, so it is re-run only where it is the file's only published value.
        cheap = [t for t in published if t.get("type") not in _FEM_TAGS]
        published = cheap or published[:1]
        if not published:
            return True, "the saved workbook reloads; the file publishes no value"
        import run_tests

        verify = _beside_companions(case.path, saved, ctx.outdir)
        broken = []
        for tag in published:
            value, error, _tol = _run_tag(run_tests, tag, verify)
            want = tag.get("expected_fs")
            if want is None:
                want = tag.get("expected_flowrate")
            if value is None:
                broken.append("%s: %s" % (tag.get("type"), error))
                continue
            _expected, tol = run_tests._expected_and_tol(tag, 0.01)
            if abs(value - float(want)) > (tol or 0.01):
                broken.append("%s: %s vs published %s"
                              % (tag.get("type"), _round(value), want))
        if broken:
            return False, "; ".join(broken[:3])
        return True, "the saved workbook reloads and holds %d published value(s)" \
                     % len(published)
    return Criterion("a saved workbook still holds its locks", check, kind="edit")


def _beside_companions(source, saved, outdir):
    """``saved``, copied where its source's sidecars can be auto-discovered.

    The sidecars are found by name beside the .xlsx, and the session writes its
    workbook into the run directory rather than beside the copy it opened — so a
    seepage model re-run on the saved file would arrive with no field to read.
    """
    import shutil

    from .faults import companions

    verify = os.path.join(outdir, "verify_saved")
    os.makedirs(verify, exist_ok=True)
    dest = os.path.join(verify, os.path.basename(source))
    shutil.copy2(saved, dest)
    for extra in companions(source):
        shutil.copy2(extra, os.path.join(verify, os.path.basename(extra)))
    return dest


# --------------------------------------------------------------------------- #
# Scenarios
# --------------------------------------------------------------------------- #
def scenario_for(case, completions=4):
    """The one-turn scenario that sweeps ``case``."""
    return Scenario(
        case.name, case.primary, case.path, [PROMPT],
        [reported_value_matches(case),
         ran_declared_method(case),
         no_snippet_errors(),
         _no_exploration(),
         _completions_at_most(completions),
         saved_model_still_holds(case)],
        what="%s — %s" % (case.kind.upper(), ", ".join(case.classes)),
        timeout_s=case.timeout_s)


def _no_exploration():
    from .scorers import no_exploration
    return no_exploration()


def _completions_at_most(n):
    from .scorers import completions_at_most
    return completions_at_most(n)


def cases(paths=None):
    """Every case, ordered so a run stopped at the budget still spans the classes.

    Round-robin over the primary-class buckets. A sweep that stopped a third of
    the way through an alphabetical list would have measured three directories;
    stopping a third of the way through this one has measured a third of every
    column, which is the thing the scorecard is grouped by.
    """
    built = [Case(p) for p in (paths or discover())]
    buckets = {}
    for case in built:
        buckets.setdefault(case.primary, []).append(case)
    for bucket in buckets.values():
        bucket.sort(key=lambda c: c.path)
    order = [c for c in CLASS_ORDER if c in buckets]
    order += [k for k in sorted(buckets) if k not in order]
    out, depth = [], 0
    while len(out) < len(built):
        for label in order:
            bucket = buckets[label]
            if depth < len(bucket):
                out.append(bucket[depth])
        depth += 1
    return out


def scenarios(paths=None, completions=4):
    """The sweep, as scenarios the runner plays. Files that do not load are left
    out — a session cannot open one, and the scorecard reports them separately."""
    return [scenario_for(case, completions) for case in cases(paths)
            if case.loads]


# --------------------------------------------------------------------------- #
# The scorecard, grouped by input class
# --------------------------------------------------------------------------- #
def render(results, meta):
    """``scorecard.md`` — the class table first, then every failing file."""
    index = {c.name: c for c in meta.get("cases") or []}
    lines = ["# Assistant corpus sweep — scorecard", ""]
    lines.append("- Run: `%s`" % meta.get("mode"))
    lines.append("- Provider / model: `%s` / `%s`" % (meta.get("provider"),
                                                      meta.get("model")))
    lines.append("- Recorded: %s" % meta.get("recorded"))
    lines.append("- Directory: `%s`" % meta.get("outdir"))
    if meta.get("mode") != "replay":
        lines.append("- Spend: **$%.2f**%s" % (
            meta.get("spend", 0.0),
            "  — **stopped at the $%.2f cap**" % meta["budget"]
            if meta.get("stopped") else ""))
    swept = len(results)
    total = meta.get("corpus_size") or swept
    passed = sum(1 for r in results if r["pass"])
    lines += ["", "**%d of %d files swept; %d pass.**" % (swept, total, passed)]
    unloadable = meta.get("unloadable") or []
    if unloadable:
        lines += ["", "%d workbook(s) do not load through `load_slope_data` and "
                      "were not run: %s" % (len(unloadable),
                                            ", ".join(sorted(unloadable)[:8]))]

    lines += ["", "## By input class", "",
              "A file is counted under every class it exercises, so the columns "
              "sum to more than the number of files.", "",
              "| Input class | Swept | Pass | Fail | Failing files |",
              "| --- | ---: | ---: | ---: | :--- |"]
    for label in _classes_present(results, index):
        rows = [r for r in results
                if label in (index[r["scenario"]].classes
                             if r["scenario"] in index else [])]
        bad = [r["scenario"] for r in rows if not r["pass"]]
        lines.append("| %s | %d | %d | %d | %s |" % (
            label, len(rows), len(rows) - len(bad), len(bad),
            ", ".join("`%s`" % n for n in bad[:6]) + (" …" if len(bad) > 6 else "")
            or "—"))

    lines += ["", "## By criterion", "",
              "| Criterion | Pass | Fail |", "| --- | ---: | ---: |"]
    for name in _criteria_names(results):
        rows = [c for r in results for c in r["criteria"] if c["name"] == name]
        lines.append("| %s | %d | %d |" % (
            name, sum(1 for c in rows if c["pass"]),
            sum(1 for c in rows if not c["pass"])))

    failing = [r for r in results if not r["pass"]]
    if failing:
        lines += ["", "## Every failure", ""]
        for row in failing:
            case = index.get(row["scenario"])
            lines += ["### `%s`" % row["scenario"], ""]
            if case:
                lines.append("- File: `%s`" % os.path.relpath(case.path,
                                                              repo()))
                lines.append("- Engine / classes: %s — %s"
                             % (case.kind.upper(), ", ".join(case.classes)))
            if row.get("error"):
                lines.append("- Session error: %s" % row["error"])
            for crit in row["criteria"]:
                if not crit["pass"]:
                    lines.append("- **%s** — %s" % (crit["name"], crit["reason"]))
            lines.append("")

    lines += ["", "## Every file", "",
              "| File | Engine | Primary class | Criteria | Result | Cost |",
              "| --- | --- | --- | ---: | :--- | ---: |"]
    for row in sorted(results, key=lambda r: r["scenario"]):
        case = index.get(row["scenario"])
        lines.append("| `%s` | %s | %s | %d/%d | %s | $%.3f |" % (
            row["scenario"], case.kind.upper() if case else "?",
            row["family"], row["passed"], row["total"],
            "pass" if row["pass"] else "FAIL", row.get("cost") or 0.0))
    return "\n".join(lines) + "\n"


def _classes_present(results, index):
    seen = set()
    for row in results:
        case = index.get(row["scenario"])
        if case:
            seen.update(case.classes)
    ordered = [c for c in CLASS_ORDER if c in seen]
    return ordered + sorted(seen - set(ordered))


def _criteria_names(results):
    names = []
    for row in results:
        for crit in row["criteria"]:
            if crit["name"] not in names:
                names.append(crit["name"])
    return names
