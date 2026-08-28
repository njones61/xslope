"""Standing checks on the scored assistant scenario suite (``tools/assistant_suite.py``).

The suite is the thing that measures the assistant, so a defect in the suite is a
measurement that reads green while the assistant is wrong. It is also the one
harness in the pipeline that cannot be exercised casually — a live run makes
billed provider calls — so everything provable without one is proved here.

What is asserted, all offscreen, no provider contacted and no network at all:

  A. THE REGISTRY IS SOUND — every scenario has a unique name, a family, at least
     one turn and at least one criterion; every workbook and every attached image
     it names exists on disk; every fault attaches to a scenario that opens a real
     model.
  B. PLANTING IS SURGICAL — each diagnosis scenario's faults apply to a COPY, the
     copy loads, the copy differs from the original in the planted fields, and the
     original file is byte-identical afterwards. A suite that quietly edited the
     repository's models would corrupt every other test that reads them.
  C. THE SCORERS DISCRIMINATE — a criterion that cannot fail measures nothing, so
     each family of criteria is run against a synthetic session built to pass it
     and one built to fail it, and both answers are required.
  D. THE TRANSCRIPT PARSER — a recorded body splits back into the assistant's
     prose, the snippets and the printed output, with the three kept apart. Every
     truth criterion rests on that separation: a number in the OUTPUT is measured,
     the same number in the PROSE is only claimed.
  E. A DRY RUN IS FREE AND COMPLETE — two cheap scenarios run end to end through
     the real window, the real dock and the real scoring path, writing a
     transcript, a dock capture, a ``session.json`` and a scorecard; and
     ``litellm.completion`` is not entered once.
  F. REPLAY IS DETERMINISTIC — re-scoring that same directory reproduces every
     criterion's answer, which is what makes ``--replay`` usable for fixing a
     scorer against a run already paid for.
  G. NOTHING LANDS IN docs/ — the run writes only under the directory it was
     given, and ``docs/tutorials/files`` gains no file.
  H. THE PRICE TABLE — one table, dated, complete (every model has both rates),
     and the cost arithmetic matches a hand computation including the cache-read
     discount.
  I. THE CORPUS SWEEP — the file list is complete and deduplicated the way the
     rule says, each file's input classes and engine are read off the loaded
     model, its criteria tell a broken snippet from a solver with no solution and
     read a factor of safety out of a table column as well as an equation, and a
     dry sweep of two files writes a scorecard grouped by input class without
     contacting a provider.
  J. AN EMPTY REINFORCEMENT LIST MEANS NONE — a sensitivity sweep addressed at
     reinforcement on a model whose ``reinforcement_lines`` was explicitly
     cleared says so, rather than falling through to the derived point lists,
     which carry no label for it to find the line by.

Skips cleanly (exit 0) when PySide6 is not installed (engine-only install — there
is no Studio layer to drive).
"""

import glob
import json
import os
import re
import shutil
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

#: The two scenarios the dry leg plays. Both are chosen because every one of their
#: criteria reads files and text only — no criterion in either one solves a model,
#: so the leg is seconds rather than the minutes a search would cost.
DRY_SCENARIOS = ("ambiguous_edit", "conceptual_search")


# --------------------------------------------------------------------------- #
# A. the registry
# --------------------------------------------------------------------------- #
def check_registry():
    from tools.assistant_scenarios import SCENARIOS

    out = []
    if len(SCENARIOS) < 25:
        out.append("only %d scenarios registered; the suite is meant to cover "
                   "every engine and every family" % len(SCENARIOS))
    seen = set()
    for scenario in SCENARIOS:
        if scenario.name in seen:
            out.append("duplicate scenario name %r" % scenario.name)
        seen.add(scenario.name)
        if not scenario.turns:
            out.append("%s: no turns" % scenario.name)
        if not scenario.criteria:
            out.append("%s: no criteria — it can neither pass nor fail"
                       % scenario.name)
        if scenario.model is not None and not os.path.exists(scenario.model):
            out.append("%s: model does not exist: %s"
                       % (scenario.name, scenario.model))
        if scenario.faults and scenario.model is None:
            out.append("%s: faults with no model to plant them in" % scenario.name)
        for turn in scenario.turns:
            if isinstance(turn, (tuple, list)):
                for image in turn[1:]:
                    for path in ([image] if isinstance(image, str) else image or []):
                        if not os.path.exists(path):
                            out.append("%s: attachment does not exist: %s"
                                       % (scenario.name, path))
    families = {s.family for s in SCENARIOS}
    for wanted in ("build", "edit", "seep", "fem", "sweep", "conceptual",
                   "diagnose", "judgment", "report"):
        if wanted not in families:
            out.append("no scenario in the %r family" % wanted)
    # The eight W-1 prompts ride along unchanged, so a scored run is comparable
    # with the recorded ones the tutorial was built from.
    w1 = {s.name for s in SCENARIOS if s.name.startswith("w1_")}
    if len(w1) != 8:
        out.append("%d W-1 scenarios carried over, expected 8 (%s)"
                   % (len(w1), ", ".join(sorted(w1))))
    return out


# --------------------------------------------------------------------------- #
# B. planting
# --------------------------------------------------------------------------- #
def check_planting():
    from tools.assistant_scenarios import SCENARIOS, load_model
    from tools.assistant_scenarios.core import digest
    from tools.assistant_scenarios.faults import plant

    out = []
    tmp = tempfile.mkdtemp(prefix="suite_plant_")
    try:
        for scenario in SCENARIOS:
            if not scenario.faults:
                continue
            before_digest = digest(scenario.model)
            sound = load_model(scenario.model)
            broken_path = plant(scenario.model, scenario.faults,
                                os.path.join(tmp, scenario.name))
            if digest(scenario.model) != before_digest:
                out.append("%s: planting modified the repository's own %s"
                           % (scenario.name, os.path.basename(scenario.model)))
            broken = load_model(broken_path)
            if broken is None:
                out.append("%s: the planted copy does not load" % scenario.name)
                continue
            if _same_model(sound, broken):
                out.append("%s: the planted copy is identical to the sound model"
                           % scenario.name)
            # And a copy with NO faults must round-trip unchanged, or a "difference"
            # above could be the save path rather than the fault.
            from tools.assistant_scenarios import Fault
            noop = Fault("noop", "no change at all", lambda sd: None, r"(?!)")
            same = plant(scenario.model, [noop],
                         os.path.join(tmp, scenario.name + "_rt"))
            if not _same_model(sound, load_model(same)):
                out.append("%s: saving the model unchanged already moves it — a "
                           "planted difference cannot be attributed to the fault"
                           % scenario.name)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return out


def _same_model(a, b):
    """Whether two loaded models agree on the fields a fault can reach."""
    if a is None or b is None:
        return a is b
    def snap(sd):
        mats = [(m.get("name"), m.get("gamma"), m.get("c"), m.get("phi"),
                 m.get("u")) for m in sd.get("materials") or []]
        loads = [[(p.get("X"), p.get("Normal")) for p in row]
                 for row in sd.get("dloads") or []]
        circles = [(c.get("Xo"), c.get("Yo"), c.get("Depth"))
                   for c in sd.get("circles") or []]
        return (mats, loads, circles, sd.get("max_depth"),
                len(sd.get("piezo_line") or []))
    return snap(a) == snap(b)


# --------------------------------------------------------------------------- #
# C. the scorers discriminate
# --------------------------------------------------------------------------- #
def _session(prose="", code="", output="", **extra):
    """A synthetic recorded session, in the shape the recorder writes."""
    body = []
    if prose:
        body.append("Assistant: " + prose)
    if code:
        body.append("Ran code:\n" + "\n".join("    " + ln
                                              for ln in code.splitlines()))
    if output:
        body.append("Output:\n" + "\n".join("    " + ln
                                            for ln in output.splitlines()))
    record = {"prompt": "p", "body": "\n\n".join(body),
              "usage": {"calls": 1, "input": 10, "output": 5}, "seconds": 1.0,
              "error": None, "workbook": None, "image": None}
    session = {"turns": [record], "usage": {"input": 10, "output": 5, "calls": 1},
               "seconds": 1.0, "error": None, "workbook": None,
               "start_model": None, "transcript": None}
    session.update(extra)
    return session


def _ask(criterion, session, scenario=None):
    from tools.assistant_scenarios import ScoreCtx, Scenario

    scenario = scenario or Scenario("t", "test", None, ["p"], [])
    return criterion(ScoreCtx(scenario, session, tempfile.gettempdir()))[0]


def check_scorers():
    from tools.assistant_scenarios import Scenario, scorers as S
    from tools.assistant_scenarios import faults as F

    out = []

    def leg(label, criterion, good, bad, scenario=None):
        if not _ask(criterion, good, scenario):
            out.append("%s: a session that should pass it does not" % label)
        if _ask(criterion, bad, scenario):
            out.append("%s: a session that should fail it passes" % label)

    leg("no_latex", S.no_latex(),
        _session(prose="tan phi / tan beta = 0.066"),
        _session(prose=r"$$\frac{\tan\varphi}{\tan\beta} = 0.066$$"))

    leg("no_exploration", S.no_exploration(),
        _session(code="res = run_lem(search=True)\nprint(res['FS'])"),
        _session(code="import inspect\nprint(inspect.getsource(xslope.solve))"))

    leg("completions_at_most", S.completions_at_most(2),
        _session(prose="ok"),
        _session(prose="ok", turns=[{"prompt": "p", "body": "Assistant: ok",
                                     "usage": {"calls": 9}, "seconds": 1.0}]))

    leg("claims_grounded", S.claims_grounded(),
        _session(prose="Spencer gives FS = 1.587.", output="FS 1.5867431"),
        _session(prose="Spencer gives FS = 1.587.", output="nothing was printed"))

    leg("cites_docs", S.cites_docs(1),
        _session(prose="See https://xslope.readthedocs.io/en/latest/lem/spencer/"),
        _session(prose="See https://xslope.readthedocs.io/en/latest/lem/invented/"))

    leg("cites_corpus", S.cites_corpus(),
        _session(prose="VP28 is the worked example.",
                 code="print(corpus_index('reliability'))",
                 output="{'title': 'VP28', 'url': 'x'}"),
        _session(prose="VP28 is the worked example."))

    leg("asks_rather_than_guesses", S.asks_rather_than_guesses(),
        _session(prose="How far up should it go?"),
        _session(prose="Raised it by 5 ft.", workbook="/nonexistent/after.xlsx"))

    leg("reports_warnings", S.reports_warnings(),
        _session(prose="The polygons were rebuilt from profile_lines, so I edited "
                       "profile_lines instead.",
                 output="WARNING: polygons were edited on a profile-line model "
                        "and have been rebuilt from profile_lines"),
        _session(prose="Done.",
                 output="WARNING: polygons were edited on a profile-line model "
                        "and have been rebuilt from profile_lines"))

    leg("mechanism_claims_tested", S.mechanism_claims_tested(),
        _session(prose="FS plateaus because the pullout caps it.",
                 code="slope_data['materials'][0]['phi'] = 30\nrun_lem(search=True)"),
        _session(prose="FS plateaus because the pullout caps it.",
                 code="print(slope_data['materials'])"))

    leg("varied_inputs", S.varied_inputs(1),
        _session(code="slope_data['max_depth'] = -10\nrun_lem(search=True)"),
        _session(code="print(slope_data['max_depth'])"))

    leg("states_the_measurement", S.states_the_measurement(),
        _session(prose="Restoring it takes FS from 0.0709 to 1.1805."),
        _session(prose="Restoring it would help a great deal."))

    leg("stated_values_grounded",
        S.stated_values_grounded(r"(?:reliability index|beta|β)", "beta"),
        _session(prose="beta = 2.943", output="beta_ln 2.9431"),
        _session(prose="beta = 2.50", output="beta_ln 2.9431"))

    leg("model_unchanged", S.model_unchanged(),
        _session(),
        _session(workbook="/nonexistent/after.xlsx"))

    # The fault scorers need a scenario carrying faults.
    planted = Scenario("t", "test", None, ["p"], [],
                       faults=[F.phi_typo("base", 3, 37)])
    leg("faults_named", S.faults_named(),
        _session(prose="The base friction angle reads 3 where 37 was meant."),
        _session(prose="Everything looks fine to me."), scenario=planted)

    leg("no_false_accusation", S.no_false_accusation(r"sliver"),
        _session(prose="The zones are as designed."),
        _session(prose="The shell zone is a degenerate sliver."))

    return out


# --------------------------------------------------------------------------- #
# D. the transcript parser
# --------------------------------------------------------------------------- #
def check_parser():
    from tools.assistant_scenarios import parse_body

    body = ("Assistant: The factor of safety is 1.587.\n\n"
            "It is the published answer.\n\n"
            "Ran code:\n    res = run_lem(search=True)\n    print(res['FS'])\n\n"
            "Output:\n    1.5867431\n    === MODEL CHECKS: clean ===\n\n"
            "[file: plot.png]")
    texts, codes, outputs, files = parse_body(body)
    out = []
    if len(texts) != 1 or "1.587" not in texts[0]:
        out.append("the assistant's prose did not come back whole: %r" % (texts,))
    if "It is the published answer." not in "\n".join(texts):
        out.append("a blank line inside the prose split it into two blocks")
    if len(codes) != 1 or "run_lem" not in codes[0]:
        out.append("the snippet did not come back: %r" % (codes,))
    if len(outputs) != 1 or "MODEL CHECKS" not in outputs[0]:
        out.append("the printed output did not come back: %r" % (outputs,))
    if "1.5867431" in "\n".join(texts):
        out.append("printed output leaked into the prose — a claim would count "
                   "as its own evidence")
    if files != ["plot.png"]:
        out.append("the file list did not come back: %r" % (files,))
    return out


# --------------------------------------------------------------------------- #
# E/F/G. a dry run, a replay, and nothing in docs/
# --------------------------------------------------------------------------- #
def check_dry_run():
    import tools.assistant_suite as suite
    from tools.assistant_scenarios import by_name

    out = []
    calls = {"n": 0}
    try:
        import litellm
        real = litellm.completion

        def counted(*args, **kwargs):
            calls["n"] += 1
            return real(*args, **kwargs)
        litellm.completion = counted
    except Exception:
        litellm = real = None

    docs_files = os.path.join(REPO_ROOT, "docs", "tutorials", "files")
    before_docs = set(os.listdir(docs_files)) if os.path.isdir(docs_files) else set()
    tmp = tempfile.mkdtemp(prefix="suite_dry_")
    try:
        scenarios = [by_name(n) for n in DRY_SCENARIOS]
        results, meta = suite.run(scenarios, tmp, dry_run=True)
        if len(results) != len(scenarios):
            out.append("%d of %d scenarios recorded" % (len(results),
                                                        len(scenarios)))
        if calls["n"]:
            out.append("a dry run entered litellm.completion %d time(s) — "
                       "--dry-run must never cost anything" % calls["n"])
        for name in DRY_SCENARIOS:
            here = os.path.join(tmp, name)
            for wanted in ("session.json", "sx_%s_transcript.md" % name):
                if not os.path.exists(os.path.join(here, wanted)):
                    out.append("%s: %s was not written" % (name, wanted))
            shot = os.path.join(here, "images", "sx_%s_1.png" % name)
            if not os.path.exists(shot) or os.path.getsize(shot) < 1000:
                out.append("%s: the dock capture is missing or empty" % name)
        for wanted in ("scorecard.md", "scorecard.json"):
            if not os.path.exists(os.path.join(tmp, wanted)):
                out.append("%s was not written" % wanted)
        with open(os.path.join(tmp, "scorecard.json"), encoding="utf-8") as fh:
            card = json.load(fh)
        if not card.get("results") or not card["results"][0].get("criteria"):
            out.append("the scorecard carries no criteria — scoring did not run")
        if meta.get("spend"):
            out.append("a dry run reported a spend of $%.4f" % meta["spend"])

        # F. replay reproduces every answer, still without a provider.
        before = {r["scenario"]: [(c["name"], c["pass"]) for c in r["criteria"]]
                  for r in results}
        again, _meta = suite.replay(tmp)
        after = {r["scenario"]: [(c["name"], c["pass"]) for c in r["criteria"]]
                 for r in again}
        if before != after:
            out.append("replay did not reproduce the scoring of the run it "
                       "re-scored")
        if calls["n"]:
            out.append("replay reached the provider")
    finally:
        if litellm is not None and real is not None:
            litellm.completion = real
        shutil.rmtree(tmp, ignore_errors=True)

    # G. the repository's docs are untouched.
    after_docs = set(os.listdir(docs_files)) if os.path.isdir(docs_files) else set()
    if after_docs - before_docs:
        out.append("the run wrote into docs/tutorials/files: %s"
                   % ", ".join(sorted(after_docs - before_docs)))
    return out


# --------------------------------------------------------------------------- #
# H. the price table
# --------------------------------------------------------------------------- #
def check_prices():
    import tools.assistant_suite as suite

    out = []
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", suite.PRICE_DATE):
        out.append("PRICE_DATE is not a date: %r" % suite.PRICE_DATE)
    for model, rates in suite.PRICES.items():
        for field in ("input", "output"):
            if not isinstance(rates.get(field), (int, float)) or rates[field] <= 0:
                out.append("%s has no %s rate" % (model, field))
    if "claude-opus-5" not in suite.PRICES:
        out.append("the model the suite runs on is not in the price table")
    # One million fresh input tokens, one million cached, one million out.
    got = suite.cost_of({"input": 2_000_000, "cached_input": 1_000_000,
                         "output": 1_000_000}, "claude-opus-5")
    want = 5.00 + 0.50 + 25.00
    if abs(got - want) > 1e-9:
        out.append("cost arithmetic: $%.4f, expected $%.4f" % (got, want))
    if suite.cost_of({"input": 10}, "a-model-that-does-not-exist"):
        out.append("an unknown model was priced rather than reported as $0")
    # A price quoted anywhere must come from the table: no other literal rate.
    source = open(os.path.join(REPO_ROOT, "tools", "assistant_suite.py"),
                  encoding="utf-8").read()
    body = source.split("PRICES = {", 1)[-1].split("}\n", 1)[-1]
    if re.search(r"\b(?:5\.00|25\.00|0\.50)\b\s*[*/]", body):
        out.append("a price literal appears outside the PRICES table")
    return out


# --------------------------------------------------------------------------- #
# I. the corpus sweep
# --------------------------------------------------------------------------- #
#: The two workbooks the dry sweep plays. Both are small, both are plain LEM, and
#: neither costs a solve to score: the stub reply states no number, so the value
#: criterion answers before it reaches the engine.
DRY_CORPUS = (os.path.join("docs", "lem", "files", "xslope_simple_embankment.xlsx"),
              os.path.join("docs", "inputs", "slope", "xslope_simple1.xlsx"))


def check_corpus_files():
    """The sweep's file list: complete, unique, and deduplicated by the rule."""
    from tools.assistant_scenarios import corpus

    out = []
    paths = corpus.discover()
    if len(paths) < 300:
        out.append("discover() found only %d workbooks; the shipped corpus is "
                   "several hundred" % len(paths))
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        out.append("discover() named %d file(s) that do not exist: %s"
                   % (len(missing), missing[0]))
    names = [corpus._case_name(p) for p in paths]
    if len(set(names)) != len(names):
        dupes = sorted({n for n in names if names.count(n) > 1})
        out.append("two workbooks map to one scenario name: %s" % ", ".join(dupes[:4]))
    for folder in corpus.CORPUS_DIRS + corpus.CORPUS_TREES:
        on_disk = glob.glob(os.path.join(REPO_ROOT, folder, "**", "*.xlsx"),
                            recursive=True)
        if on_disk and not any(p.startswith(os.path.join(REPO_ROOT, folder))
                               for p in paths):
            out.append("%s holds workbooks but none were swept" % folder)

    # The dedup rule, both ways round, on files written for the purpose.
    tmp = tempfile.mkdtemp(prefix="corpus_dedup_")
    try:
        base = os.path.join(tmp, "model.xlsx")
        same = os.path.join(tmp, "model_start.xlsx")
        other = os.path.join(tmp, "other_start.xlsx")
        with open(base, "wb") as fh:
            fh.write(b"identical bytes")
        shutil.copy2(base, same)
        with open(other, "wb") as fh:
            fh.write(b"different bytes entirely")
        if not corpus.duplicate_of_sibling(same):
            out.append("a _start copy identical to its sibling was not skipped")
        if corpus.duplicate_of_sibling(other):
            out.append("a _start file with no sibling was skipped as a duplicate")
        if corpus.duplicate_of_sibling(base):
            out.append("a plain workbook was treated as a companion copy")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return out


def check_corpus_classes():
    """Input classes and engine come off the loaded model, not off the name."""
    from tools.assistant_scenarios import corpus

    out = []
    wanted = {
        os.path.join("docs", "tutorials", "files",
                     "xslope_reinforced_slope.xlsx"):
            ("lem", {"reinforcement", "distributed loads"}),
        os.path.join("docs", "seep", "files", "xslope_earth_dam1.xlsx"):
            ("seep", {"seepage BCs"}),
        os.path.join("docs", "lem", "files", "xslope_noncircular.xlsx"):
            ("lem", {"non-circular"}),
    }
    for rel, (kind, classes) in wanted.items():
        case = corpus.Case(os.path.join(REPO_ROOT, rel))
        if not case.loads:
            out.append("%s does not load" % rel)
            continue
        if case.kind != kind:
            out.append("%s reads as a %s model, expected %s"
                       % (os.path.basename(rel), case.kind, kind))
        missing = classes - set(case.classes)
        if missing:
            out.append("%s: input class(es) not detected: %s"
                       % (os.path.basename(rel), ", ".join(sorted(missing))))
        if case.primary not in case.classes:
            out.append("%s: primary class %r is not one of its classes"
                       % (os.path.basename(rel), case.primary))
    # A file that publishes an SSRM value and also defines a failure surface is a
    # stability case: reading the kind off whichever tag came first filed the
    # LEM-8 tutorial model as FEM and scored its Spencer answer against SSRM.
    lem8 = corpus.Case(os.path.join(REPO_ROOT, "docs", "tutorials", "files",
                                    "xslope_reinforced_slope.xlsx"))
    if lem8.loads and lem8.kind != "lem":
        out.append("a model with circles and an fem_ssrm tag reads as %s"
                   % lem8.kind)
    # And every case is orderable into a bucket the scorecard can group by.
    order = set(corpus.CLASS_ORDER)
    for case in (lem8,):
        stray = [c for c in case.classes if c not in order]
        if stray:
            out.append("input class(es) outside CLASS_ORDER: %s" % ", ".join(stray))
    return out


def check_corpus_scorers():
    """The sweep's own criteria discriminate, on synthetic sessions only."""
    from tools.assistant_scenarios import corpus

    out = []
    engine_said_no = _session(
        prose="Spencer has no admissible solution on this circle.",
        code="run_lem(search=False)",
        output="RuntimeError: No solution: Spencer's method: only solutions "
               "with anomalous base tension found (96.0x cohesive capacity)")
    snippet_broke = _session(
        prose="Here is the surface.",
        code="xslope.solve.generate_failure_surface(sd)",
        output="Traceback (most recent call last):\n"
               "AttributeError: module 'xslope.solve' has no attribute "
               "'generate_failure_surface'")
    if not _ask(corpus.no_snippet_errors(), engine_said_no):
        out.append("a solver reporting no solution was counted as a broken snippet")
    if _ask(corpus.no_snippet_errors(), snippet_broke):
        out.append("a snippet calling a name that does not exist passed")

    # A factor of safety in a table COLUMN is an answer, the same as one after an
    # equals sign; reading only the second scored a correct method-by-method
    # answer as reporting no number at all.
    table = ("| Method | FS | Notes |\n|:--|--:|:--|\n"
             "| OMS | 1.398 | — |\n| Bishop | 1.237 | — |")
    if corpus._table_fs(table) != [1.398, 1.237]:
        out.append("an FS column was not read as stated values: %r"
                   % (corpus._table_fs(table),))
    if corpus._table_fs("no table at all"):
        out.append("prose with no table produced stated values")

    # Which engine ran is read off the snippets, so an SSRM answer is compared
    # with an SSRM run rather than with a limit-equilibrium one.
    from tools.assistant_scenarios import Scenario, ScoreCtx
    scenario = Scenario("t", "test", None, ["p"], [])
    for code, want in (("res = run_fem(analysis='ssrm')", "fem"),
                       ("sol = run_seep(bc=1)", "seep"),
                       ("res = run_lem(search=True)", "lem"),
                       ("print(slope_data['materials'])", None)):
        got = corpus.engine_ran(ScoreCtx(scenario, _session(code=code),
                                         tempfile.gettempdir()))
        if got != want:
            out.append("engine_ran(%r) = %r, expected %r" % (code, got, want))
    return out


def check_corpus_dry_run():
    """A dry sweep of two files: a class-grouped scorecard, and no provider."""
    import tools.assistant_suite as suite
    from tools.assistant_scenarios import corpus

    out = []
    calls = {"n": 0}
    try:
        import litellm
        real = litellm.completion

        def counted(*args, **kwargs):
            calls["n"] += 1
            return real(*args, **kwargs)
        litellm.completion = counted
    except Exception:
        litellm = real = None

    docs_files = os.path.join(REPO_ROOT, "docs", "tutorials", "files")
    before_docs = set(os.listdir(docs_files)) if os.path.isdir(docs_files) else set()
    tmp = tempfile.mkdtemp(prefix="corpus_dry_")
    try:
        cases = corpus.cases([os.path.join(REPO_ROOT, p) for p in DRY_CORPUS])
        scenarios = [corpus.scenario_for(c) for c in cases]
        if len(scenarios) != len(DRY_CORPUS):
            out.append("%d of %d corpus scenarios built"
                       % (len(scenarios), len(DRY_CORPUS)))
        results, meta = suite.run(scenarios, tmp, dry_run=True,
                                  renderer=corpus.render,
                                  meta_extra={"cases": cases,
                                              "corpus_size": len(cases),
                                              "unloadable": []})
        if calls["n"]:
            out.append("a dry sweep entered litellm.completion %d time(s)"
                       % calls["n"])
        card = os.path.join(tmp, "scorecard.md")
        if not os.path.exists(card):
            out.append("the sweep wrote no scorecard")
            return out
        text = open(card, encoding="utf-8").read()
        if "## By input class" not in text:
            out.append("the scorecard is not grouped by input class")
        for case in cases:
            if case.primary not in text:
                out.append("the scorecard never names the %r class"
                           % case.primary)
        if "## By criterion" not in text:
            out.append("the scorecard does not tally the criteria")
        # Every scenario carries the six questions the sweep asks.
        for row in results:
            names = [c["name"] for c in row["criteria"]]
            for wanted in ("reported value matches an independent run",
                           "ran the method the file declares",
                           "no exceptions in snippets",
                           "a saved workbook still holds its locks"):
                if wanted not in names:
                    out.append("%s: %r was not scored" % (row["scenario"], wanted))
        # The stub reply states no number, so the value criterion must FAIL —
        # a sweep whose truth check passes on an empty answer measures nothing.
        for row in results:
            value = next((c for c in row["criteria"]
                          if c["name"] == "reported value matches an independent run"),
                         None)
            if value and value["pass"]:
                out.append("%s: the value criterion passed on a stub reply"
                           % row["scenario"])
        # Replay re-scores a sweep by rebuilding each case from its session.
        again, _meta = suite.replay(tmp, corpus=True)
        if len(again) != len(results):
            out.append("corpus replay re-scored %d of %d sessions"
                       % (len(again), len(results)))
        before = {r["scenario"]: [(c["name"], c["pass"]) for c in r["criteria"]]
                  for r in results}
        after = {r["scenario"]: [(c["name"], c["pass"]) for c in r["criteria"]]
                 for r in again}
        if before != after:
            out.append("corpus replay did not reproduce the scoring of the run "
                       "it re-scored")
    finally:
        if litellm is not None and real is not None:
            litellm.completion = real
        shutil.rmtree(tmp, ignore_errors=True)
    after_docs = set(os.listdir(docs_files)) if os.path.isdir(docs_files) else set()
    if after_docs - before_docs:
        out.append("the sweep wrote into docs/tutorials/files: %s"
                   % ", ".join(sorted(after_docs - before_docs)))
    return out


# --------------------------------------------------------------------------- #
# J. an empty reinforcement list means none
# --------------------------------------------------------------------------- #
def check_empty_reinforcement():
    """Clearing ``reinforcement_lines`` removes the reinforcement everywhere.

    ``slice.py`` has read an explicit empty list as "no reinforcement" since the
    delete bug; the sensitivity sweep did not, and fell through to the derived
    ``reinforce_lines``. Those are point dicts with no ``label``, so a sweep
    addressed at a line by name looked for it in a list that cannot name it —
    reporting a missing line rather than a model with no reinforcement.
    """
    from xslope.sensitivity import resolve_param
    from tools.assistant_scenarios import load_model, repo

    out = []
    sd = load_model(repo("docs/tutorials/files/xslope_reinforced_slope.xlsx"))
    if sd is None:
        return ["the reinforced-slope model does not load"]
    if not sd.get("reinforcement_lines") or not sd.get("reinforce_lines"):
        return ["the fixture carries no reinforcement to clear"]
    label = sd["reinforcement_lines"][0].get("label")

    # Sound model: the line resolves.
    try:
        canonical, _setter, base = resolve_param(sd, "reinforce:%s:t_max" % label)
    except Exception as exc:
        out.append("a sound model's reinforcement line does not resolve: %s: %s"
                   % (type(exc).__name__, exc))
    else:
        if not canonical.startswith("reinforce:"):
            out.append("resolved to %r rather than a reinforcement line" % canonical)
        if base is None:
            out.append("the resolved line has no base value")

    # Source cleared, derived list left behind — what a delete leaves in memory.
    cleared = dict(sd)
    cleared["reinforcement_lines"] = []
    try:
        resolve_param(cleared, "reinforce:%s:t_max" % label)
    except ValueError as exc:
        if "no reinforcement lines" not in str(exc).lower():
            out.append("an emptied model raised the wrong message: %s" % exc)
    except Exception as exc:
        out.append("an emptied model raised %s: %s" % (type(exc).__name__, exc))
    else:
        out.append("an emptied 'reinforcement_lines' still resolved a line — the "
                   "derived point lists were used as the source")

    # No source key at all is a model assembled by hand, and the derived list is
    # still the only thing there is to read.
    legacy = {k: v for k, v in sd.items() if k != "reinforcement_lines"}
    try:
        resolve_param(legacy, "reinforce:%s:t_max" % label)
    except ValueError as exc:
        if "no reinforcement lines" in str(exc).lower():
            out.append("a model with no 'reinforcement_lines' key was told it has "
                       "no reinforcement, though 'reinforce_lines' carries it")
    except Exception:
        pass
    return out


# --------------------------------------------------------------------------- #
def run():
    try:
        import PySide6                                          # noqa: F401
    except Exception:
        print("PySide6 not installed — assistant suite checks skipped")
        return []
    failures = []
    for name, check in (("registry", check_registry),
                        ("planting", check_planting),
                        ("scorers", check_scorers),
                        ("parser", check_parser),
                        ("prices", check_prices),
                        ("corpus files", check_corpus_files),
                        ("corpus classes", check_corpus_classes),
                        ("corpus scorers", check_corpus_scorers),
                        ("empty reinforcement", check_empty_reinforcement),
                        ("dry run", check_dry_run),
                        ("corpus dry run", check_corpus_dry_run)):
        try:
            failures += ["%s: %s" % (name, msg) for msg in check()]
        except Exception as exc:
            failures.append("%s: raised %s: %s" % (name, type(exc).__name__, exc))
    return failures


if __name__ == "__main__":
    problems = run()
    for line in problems:
        print("FAIL " + line)
    print("%d problem(s)" % len(problems))
    raise SystemExit(1 if problems else 0)
