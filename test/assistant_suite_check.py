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

Skips cleanly (exit 0) when PySide6 is not installed (engine-only install — there
is no Studio layer to drive).
"""

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
                        ("dry run", check_dry_run)):
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
