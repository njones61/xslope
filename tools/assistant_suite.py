"""The scored scenario suite for the Studio assistant.

The W-1 recording measured eight conversations by hand: every number the assistant
reported was re-solved by a person, every edit was checked by opening the workbook.
That found real defects — and it does not scale, and it cannot be re-run after a
brief change to see whether the change helped. This is the same audit written down
as code.

A run plays each scenario through the real Studio window and the real assistant,
then scores the workbook it left behind and the transcript it produced against a
list of criteria that MEASURE rather than believe: the model is reloaded from disk
and re-solved here, in this process, and the assistant's numbers are compared with
that. Nothing is scored from the assistant's own account of what it did.

Three modes::

    python3 tools/assistant_suite.py --dry-run                # plumbing, no API
    python3 tools/assistant_suite.py --live --budget-usd 15   # the real thing
    python3 tools/assistant_suite.py --replay <dir>           # re-score a run

and a fourth axis across them. ``--corpus`` replaces the thirty hand-written
scenarios with EVERY WORKBOOK the repository ships — one session and one prompt
per file, scored against the same tag runner ``run_tests.py`` drives, and
reported by the input columns each file exercises rather than by scenario. Where
the scenarios cover the kinds of thing a person asks for, the sweep covers the
kinds of MODEL they ask about; see :mod:`tools.assistant_scenarios.corpus`::

    python3 tools/assistant_suite.py --corpus --live --budget-usd 10
    python3 tools/assistant_suite.py --corpus --replay <dir>

Sweeping every file with one prompt measures one verb. ``--sample`` narrows the
files and widens the questions: a stratified draw of the corpus — every input
class filled where the corpus can fill it — with a MENU of tasks asked of each
file, drawn at random from the ones that file can support. Describing a model,
editing it, explaining which slices carry the surface, comparing methods,
sweeping a parameter, finding a planted fault, writing the report; and, with
``--builds``, building a model from a tutorial's problem drawing or from the
prose it is described in, scored against the workbook that tutorial ships. See
:mod:`tools.assistant_scenarios.tasks`::

    python3 tools/assistant_suite.py --corpus --sample 60 --tasks 2 --seed 1 \\
        --live --budget-usd 25
    python3 tools/assistant_suite.py --corpus --sample 12 --tasks 2 --seed 1 --dry-run
    python3 tools/assistant_suite.py --corpus --sample 60 --builds --list

The draw is seeded, so a seed names a run and ``--replay`` re-scores exactly the
set that was played. Files the first live round already swept are skipped by
default (``--include-swept`` puts them back).

``--dry-run`` answers every turn with a stub and calls no provider at all; it is
what ``run_tests.py`` runs (through ``test/assistant_suite_check.py``) so the
harness itself stays regression-tested without spending anything. ``--replay``
re-scores saved sessions, so a scorer can be fixed and the same recorded run
re-graded without buying it twice.

Everything a run writes goes to a scratch directory — transcripts, dock captures,
the workbooks each session left behind, ``scorecard.md`` and ``scorecard.json``.
Nothing is ever written into ``docs/``.
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
import sys
import time

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from tools.assistant_scenarios import SCENARIOS, ScoreCtx, use_solve_cache  # noqa: E402
from tools.assistant_scenarios.faults import plant                          # noqa: E402


# --------------------------------------------------------------------------- #
# What a run costs
# --------------------------------------------------------------------------- #
#: Anthropic list prices in US dollars per million tokens, as published on
#: PRICE_DATE. This is the ONLY place a price appears: a number quoted anywhere in
#: a scorecard came from here. Prices move, so the date travels with the table and
#: is printed on every scorecard — a stale figure is then visible as a stale date
#: rather than hiding inside a total.
PRICE_DATE = "2026-08-27"
PRICES = {
    "claude-opus-5":    {"input": 5.00, "output": 25.00},
    "claude-opus-4-8":  {"input": 5.00, "output": 25.00},
    "claude-fable-5":   {"input": 10.00, "output": 50.00},
    "claude-sonnet-5":  {"input": 2.00, "output": 10.00},
    "claude-haiku-4-5": {"input": 1.00, "output": 5.00},
}
#: A prompt-cache read is a tenth of an input token; a five-minute cache WRITE is
#: 1.25x. The assistant's accumulator does not separate cache writes from ordinary
#: input, so a total here is a floor: it prices the writes as plain input.
CACHE_READ_FACTOR = 0.10

#: What one completion has actually cost, measured over the 82 sessions of the
#: first live corpus sweep rather than guessed: 17,887 input tokens of which 83%
#: were cache reads, and 681 output tokens. An ESTIMATE printed before a run is
#: built from these, so the number a budget is set against is the number the last
#: run produced. Re-measure them when the brief grows.
TOKENS_PER_COMPLETION = {"input": 17_887, "output": 681, "cached_fraction": 0.83}


def estimate_usd(completions, model="claude-opus-5"):
    """What ``completions`` turns of work costs, at the measured token rates."""
    rates = PRICES.get(model)
    if not rates or not completions:
        return 0.0
    per = TOKENS_PER_COMPLETION
    cached = per["input"] * per["cached_fraction"]
    fresh = per["input"] - cached
    return completions * (fresh * rates["input"]
                          + cached * rates["input"] * CACHE_READ_FACTOR
                          + per["output"] * rates["output"]) / 1e6


def cost_of(usage, model):
    """Dollars for one session's usage, at the list prices above."""
    rates = PRICES.get(model)
    if not rates:
        return 0.0
    total_in = int(usage.get("input") or 0)
    cached = min(int(usage.get("cached_input") or 0), total_in)
    fresh = total_in - cached
    out = int(usage.get("output") or 0)
    return (fresh * rates["input"]
            + cached * rates["input"] * CACHE_READ_FACTOR
            + out * rates["output"]) / 1e6


# --------------------------------------------------------------------------- #
# Where a run writes
# --------------------------------------------------------------------------- #
def default_root():
    """The scratchpad this session was given, or a temp directory beside it."""
    scratch = os.environ.get("XSLOPE_SCRATCHPAD")
    if not scratch:
        scratch = os.path.join(
            "/private/tmp/claude-501/-Users-njones-python-projects-xslope",
            "03b4c37c-d334-4dd8-a131-a66aef1edf42", "scratchpad")
    if not os.path.isdir(scratch):
        import tempfile
        scratch = tempfile.mkdtemp(prefix="assistant_suite_")
    return os.path.join(scratch, "assistant_suite")


def stamped_dir(root):
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    path = os.path.join(root, stamp)
    os.makedirs(path, exist_ok=True)
    return path


# --------------------------------------------------------------------------- #
# One scenario
# --------------------------------------------------------------------------- #
def prepare(scenario, outdir):
    """The workbook the session will open — always a copy, under this run.

    Nothing opens a repository file directly. A session can write beside the
    project it has open (``generate_report`` puts the .docx there), so the project
    it has open is a copy inside the run's own directory; and a scenario with
    faults gets those written into that copy, leaving the sound original for every
    other test that reads it.
    """
    if scenario.model is None:
        return None
    return plant(scenario.model, scenario.faults, os.path.join(outdir, "start"))


def play(scenario, outdir, provider, model, dry_run):
    """Record one scenario, and return the session dict the scorer reads."""
    from tools.assistant_sessions import run_assistant_session

    os.makedirs(outdir, exist_ok=True)
    start = prepare(scenario, outdir)
    result = run_assistant_session(
        scenario.name, start, scenario.turns,
        provider=provider, model=model, images=scenario.images,
        out_dir=os.path.join(outdir, "images"), files_dir=outdir,
        timeout_s=scenario.timeout_s, dry_run=dry_run,
        save_after=scenario.save_after, save_each_turn=scenario.save_each_turn,
        settings=scenario.settings, prefix="sx")
    session = dict(result)
    session["start_model"] = start
    session["lock_model"] = scenario.model      # where a published lock is indexed
    session["scenario"] = scenario.name
    session["family"] = scenario.family
    session["provider"] = provider
    session["model"] = model
    session["dry_run"] = bool(dry_run)
    session["recorded"] = datetime.datetime.now().isoformat(timespec="seconds")
    with open(os.path.join(outdir, "session.json"), "w", encoding="utf-8") as fh:
        json.dump(session, fh, indent=1, default=str)
    return session


def score(scenario, session, outdir):
    """Every criterion, against the finished session. Returns a result dict."""
    ctx = ScoreCtx(scenario, session, outdir)
    rows = []
    for criterion in scenario.criteria:
        started = time.time()
        ok, reason = criterion(ctx)
        rows.append({"name": criterion.name, "kind": criterion.kind,
                     "pass": ok, "reason": reason,
                     "seconds": round(time.time() - started, 1)})
    passed = sum(1 for r in rows if r["pass"])
    return {"scenario": scenario.name, "family": scenario.family,
            "what": scenario.what, "criteria": rows,
            "passed": passed, "total": len(rows),
            "pass": passed == len(rows),
            "usage": dict(session.get("usage") or {}),
            "seconds": session.get("seconds") or 0.0,
            "error": session.get("error"),
            "start_model": session.get("start_model"),
            "faults": [f.describe for f in scenario.faults]}


# --------------------------------------------------------------------------- #
# The scorecard
# --------------------------------------------------------------------------- #
def render(results, meta):
    """``scorecard.md`` — the table first, then every criterion that failed."""
    lines = ["# Assistant scenario suite — scorecard", ""]
    lines.append("- Run: `%s`" % meta.get("mode"))
    lines.append("- Provider / model: `%s` / `%s`" % (meta.get("provider"),
                                                      meta.get("model")))
    lines.append("- Recorded: %s" % meta.get("recorded"))
    lines.append("- Directory: `%s`" % meta.get("outdir"))
    if meta.get("mode") != "replay":
        lines.append("- Spend: **$%.2f** at %s list prices (%s); cache writes are "
                     "priced as plain input, so this is a floor."
                     % (meta.get("spend", 0.0), meta.get("model"), PRICE_DATE))
        if meta.get("budget"):
            lines.append("- Budget: $%.2f%s" % (
                meta["budget"],
                "  — **stopped at the cap**" if meta.get("stopped") else ""))
    passed = sum(1 for r in results if r["pass"])
    lines += ["", "**%d of %d scenarios pass.**" % (passed, len(results)), "",
              "| Scenario | Family | Criteria | Result | Wall | Cost |",
              "| --- | --- | ---: | :--- | ---: | ---: |"]
    for row in results:
        lines.append("| `%s` | %s | %d/%d | %s | %.0f s | $%.3f |" % (
            row["scenario"], row["family"], row["passed"], row["total"],
            "pass" if row["pass"] else "FAIL", row.get("seconds") or 0.0,
            row.get("cost") or 0.0))
    failing = [r for r in results if not r["pass"]]
    if failing:
        lines += ["", "## What failed", ""]
        for row in failing:
            lines += ["### `%s` — %s" % (row["scenario"], row["what"]), ""]
            if row.get("faults"):
                lines.append("Planted: %s." % "; ".join(row["faults"]))
                lines.append("")
            for crit in row["criteria"]:
                if not crit["pass"]:
                    lines.append("- **%s** — %s" % (crit["name"], crit["reason"]))
            lines.append("")
    lines += ["", "## Every criterion", ""]
    for row in results:
        lines += ["### `%s` — %s" % (row["scenario"], row["what"]), ""]
        if row.get("error"):
            lines += ["Session ended on an error: %s" % row["error"], ""]
        for crit in row["criteria"]:
            lines.append("- %s **%s** — %s" % ("PASS" if crit["pass"] else "FAIL",
                                               crit["name"], crit["reason"]))
        lines.append("")
    return "\n".join(lines) + "\n"


def summarize(results, meta):
    """The same thing, short, for the console."""
    out = ["", "=" * 72,
           "%d of %d scenarios pass  (%s)" % (
               sum(1 for r in results if r["pass"]), len(results), meta["mode"])]
    for row in results:
        flag = "pass" if row["pass"] else "FAIL"
        out.append("  %-4s %-24s %2d/%-2d  %5.0fs  $%.3f"
                   % (flag, row["scenario"], row["passed"], row["total"],
                      row.get("seconds") or 0.0, row.get("cost") or 0.0))
        for crit in row["criteria"]:
            if not crit["pass"]:
                out.append("        - %s: %s" % (crit["name"], crit["reason"]))
    if meta.get("mode") != "replay":
        out.append("  total spend $%.2f (%s list prices, %s)"
                   % (meta.get("spend", 0.0), meta.get("model"), PRICE_DATE))
    out.append("=" * 72)
    return "\n".join(out)


# --------------------------------------------------------------------------- #
# The run
# --------------------------------------------------------------------------- #
def select(names=None, families=None):
    chosen = SCENARIOS
    if families:
        wanted = {f.strip() for f in families}
        chosen = [s for s in chosen if s.family in wanted]
    if names:
        wanted = [n.strip() for n in names]
        chosen = [s for s in chosen if s.name in wanted]
        missing = [n for n in wanted if not any(s.name == n for s in chosen)]
        if missing:
            raise SystemExit("no scenario named: %s" % ", ".join(missing))
    return chosen


def run(scenarios, outdir, provider="anthropic", model="claude-opus-5",
        dry_run=False, budget=None, do_score=True, renderer=None,
        meta_extra=None):
    """Play and score ``scenarios``, stopping at ``budget``.

    ``renderer`` writes ``scorecard.md`` — the corpus sweep groups its files by
    input class rather than listing scenarios, and everything else about a run is
    the same. ``meta_extra`` rides into the meta the renderer reads.
    """
    os.makedirs(outdir, exist_ok=True)
    use_solve_cache(os.path.join(os.path.dirname(outdir), "solve_cache.json"))
    results, spend, stopped = [], 0.0, False
    estimated_total = sum(getattr(s, "estimate_usd", 0.0) or 0.0
                          for s in scenarios)
    if estimated_total:
        print("   estimated %.2f USD for %d scenario(s) at %s list prices (%s)"
              % (estimated_total, len(scenarios), model, PRICE_DATE))
    for scenario in scenarios:
        if budget is not None and spend >= budget:
            stopped = True
            print("! budget of $%.2f reached — %s and the rest were not run"
                  % (budget, scenario.name))
            break
        here = os.path.join(outdir, scenario.name)
        try:
            session = play(scenario, here, provider, model, dry_run)
        except Exception as exc:
            # A run is paid for scenario by scenario; one that cannot even be set
            # up must not take the ones after it down with it.
            print("!! %s could not be played: %s: %s"
                  % (scenario.name, type(exc).__name__, exc))
            results.append({"scenario": scenario.name, "family": scenario.family,
                            "what": scenario.what, "criteria": [], "passed": 0,
                            "total": 1, "pass": False, "cost": 0.0,
                            "seconds": 0.0, "faults": [],
                            "error": "%s: %s" % (type(exc).__name__, exc)})
            continue
        cost = 0.0 if dry_run else cost_of(session.get("usage") or {}, model)
        spend += cost
        row = (_scored(scenario, session, here) if do_score
               else {"scenario": scenario.name, "family": scenario.family,
                     "what": scenario.what, "criteria": [], "passed": 0,
                     "total": 0, "pass": True,
                     "seconds": session.get("seconds"),
                     "error": session.get("error"), "faults": []})
        row["cost"] = cost
        row["estimate"] = getattr(scenario, "estimate_usd", None)
        results.append(row)
        print("   %s: %d/%d criteria, $%.3f%s"
              % (scenario.name, row["passed"], row["total"], cost,
                 " (estimated $%.3f; run so far $%.2f of $%.2f estimated)"
                 % (row["estimate"], spend, estimated_total)
                 if row["estimate"] else ""))
    meta = {"mode": "dry run" if dry_run else "live", "provider": provider,
            "model": model, "outdir": outdir, "spend": spend, "budget": budget,
            "stopped": stopped,
            "recorded": datetime.datetime.now().isoformat(timespec="seconds")}
    meta.update(meta_extra or {})
    write(results, meta, outdir, renderer)
    return results, meta


def _scored(scenario, session, here):
    """:func:`score`, with a scorer crash reported rather than raised."""
    try:
        return score(scenario, session, here)
    except Exception as exc:
        return {"scenario": scenario.name, "family": scenario.family,
                "what": scenario.what, "criteria": [], "passed": 0, "total": 1,
                "pass": False, "seconds": session.get("seconds") or 0.0,
                "faults": [], "error": "scoring raised %s: %s"
                                       % (type(exc).__name__, exc)}


def replay(outdir, do_score=True, corpus=False, sampled=False):
    """Re-score sessions already recorded under ``outdir``. No provider, no cost.

    ``corpus`` re-scores a sweep instead of the registry: a sweep's scenarios are
    not written down anywhere, so each one is rebuilt from the workbook its
    session recorded opening — which is the whole of what defines it. ``sampled``
    does the same for a round-two run, where the scenario's name also carries the
    TASK that was asked, and the spec inside that task is seeded on the file and
    the task alone so it rebuilds identically with nothing stored.
    """
    from tools.assistant_scenarios import by_name

    resolve, renderer, cases = by_name, None, []
    if sampled:
        from tools.assistant_scenarios import tasks as menu

        renderer = menu.render

        def resolve(name, session=None):
            row = menu.row_for(name, session)
            if row is not None:
                cases.append(row)
            return menu.resolve(name, session or {})
    elif corpus:
        from tools.assistant_scenarios import corpus as sweep

        renderer = sweep.render

        def resolve(_name, session=None):
            case = sweep.Case(session["lock_model"])
            cases.append(case)
            return sweep.scenario_for(case)

    use_solve_cache(os.path.join(os.path.dirname(outdir.rstrip("/")),
                                 "solve_cache.json"))
    results = []
    for name in sorted(os.listdir(outdir)):
        here = os.path.join(outdir, name)
        record = os.path.join(here, "session.json")
        if not os.path.exists(record):
            continue
        with open(record, encoding="utf-8") as fh:
            session = json.load(fh)
        scenario = (resolve(session.get("scenario") or name, session)
                    if (corpus or sampled)
                    else resolve(session.get("scenario") or name))
        if scenario is None:
            print("! no scenario named %r in the registry — skipped" % name)
            continue
        row = _scored(scenario, session, here)
        row["cost"] = 0.0
        results.append(row)
        print("   %s: %d/%d criteria" % (row["scenario"], row["passed"],
                                         row["total"]))
    meta = {"mode": "replay", "provider": "-", "model": "-", "outdir": outdir,
            "recorded": datetime.datetime.now().isoformat(timespec="seconds")}
    if corpus or sampled:
        meta["cases"] = cases
    if sampled:
        meta.update(_stored_draw(outdir))
    write(results, meta, outdir, renderer)
    return results, meta


# --------------------------------------------------------------------------- #
# The corpus sweep
# --------------------------------------------------------------------------- #
def recorded(outdir):
    """The cases ``outdir`` already holds a finished session for.

    A sweep is stopped by a budget rather than by finishing, so the normal way
    to run the rest of it is to point a second run at the same directory. A
    directory with a ``session.json`` in it is a session that was played and
    paid for; a directory without one is a session that was interrupted, and it
    is run again.
    """
    if not os.path.isdir(outdir):
        return set()
    return {name for name in os.listdir(outdir)
            if os.path.exists(os.path.join(outdir, name, "session.json"))}


def sweep(outdir, provider="anthropic", model="claude-opus-5", dry_run=False,
          budget=None, limit=None, only=None, skip=None, resume=False):
    """Every shipped workbook, one session each, scored by input class."""
    from tools.assistant_scenarios import corpus

    every = corpus.cases()
    if only:
        wanted = {n.strip() for n in only}
        every = [c for c in every if c.name in wanted]
    runnable = [c for c in every if c.loads]
    unloadable = [c.name for c in every if not c.loads]
    done = recorded(outdir) if resume else set()
    passed_over = {n.strip() for n in (skip or [])}
    if done or passed_over:
        runnable = [c for c in runnable
                    if c.name not in done and c.name not in passed_over]
        print("   %d already recorded here, %d named on --corpus-skip; %d left"
              % (len(done), len(passed_over), len(runnable)))
    if limit:
        runnable = runnable[:int(limit)]
    scenarios = [corpus.scenario_for(c) for c in runnable]
    print("corpus sweep: %d workbook(s), %d runnable -> %s"
          % (len(every), len(runnable), outdir))
    if unloadable:
        print("  %d do not load and were not run: %s"
              % (len(unloadable), ", ".join(sorted(unloadable)[:6])))
    return run(scenarios, outdir, provider=provider, model=model,
               dry_run=dry_run, budget=budget, renderer=corpus.render,
               meta_extra={"cases": runnable, "corpus_size": len(every),
                           "unloadable": unloadable})


# --------------------------------------------------------------------------- #
# The sampled sweep — a stratified draw, and a menu of tasks
# --------------------------------------------------------------------------- #
#: The run whose files a second round skips by default: sweeping a file the first
#: round already answered measures the same thing twice. ``--include-swept`` puts
#: them back, and ``--swept-from`` points at a different first round.
FIRST_ROUND = "live1"


def first_round():
    """Where the first live corpus round wrote its sessions."""
    return os.path.join(default_root(), "corpus", FIRST_ROUND)


def sampled_sweep(outdir, provider="anthropic", model="claude-opus-5",
                  dry_run=False, budget=None, n=60, per_file=2, seed=1,
                  swept_from=None, include_swept=False, resume=False,
                  builds=None, min_per_class=2):
    """A stratified sample of the corpus, with a menu of tasks asked of each file.

    The draw is written to the run directory before anything is played, so a run
    stopped at its budget can be finished from the same directory and a replay
    re-scores the same set.
    """
    from tools.assistant_scenarios import corpus
    from tools.assistant_scenarios import tasks as menu

    # ``--builds`` on its own is the build family and nothing else: there is no
    # sample to draw, and drawing one anyway would buy a sweep nobody asked for.
    every = corpus.cases() if n else []
    exclude = set()
    if n and not include_swept:
        source = swept_from or first_round()
        exclude = menu.swept_in(source)
        print("   %d file(s) already swept in %s are excluded"
              % (len(exclude), os.path.basename(source.rstrip("/")) or source))
    chosen, draw = (menu.sample(every, n, seed, exclude=exclude,
                                min_per_class=min_per_class) if n
                    else ([], {"pool": 0, "quota": 0, "asked": 0, "drawn": 0,
                               "min_per_class": min_per_class}))
    pairs = menu.plan(chosen, per_file, seed)
    scenarios, rows, skipped = [], [], []
    for case, task in pairs:
        built = menu.scenario_for(case, task)
        if built is None:
            skipped.append("%s/%s" % (case.name, task.name))
            continue
        built.estimate_usd = estimate_usd(task.cost, model)
        scenarios.append(built)
        rows.append(menu.Row(built.name, case.classes, case.primary, case.kind,
                              case.path, task.name))
    for build in _build_cases(builds, seed, menu):
        built = menu.build_scenario(build)
        built.estimate_usd = estimate_usd(5, model)
        scenarios.append(built)
        rows.append(menu.Row(built.name, build.classes, build.primary, "build",
                              build.path, "build_%s" % build.kind))
    done = recorded(outdir) if resume else set()
    if done:
        keep = [(s, r) for s, r in zip(scenarios, rows) if s.name not in done]
        scenarios = [s for s, _r in keep]
        print("   %d already recorded here; %d left" % (len(done), len(scenarios)))
    draw["seed"] = seed
    draw["tasks_per_file"] = per_file
    print("sampled sweep: %d file(s) drawn from %d eligible (%d filled the class "
          "quota), %d scenario(s) -> %s"
          % (draw["drawn"], draw["pool"], draw["quota"], len(scenarios), outdir))
    if skipped:
        print("   %d (file, task) pair(s) the file could not support: %s"
              % (len(skipped), ", ".join(skipped[:4])))
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "sample.json"), "w", encoding="utf-8") as fh:
        json.dump({"draw": draw, "seed": seed,
                   "files": [c.name for c in chosen],
                   "scenarios": [s.name for s in scenarios]}, fh, indent=1)
    return run(scenarios, outdir, provider=provider, model=model,
               dry_run=dry_run, budget=budget, renderer=menu.render,
               meta_extra={"cases": rows, "corpus_size": len(every),
                           "unloadable": [c.name for c in every if not c.loads],
                           "draw": draw, "seed": seed})


def _build_cases(builds, seed, menu):
    """The build-from-nothing cases this run plays: none, ``all``, or ``N`` drawn."""
    if builds in (None, "", 0, "0"):
        return []
    every = menu.builds()
    if str(builds).lower() in ("all", "-1"):
        return every
    import random

    rng = random.Random("%s|builds" % seed)
    return rng.sample(every, min(int(builds), len(every)))


def _stored_draw(outdir):
    """The draw a sampled run wrote, for a replay's scorecard."""
    path = os.path.join(outdir, "sample.json")
    if not os.path.exists(path):
        return {}
    try:
        with open(path, encoding="utf-8") as fh:
            stored = json.load(fh)
    except Exception:
        return {}
    return {"draw": stored.get("draw"), "seed": stored.get("seed")}


def write(results, meta, outdir, renderer=None):
    with open(os.path.join(outdir, "scorecard.md"), "w", encoding="utf-8") as fh:
        fh.write((renderer or render)(results, meta))
    # ``cases`` carries loaded models the renderer needs and JSON cannot hold.
    stored = {k: v for k, v in meta.items() if k != "cases"}
    with open(os.path.join(outdir, "scorecard.json"), "w", encoding="utf-8") as fh:
        json.dump({"meta": stored, "results": results}, fh, indent=1, default=str)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--dry-run", action="store_true",
                      help="stub replies; exercises everything but the provider")
    mode.add_argument("--live", action="store_true",
                      help="real turns against the provider — this spends money")
    mode.add_argument("--replay", metavar="DIR",
                      help="re-score the sessions already recorded in DIR")
    parser.add_argument("--provider", default="anthropic")
    parser.add_argument("--model", default="claude-opus-5")
    parser.add_argument("--budget-usd", type=float, default=None,
                        help="stop before the scenario that would pass this cap")
    parser.add_argument("--only", default=None,
                        help="comma-separated scenario names")
    parser.add_argument("--family", default=None,
                        help="comma-separated families")
    parser.add_argument("--out", default=None, help="where to write this run")
    parser.add_argument("--corpus", action="store_true",
                        help="sweep every shipped workbook instead of the "
                             "registry: one session and one prompt per file, "
                             "scored by input class")
    parser.add_argument("--corpus-limit", type=int, default=None,
                        help="with --corpus, sweep only the first N files")
    parser.add_argument("--corpus-skip", default=None,
                        help="with --corpus, comma-separated case names not to "
                             "run at all")
    parser.add_argument("--resume", action="store_true",
                        help="with --corpus and --out, skip the files this "
                             "directory already holds a session for — how the "
                             "rest of a budget-stopped sweep is run")
    parser.add_argument("--sample", type=int, default=None,
                        help="with --corpus, draw this many files rather than "
                             "sweeping all of them — stratified, so every input "
                             "class the corpus can fill is filled; a floor, not "
                             "a cap, since a sample too small to be stratified "
                             "comes back larger than asked")
    parser.add_argument("--tasks", type=int, default=2,
                        help="with --sample, how many tasks to ask of each file "
                             "(run_declared is always one of them)")
    parser.add_argument("--seed", type=int, default=1,
                        help="with --sample, the seed the draw is made from")
    parser.add_argument("--min-per-class", type=int, default=2,
                        help="with --sample, how many files each input class is "
                             "represented by where the corpus allows")
    parser.add_argument("--include-swept", action="store_true",
                        help="with --sample, do not skip the files the first "
                             "live round already swept")
    parser.add_argument("--swept-from", default=None,
                        help="with --sample, the finished run whose files are "
                             "skipped (default: the first live corpus round)")
    parser.add_argument("--builds", nargs="?", const="all", default=None,
                        help="with --sample, also play the build-from-nothing "
                             "family: `--builds` for every tutorial drawing and "
                             "description, `--builds N` for N of them")
    parser.add_argument("--no-score", action="store_true",
                        help="record only; score later with --replay")
    parser.add_argument("--list", action="store_true",
                        help="print the scenarios and exit")
    args = parser.parse_args(argv)

    if args.list:
        if args.corpus and args.builds is not None:
            from tools.assistant_scenarios import tasks as menu
            for build in menu.builds():
                print("%-40s %-12s %s" % (build.name, build.kind,
                                          os.path.basename(build.path)))
            return 0
        if args.corpus and args.sample:
            from tools.assistant_scenarios import corpus
            from tools.assistant_scenarios import tasks as menu
            exclude = (set() if args.include_swept
                       else menu.swept_in(args.swept_from or first_round()))
            chosen, draw = menu.sample(corpus.cases(), args.sample, args.seed,
                                       exclude=exclude,
                                       min_per_class=args.min_per_class)
            total = 0
            for case, task in menu.plan(chosen, args.tasks, args.seed):
                total += task.cost
                print("%-52s %-16s %-22s %s"
                      % (case.name, task.name, case.primary,
                         ", ".join(case.classes)))
            print("%d file(s) of %d eligible, %d scenario(s), about %d "
                  "completion(s) — an estimated $%.2f at %s list prices (%s)"
                  % (draw["drawn"], draw["pool"],
                     len(menu.plan(chosen, args.tasks, args.seed)), total,
                     estimate_usd(total, args.model), args.model, PRICE_DATE))
            return 0
        if args.corpus:
            from tools.assistant_scenarios import corpus
            for case in corpus.cases():
                print("%-52s %-5s %-22s %s"
                      % (case.name, case.kind if case.loads else "-",
                         case.primary, ", ".join(case.classes)))
            return 0
        for scenario in SCENARIOS:
            print("%-24s %-12s %d turn(s) %2d criteria  %s"
                  % (scenario.name, scenario.family, len(scenario.turns),
                     len(scenario.criteria), scenario.what))
        return 0

    if args.replay:
        sampled = bool(args.sample or args.builds is not None
                       or os.path.exists(os.path.join(args.replay,
                                                      "sample.json")))
        results, meta = replay(args.replay, do_score=not args.no_score,
                               corpus=args.corpus and not sampled,
                               sampled=args.corpus and sampled)
        print(summarize(results, meta))
        return 0 if all(r["pass"] for r in results) else 1

    if not args.dry_run and not args.live:
        parser.error("choose --dry-run, --live or --replay")

    if args.corpus and (args.sample or args.builds is not None):
        outdir = args.out or stamped_dir(os.path.join(default_root(), "corpus"))
        results, meta = sampled_sweep(
            outdir, provider=args.provider, model=args.model,
            dry_run=args.dry_run, budget=args.budget_usd, n=args.sample or 0,
            per_file=args.tasks, seed=args.seed, swept_from=args.swept_from,
            include_swept=args.include_swept, resume=args.resume,
            builds=args.builds, min_per_class=args.min_per_class)
        print(summarize(results, meta))
        print("scorecard: %s" % os.path.join(outdir, "scorecard.md"))
        return 0 if all(r["pass"] for r in results) else 1

    if args.corpus:
        outdir = args.out or stamped_dir(os.path.join(default_root(), "corpus"))
        results, meta = sweep(outdir, provider=args.provider, model=args.model,
                              dry_run=args.dry_run, budget=args.budget_usd,
                              limit=args.corpus_limit, resume=args.resume,
                              only=args.only.split(",") if args.only else None,
                              skip=(args.corpus_skip.split(",")
                                    if args.corpus_skip else None))
        print(summarize(results, meta))
        print("scorecard: %s" % os.path.join(outdir, "scorecard.md"))
        return 0 if all(r["pass"] for r in results) else 1

    scenarios = select(args.only.split(",") if args.only else None,
                       args.family.split(",") if args.family else None)
    outdir = args.out or stamped_dir(default_root())
    print("assistant suite: %d scenario(s) -> %s" % (len(scenarios), outdir))
    results, meta = run(scenarios, outdir, provider=args.provider,
                        model=args.model, dry_run=args.dry_run,
                        budget=args.budget_usd, do_score=not args.no_score)
    print(summarize(results, meta))
    print("scorecard: %s" % os.path.join(outdir, "scorecard.md"))
    return 0 if all(r["pass"] for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
