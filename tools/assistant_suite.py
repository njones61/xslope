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
        dry_run=False, budget=None, do_score=True):
    """Play and score ``scenarios``, stopping at ``budget``."""
    os.makedirs(outdir, exist_ok=True)
    use_solve_cache(os.path.join(os.path.dirname(outdir), "solve_cache.json"))
    results, spend, stopped = [], 0.0, False
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
        results.append(row)
        print("   %s: %d/%d criteria, $%.3f"
              % (scenario.name, row["passed"], row["total"], cost))
    meta = {"mode": "dry run" if dry_run else "live", "provider": provider,
            "model": model, "outdir": outdir, "spend": spend, "budget": budget,
            "stopped": stopped,
            "recorded": datetime.datetime.now().isoformat(timespec="seconds")}
    write(results, meta, outdir)
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


def replay(outdir, do_score=True):
    """Re-score sessions already recorded under ``outdir``. No provider, no cost."""
    from tools.assistant_scenarios import by_name

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
        scenario = by_name(session.get("scenario") or name)
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
    write(results, meta, outdir)
    return results, meta


def write(results, meta, outdir):
    with open(os.path.join(outdir, "scorecard.md"), "w", encoding="utf-8") as fh:
        fh.write(render(results, meta))
    with open(os.path.join(outdir, "scorecard.json"), "w", encoding="utf-8") as fh:
        json.dump({"meta": meta, "results": results}, fh, indent=1, default=str)


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
    parser.add_argument("--no-score", action="store_true",
                        help="record only; score later with --replay")
    parser.add_argument("--list", action="store_true",
                        help="print the scenarios and exit")
    args = parser.parse_args(argv)

    if args.list:
        for scenario in SCENARIOS:
            print("%-24s %-12s %d turn(s) %2d criteria  %s"
                  % (scenario.name, scenario.family, len(scenario.turns),
                     len(scenario.criteria), scenario.what))
        return 0

    if args.replay:
        results, meta = replay(args.replay, do_score=not args.no_score)
        print(summarize(results, meta))
        return 0 if all(r["pass"] for r in results) else 1

    if not args.dry_run and not args.live:
        parser.error("choose --dry-run, --live or --replay")

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
