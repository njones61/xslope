"""Capability questions, asked of the real assistant and recorded.

The assistant was asked whether xslope supports line loads. It said no, and
offered a short distributed load instead, on an input the template, the loader,
the FEM and six documentation pages all carry. Nothing in the harness could have
caught that: it is not a wrong number or a broken snippet, it is an answer given
from memory because the model had no way to read the documentation of the version
it was running.

So the fix — the shipped documentation index, the ``docs()`` helper, and the rule
in §7 of the brief that says to call it first — is proved the only way it can be:
by asking the questions and reading the answers. Each question here is a
capability question of the shape that failed, paired with the pages that document
the capability. One session each, one turn each, fresh conversation every time, so
no answer is propped up by an earlier turn.

The turns cost real money against the key in the OS keychain, so nothing here runs
by accident. ``--dry-run`` exercises the whole path with a stub reply and no
provider call.

Run::

    python3 tools/assistant_docs_questions.py --dry-run     # plumbing only
    python3 tools/assistant_docs_questions.py               # real turns
    python3 tools/assistant_docs_questions.py line_loads    # just one

The transcripts land in ``test/fixtures/docs_questions/`` and
``test/assistant_docs_questions_check.py`` reads them back offline: every recorded
answer must have called ``docs()`` and must cite one of its question's pages.
"""

from __future__ import annotations

import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

#: Where the recorded transcripts live. Under test/ rather than under docs/,
#: because these are evidence for a check, not figures for a page.
FIXTURES = os.path.join(REPO_ROOT, "test", "fixtures", "docs_questions")

#: Dock grabs are written too (the harness always makes one) but nothing reads
#: them, so they go to a scratch directory rather than into the repository.
GRABS = os.path.join(REPO_ROOT, "test", "fixtures", "docs_questions", "_grabs")

#: The filename prefix every artifact of this set carries.
PREFIX = "docs"

#: The model the questions are asked over: LEM-3's layered slope, the same one
#: the W-1 tutorial records against. The questions are about xslope rather than
#: about this model, but a project has to be open for the dock to be usable, and
#: an ordinary model is the honest setting — the one a person would be looking at
#: when the question occurs to them.
MODEL = os.path.join(REPO_ROOT, "docs/lem/files/xslope_simple_mult_layers.xlsx")

#: Each question, and the pages that document the capability it asks about. The
#: answer must cite one of them: not "a page", and not a page that merely mentions
#: the word — the page a reader sent there would find the feature described on.
QUESTIONS = [
    {
        "name": "line_loads",
        "prompt": "Does xslope support line loads?",
        "expect": ("usage/input_template/", "lem/overview/", "fem/overview/",
                   "lem/oms/", "lem/bishop/", "lem/janbu/", "lem/spencer/",
                   "lem/mprice/", "lem/force_eq/", "lem/reinforcement/",
                   "studio/editing/", "tutorials/lem02_loads_on_the_crest/"),
    },
    {
        "name": "rapid_drawdown",
        "prompt": "Can xslope do rapid drawdown?",
        "expect": ("lem/rapid/", "lem/overview/",
                   "tutorials/combo02_rapid_drawdown/"),
    },
    {
        "name": "tension_crack",
        "prompt": "Does xslope have a tension crack?",
        "expect": ("lem/overview/", "lem/samples/", "usage/input_template/",
                   "usage/geostudio/"),
    },
    {
        "name": "seismic",
        "prompt": "Can it do seismic loading?",
        "expect": ("lem/overview/", "fem/overview/", "usage/input_template/",
                   "lem/oms/", "lem/bishop/", "lem/janbu/", "lem/spencer/",
                   "lem/mprice/", "lem/force_eq/"),
    },
    {
        "name": "piezo_line",
        "prompt": "How do I add a piezometric line?",
        "expect": ("usage/input_template/", "studio/editing/",
                   "seep/seep_slope/", "lem/overview/",
                   "tutorials/lem04_water_in_the_slope/"),
    },
    {
        "name": "anchors",
        "prompt": "Does xslope support anchors or tiebacks?",
        "expect": ("lem/reinforcement/", "fem/reinforcement/",
                   "usage/input_template/", "studio/editing/",
                   "tutorials/lem09_tieback_wall/"),
    },
]

#: Transcript path for one question, the single place the name -> file rule lives.
def transcript_path(name):
    return os.path.join(FIXTURES, "%s_%s_transcript.md" % (PREFIX, name))


def run_one(question, dry_run=False):
    """Ask one question and record the turn. Returns the harness's result dict."""
    from tools.assistant_sessions import run_assistant_session

    return run_assistant_session(
        question["name"], MODEL, [question["prompt"]],
        prefix=PREFIX, out_dir=GRABS, files_dir=FIXTURES,
        # A question changes nothing, so nothing should be written back. Saying
        # so explicitly means a session that DID change the model leaves the
        # workbook behind as evidence rather than being tidied away.
        save_after=None, timeout_s=300, dry_run=dry_run)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    dry_run = "--dry-run" in argv
    argv = [a for a in argv if not a.startswith("--")]
    wanted = [q for q in QUESTIONS if not argv or q["name"] in argv]
    unknown = set(argv) - {q["name"] for q in QUESTIONS}
    if unknown:
        raise SystemExit("unknown question(s): %s — known: %s"
                         % (", ".join(sorted(unknown)),
                            ", ".join(q["name"] for q in QUESTIONS)))
    os.makedirs(FIXTURES, exist_ok=True)
    results = [run_one(q, dry_run=dry_run) for q in wanted]

    print("\nrecorded %d question(s):" % len(results))
    totals = {}
    for res in results:
        usage = res.get("usage") or {}
        for key in ("input", "cached_input", "output"):
            totals[key] = totals.get(key, 0) + int(usage.get(key) or 0)
        print("  %-16s %-70s %s"
              % (res["name"], os.path.relpath(res["transcript"], REPO_ROOT),
                 "ERROR: " + res["error"] if res.get("error") else "ok"))
    print("  tokens: in %(input)s (cached %(cached_input)s), out %(output)s"
          % {k: totals.get(k, 0) for k in ("input", "cached_input", "output")})
    return 1 if any(r.get("error") for r in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
