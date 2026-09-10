"""The recorded capability answers: the search was used, and the citation is real.

``tools/assistant_docs_questions.py`` asks the live assistant six capability
questions — the shape of question it got wrong, when it denied line loads — and
records each turn. This check reads those transcripts back, offline and free, and
holds each answer to what the fix promises:

* ``docs()`` was actually called. The rule in §7 of the brief is "call it FIRST",
  and an answer that skipped it is an answer from memory whether or not it
  happened to be right this time.
* the answer cites one of the pages that document the capability, and cites it as
  a real address on the documentation site.
* nothing was denied: no recorded answer says xslope cannot do the thing it was
  asked about.

A transcript is evidence of one live run, so this check does not prove today's
model would answer the same way — it proves what was measured, and it fails the
moment a recorded answer stops meeting the bar (a re-recording that regressed, or
a fixture edited to say something the run did not).

Re-record with::

    python3 tools/assistant_docs_questions.py            # all six, billed
    python3 tools/assistant_docs_questions.py line_loads # just one
"""

from __future__ import annotations

import os
import re
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

#: Phrases that deny a capability. Matched against the assistant's own prose only
#: (never the snippets or their output, where "not supported" may be quoted from a
#: page about something else).
_DENIALS = (
    "does not support", "doesn't support", "not supported",
    "xslope cannot", "xslope can't", "there is no support",
    "is not available in xslope", "no built-in",
)


def _split_turn(text):
    """(prose, code) for one transcript — what the assistant SAID and what it RAN.

    The harness writes a turn as a single fenced block of labelled lines: 'You:',
    'Assistant:', 'Ran code:' with the snippet indented under it, and 'Output:'
    likewise. The two are separated here because the rules differ: a citation is
    only a citation when the assistant wrote it to the user, and a `docs(` call is
    only a call when it is in a snippet that ran.
    """
    prose, code = [], []
    bucket = None
    for line in text.splitlines():
        if line.startswith("Assistant:"):
            bucket = prose
            prose.append(line[len("Assistant:"):])
            continue
        if line.startswith("Ran code:") or line.startswith("Output:"):
            bucket = code
            continue
        # The transcript's OWN headings, not the assistant's: the answers use
        # markdown headings themselves, and treating every '## ' as a transcript
        # heading threw away everything an answer said after its first one.
        if (line.startswith("You:") or line.startswith("Tokens:")
                or line.startswith("## Turn ") or line.startswith("## Session total")):
            bucket = None
            continue
        if bucket is not None:
            bucket.append(line)
    return "\n".join(prose), "\n".join(code)


def run():
    """Failures as a list of strings; empty when every recorded answer holds."""
    from studio.ai.kernel import _docs_search
    from tools import assistant_docs_questions as adq

    search = _docs_search()
    if search is None:
        return ["the shipped documentation index could not be loaded, so a "
                "citation cannot be checked against the real page list"]
    #: Every page the documentation has, as the index writes it (`lem/rapid/`).
    pages = {doc["url"] for doc in search.docs}

    out = []
    missing = [q["name"] for q in adq.QUESTIONS
               if not os.path.exists(adq.transcript_path(q["name"]))]
    if len(missing) == len(adq.QUESTIONS):
        return ["no recorded capability answers under test/fixtures/docs_questions "
                "— record them with: python3 tools/assistant_docs_questions.py"]
    for name in missing:
        out.append(f"no recorded answer for {name!r} — re-record it with: "
                   f"python3 tools/assistant_docs_questions.py {name}")

    for question in adq.QUESTIONS:
        name = question["name"]
        path = adq.transcript_path(name)
        if not os.path.exists(path):
            continue
        text = open(path, encoding="utf-8").read()
        if "(dry run — stub reply)" in text:
            out.append(f"{name}: the recorded transcript is a --dry-run stub, not "
                       "a live answer")
            continue
        prose, code = _split_turn(text)
        if not prose.strip():
            out.append(f"{name}: the transcript carries no answer")
            continue

        if not re.search(r"\bdocs\s*\(", code):
            out.append(f"{name}: the answer never called docs() — it was answered "
                       "from memory, which is the failure this path exists to stop")

        # A citation counts in either form the model has been taught: the full
        # address docs() hands back, or the site-relative path the brief's page
        # table is written in (`lem/rapid/`). Both name the same page, and the
        # recorded answers use both. Matched against the shipped index's own page
        # list, so a path that is not a real page is not read as a citation.
        cited = sorted(u for u in pages if u and u in prose)
        if not cited:
            out.append(f"{name}: the answer cites no documentation page")
        elif not any(c.startswith(e) for c in cited for e in question["expect"]):
            out.append(f"{name}: cites {cited[:3]}, none of which is "
                       f"one of {list(question['expect'][:3])}…")

        low = prose.lower()
        for phrase in _DENIALS:
            if phrase in low:
                out.append(f"{name}: the answer says {phrase!r} about a shipped "
                           "capability")
                break
    return out


def main():
    print("recorded capability answers:")
    failures = run()
    for f in failures:
        print("  - " + f)
    if failures:
        raise SystemExit(1)
    print("  every recorded answer searched the documentation and cited a real page.")


if __name__ == "__main__":
    main()
