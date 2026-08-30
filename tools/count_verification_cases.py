#!/usr/bin/env python3
"""Census of xslope's verification corpus, counted from the source of truth.

Every verification comparison the regression suite executes is declared by a
``<!-- test: ... -->`` tag in the documentation itself (``run_tests.py`` builds
its rows by parsing exactly these tags out of ``docs/**/*.md``).  This script
re-derives the published counts from those tags so any claim of the form
"N verification cases" can be recounted from a clean checkout in one command:

    python3 tools/count_verification_cases.py

It reads the tags only; it solves nothing and imports no solver, so it is
deterministic and runs in well under a second.

Counting definitions
--------------------
tags      Verification-bearing tags: a tag that holds at least one published or
          locked NUMBER.  Excludes pass/fail equivalence tags that hold no value
          (``mp_spencer``, ``gsat_pair``, ``roundtrip``) and prose that merely
          quotes the tag syntax.  Infrastructure suite rows (template_sync,
          deps_declared, preflight, updater, the dialog and round-trip guards,
          ...) are registered in ``run_tests.main`` in code, never as document
          tags, so they are outside this census by construction.

cases     Locked comparisons the suite executes, reported three ways.  The
          strictest is the suite ROW -- what ``run_tests.py`` prints one
          pass/fail line for; the same number comes out of
          ``run_tests.parse_test_tags`` run over ``docs/**/*.md``, which is the
          independent way to recount it.  A LEM tag naming seven methods
          is seven cases (each method's factor of safety is checked separately);
          an element-type coverage tag is one case per element type, because
          each leg is an independent mesh + solve; a head-field or FS-vs-time
          list is reported under both conventions --

            field : one solved field / one curve = 1 case
            probe : one case per locked station or march step

values    Individual numbers held: methods x models, list elements, point
          probes, and the scalar side-locks (critical time, minimum FS,
          probability of failure, sensitivity base/low/high).  Element-type legs
          do NOT multiply this: one expected value covers every leg.

models    Distinct input files behind the tags.

The certified verification pages carry a second, independent counter:
``tools/verification_checks/tags.py`` reports "values checked" for the subset of
keys it audits against the printed text.  This script prints that subset as a
cross-check line so the two counters can be reconciled.
"""
from __future__ import annotations

import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DOCS = REPO / "docs"

TAG = re.compile(r"<!--\s*test:\s*(.*?)\s*-->", re.S)

#: Tag types that name one factor of safety per method (``fs_spencer=...``) or a
#: single ``method=``/``expected_fs=`` pair.
LEM_TYPES = {"single_circle", "circular_search", "single_noncirc",
             "noncircular_search", "design_search"}

#: Types whose one expected value is checked once per element type.  The default
#: leg sets mirror ``run_tests.run_seep_elements_test`` / ``run_fem_elements_test``.
ELEMENT_LEGS = {"seep_elements": ("tri3", "tri6", "quad4", "quad8", "quad9"),
                "fem_elements": ("tri6", "quad8", "quad9")}

#: Tags that assert an equivalence rather than a published number (no locked
#: value of their own).  Counted and reported, but excluded from every number.
PASSFAIL_TYPES = {"mp_spencer", "gsat_pair", "roundtrip"}

#: Source page -> verification family.  Families are pages because each page is
#: one source: one vendor manual, one code, or the sample/textbook set.
FAMILIES = [
    ("docs/verification/rocscience.md", "Rocscience Slide2 (VP manual, LEM)"),
    ("docs/verification/rocscience_groundwater.md",
     "Rocscience Slide2 (VP manual, groundwater)"),
    ("docs/verification/rs2.md", "Rocscience RS2 (FEM/SSRM)"),
    ("docs/verification/geostudio.md", "GeoStudio SLOPE/W + SEEP/W"),
    ("docs/verification/ssrm.md", "Analytical / classical literature (FEM)"),
    ("docs/verification/seep.md", "Analytical / classical literature (seepage)"),
    ("docs/verification/published.md", "Published worked problems (design manuals)"),
]
SAMPLES_FAMILY = "Samples, design and parametric pages"

#: The certified verification pages, and the tag keys their independent
#: text-vs-tag auditor reads (``tools/verification_checks/config.py``).  Counting
#: this subset here reproduces that auditor's own "values checked" total, so the
#: two counters can be reconciled without running the solver.
CERTIFIED_PAGES = ["geostudio", "published", "rocscience",
                   "rocscience_groundwater", "rs2", "seep", "ssrm"]
CHECKER_KEYS = ["expected_fs*", "fs_*", "expected_beta", "expected_kc",
                "expected_flowrate*", "expected_head*", "expected_pullout",
                "expected_envelope", "points", "expected"]
CHECKER_TOTAL = 803     # as reported by tools/verification_checks/tags.py


def parse_tags(md: Path):
    """The tags of one markdown file, as dicts, with the file path resolved.

    Mirrors ``run_tests.parse_test_tags`` (same regex, same comma/equals split,
    same relative-path resolution) minus the numeric coercion, which no count
    needs.  Tags with no ``file=`` key are prose examples of the syntax and are
    dropped, exactly as the suite drops them.
    """
    out = []
    for m in TAG.finditer(md.read_text()):
        kv = {}
        for pair in m.group(1).split(","):
            k, _, v = pair.partition("=")
            kv[k.strip()] = v.strip()
        if not kv.get("file"):
            continue
        kv["_file"] = str((md.parent / kv["file"]).resolve())
        kv["_src"] = str(md.relative_to(REPO))
        kv["_type"] = kv.get("type", "single_circle")
        out.append(kv)
    return out


def n_methods(kv):
    """Locked factors of safety in one LEM tag: one per ``fs_<method>`` key, or
    one for the classic ``method=`` + ``expected_fs=`` form."""
    fs = [k for k in kv if k.startswith("fs_")]
    return len(fs) if fs else (1 if "expected_fs" in kv else 0)


def n_legs(kv):
    """Element-type legs: independent mesh + solve runs the tag checks."""
    default = ELEMENT_LEGS.get(kv["_type"])
    if default is None:
        return 1
    override = kv.get("element_types")
    return len([s for s in override.split(";") if s.strip()]) if override else len(default)


def tally(kv):
    """(values, comparisons_field, comparisons_probe) for one tag.

    ``values``   numbers the tag holds.
    ``field``    comparisons counting a list lock as one field / one curve.
    ``probe``    comparisons counting every list element separately.
    Element legs multiply the comparisons, never the values.
    """
    t = kv["_type"]
    if t in PASSFAIL_TYPES:
        return 0, 0, 0
    if t in LEM_TYPES:
        n = n_methods(kv)
        return n, n, n
    if t in ("fem_ssrm", "fem_elements", "fem_reliability", "seep",
             "seep_elements", "critical_kc", "reliability"):
        return 1, 1, 1
    if t in ("seep_head", "tseep_head"):
        n = len([p for p in kv.get("points", "").split(";") if p.strip()])
        return n, 1, n
    if t == "fs_vs_time":
        n = len([v for v in kv.get("expected", "").split(";") if v.strip()])
        n += sum(1 for k in ("critical_time", "min_fs") if k in kv)
        return n, 1, n
    if t == "pullout_envelope":
        # A pullout table: one element per reinforcement layer in each list, and
        # one envelope evaluation per element.  The two lists are two quantities
        # read at the same stations, so each is one "field" comparison.
        keys = [k for k in ("expected_pullout", "expected_envelope") if k in kv]
        n = sum(len([e for e in kv[k].split(";") if e.strip()]) for k in keys)
        return n, len(keys), n
    if t == "reliability_mc":
        n = sum(1 for k in ("expected_beta", "expected_pf") if k in kv)
        return n, n, n
    if t == "sensitivity":
        n = sum(1 for k in ("expected_base", "expected_low", "expected_high")
                if k in kv)
        return n, n, n
    raise SystemExit(f"unknown tag type {t!r} in {kv['_src']} -- "
                     "add it to tools/count_verification_cases.py")


def checker_values(kv):
    """Values of one tag under the certified-page auditor's key subset."""
    n = 0
    for k, v in kv.items():
        if k.startswith("_"):
            continue
        hit = any(k.startswith(p[:-1]) if p.endswith("*") else k == p
                  for p in CHECKER_KEYS)
        if hit:
            n += len([e for e in v.split(";") if e.strip()])
    return n


def family_of(src):
    for page, name in FAMILIES:
        if src == page:
            return name
    return SAMPLES_FAMILY


def main():
    tags = []
    for md in sorted(DOCS.rglob("*.md")):
        tags.extend(parse_tags(md))

    passfail = [t for t in tags if t["_type"] in PASSFAIL_TYPES]
    bearing = [t for t in tags if t["_type"] not in PASSFAIL_TYPES]

    tot_v = tot_f = tot_p = tot_rows = 0
    fam = defaultdict(lambda: [0, 0, 0, 0, set()])   # tags, field, probe, values, files
    by_type = Counter()
    for kv in bearing:
        v, f, p = tally(kv)
        legs = n_legs(kv)
        f, p = f * legs, p * legs
        tot_v += v
        tot_f += f
        tot_p += p
        # One suite ROW: what run_tests reports pass/fail on.  A LEM tag becomes
        # one row per method; every other tag is a single row, whatever it holds.
        tot_rows += n_methods(kv) if kv["_type"] in LEM_TYPES else 1
        row = fam[family_of(kv["_src"])]
        row[0] += 1
        row[1] += f
        row[2] += p
        row[3] += v
        row[4].add(kv["_file"])
        by_type[kv["_type"]] += 1

    models = {t["_file"] for t in bearing}
    benchmarks = {t["benchmark"] for t in bearing if t.get("benchmark")}

    print("xslope verification census -- counted from the <!-- test: ... --> tags")
    print(f"source: {DOCS.relative_to(REPO)}/**/*.md in {REPO}")
    print()
    print("  (a) verification-bearing tags          %5d" % len(bearing))
    print("      pass/fail equivalence tags (no value, excluded)  %d"
          % len(passfail))
    print("  (b) CASES -- locked comparisons executed")
    print("        suite rows (strictest: one reported")
    print("        pass/fail each)                  %5d" % tot_rows)
    print("        field convention (a solved field or an FS-vs-time")
    print("        curve counts once)               %5d" % tot_f)
    print("        probe convention (one per locked")
    print("        station / march step)            %5d" % tot_p)
    print("  (c) locked VALUES held                 %5d" % tot_v)
    print("  (d) distinct models (input files)      %5d" % len(models))
    print("      distinct benchmark ids             %5d" % len(benchmarks))
    print()
    print("By family (cases, both conventions):")
    print("  %-46s %6s %6s %7s %7s %7s" %
          ("family", "tags", "field", "probe", "values", "models"))
    for name in [n for _, n in FAMILIES] + [SAMPLES_FAMILY]:
        if name not in fam:
            continue
        t, f, p, v, files = fam[name]
        print("  %-46s %6d %6d %7d %7d %7d" % (name, t, f, p, v, len(files)))
    print("  %-46s %6d %6d %7d %7d %7d" %
          ("TOTAL", len(bearing), tot_f, tot_p, tot_v, len(models)))
    print()
    print("By tag type:")
    for t, n in sorted(by_type.items(), key=lambda kv: -kv[1]):
        print("  %-24s %5d" % (t, n))
    print()
    cert_pages = {f"docs/verification/{p}.md" for p in CERTIFIED_PAGES}
    cert = [t for t in bearing if t["_src"] in cert_pages]
    cert_all = sum(tally(t)[0] for t in cert)
    cert_sub = sum(checker_values(t) for t in cert)
    print("Cross-check against the independent text-vs-tag auditor")
    print("(tools/verification_checks/tags.py, the six certified pages):")
    print("  values held on those pages                    %5d" % cert_all)
    print("  of them, under the auditor's key subset       %5d  (it reports %d)"
          % (cert_sub, CHECKER_TOTAL))
    print("  difference = keys the auditor does not read: expected_pf,")
    print("  expected_base/low/high, critical_time, min_fs %5d"
          % (cert_all - cert_sub))
    print("  recount it with: for p in %s; do \\" % " ".join(CERTIFIED_PAGES))
    print("      python3 -m tools.verification_checks.tags"
          " docs/verification/$p.md; done")
    return 0 if cert_sub == CHECKER_TOTAL else 1


if __name__ == "__main__":
    sys.exit(main())
