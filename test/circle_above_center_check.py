"""A circle that meets the ground ABOVE its own center.

A failure surface is built as the bottom semicircle (``generate_failure_surface``
sweeps theta from pi to 2pi), so the two daylight points have to bound an arc no
longer than half the circle. ``circle_polyline_intersections`` enforces that by
keeping only crossings below the center: a crossing above it bounds a MAJOR arc,
whose end slices would overhang — a base steeper than vertical, which no method
here can carry. Dropping the point is the right answer, and removing the filter
costs vp023 bishop 1.130 -> 0.820 and five other locked searches.

What was missing was the WORDING. A circle drawn plainly across the section came
back as "Expected at least 2 intersection points, but got 1", which reads as a
circle that misses the ground and sends the reader to the ground surface. A search
handed such a seed said nothing at all: it scored the seed's grid instead and
reported an answer for a surface the user never specified. And nine shipped sample
circles were dead this way without anything saying so.

What is pinned here:

  1. the slicer NAMES the cause — the message keeps its original opening (a
     tutorial quotes that sentence) and adds the offending crossing, the center it
     sits above, and the remedy;
  2. a genuinely one-crossing circle keeps the plain message, with no cause
     appended: the diagnosis must not fire on the case it is not about;
  3. the preflight rule's severity split — ERROR where the circle IS the run's
     surface, WARNING for a search seed, and silent on a sound circle;
  4. a stated tension crack silences the rule, because the crack resolves the
     uphill end and the arc exits at it (vp030a/b solve this way);
  5. the search discloses a seed it could not build, and still answers;
  6. every starting circle in the files that carried this defect now builds;
  7. a mutation: with the crossings-above-center helper stubbed empty, leg 1 loses
     its cause and leg 3 goes silent — so both read the real rule.

Run directly:  PYTHONPATH=. python3 test/circle_above_center_check.py
"""

import contextlib
import copy
import io
import warnings

from xslope import slice as sl
from xslope.fileio import load_slope_data
from xslope.preflight import ERROR, WARNING, preflight
from xslope.search import circular_search
from xslope.slice import crossings_above_center, generate_slices

warnings.filterwarnings("ignore")

#: The levee section vp020 is built on: a left-facing slope whose crest sits at
#: y = 70, so a center at y = 60 is below the ground it is drawn under.
LEVEE = "docs/verification/files/rocscience/vp020.xlsx"

#: The circle vp020 shipped with, and the crossing that killed it.
BURIED = {"Xo": 90.0, "Yo": 60.0, "Depth": 15.0, "R": 45.0}
BURIED_CROSSING = (134.88, 63.25)

#: The circle it carries now: the same Xo and the same tangent elevation, with the
#: center raised clear of the ground so both crossings sit below it.
SOUND = {"Xo": 90.0, "Yo": 100.0, "Depth": 15.0, "R": 85.0}

#: A circle small enough and far enough left that it crosses the ground once. The
#: control for leg 2: one crossing, and nothing to do with the equator filter.
TOO_SMALL = {"Xo": 20.0, "Yo": 60.0, "Depth": 55.0, "R": 5.0}

#: The opening sentence a tutorial quotes verbatim (lem05). It must survive.
PLAIN = "Expected at least 2 intersection points, but got 1."

RULE = "surface.circle_daylights_above_center"

#: The files whose shipped starting circle was dead this way.
REPAIRED = [
    "docs/fem/files/xslope_griffiths3_r1.xlsx",
    "docs/fem/files/xslope_griffiths3_r0p8.xlsx",
    "docs/fem/files/xslope_griffiths3_r0p6.xlsx",
    "docs/fem/files/xslope_griffiths3_r0p5.xlsx",
    "docs/fem/files/xslope_griffiths3_r0p4.xlsx",
    "docs/fem/files/xslope_griffiths3_r0p2.xlsx",
    "docs/fem/files/xslope_griffiths3_r0p2_thin.xlsx",
    "docs/fem/files/xslope_griffiths4_r1.xlsx",
    "docs/fem/files/xslope_griffiths4_r2.xlsx",
    LEVEE,
]


def _slice_message(sd, circle):
    """The refusal text ``generate_slices`` returns for one circle."""
    ok, res = generate_slices(sd, circle=dict(circle), num_slices=50,
                              check_inputs=False)
    return ok, str(res)


def _findings(sd, selection):
    return [f for f in preflight(sd, "lem", selection).findings if f.rule_id == RULE]


def check_slicer_names_the_cause(sd):
    failures = []
    ok, msg = _slice_message(sd, BURIED)
    if ok:
        return ["the buried circle sliced; it should be refused"]
    if PLAIN not in msg:
        failures.append(f"the quoted opening sentence is gone: {msg[:120]!r}")
    for want in ("above its own center", "Lower the center or enlarge the radius"):
        if want not in msg:
            failures.append(f"refusal does not say {want!r}: {msg[:160]!r}")
    x, y = BURIED_CROSSING
    if f"{x:.4g}" not in msg or f"{y:.4g}" not in msg:
        failures.append(f"refusal does not quote the crossing {BURIED_CROSSING}: "
                        f"{msg[:160]!r}")
    if "Xo=90" not in msg:
        failures.append(f"refusal does not quote the circle: {msg[:160]!r}")
    return failures


def check_plain_case_unchanged(sd):
    ok, msg = _slice_message(sd, TOO_SMALL)
    if ok:
        return ["the one-crossing circle sliced; it should be refused"]
    if "above its own center" in msg:
        return [f"the cause was appended to a circle it is not about: {msg[:160]!r}"]
    if "intersection points" not in msg:
        return [f"unexpected refusal for the one-crossing circle: {msg[:160]!r}"]
    return []


def check_severity_split(sd):
    failures = []
    buried = copy.deepcopy(sd)
    buried["circles"] = [dict(BURIED)]

    single = _findings(buried, {"surface": "circular"})
    if len(single) != 1:
        failures.append(f"single-surface run: expected 1 finding, got {len(single)}")
    elif single[0].severity != ERROR:
        failures.append(f"single-surface run: severity is {single[0].severity}, "
                        f"expected {ERROR}")
    elif "no other" not in single[0].message:
        failures.append("single-surface finding does not say the run has no other "
                        "circle")

    search = _findings(buried, {"surface": "circular", "search": True})
    if len(search) != 1:
        failures.append(f"search run: expected 1 finding, got {len(search)}")
    elif search[0].severity != WARNING:
        failures.append(f"search run: severity is {search[0].severity}, expected "
                        f"{WARNING}")
    elif "grid around it" not in search[0].message:
        failures.append("search finding does not say where the answer comes from")

    sound = copy.deepcopy(sd)
    sound["circles"] = [dict(SOUND)]
    for sel in ({"surface": "circular"}, {"surface": "circular", "search": True}):
        if _findings(sound, sel):
            failures.append(f"fired on the sound circle under {sel}")
    return failures


def check_tension_crack_silences(sd):
    cracked = copy.deepcopy(sd)
    cracked["circles"] = [dict(BURIED)]
    # Deep enough to put the uphill exit at or below the equator: the crossing sits
    # 3.25 above the center, so any crack past that recovers the end.
    cracked["tcrack_depth"] = 10.0
    if _findings(cracked, {"surface": "circular"}):
        return ["fired on a model that states a tension crack"]
    return []


def check_search_discloses(sd):
    buried = copy.deepcopy(sd)
    buried["circles"] = [dict(BURIED)]
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        fs_cache, _conv, _path, _cc = circular_search(buried, "bishop",
                                                      num_slices=50, max_iter=1)
    text = out.getvalue()
    failures = []
    for want in ("cannot be built", "crossing above center",
                 "searching from the launch grid around it"):
        if want not in text:
            failures.append(f"the search did not say {want!r}")
    if not fs_cache or fs_cache[0]["FS"] >= 9999:
        failures.append("the search reported no answer from the launch grid")

    sound = copy.deepcopy(sd)
    sound["circles"] = [dict(SOUND)]
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        circular_search(sound, "bishop", num_slices=50, max_iter=1)
    if "cannot be built" in out.getvalue():
        failures.append("the search announced a dead seed on a sound circle")
    return failures


def check_repaired_files_build():
    failures = []
    for path in REPAIRED:
        sd = load_slope_data(path)
        ground = sd["ground_surface"]
        for i, c in enumerate(sd.get("circles") or []):
            above = crossings_above_center(c["Xo"], c["Yo"], c["R"], ground)
            if above:
                p = max(above, key=lambda q: q.y)
                failures.append(f"{path}: circle {i + 1} still meets the ground at "
                                f"({p.x:.4g}, {p.y:.4g}), above its center "
                                f"y = {c['Yo']:g}")
            ok, _res = generate_slices(sd, circle=dict(c), num_slices=40,
                                       check_inputs=False)
            if not ok:
                failures.append(f"{path}: circle {i + 1} does not slice")
    return failures


def check_mutation(sd):
    """Stub the helper empty: the cause must vanish from both readers."""
    failures = []
    original = sl.crossings_above_center
    sl.crossings_above_center = lambda *a, **k: []
    try:
        _ok, msg = _slice_message(sd, BURIED)
        if "above its own center" in msg:
            failures.append("the slicer still named the cause with the helper "
                            "stubbed empty — it is not reading it")
        buried = copy.deepcopy(sd)
        buried["circles"] = [dict(BURIED)]
        if _findings(buried, {"surface": "circular"}):
            failures.append("the preflight rule still fired with the helper "
                            "stubbed empty — it is not reading it")
    finally:
        sl.crossings_above_center = original
    return failures


def run():
    sd = load_slope_data(LEVEE)
    failures = []

    print("  the slicer names the cause:")
    failures += check_slicer_names_the_cause(sd)

    print("  a genuinely one-crossing circle keeps the plain message:")
    failures += check_plain_case_unchanged(sd)

    print("  preflight severity: error on the run's own surface, warning on a seed:")
    failures += check_severity_split(sd)

    print("  a stated tension crack silences the rule:")
    failures += check_tension_crack_silences(sd)

    print("  the search discloses a seed it could not build:")
    failures += check_search_discloses(sd)

    print("  every repaired file's starting circle builds:")
    failures += check_repaired_files_build()

    print("  mutation — the cause vanishes with the helper stubbed empty:")
    failures += check_mutation(sd)

    return failures


def main():
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll above-center circle checks passed.")


if __name__ == "__main__":
    main()
