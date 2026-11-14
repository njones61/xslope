# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""A rapid drawdown surface that crosses no material with d and ψ is analyzed,
not refused.

Until 2026-11 ``advanced.rapid_drawdown`` returned a failure for any surface
whose slices all had d = ψ = 0 ("Rapid drawdown requires d and psi parameters
for low-K materials. All values are zero"), so a search scored every such
circle as an unsolved trial and threw it away. Such a circle is a real surface:
every material it cuts drains during the drawdown, so its answer is the drained
analysis, c′ and φ′ everywhere, under the water before and after drawdown. On a
dam whose free-draining shell is weak, it is the critical surface, and the
search used to miss it.

  A. THE SEARCH FINDS THE SHELL CIRCLE. On the weak-shell input (COMBO-2's
     starter with the shell at c′ = 25 psf, φ′ = 28°, built by
     ``tools/build_johnson_rapid.py weak_shell``) the rapid drawdown search
     from the shell's starting circle settles at F = 1.144 on a circle that
     crosses none of the core: no slice is low-permeability, stage 3 did not
     run, and F is stage 2's. Under the old refusal that circle could not be
     scored at all.
  B. THE RESULT SAYS SO, EVERYWHERE IT IS SHOWN. The Log prints one plain
     sentence after the stage line; the Analysis Report's stage table gives
     c′, φ′ in all materials for stage 2 and the section says the circle crosses
     no material with d and ψ; Studio's LEM Solution tab shows the sentence above
     the plot and the status bar carries its first clause.
  C. A SINGLE CIRCLE THAT MISSES THE CORE SOLVES. On the unaltered starter, a
     circle in the upstream shell is solved as a single surface, with the
     Log's sentence naming "the circle".
  D. THE FULL-POOL RULE STILL HOLDS FOR THESE CIRCLES. A shell circle whose
     stage 1 (full pool) factor of safety is below 1 is refused with the same
     message as any other surface, so the three-stage procedure behaves the
     same whether or not the surface reaches a d/ψ material.
  E. A SURFACE THROUGH THE CORE IS UNCHANGED. The tutorial's critical circle on
     the starter still reads stage 3 not required at 1.196, and records its 13
     core slices as low-permeability.

The search in A is the only one here (about 35 s); everything else is a
single-surface solve.
"""

import contextlib
import io
import os
import sys
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

WEAK_SHELL_XLSX = os.path.join(_REPO, "test", "fixtures",
                               "xslope_johnson_rapid_weak_shell.xlsx")
STARTER_XLSX = os.path.join(_REPO, "docs", "tutorials", "files",
                            "xslope_johnson_rapid_start.xlsx")

#: What the weak-shell search settles at, and how far it may move.
WEAK_SHELL_FS, FS_TOL = 1.144, 0.005

#: A circle in the starter's upstream shell, clear of the core: tangent at
#: elevation 120, 15 ft below the face at x = 270.
SHELL_CIRCLE = {"Xo": 270.0, "Yo": 220.0, "Depth": 120.0, "R": 100.0}

#: The tutorial's critical circle through the core (COMBO-2 Part 1).
CORE_CIRCLE = {"Xo": 246.3773193359375, "Yo": 242.11181640625,
               "Depth": 82.05196142196655, "R": 160.05985498428345}

#: The sentence, as the Log, the report and Studio print it.
SENTENCE_TAIL = ("crosses no material with d and ψ, so stages 1 and 2 used "
                 "drained strengths (c′, φ′) and stage 3 was not required.")


def _load(path):
    from xslope.fileio import load_slope_data
    with contextlib.redirect_stdout(io.StringIO()):
        return load_slope_data(path)


def _run(slope_data, **kw):
    """``(bundle, log)`` of one run through the function every caller uses."""
    from xslope.search import run_lem_analysis
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        bundle = run_lem_analysis(slope_data, "spencer", num_slices=40,
                                  rapid=True, **kw)
    return bundle, out.getvalue()


def _drained(results, where):
    """The failures that say a result is not the all-drained analysis."""
    fails = []
    low_k = results.get("rapid_low_k")
    if low_k is None or any(low_k):
        fails.append(f"{where}: {sum(low_k or [])} slices are recorded as "
                     f"low-permeability on a surface that misses the core")
    if results.get("stage3_run"):
        fails.append(f"{where}: stage 3 ran on a surface with no d/ψ slice")
    if results.get("FS") != results.get("stage2_FS"):
        fails.append(f"{where}: F = {results.get('FS')} is not stage 2's "
                     f"{results.get('stage2_FS')}")
    return fails


def _report_section(slope_data, bundle):
    """The Rapid Drawdown section the report builds for this bundle."""
    from xslope.report import build_report
    opts = {"title": "rapid drained surface", "method": "spencer",
            "pd_figure": False, "lem_inputs_figure": False,
            "lem_search_figure": False, "lem_solution_figure": False,
            "lem_calculations": False}
    with contextlib.redirect_stdout(io.StringIO()):
        report = build_report(slope_data, {"lem": [bundle]}, opts,
                              tempfile.mkdtemp(prefix="xslope_rapid_"))
    for top in report.sections:
        for _lvl, sec in top.walk():
            if sec.title == "Rapid Drawdown":
                return sec
    return None


def _studio(bundle, where):
    """The note the LEM Solution tab shows and the status-bar line."""
    try:
        from PySide6.QtWidgets import QApplication
    except ImportError:
        print("  PySide6 not installed: Studio legs skipped")
        return []
    QApplication.instance() or QApplication([])
    from studio.main_window import (SolutionView, _drained_rapid_note,
                                    lem_done_message)
    fails = []
    note = _drained_rapid_note(bundle)
    if not note.endswith(SENTENCE_TAIL):
        fails.append(f"{where}: Studio's note reads {note!r}")
    view = SolutionView()
    view.set_note("Rapid drawdown", note)
    if view._note.isHidden() or SENTENCE_TAIL[:30] not in view._note.text():
        fails.append(f"{where}: the LEM Solution tab does not show the note "
                     f"({view._note.text()!r})")
    view.set_note("Rapid drawdown", "")
    if not view._note.isHidden():
        fails.append(f"{where}: the note stays up with nothing to say")
    status = lem_done_message(bundle)
    if "crosses no material with d and ψ" not in status:
        fails.append(f"{where}: the status bar reads {status!r}")
    return fails


def run():
    fails = []
    from xslope.advanced import NO_LOW_K_SENTENCE, crosses_no_low_k

    # --- A and B: the weak-shell search -----------------------------------
    sd = _load(WEAK_SHELL_XLSX)
    shell = next(m for m in sd["materials"] if m["name"] == "shell")
    if (shell["c"], shell["phi"]) != (25.0, 28.0) or shell.get("d") or shell.get("psi"):
        fails.append(f"the weak-shell input's shell is c={shell['c']}, "
                     f"phi={shell['phi']}, d={shell.get('d')}, psi={shell.get('psi')}; "
                     f"rebuild it with tools/build_johnson_rapid.py weak_shell")
    bundle, log = _run(sd, analysis="auto_search")
    res = bundle["results"]
    where = "weak-shell search"
    if abs(res["FS"] - WEAK_SHELL_FS) > FS_TOL:
        fails.append(f"{where}: F = {res['FS']:.4f}, expected {WEAK_SHELL_FS} "
                     f"± {FS_TOL}")
    fails += _drained(res, where)
    core = next(i for i, m in enumerate(sd["materials"], 1) if m["name"] == "core")
    if any(int(m) == core for m in bundle["slice_df"]["mat"]):
        fails.append(f"{where}: the critical circle cuts the core")
    if not crosses_no_low_k(res):
        fails.append(f"{where}: crosses_no_low_k is false on the critical circle")
    want = NO_LOW_K_SENTENCE.format(surface="critical circle")
    lines = log.splitlines()
    stage = [i for i, l in enumerate(lines)
             if l.startswith("Rapid drawdown stages on the critical surface")]
    if not stage or stage[0] + 1 >= len(lines) or lines[stage[0] + 1] != want:
        fails.append(f"{where}: the Log does not follow the stage line with "
                     f"{want!r}")
    if "Stage 3 FS = not required" not in log:
        fails.append(f"{where}: the Log's stage line does not say stage 3 was "
                     f"not required")

    sec = _report_section(sd, bundle)
    if sec is None:
        fails.append(f"{where}: the report has no Rapid Drawdown section")
    else:
        table = next((b for b in sec.blocks if b.kind == "table"), None)
        prose = [b.text for b in sec.blocks if b.kind == "prose"]
        if table is None or table.rows[1][2] != "c′, φ′ in all materials":
            fails.append(f"{where}: the report's stage 2 strengths read "
                         f"{table.rows[1][2] if table else None!r}")
        elif table.rows[2][2:] != ["no slice is in a material with d and ψ",
                                   "not required"]:
            fails.append(f"{where}: the report's stage 3 row reads "
                         f"{table.rows[2][2:]!r}")
        if NO_LOW_K_SENTENCE.format(surface="critical circle") not in prose:
            fails.append(f"{where}: the report does not say the critical circle "
                         f"crosses no material with d and ψ")
    fails += _studio(bundle, where)

    # --- C: a single circle that misses the core --------------------------
    sd = _load(STARTER_XLSX)
    sd["circles"] = [dict(SHELL_CIRCLE)]
    where = "single shell circle"
    try:
        bundle, log = _run(sd, analysis="single_surface")
    except Exception as exc:                                  # noqa: BLE001
        fails.append(f"{where}: refused: {exc}")
    else:
        if not isinstance(bundle.get("results"), dict):
            fails.append(f"{where}: no result ({bundle.get('failure')})")
        else:
            fails += _drained(bundle["results"], where)
            want = NO_LOW_K_SENTENCE.format(surface="circle")
            if want not in log.splitlines():
                fails.append(f"{where}: the Log does not print {want!r}")
            sec = _report_section(sd, bundle)
            prose = [b.text for b in sec.blocks if b.kind == "prose"] if sec else []
            if NO_LOW_K_SENTENCE.format(surface="specified circle") not in prose:
                fails.append(f"{where}: the report does not say the specified "
                             f"circle crosses no material with d and ψ")

    # --- D: the full-pool rule applies to these circles too ---------------
    from xslope.advanced import rapid_drawdown
    from xslope.slice import generate_slices
    sd = _load(STARTER_XLSX)
    shell = next(m for m in sd["materials"] if m["name"] == "shell")
    shell["c"], shell["phi"] = 0.0, 20.0
    with contextlib.redirect_stdout(io.StringIO()):
        ok, made = generate_slices(sd, circle=dict(SHELL_CIRCLE), num_slices=40)
    if not ok:
        fails.append(f"stage 1 rule: the shell circle would not slice: {made}")
    else:
        ok, msg = rapid_drawdown(made[0], "spencer", debug_level=0)
        if ok:
            fails.append(f"stage 1 rule: a shell circle at stage 1 F < 1 was "
                         f"solved (F = {msg['FS']:.3f}, stage 1 "
                         f"{msg['stage1_FS']:.3f})")
        elif "stable before drawdown" not in str(msg):
            fails.append(f"stage 1 rule: refused for another reason: {msg}")

    # --- E: the core circle is unchanged ----------------------------------
    sd = _load(STARTER_XLSX)
    with contextlib.redirect_stdout(io.StringIO()):
        ok, made = generate_slices(sd, circle=dict(CORE_CIRCLE), num_slices=40)
        ok, res = rapid_drawdown(made[0], "spencer", debug_level=0)
    if not ok:
        fails.append(f"core circle: refused: {res}")
    else:
        if f"{res['FS']:.3f}" != "1.196" or res["stage3_run"]:
            fails.append(f"core circle: F = {res['FS']:.4f}, stage 3 run "
                         f"{res['stage3_run']}; expected 1.196, not run")
        if sum(res["rapid_low_k"]) != 13 or crosses_no_low_k(res):
            fails.append(f"core circle: {sum(res['rapid_low_k'])} low-permeability "
                         f"slices recorded, expected the 13 in the core")
    return fails


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    f = run()
    print("\n".join(["FAILED:"] + f) if f else
          "\nAll rapid drawdown drained-surface checks passed.")
    sys.exit(1 if f else 0)
