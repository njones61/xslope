"""Build the Johnson Reservoir rapid-drawdown workbooks from one parameter block.

Three files come out of the single set of constants below, so the rapid-drawdown
family cannot drift apart:

  ``docs/lem/files/xslope_johnson_res_rapid.xlsx``
      The worked example of ``docs/lem/rapid.md`` ("Rapid Drawdown from a Transient
      Solution"): the transient route alone, no second boundary set, no piezometric
      lines.

  ``docs/tutorials/files/xslope_johnson_rapid_start.xlsx``
      The starter file of tutorial COMBO-2. Strengths, the undrained core's d/psi,
      ``Water loads: auto`` and the coarse **piezometric pair** — and nothing else:
      no seepage boundary sets and no transient schedule, because building those is
      what the tutorial does. Every material takes ``u = piezo``.

  ``docs/tutorials/files/xslope_johnson_rapid.xlsx``
      The completed COMBO-2 model: the same section with **both** seepage boundary
      sets (full pool and drawn down), the transient pool schedule with its
      ``stage_1`` / ``stage_2`` times, and ``u = seep`` on every material. The
      piezometric pair is left on the file, so switching the ``u`` column back to
      ``piezo`` reproduces the tutorial's first run. No ``_seep.csv`` /
      ``_seep2.csv`` companions ship with it: a steady run writes them, and a
      transient run stages its two frames in memory.

Everything is built deterministically from the committed transient **seepage**
sample ``docs/seep/files/xslope_johnson_res_tseep.xlsx``, so the cross-section,
zones, conductivities, storage properties, saved-frame schedule and stage times
can never diverge from the seepage side of the same dam. That file is never
modified. The one value the tutorial pair overrides is the level the reservoir is
lowered to: the tutorial's drawdown stops at a residual pool 10 ft deep, while the
worked example of ``docs/lem/rapid.md`` keeps the base file's total drawdown to
the tailwater datum.

Run:  PYTHONPATH=. python3 tools/build_johnson_rapid.py
"""

from __future__ import annotations

import copy
import os

from xslope.fileio import (load_slope_data, save_slope_data_to_xlsx,
                           default_template_path)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEEP_FILES = os.path.join(REPO_ROOT, "docs", "seep", "files")
LEM_FILES = os.path.join(REPO_ROOT, "docs", "lem", "files")
TUT_FILES = os.path.join(REPO_ROOT, "docs", "tutorials", "files")

BASE = os.path.join(SEEP_FILES, "xslope_johnson_res_tseep.xlsx")

# --- the one parameter block ------------------------------------------------

#: Undrained core, Kc = 1 (CU) envelope, below the drained c' = 400 / phi' = 18.
#: The free-draining sand shell and the silty-sand foundation carry no d/psi and
#: are analyzed drained through the drawdown.
CORE_D, CORE_PSI = 250.0, 14.0

#: The critical upstream circle, located by a rapid-drawdown circular search over
#: the transient-staged fields. Yo - Depth = R keeps the tangent elevation
#: consistent. It toes near the upstream base, daylights just past the crest, and
#: cuts twelve of forty slices through the core.
CIRCLE = dict(Xo=275.0, Yo=235.0, R=160.0, Depth=75.0)

#: The pool pair. The full pool is the reservoir at elevation 160 that the base
#: file's `pool` series starts from. The drawdown stops 10 ft above the
#: elevation-100 tailwater datum: a residual pool at elevation 110 stands against
#: the upstream toe at the end of it, which is the level the second seepage
#: boundary set states as a steady problem and the level the tutorial's `pool`
#: series ramps down to.
POOL_FULL, POOL_DOWN = 160.0, 110.0

#: The drawn-down steady boundary set: the residual pool along the upstream
#: foreshore and up the face to elevation 110, the tailwater at 100, and the
#: downstream slope as an exit face. Ten feet of head still cross the section, so
#: this is an unconfined flow problem with a phreatic surface of its own, stated
#: exactly the way set 1 states the full-pool one.
BC2 = {
    "specified_heads": [
        {"head": POOL_DOWN,
         "coords": [(0.0, 100.0), (200.0, 100.0), (220.0, 110.0)], "kind": "head"},
        {"head": 100.0, "coords": [(550.0, 100.0), (750.0, 100.0)], "kind": "head"},
    ],
    "specified_fluxes": [],
    "exit_face": [(380.0, 180.0), (550.0, 100.0)],
}

#: The tutorial's reservoir schedule: full pool held for five days, then lowered
#: to the residual pool over the following 45. Only the last breakpoint differs
#: from the base seepage sample's series, whose times and stage instants are kept.
POOL_SCHEDULE = [POOL_FULL, POOL_FULL, POOL_DOWN]

#: The coarse piezometric pair. These are five- and six-point polylines a designer
#: would sketch from piezometer readings, not traces of a solved field: Line 1
#: mimics the full-pool phreatic surface of boundary set 1, Line 2 the drawn-down
#: one of boundary set 2 — the residual pool at 110 upstream, falling through the
#: core and running out to tailwater at 100. Both were read off the solved fields
#: and rounded. Line 1 leaves the upstream shell AT the reservoir level and Line 2
#: at the residual pool level, so the pools they describe are the ones the two
#: seepage boundary sets state; a line that came off the shell a foot low would put
#: two different pools on one model. Both tails sit on the ground surface from
#: x = 550 out, so neither invents a pond on the downstream foreshore.
PIEZO_1 = [(0.0, 160.0), (360.0, 160.0), (410.0, 120.0), (550.0, 100.0),
           (750.0, 100.0)]
PIEZO_2 = [(0.0, 110.0), (220.0, 110.0), (360.0, 107.0), (410.0, 103.0),
           (550.0, 100.0), (750.0, 100.0)]


def _tutorial_schedule(sd):
    """The base file's transient schedule with the pool ramped to the residual
    level instead of to the tailwater datum. The saved-frame schedule and both
    stage times are the base file's own."""
    ts = copy.deepcopy(sd["tseep"])
    ts["series"]["pool"] = list(POOL_SCHEDULE)
    return ts


def _base():
    """The seepage sample with the two rapid-drawdown additions every file shares:
    the undrained core and the critical upstream circle."""
    sd = load_slope_data(BASE)
    sd["materials"][1]["d"] = CORE_D
    sd["materials"][1]["psi"] = CORE_PSI
    # max_depth is left at the base file's value: the fixed circle carries its own
    # tangent Depth, and overriding max_depth would clip the domain through the core.
    sd["circular"] = True
    sd["circles"] = [dict(CIRCLE)]
    return sd


def _save(sd, path):
    save_slope_data_to_xlsx(sd, path, template=default_template_path())
    print(f"wrote {os.path.relpath(path, REPO_ROOT)}")
    return path


def build_lem_worked_example():
    """``docs/lem/rapid.md``'s transient worked example — the base file plus the
    undrained core and the circle, and nothing else."""
    return _save(_base(), os.path.join(LEM_FILES, "xslope_johnson_res_rapid.xlsx"))


def build_tutorial_completed():
    """COMBO-2's completed model: both boundary sets, the schedule, u = seep."""
    sd = _base()
    sd["tseep"] = _tutorial_schedule(sd)
    sd["seepage_bc2"] = copy.deepcopy(BC2)
    sd["has_seepage_bc2"] = True
    sd["piezo_line"] = list(PIEZO_1)
    sd["piezo_line2"] = list(PIEZO_2)
    for m in sd["materials"]:
        m["u"] = "seep"
    return _save(sd, os.path.join(TUT_FILES, "xslope_johnson_rapid.xlsx"))


def build_tutorial_starter():
    """COMBO-2's starter: the piezometric pair only — no boundary sets, no
    schedule, u = piezo."""
    sd = _base()
    sd["piezo_line"] = list(PIEZO_1)
    sd["piezo_line2"] = list(PIEZO_2)
    for m in sd["materials"]:
        m["u"] = "piezo"
    empty = {"specified_heads": [], "specified_fluxes": [], "exit_face": []}
    sd["seepage_bc"] = copy.deepcopy(empty)
    sd["seepage_bc2"] = copy.deepcopy(empty)
    sd["has_seepage_bc2"] = False
    sd["tseep"] = None
    return _save(sd, os.path.join(TUT_FILES, "xslope_johnson_rapid_start.xlsx"))


def build():
    return [build_lem_worked_example(),
            build_tutorial_starter(),
            build_tutorial_completed()]


if __name__ == "__main__":
    build()
