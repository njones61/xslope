"""Build ``xslope_earth_dam_fs_time.xlsx``, Tutorial COMBO-3's model — the cored
earth dam of Tutorial SEEP-3 given the strengths a stability run needs.

COMBO-3 puts a factor-of-safety-versus-time curve on the drawdown SEEP-3 solved,
and it needs no new seepage input to do it. The build starts from
``build_earth_dam_tseep.seep03_model()`` — the same section, the same two zones,
the same storage, the same boundary set and the same schedule — and adds only
what a limit equilibrium run reads: a strength band on the materials table, a
starting circle on the upstream face and the search window that holds the curve
on one mechanism.

One file rather than a starter/completed pair, on COMBO-1's pattern: the page
opens a finished model and spends its length on the runs. It is sidecar-free —
meshing, the steady solve and the march are what the page teaches, so nothing may
arrive already solved.

The strengths are typical values for a granular shell and a compacted clay core,
chosen for the exercise rather than measured on this dam; SEEP-3's file carries
none, because a seepage analysis reads none. The tutorial page states this in the
same words.

Run:  PYTHONPATH=. python3 tools/build_earth_dam_fs_time.py
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from build_earth_dam_tseep import TUTORIAL_FILES, _write, seep03_model  # noqa: E402

OUT = "xslope_earth_dam_fs_time.xlsx"

#: The limit equilibrium band, indexed to the base file's material order. Drained
#: effective-stress strengths on both zones, and ``u = seep`` on both so every
#: slice base reads its pore pressure from the solved field rather than from a
#: sketched line — which is the whole point of a curve driven by a march.
#: ``gamma`` is the moist unit weight above the water table and ``gamma_sat`` the
#: saturated one below it; the slicer splits each slice's weight at the water
#: table, so a drawdown that drains the shell lightens it.
STRENGTHS = [
    dict(option="mc", gamma=20.0, gamma_sat=21.0, c=0.0, phi=35.0, u="seep"),
    dict(option="mc", gamma=19.0, gamma_sat=20.0, c=10.0, phi=25.0, u="seep"),
]

#: The starting circle, by the usual construction on the UPSTREAM face — the side
#: a drawdown attacks. Center at mid-slope in x (the face runs from the heel at
#: x = 0 to the upstream crest edge at x = 51) and two dam heights above the toe
#: in y; ``Depth`` is an ELEVATION, and elevation 0 is the rock the dam sits on.
CIRCLE = dict(Xo=25.0, Yo=44.0, Depth=0.0, R=44.0)

#: The search window, in the circles sheet's own keys. Entry is the crest-side
#: end of the trace and exit the toe-side end, so the pair confines the search to
#: surfaces that break above the full-pool waterline — the reservoir meets the
#: face at (42, 18) — and daylight below it, on the upstream slope. Without them
#: the critical surface at full pool is a downstream toe circle, which the
#: reservoir never touches and a drawdown curve says nothing about.
#: ``min_slip_depth`` puts a floor on the SIZE of the mechanism. The shell is
#: cohesionless, so nothing in it sets a scale: the shallower a trial surface is
#: drawn, the lower it scores, and an unfiltered search walks down to a surficial
#: wedge whose factor of safety is the infinite-slope value and whose failure is
#: raveling rather than a slide. 8 m below the ground surface, on a 22 m dam,
#: holds the search on slides of the embankment. Without it the reported minimum
#: is set by the smallest surface the window happens to allow.
SEARCH_WINDOW = dict(entry_x_min=42.0, entry_x_max=59.0,
                     exit_x_min=0.0, exit_x_max=42.0,
                     min_slip_depth=8.0)


def build():
    """SEEP-3's completed model with the strengths, the surface and the window."""
    sd = seep03_model()
    for mat, upd in zip(sd["materials"], STRENGTHS):
        mat.update(upd)
    sd["circles"] = [dict(CIRCLE)]
    sd["circular"] = True
    sd["search_window"] = dict(SEARCH_WINDOW)
    return _write(sd, TUTORIAL_FILES, OUT)


if __name__ == "__main__":
    build()
