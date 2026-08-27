"""Build ``xslope_earth_dam_fs_time.xlsx``, Tutorial COMBO-3's model — the cored
earth dam of Tutorial SEEP-3 given the strengths a stability run needs.

COMBO-3 puts a factor-of-safety-versus-time curve on the drawdown SEEP-3 solved,
and it needs no new seepage input to do it. The build starts from
``build_earth_dam_tseep.seep03_model()`` — the same section, the same two zones,
the same storage, the same boundary set and the same schedule — and adds only
what a limit equilibrium run reads: a strength band on the materials table, a
starting circle on each face, and a minimum slip depth that keeps the search off
surficial slivers.

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
    dict(option="mc", gamma=20.0, gamma_sat=21.0, c=0.0, phi=32.0, u="seep"),
    dict(option="mc", gamma=19.0, gamma_sat=20.0, c=10.0, phi=25.0, u="seep"),
]

#: Two starting circles, one per face, each placed on the deep mechanism its
#: slope can make: entering near the crest, tangent to the rock (``Depth`` is an
#: ELEVATION, and 0 is the rock the dam sits on), and daylighting near the toe.
#: The upstream face runs from the heel at x = 0 to the upstream crest edge at
#: x = 51, the downstream face from the downstream crest edge at x = 59 to the toe
#: at x = 110, so a circle spanning either face has its center out beyond the toe
#: and a radius near the span. Both faces are offered because both govern over a
#: run this long: a drawdown attacks the upstream slope, but a full reservoir
#: loads it, and under a full pool the weaker mechanism is the downstream one.
#: Which one wins at a given instant is left to the search.
#:
#: Both centers are one metre off the round value a hand would reach for, and that
#: is deliberate on the upstream side: a circle centered at (7, 55) with R 55
#: passes exactly through the crest-edge vertex at (51, 22), and a starting circle
#: that meets the ground surface at one of its own vertices has a degenerate trace.
CIRCLES = [dict(Xo=7.0, Yo=56.0, Depth=0.0, R=56.0),      # upstream face
           dict(Xo=103.0, Yo=59.0, Depth=0.0, R=59.0)]    # downstream face

#: The one search limit the file keeps. ``min_slip_depth`` puts a floor on the
#: SIZE of the mechanism. The shell is cohesionless, so nothing in it sets a
#: scale: the shallower a trial surface is drawn, the lower it scores, and an
#: unfiltered search walks down to a surficial wedge whose factor of safety is the
#: infinite-slope value and whose failure is raveling rather than a slide. 8 m
#: below the ground surface, on a 22 m dam, holds the search on slides of the
#: embankment rather than on crest slivers. No entry or exit range is declared:
#: which face governs is what the curve measures, so confining the trace to one of
#: them would answer the question in the input.
SEARCH_WINDOW = dict(min_slip_depth=8.0)

#: The saved-frame schedule this page sweeps, replacing SEEP-3's coarser one. A
#: curve is only as fine as the frames under it — a time that is not saved has no
#: field to read and is never interpolated — so the instants are packed where the
#: answer moves and spread where it does not. Every five days through the
#: drawdown, which ends on day 47; then widening steps through the recovery, when
#: the core is draining and each frame gains less than the one before. SEEP-3's
#: own file keeps its schedule: this is a stability question asked of the same
#: march, not a change to the seepage tutorial.
#:
#: ``save_interval`` is set to the full duration so the interval grid contributes
#: only the last frame and the list below states the schedule outright. The
#: stepper lands exactly on every saved time, and the widest gap between two of
#: them is 60 days, so the interval never caps a step.
DURATION = 300.0
SAVE_TIMES = [2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 47.0,
              55.0, 65.0, 80.0, 100.0, 130.0, 180.0, 240.0, 300.0]


def build():
    """SEEP-3's completed model with the strengths, both surfaces, the floor and
    the denser saved-frame schedule."""
    sd = seep03_model()
    for mat, upd in zip(sd["materials"], STRENGTHS):
        mat.update(upd)
    sd["tseep"] = dict(sd["tseep"], duration=DURATION,
                       save_interval=DURATION, save_times=list(SAVE_TIMES))
    sd["circles"] = [dict(c) for c in CIRCLES]
    sd["circular"] = True
    sd["search_window"] = dict(SEARCH_WINDOW)
    return _write(sd, TUTORIAL_FILES, OUT)


if __name__ == "__main__":
    build()
