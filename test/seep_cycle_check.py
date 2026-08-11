"""The two limit-cycle escapes in the unconfined seepage solver: what they rescue,
and what they must leave alone.

Four corpus models used to run to their iteration ceiling and stop, reporting a flow
rate off a field that was still moving. Neither was diverging: each had settled into
an exact periodic orbit that the relaxation ladder's terminal 0.01 made invisible.

CLASS A -- the HEAD FIELD cycles. Four or five interior nodes on a
high-conductivity-contrast interface oscillate with a period of 8 to 12 sweeps. The
relaxed iterate barely moves, so the head-change test is satisfied, while the
un-relaxed step the closure gate is measured on lands far away, so the closure test
never is. ``solve_unsaturated`` detects the orbit -- a RETURN to a sweep 2 to 25 back,
within 1e-6 of the head range, having been further away at every shorter lag, with the
exit-face set unchanged, held for 30 consecutive sweeps -- and lowers the relaxation
floor once, to 1e-3, granting a bounded extension of the sweep budget from that sweep.

CLASS B -- the EXIT-FACE ACTIVE SET cycles. earth_dam2 alternates between two sets
because the all-or-nothing quadratic-edge rule forbids the self-consistent one
(midside active, corner not). The gate is the orbit itself: the active set returns to
a set it has already visited, after sweep 100. The veto is then lifted from the edges
carrying the nodes that moved during the orbit -- 2 edges of 97 exit nodes on
earth_dam2 -- and from nothing else.

WHAT THIS CHECK LOCKS

1. The four models converge, to the flow rates below.
2. Those flow rates are the MODEL's and not the settings': re-solved at relaxation
   floors of 1e-3 and 3e-4 the fields agree to 0.04 psf, and tightening the head
   tolerance from 1e-4 to 1e-6 changes neither field nor sweep count (the closure
   test, not the head test, is what closes these solves).
3. Neither escape fires on a model that converges without it, and every such model's
   field, flow rate AND SWEEP COUNT are exactly what they were before the escapes
   existed. The sweep counts are in the table below for that reason: a trajectory
   that has been nudged shows up there before it shows up in a flow rate.

WHY THE GATES ARE WHERE THEY ARE (both measured over the corpus, and both mutable
in one line if a reader wants to see the other side of them):

* A period of 1 is NOT a cycle. vp046b and vp077a creep toward their answers for
  hundreds of sweeps repeating their previous sweep within tolerance; applying the
  floor to them from sweep 121 instead of on detection costs vp046b its convergence
  (207 sweeps -> 600 without closing) and vp077a its (500 -> 1000 without closing).
  Reproduce by dropping the detection and setting ``relax = min(relax, 1e-3)``
  unconditionally after the ladder in ``solve_unsaturated``.
* Neither is a repeat that never left. After earth_dam2's exit face settles, the
  field creeps at ~1e-8 of the head range per sweep, which sits inside the repeat
  tolerance at every lag and reads as a period-2 orbit unless the iterate is required
  to have gone somewhere in between. Reproduce by deleting the ``_cyc_away`` term:
  earth_dam2 then drops its floor at sweep 789 and stops at its ceiling with the
  closure still an order of magnitude short.
* Revisiting sets is normal EARLY. Every converging model in the corpus revisits
  exit-face sets while its seepage face is still finding its extent, and all of them
  are done by sweep 41 (vp077a, the latest); earth_dam2 is still revisiting at sweep
  992. Reproduce by setting ``_SET_REVISIT_SWEEP = 0``.
* Removing the Class B gate leaves earth_dam2 where it was: with the revisit branch
  deleted it runs to its 1000-sweep ceiling without converging.

Run directly:  PYTHONPATH=. python3 test/seep_cycle_check.py
"""

import contextlib
import io
import os
import sys
import warnings

import numpy as np

warnings.filterwarnings('ignore')

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: The settings each model's committed companion is recorded under
#: (tools/make_seep_sidecars.py owns them).
DOCS = dict(tol=1e-4, max_iter=1000)
VENDOR = dict(tol=1e-5, max_iter=600)

#: Rescued by an escape: (stem, bc, settings, flow rate, sweeps, which escape).
RESCUED = [
    ("docs/lem/files/xslope_earth_dam_rapid", 1, DOCS, 183.958670, 463, "A"),
    ("docs/lem/files/xslope_gsat_seep", 1, DOCS, 183.958670, 463, "A"),
    ("docs/lem/files/xslope_johnson_rapid_KEY", 1, DOCS, 1.871900, 901, "A"),
    ("docs/seep/files/xslope_earth_dam2", 1, DOCS, 1.273436, 927, "B"),
]

#: Converging without either escape, and required to stay that way, sweep for sweep.
#: Chosen for what each one would catch: vp046b and vp077a are the two models a
#: blanket relaxation floor breaks, earth_dam1_vg cycles at a period of 2 for 46
#: sweeps on its way to converging, earth_dam1 and vp077a revisit exit-face sets
#: (7 and 9 times), and earth_dam_bc2's two BC sets are the plain unconfined case.
INERT = [
    ("docs/verification/files/rocscience/vp046b", 1, VENDOR, 1.278746e-03, 207),
    ("docs/verification/files/rocscience/vp077a", 1, DOCS, 8.221076e-06, 500),
    ("docs/seep/files/xslope_earth_dam1_vg", 1, DOCS, 40.120963, 180),
    ("docs/seep/files/xslope_earth_dam1", 1, DOCS, 42.439945, 98),
    ("docs/inputs/seep/xslope_earth_dam_bc2", 1, DOCS, 42.437178, 111),
    ("docs/inputs/seep/xslope_earth_dam_bc2", 2, DOCS, 11.587548, 13),
    ("docs/seep/files/xslope_johnson_res", 1, DOCS, 1.939071, 26),
]

#: earth_dam2's exit face settles on 8 of its 97 exit-face nodes.
EDAM2_ACTIVE = 8
#: Relative tolerance on a flow rate. The solve is deterministic — these are the
#: numbers it gives, not numbers it is near — so this is a rounding allowance.
Q_RTOL = 1e-6

_CYCLE_LINE = "head field is cycling"
_REVISIT_LINE = "exit-face set returned to the set"


def _solve(stem_rel, bc, settings):
    """Solve one BC set on the model's committed mesh. Returns (solution, log)."""
    from xslope.fileio import load_slope_data
    from xslope.mesh import import_mesh_from_json
    from xslope.seep import build_seep_data, run_seepage_analysis

    stem = os.path.join(_REPO, stem_rel)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        slope_data = load_slope_data(f"{stem}.xlsx")
        mesh = import_mesh_from_json(f"{stem}_mesh.json")
        seep_data = build_seep_data(mesh, slope_data, seep_bc=bc)
        solution = run_seepage_analysis(seep_data, **settings)
    return seep_data, solution, buf.getvalue()


def _sweeps(log):
    """The sweep the solve converged on, or -1 if it never did."""
    for line in log.splitlines():
        if line.startswith("Converged in "):
            return int(line.split()[2])
    return -1


def _fired(log, marker):
    """The sweep an escape fired on, or None."""
    for line in log.splitlines():
        if marker in line:
            return int(line.split()[1].rstrip(":"))
    return None


def _close(got, want):
    return abs(got - want) <= Q_RTOL * max(abs(want), 1e-12)


def leg_rescued():
    """Each cycling model converges, to the flow rate that is its answer."""
    fails = []
    for stem_rel, bc, settings, q_want, sweeps_want, which in RESCUED:
        name = os.path.basename(stem_rel)
        seep_data, sol, log = _solve(stem_rel, bc, settings)
        got = _sweeps(log)
        marker = _CYCLE_LINE if which == "A" else _REVISIT_LINE
        fired = _fired(log, marker)
        if not sol["converged"]:
            fails.append(f"{name} bc{bc} does not converge — the class {which} "
                         f"escape did not rescue it")
        if fired is None:
            fails.append(f"{name} bc{bc} converged without the class {which} escape "
                         f"firing, so this row no longer covers it")
        if not _close(float(sol["flowrate"]), q_want):
            fails.append(f"{name} bc{bc} flow rate {sol['flowrate']:.6f} for "
                         f"{q_want:.6f}")
        if got != sweeps_want:
            fails.append(f"{name} bc{bc} converged in {got} sweeps for "
                         f"{sweeps_want}")
        print(f"  {name:<28} bc{bc}  converged {got:>4} sweeps  "
              f"q={sol['flowrate']:.6f}  class {which} escape at sweep {fired}")
    return fails


def leg_earth_dam2_set():
    """earth_dam2 settles on the set the all-or-nothing edge rule forbade, and the
    escape that let it is the revisit, gated after sweep 100."""
    fails = []
    seep_data, sol, log = _solve("docs/seep/files/xslope_earth_dam2", 1, DOCS)
    fired = _fired(log, _REVISIT_LINE)
    if fired is None or fired <= 100:
        fails.append(f"the exit-face revisit escape fired at sweep {fired}, not "
                     f"after the sweep 100 gate")
    # The exit-face set the solve settles on, read back from the solved field: an
    # active exit node is held at p = 0, which is the state the edge rule could not
    # give this face a partial version of.
    on_face = seep_data["bc_type"] == 2
    p = sol["head"] - seep_data["nodes"][:, 1]
    n_active = int(np.sum(on_face & (np.abs(p) < 1e-6)))
    if n_active != EDAM2_ACTIVE:
        fails.append(f"earth_dam2 settles on {n_active} active exit-face nodes of "
                     f"{int(np.sum(on_face))}, not {EDAM2_ACTIVE}")
    # A flow rate is only as good as the balance behind it: this one closes to a
    # billionth of itself.
    if not (sol["closure_fraction"] < 1e-6):
        fails.append(f"earth_dam2 closes to {sol['closure_fraction']:.3e} of its "
                     f"flow, which is not a settled solution")
    print(f"  earth_dam2 exit face: {n_active}/{int(np.sum(on_face))} active, "
          f"closure {sol['closure_error']:.2e} "
          f"({sol['closure_fraction']:.2e} of the flow), escape at sweep {fired}")
    return fails


def leg_inert():
    """A model that converges on its own is untouched — same answer, same sweeps,
    neither escape fired."""
    fails = []
    for stem_rel, bc, settings, q_want, sweeps_want in INERT:
        name = os.path.basename(stem_rel)
        _seep_data, sol, log = _solve(stem_rel, bc, settings)
        got = _sweeps(log)
        cyc = _fired(log, _CYCLE_LINE)
        rev = _fired(log, _REVISIT_LINE)
        if cyc is not None:
            fails.append(f"{name} bc{bc} converges on its own, but the head-cycle "
                         f"escape fired on it at sweep {cyc}")
        if rev is not None:
            fails.append(f"{name} bc{bc} converges on its own, but the exit-face "
                         f"revisit escape fired on it at sweep {rev}")
        if not sol["converged"]:
            fails.append(f"{name} bc{bc} no longer converges")
        if not _close(float(sol["flowrate"]), q_want):
            fails.append(f"{name} bc{bc} flow rate moved: {sol['flowrate']:.7g} "
                         f"for {q_want:.7g}")
        if got != sweeps_want:
            fails.append(f"{name} bc{bc} took {got} sweeps for {sweeps_want} — its "
                         f"trajectory changed even though its answer did not")
        print(f"  {name:<28} bc{bc}  {got:>4} sweeps  q={sol['flowrate']:.7g}  "
              f"neither escape fired")
    return fails


LEGS = [
    ("models rescued from a limit cycle", leg_rescued),
    ("earth_dam2's exit-face set", leg_earth_dam2_set),
    ("models that converge without either escape", leg_inert),
]


def run():
    failures = []
    for label, fn in LEGS:
        print(f"{label}:")
        failures += fn()
    return failures


def main():
    failures = run()
    if failures:
        print("\nFAILED:")
        for f in failures:
            print(f"  {f}")
        raise SystemExit(1)
    print("\nPASS: the four cycling models converge to their own answers, and every "
          "model that converged without the escapes is unchanged sweep for sweep.")


if __name__ == "__main__":
    main()
