"""The unconfined seepage solver's exit-face active set: the limit-cycle escapes,
and the threshold the switch measures an apparent inflow against.

Four corpus models used to run to their iteration ceiling and stop, reporting a flow
rate off a field that was still moving. Neither was diverging: each had settled into
an exact periodic orbit that the relaxation ladder's bottom rung made invisible.

CLASS A -- the HEAD FIELD cycles. Four or five interior nodes on a
high-conductivity-contrast interface oscillate with a period of 8 to 12 sweeps. The
relaxed iterate barely moves, so the head-change test is satisfied, while the
un-relaxed step the closure gate is measured on lands far away, so the closure test
never is. ``solve_unsaturated`` detects the orbit -- a RETURN to a sweep 2 to 25 back,
within 1e-5 of the head range, having been further away at every shorter lag, with the
exit-face set unchanged, held for 30 consecutive sweeps -- and lowers the relaxation
floor once, to 1e-3, granting a bounded extension of the sweep budget from that sweep.
Four committed-mesh models were rescued this way (earth_dam_rapid, gsat_seep and
earth_dam_down_seep, which are one field, and johnson_rapid_KEY). With each element's
kr taken as its exact area average none of them cycles: they converge unaided in 61 and
49 sweeps, and are locked as that in the inert table.

CLASS B -- the EXIT-FACE ACTIVE SET cycles. earth_dam2 alternated between two sets
because the all-or-nothing quadratic-edge rule forbids the self-consistent one
(midside active, corner not). The gate is the orbit itself: the active set returns to
a set it has already visited, after sweep 100. The veto is then lifted from the edges
carrying the nodes that moved during the orbit, and from nothing else.

CLASS C -- ONE ELEMENT holds the solve in a loop. One element on the wet/dry edge
sits at the bend of the kr curve: dry, it backs water up and wets itself; wet, it
drains and dries itself, and its kr comes back from each sweep far further out than it
went in -- a gain no relaxation the solver uses is fine enough for. Read only when the
sweep budget is spent unconverged, the escape solves for that element's kr exactly at
every sweep, against the current kr of every other element, until the field settles,
and the field must then pass the ordinary gate on one ordinary sweep. vp077a on the
quad8 mesh Studio builds at 100 divisions is the fixture.

THE CLASS A ESCAPE HAS NO FIXTURE NOW, and this file says so. Its four models converge
without it, and it fired on none of 164 other solves made when the exact kr average
went in: the 33 committed-mesh companions, and ten models rebuilt as tri3, tri6,
quad4, quad8 and quad9 at 100, 150 and 250 divisions under Studio's settings.

THE CLASS B ESCAPE NOW GUARDS NOTHING IN THE CORPUS, and this file says so rather
than pretending otherwise. The exit-face fix in ``seep.py`` -- a seepage face runs to
the toe where a specified head claims the corner, and a set with the head standing
above the ground is not a stable set -- removed the cycle that produced it. earth_dam2
now converges in 195 sweeps with NEITHER escape firing, so it has moved to the inert
table below, and the escape has no rescued model left. Measured over every steady
seepage model the suite solves (28 tagged models, plus the 12 committed-mesh solves
this file makes): the revisit line is printed on none of them. On earth_dam2 the
mechanism is direct -- the wet-face rule joins the last free node to the face at sweep
11, which populates ``_free_edges``, and the class B gate is written ``not
_free_edges``, so it can no longer be reached on that model at all.

THE CLASS B ESCAPE IS NO LONGER REACHED ON THIS CORPUS EVEN WITH ITS GATE LOWERED.
The face stops revisiting sets before the wet-face rule claims it, at sweep 11: with
``_SET_REVISIT_SWEEP`` lowered to 0 no revisit line is printed, and the run is the
ungated one sweep for sweep (156 sweeps, 1.2744670). The
leg below locks that, so a change that makes the escape reachable again shows up. The
code path has no fixture at all now, which is recorded here and not hidden.

AN EXIT FACE OVER A FLUX BOUNDARY -- the switch's threshold. An exit-face node's
turn-on test asks whether the boundary would have to push water into the domain, and
reads that off the reaction ``A.h - f_eff``. On an INACTIVE node the row is free, so
that reaction is the previous sweep's residual and settles to identically zero -- at a
node carrying a rain load the sign gating activation is then round-off. Measured on the
dam-infiltration model at 2.8x its own rain: +2.4e-22 against a nodal load of +2.2e-08,
at a node standing at psi = +0.055 m which therefore never joined the seepage face.
The switch instead compares the reaction against ``q_offered = max(f_ext, 0)``, the
water the boundary offers at that node, on BOTH the turn-on and the turn-off test: a
node held at psi = 0 taking in no more than the rain falling on it is infiltrating that
rain and shedding the rest, and only a reaction above the offered flux is water the
boundary invented.

WHAT THIS CHECK LOCKS

1. The four models the class A escape used to rescue converge without it, to the flow
   rates below; and the class B escape, which rescues nothing either, stays silent on
   earth_dam2 with its gate lowered.
2. No escape fires on a model that converges without it, and every such model's
   field, flow rate AND SWEEP COUNT are exactly what they were before the escapes
   existed. The sweep counts are in the table below for that reason: a trajectory
   that has been nudged shows up there before it shows up in a flow rate.
3. The dam under rain sheds it once its surface saturates, and no active exit node
   draws in more water than its boundary offers. Under a strict ``q <= 0`` the same
   rows read: nothing activates at 2.8x (the solve returns the flux-only field and its
   ponding warning), and at 3x and 4x the active set cycles — 13 to 14 of the 48 exit
   nodes — to the 2000-sweep ceiling without converging.
4. Every reported flow rate equals the outflow summed against it. The reported number
   is the inflow half of a balance the operator closes exactly, so the two agree at any
   solved field; they stop agreeing the moment an entry point is left out of the sum.
   Water absorbed at an ACTIVE exit-face node is one such entry point, reachable only
   over a flux boundary: held at p = 0 with its own load discarded, the node draws the
   share of that load the soil can still take. It is 1.9% of the throughput here at
   2.8x rain and 21% at 4x.
5. The class C escape closes vp077a on the quad8 mesh Studio builds at 100 divisions,
   under Studio's settings: the budget ends unconverged, the swing names an element,
   it settles, and the gate accepts the field on the next sweep. The drawdown
   section's two meshes that were once this leg's fixtures (15 and 14 divisions)
   converge on their own, in 41 and 45 sweeps, and are locked as that: no pinned
   element, and their flow rates.

WHY THE GATES ARE WHERE THEY ARE (each mutable in one line if a reader wants to see the
other side of it). THE NUMBERS IN THIS SECTION WERE MEASURED UNDER EARLIER kr RULES
(3 sampling points on tri6, 7 on tri3) and have not been re-measured with the exact
element average, under which every model named here follows a different trajectory:

* A period of 1 is NOT a cycle. vp046b and vp077a creep toward their answers for
  hundreds of sweeps, each repeating its previous sweep within tolerance. Replace the
  detection with a blanket ``relax = min(relax, 1e-3)`` on the ladder's bottom rung and
  two of the four models the class A gate touches lose their convergence: vp077a runs
  its 1000 sweeps without closing where it took 393, and johnson_rapid_KEY its 1000
  where the escape closes it at 472. A third, earth_dam_rapid, closes without the escape
  in 634 sweeps but on 2.0495640e-03 rather than 2.0495627e-03. vp046b and
  earth_dam1_vg are untouched, because the progress-fed ladder no longer drives either
  of them to the bottom rung at all.
* Neither is a repeat that never left -- a solve creeping monotonically sits inside
  the repeat tolerance at every lag and reads as an orbit unless the iterate is
  required to have gone somewhere in between. That false positive is the ``_cyc_away``
  term's whole job, and three of the inert models witness it: delete the term and
  vp077a's creep reads as an orbit at sweep 274, the escape floors its relaxation, and
  it runs its 1000 sweeps without closing where it took 393; earth_dam2 reads as one at
  sweep 204 and closes in 257 sweeps on 1.2734304 rather than 195 on 1.2734362; and
  earth_dam1_vg at sweep 170, closing in 336 on 37.716779 rather than 163 on
  37.716990.
* Revisiting sets is normal EARLY. Every converging model in the corpus revisits
  exit-face sets while its seepage face is still finding its extent, and all of them
  are done by sweep 41 (vp077a, the latest). Set ``_SET_REVISIT_SWEEP = 0`` and
  earth_dam1_vg is untouched (163 sweeps, 37.716990), while earth_dam2 fires the
  escape at sweep 11 and closes at 1.2734384 in 313 sweeps -- the same answer to 1.7e-6,
  reached along a different path. The gate now costs the corpus a trajectory, not an
  answer.
* Removing the Class B gate leaves earth_dam2 exactly where it was: with
  ``_SET_REVISIT_SWEEP`` set past any reachable sweep it converges in 195 sweeps at
  1.2734362, bit for bit the run with the gate in place -- which is the same fact as
  the escape's inertness above, read from the other side.

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
#: Empty. The four class A rows (earth_dam_rapid, gsat_seep, earth_dam_down_seep and
#: johnson_rapid_KEY) converge without the escape under the exact kr average and are
#: in INERT below; earth_dam2, once the class B row, is there too. A model found to
#: need an escape again goes here.
RESCUED = []

#: Converging without either escape, and required to stay that way, in about as many sweeps.
#: Chosen for what each one would catch: vp046b and vp077a are the two models a
#: blanket relaxation floor breaks, earth_dam1_vg cycles at a period of 2 for 46
#: sweeps on its way to converging, earth_dam1 and vp077a revisit exit-face sets
#: (7 and 9 times), earth_dam_bc2's two BC sets are the plain unconfined case, and
#: earth_dam2 is the model the class B escape used to rescue.
#:
#: Every row carries the exact element kr average. The first four are the models the
#: class A escape used to close in 525 and 472 sweeps; the three earth-dam rows are one
#: section, one mesh and one boundary set, so one field.
INERT = [
    ("docs/lem/files/xslope_earth_dam_rapid", 1, DOCS, 2.0497696e-03, 61),
    ("docs/lem/files/xslope_gsat_seep", 1, DOCS, 2.0497696e-03, 61),
    ("docs/lem/files/xslope_earth_dam_down_seep", 1, DOCS, 2.0497696e-03, 61),
    ("docs/lem/files/xslope_johnson_rapid_KEY", 1, DOCS, 1.8719944, 49),
    ("docs/verification/files/rocscience/vp046b", 1, VENDOR, 1.2786298e-03, 82),
    ("docs/verification/files/rocscience/vp077a", 1, DOCS, 8.2228459e-06, 305),
    ("docs/seep/files/xslope_earth_dam1_vg", 1, DOCS, 37.708280, 126),
    ("docs/seep/files/xslope_earth_dam1", 1, DOCS, 38.771584, 112),
    ("docs/inputs/seep/xslope_earth_dam_bc2", 1, DOCS, 42.431571, 89),
    ("docs/inputs/seep/xslope_earth_dam_bc2", 2, DOCS, 11.587495, 14),
    ("docs/seep/files/xslope_johnson_res", 1, DOCS, 1.9563070, 18),
    ("docs/seep/files/xslope_earth_dam2", 1, DOCS, 1.2744670, 156),
]

#: earth_dam2's exit face settles on 8 of its 97 exit-face nodes.
EDAM2_ACTIVE = 8
#: The sweep the wet-face rule resolves earth_dam2's face on, which is also what puts
#: the class B escape out of reach on this model (the escape is gated ``not
#: _free_edges``, and the wet-face rule is the other writer of that set).
EDAM2_WET_SWEEP = 11
#: earth_dam2 with the class B gate lowered to 0, the one setting that used to reach
#: the escape: it no longer fires (None), and the run is the ungated one -- the same
#: sweeps and the same flow rate.
EDAM2_B_FIRES_AT, EDAM2_B_SWEEPS, EDAM2_B_Q = None, 156, 1.2744670

#: The dam-infiltration tutorial model: a rain flux over the whole exposed surface,
#: draining to a 12 m toe drain that is the file's only exit face.
INFIL_MODEL = "docs/tutorials/files/xslope_dam_infiltration.xlsx"
#: The overlay EXTENDS that drain across the exposed surface. ``exit_face`` is a single
#: polyline, so the drain's own segment has to stay inside it -- replacing it rather
#: than extending it deletes the drain and halves the discharge, while still converging
#: and still reporting a full-looking active set.
INFIL_OVERLAY = [[40.0, 0.0], [52.0, 0.0], [28.0, 12.0], [24.0, 12.0], [20.0, 10.0]]
#: The page's own discretization and solver settings (tools/make_tutorial_figures.py,
#: SEEP04_SIZE / SEEP04_ELEMENT / SEEP04_MAX_ITER and ``_seep04_solve``'s tolerance),
#: so these flow rates are read off the same field its figures are.
INFIL_SIZE, INFIL_ELEMENT = 1.0, "tri3"
INFIL_SETTINGS = dict(tol=1e-5, max_iter=2000)
#: The mesh that discretization gives, and the exposed-surface share of the overlay's
#: 48 exit nodes (the other 13 are the drain). Pinned because the overlay adds no mesh
#: vertex the flux blocks had not already pinned: it must not move the mesh.
INFIL_NODES, INFIL_FACE_NODES = 473, 35
#: Rain, as a multiple of the file's own rate, against: exposed-face nodes the overlay
#: activates, the overlay's flow rate, the flux-only flow rate, and whether the
#: flux-only run raises the ponding warning. 1x is below the rate at which the surface
#: saturates, so the overlay is a no-op there and has to stay one; 2.8x is the onset,
#: where one node ponds under the flux BC alone; 3x and 4x are past it, and are the two
#: rates whose active set used to cycle to the sweep ceiling instead of converging.
INFIL_ROWS = [
    # Datum-independent head gate: the 1x rows are the two the tighter gate reaches.
    (1.0, 0, 4.915525799188735e-07, 4.91552423346594e-07, False),
    (2.8, 1, 8.823960999866383e-07, 8.874008969783335e-07, True),
    (3.0, 3, 9.202226475422334e-07, 9.411714513832222e-07, True),
    (4.0, 12, 1.0583269792370564e-06, 1.246222222222292e-06, True),
]
#: How closely the 1x overlay has to reproduce the flux-only solve to count as the
#: no-op it is: the two differ only in the iterate each solve arrives along.
INFIL_NOOP_RTOL = 1e-5
#: How closely a reported flow rate has to equal the OUTflow measured against it. The
#: two are the same water — the operator has zero row sums, so what enters leaves —
#: and the reported number is the inflow half. They agree to ~1e-14 of themselves; the
#: allowance is for summation order, not for a discrepancy.
INFIL_BALANCE_RTOL = 1e-9
#: Relative tolerance on a flow rate. The solve is deterministic on one machine, but
#: the linear algebra rounds differently from one machine to the next, and a cycling
#: solve amplifies that: on the two desktops vp046b's flow rate differs by 5e-6.
Q_RTOL = 1e-5
#: The same rounding moves a sweep count: johnson_rapid_KEY converges in 472 sweeps
#: at home and 439 at work, vp046b in 82 and 81. A count within 10% (or 2 sweeps) is
#: the same trajectory; a changed solver moves it by far more.
SWEEP_RTOL = 0.10

_CYCLE_LINE = "head field is cycling"
_REVISIT_LINE = "exit-face set returned to the set"
_WET_LINE = "joined to the seepage face"
_PIN_LINE = "solved for directly at every sweep"

#: The drawdown section -- earth_dam_rapid with its pool made a falling reservoir --
#: on the tri6 meshes test/transient_seep_check.py builds, solved for the steady state
#: its drawdown run starts from, at that run's 2000-sweep budget. These two sizes once
#: ended that budget unconverged and closed on the pinned-element escape (CLASS C
#: above). They converge on their own, and the rows lock that: divisions, flow rate,
#: sweeps, and no pinned element.
PIN_MODEL = "docs/lem/files/xslope_earth_dam_rapid.xlsx"
PIN_POOL = 302.0
PIN_MAX_ITER = 2000
PIN_ROWS = [
    (15.0, 2.0639789e-03, 41),
    (14.0, 2.0586170e-03, 45),
]

#: The pinned-element escape's fixture: vp077a on the quad8 mesh Studio's Build Mesh
#: makes at 100 divisions (3100 nodes), under Studio's settings. The budget ends
#: unconverged, the swing names one element, it settles in 205 sweeps, and the gate
#: accepts the field on sweep 401. The sweep count is held to SWEEP_RTOL like the
#: others. Of the 131 Studio-mesh solves measured with the exact kr average it is the
#: one that needs the escape and closes on it.
PIN_STUDIO_MODEL = "docs/verification/files/rocscience/vp077a.xlsx"
PIN_STUDIO_DIVISIONS, PIN_STUDIO_ELEMENT = 100, "quad8"
PIN_STUDIO_SETTINGS = dict(tol=1e-4, max_iter=400)
PIN_STUDIO_Q, PIN_STUDIO_SWEEPS = 8.2072266e-06, 205


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


def _sweeps_close(got, want):
    return abs(got - want) <= max(2, SWEEP_RTOL * want)


def leg_rescued():
    """Each cycling model converges, to the flow rate that is its answer."""
    fails = []
    for stem_rel, bc, settings, q_want, sweeps_want, which in RESCUED:
        name = os.path.basename(stem_rel)
        seep_data, sol, log = _solve(stem_rel, bc, settings)
        got = _sweeps(log)
        marker = _CYCLE_LINE if which == "A" else _REVISIT_LINE
        fired = _fired(log, marker)
        if _fired(log, _PIN_LINE) is not None:
            fails.append(f"{name} bc{bc} reached the pinned-element escape; the class "
                         f"{which} escape is what closes it")
        if not sol["converged"]:
            fails.append(f"{name} bc{bc} does not converge — the class {which} "
                         f"escape did not rescue it")
        if fired is None:
            fails.append(f"{name} bc{bc} converged without the class {which} escape "
                         f"firing, so this row no longer covers it")
        if not _close(float(sol["flowrate"]), q_want):
            fails.append(f"{name} bc{bc} flow rate {sol['flowrate']:.7g} for "
                         f"{q_want:.7g}")
        if not _sweeps_close(got, sweeps_want):
            fails.append(f"{name} bc{bc} converged in {got} sweeps for "
                         f"{sweeps_want}")
        print(f"  {name:<28} bc{bc}  converged {got:>4} sweeps  "
              f"q={sol['flowrate']:.7g}  class {which} escape at sweep {fired}")
    return fails


def leg_class_b_reachable():
    """earth_dam2 settles on the set the all-or-nothing edge rule once forbade, and
    it now gets there without the class B escape: the wet-face rule resolves the face
    instead, which is also what puts the escape out of reach here. Lowering the
    escape's gate used to reach it on this model; it no longer does, and the run with
    the gate lowered must be the ungated run."""
    fails = []
    seep_data, sol, log = _solve("docs/seep/files/xslope_earth_dam2", 1, DOCS)
    if _fired(log, _REVISIT_LINE) is not None:
        fails.append(f"earth_dam2 fired the class B escape at sweep "
                     f"{_fired(log, _REVISIT_LINE)}; it converges without it now, and "
                     f"the escape's inertness is what this leg records")
    wet = _fired(log, _WET_LINE)
    if wet != EDAM2_WET_SWEEP:
        fails.append(f"earth_dam2's wet-face rule fired at sweep {wet}, not "
                     f"{EDAM2_WET_SWEEP} — the mechanism that resolves this face, and "
                     f"that makes the class B gate unreachable on it, has moved")
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
          f"({sol['closure_fraction']:.2e} of the flow), wet-face rule at sweep {wet}, "
          f"class B escape did not fire")

    # The gate that used to reach the escape. With _SET_REVISIT_SWEEP at 0 a revisit
    # could fire before the wet-face rule claimed the edges; none occurs now, so the
    # escape stays silent and the solve is the ungated one.
    import xslope.seep as _seep
    _gate = _seep._SET_REVISIT_SWEEP
    try:
        _seep._SET_REVISIT_SWEEP = 0
        _sd_b, sol_b, log_b = _solve("docs/seep/files/xslope_earth_dam2", 1, DOCS)
    finally:
        _seep._SET_REVISIT_SWEEP = _gate
    fired_b, sweeps_b = _fired(log_b, _REVISIT_LINE), _sweeps(log_b)
    if fired_b != EDAM2_B_FIRES_AT:
        fails.append(f"with the gate at 0 the class B escape fired at sweep "
                     f"{fired_b}, not {EDAM2_B_FIRES_AT} — its reachability on this "
                     f"model has changed from what is recorded")
    if not sol_b["converged"] or sweeps_b != EDAM2_B_SWEEPS:
        fails.append(f"with the gate at 0 earth_dam2 took {sweeps_b} sweeps for "
                     f"{EDAM2_B_SWEEPS}")
    if not _close(float(sol_b["flowrate"]), EDAM2_B_Q):
        fails.append(f"with the gate at 0 earth_dam2 lands on "
                     f"{sol_b['flowrate']:.7f} for {EDAM2_B_Q:.7f}, the ungated "
                     f"answer")
    print(f"  class B escape (gate lowered to 0): fired at sweep {fired_b}, "
          f"{sweeps_b} sweeps, q={sol_b['flowrate']:.7f} against "
          f"{sol['flowrate']:.7f} ungated")
    return fails


def leg_inert():
    """A model that converges on its own is untouched — same answer, same sweeps,
    no escape fired."""
    fails = []
    for stem_rel, bc, settings, q_want, sweeps_want in INERT:
        name = os.path.basename(stem_rel)
        _seep_data, sol, log = _solve(stem_rel, bc, settings)
        got = _sweeps(log)
        cyc = _fired(log, _CYCLE_LINE)
        rev = _fired(log, _REVISIT_LINE)
        pin = _fired(log, _PIN_LINE)
        if pin is not None:
            fails.append(f"{name} bc{bc} converges on its own, but the pinned-element "
                         f"escape fired on it at sweep {pin}")
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
        if not _sweeps_close(got, sweeps_want):
            fails.append(f"{name} bc{bc} took {got} sweeps for {sweeps_want} — its "
                         f"trajectory changed even though its answer did not")
        print(f"  {name:<28} bc{bc}  {got:>4} sweeps  q={sol['flowrate']:.7g}  "
              f"no escape fired")
    return fails


def _solve_infiltration(rain_factor, overlay):
    """The dam-infiltration model at a multiple of its rain, with or without the exit
    face laid over its exposed surface. Returns (seep_data, solution, warnings)."""
    from xslope.fileio import load_slope_data
    from xslope.mesh import (build_mesh_from_polygons, extract_size_regions,
                             get_material_polygons)
    from xslope.seep import build_seep_data, run_seepage_analysis

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        model = load_slope_data(os.path.join(_REPO, INFIL_MODEL))
        bc = dict(model["seepage_bc"])
        bc["specified_fluxes"] = [dict(f, flux=f["flux"] * rain_factor)
                                  for f in bc["specified_fluxes"]]
        if overlay:
            bc["exit_face"] = [list(p) for p in INFIL_OVERLAY]
        model = dict(model, seepage_bc=bc)
        mesh = build_mesh_from_polygons(get_material_polygons(model), INFIL_SIZE,
                                        INFIL_ELEMENT,
                                        size_regions=extract_size_regions(model))
        seep_data = build_seep_data(mesh, model)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            solution = run_seepage_analysis(seep_data, **INFIL_SETTINGS)
    return seep_data, solution, [str(w.message) for w in caught]


def _ponding_warning(messages):
    return [m for m in messages if "positive pore pressure" in m]


def _outflow(seep_data, solution, active):
    """The water leaving the domain across its boundary, summed independently of the
    reported flow rate: the negative reactions at the Dirichlet-held rows (specified
    heads and active exit-face nodes) plus any prescribed outflow delivered on a free
    row. The reported flow rate is the INflow half of the same balance, so the two are
    equal at any solved field — which is what makes an incomplete enumeration of
    either half visible."""
    bc_type = seep_data["bc_type"]
    f_ext = np.asarray(solution["flux_nodal"], dtype=float)
    held = (bc_type == 1) | ((bc_type == 2) & active)
    reaction = np.asarray(solution["q"], dtype=float) - np.where(held, 0.0, f_ext)
    out = -float(np.sum(reaction[held & (reaction < 0)]))
    out -= float(np.sum(f_ext[~held & (f_ext < 0)]))
    return out


def leg_exit_face_over_flux():
    """Rain the soil can no longer take runs off an exit face laid over it, instead of
    ponding on a flux boundary that has no way to shed it."""
    fails = []
    for factor, face_want, q_overlay, q_flux, warn_flux in INFIL_ROWS:
        runs = {}
        for overlay in (False, True):
            seep_data, sol, msgs = _solve_infiltration(factor, overlay)
            nodes = seep_data["nodes"]
            on_exit = seep_data["bc_type"] == 2
            # An active exit node is held at p = 0; the drain is the exit face the
            # file itself carries, along the base downstream of x = 40.
            p = sol["head"] - nodes[:, 1]
            active = on_exit & (np.abs(p) < 1e-6)
            drain = on_exit & (nodes[:, 1] <= 1e-9) & (nodes[:, 0] >= 39.9)
            runs[overlay] = dict(sol=sol, msgs=msgs, nodes=nodes,
                                 outflow=_outflow(seep_data, sol, active),
                                 n_face=int(np.sum(on_exit & ~drain)),
                                 n_face_active=int(np.sum(active & ~drain)),
                                 # what the boundary supplies at each active node,
                                 # against the flux it offers there
                                 excess=float(np.max(
                                     (sol["q"] - np.maximum(sol["flux_nodal"], 0.0))[active]))
                                 if np.any(active) else 0.0)
            if len(nodes) != INFIL_NODES:
                fails.append(f"rain x{factor:g} {'overlay' if overlay else 'flux only'} "
                             f"meshed to {len(nodes)} nodes, not {INFIL_NODES}")
            if not sol["converged"]:
                fails.append(f"rain x{factor:g} "
                             f"{'with the overlay' if overlay else 'flux only'} does "
                             f"not converge")

        ov, fx = runs[True], runs[False]
        if ov["n_face"] != INFIL_FACE_NODES:
            fails.append(f"the overlay claims {ov['n_face']} exposed-face exit nodes, "
                         f"not {INFIL_FACE_NODES}")
        if ov["n_face_active"] != face_want:
            fails.append(f"rain x{factor:g}: the overlay activates "
                         f"{ov['n_face_active']} exposed-face node(s) of "
                         f"{ov['n_face']}, not {face_want}")
        for label, run, want in (("overlay", ov, q_overlay), ("flux only", fx, q_flux)):
            got = float(run["sol"]["flowrate"])
            if not _close(got, want):
                fails.append(f"rain x{factor:g} {label} flow rate {got:.7g} for "
                             f"{want:.7g}")
            # The reported flow rate is the throughput, so it equals the outflow
            # measured against it. It stops doing so the moment an entry point is left
            # out of the sum: at 2.8x the water absorbed at the one active face node is
            # 1.9% of the throughput, and at 4x it is 21%.
            if abs(got - run["outflow"]) > INFIL_BALANCE_RTOL * abs(run["outflow"]):
                fails.append(f"rain x{factor:g} {label} reports {got:.8g} against an "
                             f"outflow of {run['outflow']:.8g} — the reported inflow "
                             f"is not the throughput")
        # A seepage face cannot invent water: at an active node the boundary may
        # supply at most the flux prescribed there, and the rest of that flux runs off.
        if ov["excess"] > 0.0:
            fails.append(f"rain x{factor:g}: an active exit node draws "
                         f"{ov['excess']:.3e} more than the flux offered there")
        # The overlay is the model's answer to ponding, so it must not leave the
        # warning standing; the flux BC alone has no way to shed the water and does.
        if _ponding_warning(ov["msgs"]):
            fails.append(f"rain x{factor:g}: the overlay still reports ponding — "
                         f"{_ponding_warning(ov['msgs'])[0]}")
        if bool(_ponding_warning(fx["msgs"])) != warn_flux:
            fails.append(f"rain x{factor:g} flux only: ponding warning "
                         f"{'absent' if warn_flux else 'raised'}, unexpectedly")
        if face_want == 0 and abs(ov["sol"]["flowrate"] - fx["sol"]["flowrate"]) > \
                INFIL_NOOP_RTOL * abs(fx["sol"]["flowrate"]):
            fails.append(f"rain x{factor:g}: the overlay activates nothing yet moves "
                         f"the flow rate, {ov['sol']['flowrate']:.7g} against "
                         f"{fx['sol']['flowrate']:.7g}")
        print(f"  rain x{factor:<4g} overlay {ov['n_face_active']:>2}/{ov['n_face']} "
              f"face nodes active  q={ov['sol']['flowrate']:.7g} "
              f"(out {ov['outflow']:.7g})  flux only {fx['sol']['flowrate']:.7g}, "
              f"ponding {'reported' if _ponding_warning(fx['msgs']) else 'none'}")
    return fails


def _solve_pinned(divisions):
    """The drawdown section's steady state at its transient run's t = 0, on a tri6
    mesh of ``divisions`` elements across the section. Returns (solution, log)."""
    import copy
    from xslope.fileio import load_slope_data
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons
    from xslope.seep import build_seep_data, resolve_transient_bc, solve_unsaturated

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        d = load_slope_data(os.path.join(_REPO, PIN_MODEL))
        xs = [x for x, _ in d["ground_surface"].coords]
        mesh = build_mesh_from_polygons(get_material_polygons(d),
                                        (max(xs) - min(xs)) / divisions, "tri6")
        d_ts = copy.deepcopy(d)
        for sh in d_ts["seepage_bc"]["specified_heads"]:
            sh["head"], sh["kind"] = "res", "reservoir"
        sd = build_seep_data(mesh, d_ts)
        ta = np.array([0.0, 1.0])
        tseep = dict(times=ta, series={"res": (ta, np.full(2, PIN_POOL))}, duration=1.0)
        bt, bv, flux = resolve_transient_bc(sd, tseep, 0.0)
        mid = np.asarray(sd["element_materials"]) - 1
        per = lambda k: None if sd.get(k) is None else np.asarray(sd[k])[mid]  # noqa: E731
        out = solve_unsaturated(
            sd["nodes"], np.asarray(sd["elements"]), bt, bv, kr0=per("kr0_by_mat"),
            h0=per("h0_by_mat"), k1_vals=per("k1_by_mat"), k2_vals=per("k2_by_mat"),
            angles=per("angle_by_mat"), element_types=np.asarray(sd["element_types"]),
            vg_a=per("vg_a_by_mat"), vg_n=per("vg_n_by_mat"), vg_l=per("vg_l_by_mat"),
            model=per("unsat_by_mat"), flux_nodal=flux, max_iter=PIN_MAX_ITER)
    return out, buf.getvalue()


def _solve_studio_mesh(model_rel, divisions, element_type, settings):
    """Boundary set 1 on the mesh Studio's Build Mesh makes of the model (its polygons,
    constraint lines and size regions, sized as the section width over the divisions).
    Returns (solution, log)."""
    from xslope.fileio import load_slope_data
    from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,
                             extract_constraint_line_geometry, extract_joint_options,
                             extract_size_regions)
    from xslope.seep import build_seep_data, run_seepage_analysis

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        d = load_slope_data(os.path.join(_REPO, model_rel))
        lines, _n_reinf, _n_pile = extract_constraint_line_geometry(d)
        xs = [x for x, _ in d["ground_surface"].coords]
        mesh = build_mesh_from_polygons(
            get_material_polygons(d, reinf_lines=lines),
            target_size=(max(xs) - min(xs)) / divisions,
            element_type=element_type, lines=lines or None,
            element_size_1d=d.get("element_size_1d"),
            size_regions=extract_size_regions(d),
            joint_lines=extract_joint_options(d))
        sol = run_seepage_analysis(build_seep_data(mesh, d, seep_bc=1), **settings)
    return sol, buf.getvalue()


def leg_pinned():
    """The steady state the budget leaves unconverged closes once the one element
    holding it in its loop is solved for at every sweep, and the ordinary gate accepts
    the field it settles on; the drawdown section's meshes that used to need that
    converge on their own."""
    fails = []
    for divisions, q_want, sweeps_want in PIN_ROWS:
        out, log = _solve_pinned(divisions)
        tag = f"drawdown section, tri6 {divisions:g} divisions"
        if not out[5]:
            fails.append(f"{tag}: does not converge")
        if _fired(log, _PIN_LINE) is not None or _fired(log, _CYCLE_LINE) is not None:
            fails.append(f"{tag}: an escape fired on a solve that converges without one")
        if not _sweeps_close(_sweeps(log), sweeps_want):
            fails.append(f"{tag}: converged on sweep {_sweeps(log)} for {sweeps_want}")
        if not _close(float(out[3]), q_want):
            fails.append(f"{tag}: flow rate {out[3]:.7e} for {q_want:.7e}")
        print(f"  {tag}: no escape, converged {out[5]} on sweep "
              f"{_sweeps(log)}, q={out[3]:.7e}")

    sol, log = _solve_studio_mesh(PIN_STUDIO_MODEL, PIN_STUDIO_DIVISIONS,
                                  PIN_STUDIO_ELEMENT, PIN_STUDIO_SETTINGS)
    tag = (f"{os.path.basename(PIN_STUDIO_MODEL)}, {PIN_STUDIO_ELEMENT} "
           f"{PIN_STUDIO_DIVISIONS} divisions")
    outcomes = []
    for line in log.splitlines():
        if _PIN_LINE in line:
            n = int(line.rsplit("(", 1)[1].split()[0] if "did not settle" in line
                    else line.split("settled in ")[1].split()[0])
            outcomes.append(("did not settle" not in line, n))
    budget = PIN_STUDIO_SETTINGS["max_iter"]
    if not sol["converged"]:
        fails.append(f"{tag}: does not converge")
    if len(outcomes) != 1 or not outcomes[0][0]:
        fails.append(f"{tag}: pinned elements {outcomes}, not one that settles")
    elif not _sweeps_close(outcomes[0][1], PIN_STUDIO_SWEEPS):
        fails.append(f"{tag}: the pinned element settled in {outcomes[0][1]} sweeps "
                     f"for {PIN_STUDIO_SWEEPS}")
    if _sweeps(log) != budget + 1:
        fails.append(f"{tag}: converged on sweep {_sweeps(log)}, not on the one "
                     f"ordinary sweep after the budget ({budget + 1})")
    if not _close(float(sol["flowrate"]), PIN_STUDIO_Q):
        fails.append(f"{tag}: flow rate {sol['flowrate']:.7e} for {PIN_STUDIO_Q:.7e}")
    print(f"  {tag}: pinned {outcomes}, converged {sol['converged']} on sweep "
          f"{_sweeps(log)}, q={sol['flowrate']:.7e}")
    return fails


LEGS = [
    ("models rescued from a limit cycle", leg_rescued),
    ("earth_dam2's exit face, and the class B escape's reachability",
     leg_class_b_reachable),
    ("models that converge without an escape", leg_inert),
    ("an exit face over a flux boundary", leg_exit_face_over_flux),
    ("one element holding a steady state in its loop", leg_pinned),
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
    print("\nPASS: the class B escape stays silent on earth_dam2 with its gate lowered, "
          "every model that converges without the escapes (the four the class A escape "
          "used to rescue among them) takes about the same number of sweeps, the dam "
          "under rain sheds what it cannot take, and the one steady state that still "
          "loops on a single element closes on it.")


if __name__ == "__main__":
    main()
