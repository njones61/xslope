"""The unconfined seepage solver's exit-face active set: the two limit-cycle escapes,
and the threshold the switch measures an apparent inflow against.

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

1. The four models converge, to the flow rates below.
2. Neither escape fires on a model that converges without it, and every such model's
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

MEASURED DURING THE ROUND, NOT LOCKED HERE

The rescued flow rates are the MODELs' and not the settings': re-solved at relaxation
floors of 1e-3 and 3e-4 the fields agree to 0.04 psf (0.012 psf on earth_dam_rapid),
and tightening the head tolerance from 1e-4 to 1e-6 changes neither field nor sweep
count on any of the four -- the closure test, not the head test, is what closes these
solves. No leg re-solves them at a second setting to prove it: that would double the
runtime for a property the trajectory pins above already constrain, since a solve
whose answer had become a function of its settings would have to reach it along a
different path, and the escape sweeps and sweep counts locked below are that path.

WHY THE GATES ARE WHERE THEY ARE (both measured over the corpus, and both mutable
in one line if a reader wants to see the other side of them):

* A period of 1 is NOT a cycle. vp046b and vp077a creep toward their answers for
  hundreds of sweeps, each repeating its previous sweep within tolerance. Replace the
  detection with ``relax = min(relax, 1e-3)`` for every sweep past 120 and both lose
  their convergence (207 sweeps -> 600 without closing, and 500 -> 1000), earth_dam1_vg
  takes 687 sweeps for 180, and johnson_rapid_KEY does not converge in 1000 either.
* Neither is a repeat that never left. After earth_dam2's exit face settles, the
  field creeps at ~1e-8 of the head range per sweep, which sits inside the repeat
  tolerance at every lag and reads as a period-2 orbit unless the iterate is required
  to have gone somewhere in between. Delete the ``_cyc_away`` term and earth_dam2
  drops its floor at sweep 789 and stops at its ceiling with its closure test still
  short of closing, while earth_dam1_vg fires at 164 and takes 305 sweeps for 180.
* Revisiting sets is normal EARLY. Every converging model in the corpus revisits
  exit-face sets while its seepage face is still finding its extent, and all of them
  are done by sweep 41 (vp077a, the latest); earth_dam2 is still revisiting at sweep
  992. Set ``_SET_REVISIT_SWEEP = 0`` and earth_dam1_vg's committed field moves --
  its flow rate goes to 40.271920 for 40.120963 -- while earth_dam2 converges to
  1.272890, a different answer from the one it reaches when the orbit is what frees
  the edge.
* Removing the Class B gate leaves earth_dam2 where it was: with
  ``_SET_REVISIT_SWEEP`` set past any reachable sweep it runs to its 1000-sweep
  ceiling and reports 1.275147, the flow rate off the moving field.

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
    ("docs/seep/files/xslope_earth_dam1_vg", 1, DOCS, 37.711830, 346),
    ("docs/seep/files/xslope_earth_dam1", 1, DOCS, 38.781841, 115),
    ("docs/inputs/seep/xslope_earth_dam_bc2", 1, DOCS, 42.437178, 111),
    ("docs/inputs/seep/xslope_earth_dam_bc2", 2, DOCS, 11.587548, 13),
    ("docs/seep/files/xslope_johnson_res", 1, DOCS, 1.939071, 26),
]

#: earth_dam2's exit face settles on 8 of its 97 exit-face nodes.
EDAM2_ACTIVE = 8

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
    (1.0, 0, 4.915513796099751e-07, 4.915514596217149e-07, False),
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


LEGS = [
    ("models rescued from a limit cycle", leg_rescued),
    ("earth_dam2's exit-face set", leg_earth_dam2_set),
    ("models that converge without either escape", leg_inert),
    ("an exit face over a flux boundary", leg_exit_face_over_flux),
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
    print("\nPASS: the four cycling models converge to their own answers, every model "
          "that converged without the escapes is unchanged sweep for sweep, and the "
          "dam under rain sheds what it cannot take.")


if __name__ == "__main__":
    main()
