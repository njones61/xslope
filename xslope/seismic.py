"""Seismic slope stability beyond the pseudostatic factor of safety.

The pseudostatic force itself lives in every limit equilibrium method (a
horizontal force kW through each slice's center of gravity, directed by the
failure surface) and the two-stage analysis in :mod:`xslope.advanced`. This module
holds what is computed FROM them:

- :func:`yield_coefficient`: the seismic yield coefficient k_y, the seismic
  coefficient at which the factor of safety is exactly 1 (Duncan, Wright, and
  Brandon 2014, Sec. 10.1.3), on a given surface or by search.
- :func:`fem_yield_coefficient`: k_y by the finite element method, the largest
  seismic coefficient under which the slope stands at full strength.
- :func:`newmark_displacement`: the permanent displacement of the sliding mass
  under an acceleration record, as a rigid block that slides whenever the
  ground acceleration exceeds k_y·g (Newmark 1965; Duncan, Wright, and Brandon
  2014, Sec. 10.1.3), with :func:`read_record` to load the record.

References: Duncan, J.M., Wright, S.G., and Brandon, T.L. (2014), *Soil Strength
and Slope Stability*, 2nd ed., Wiley, Chapter 10. Newmark, N.M. (1965), "Effects
of earthquakes on dams and embankments," *Géotechnique* 15(2), 139–160.
"""
import os
import re

import numpy as np

from . import solve
from .advanced import seismic2_strength, seismic_two_stage
from .slice import generate_slices

#: The largest seismic coefficient the bracket search tries before it reports
#: that the slope does not reach F = 1.
K_MAX = 2.0

#: The root on one surface is taken when |F - 1| is this small, or when the
#: bracket on k is narrower than K_TOL.
F_TOL = 1e-5
K_TOL = 1e-6

#: A confirming search that finds a surface this far below F = 1 at the current
#: k_y moves the answer to that surface; one that finds nothing lower ends the
#: search for the minimum k_y.
CONFIRM_TOL = 0.001

#: The most confirming searches the minimum-k_y search makes.
MAX_ROUNDS = 6


class _OneSurface:
    """F(k) on one fixed surface.

    The slices are cut once. The slicer computes the seismic force of each slice
    as |k| times its soil weight (:mod:`xslope.slice`), so F at any k is a solve
    of the same slices with that column scaled: no re-slicing, and the pore
    pressures, loads, and strengths stay exactly as cut. Every solve starts from
    a fresh copy, because the solvers and the two-stage analysis write into the
    table they are given. A vertical seismic force scales with k the same way, at
    the model's own ratio k_v/k (``kv_ratio``). ``trace`` keeps every (k, F) the root search passed
    through, which is what the Log and the report show.
    """

    def __init__(self, slope_data, method, circle, non_circ, num_slices,
                 composite, seismic2, kv_ratio=0.0):
        from .search import AnalysisError
        sd = dict(slope_data)
        sd["k_seismic"] = 1.0
        sd["k_vertical"] = kv_ratio
        ok, out = generate_slices(sd, circle=circle, non_circ=non_circ,
                                  num_slices=num_slices, composite=composite,
                                  debug=False)
        if not ok:
            raise AnalysisError(str(out))
        self.df, self.failure_surface = out
        self.kw1 = self.df["kw"].to_numpy(copy=True)
        # The vertical force at k = 1: k_v/k times the soil weight, signed. It
        # scales with k as the horizontal one does (see yield_coefficient).
        self.kvw1 = self.df["kvw"].to_numpy(copy=True)
        self.method = method
        self.strength = seismic2_strength(seismic2)
        self.func = getattr(solve, method)
        self.trace = []

    def __call__(self, k):
        """(F, None) at coefficient k, or (None, the solver's reason)."""
        df = self.df.copy()
        df["kw"] = k * self.kw1
        df["kvw"] = k * self.kvw1
        if self.strength:
            ok, res = seismic_two_stage(df, self.method, debug_level=0,
                                        strength=self.strength)
        else:
            ok, res = self.func(df)
        if not ok:
            self.trace.append((float(k), None))
            return None, str(res)
        f = float(res["FS"])
        self.trace.append((float(k), f))
        return f, None


def _k_at_unity(F, where, hint=None, cancel_check=None):
    """The k at which ``F(k) = 1`` on one surface, and F there with no seismic
    force.

    F falls as k rises (the seismic force adds to the driving side and, on a
    frictional base, takes normal force off it), so F = 1 is a single crossing.
    k = 0 is checked first: a surface already at or below F = 1 with no seismic
    force has no yield coefficient. The bracket then grows from ``hint`` (the
    previous round's k_y) or 0.1, doubling up to :data:`K_MAX` and backing off
    from any k where the method finds no solution, and the root is
    closed by regula falsi with the Illinois modification, which keeps the
    bracket and needs a handful of solves where bisection needs twenty.
    """
    from .search import AnalysisError, _check_cancel

    def solve_at(k):
        _check_cancel(cancel_check)
        f, err = F(k)
        if f is None:
            raise AnalysisError(f"The method found no solution on {where} at "
                                f"k = {k:.4f}: {err}")
        return f

    f0 = solve_at(0.0)
    if f0 <= 1.0:
        raise AnalysisError(
            f"F = {f0:.3f} on {where} with no seismic force. A slope that is not "
            f"stable before the earthquake has no yield coefficient.")
    # The bracket. A method can lose its solution at a large k well past the
    # crossing (Spencer's and Morgenstern-Price's moment equilibrium often does),
    # so a trial with no solution is not the end: the next trial halves the gap
    # back toward the last k that solved, until one solves below F = 1 or the gap
    # closes with F still above 1.
    a, fa = 0.0, f0 - 1.0
    k = hint if hint and hint > 0 else 0.1
    k_fail, why = None, None
    while True:
        _check_cancel(cancel_check)
        f, err = F(k)
        if f is None:
            k_fail, why = k, err
            if k_fail - a < 1e-4:
                raise AnalysisError(
                    f"The method found no solution on {where} above k = {a:.4f}, "
                    f"where F = {fa + 1.0:.3f}, so F = 1 cannot be reached: {why}")
            k = 0.5 * (a + k_fail)
            continue
        fk = f - 1.0
        if fk <= 0.0:
            b, fb = k, fk
            break
        a, fa = k, fk
        if k >= K_MAX:
            raise AnalysisError(
                f"F = {f:.3f} on {where} at k = {K_MAX:g}, the largest "
                f"coefficient tried, so its yield coefficient is above {K_MAX:g}.")
        k = 0.5 * (a + k_fail) if k_fail is not None else min(2.0 * k, K_MAX)
    if fb == 0.0:
        return b, f0

    last = None
    c = b
    for _ in range(100):
        c = b - fb * (b - a) / (fb - fa)
        if not (a < c < b):
            c = 0.5 * (a + b)
        _check_cancel(cancel_check)
        f, _err = F(c)
        if f is None and c != 0.5 * (a + b):
            c = 0.5 * (a + b)            # one retry at the bracket's midpoint
        fc = solve_at(c) - 1.0 if f is None else f - 1.0
        if abs(fc) <= F_TOL or (b - a) <= K_TOL:
            break
        if fc > 0.0:                      # still stable: c is the new lower end
            a, fa = c, fc
            if last == "a":
                fb *= 0.5
            last = "a"
        else:
            b, fb = c, fc
            if last == "b":
                fa *= 0.5
            last = "b"
    return c, f0


def _surface_of(entry, circular, movement):
    """The surface a search's cache entry describes, in the form the slicer and
    the model take it: a circle dict, or non-circular points with the model's
    own Movement flags."""
    if circular:
        return {"Xo": float(entry["Xo"]), "Yo": float(entry["Yo"]),
                "Depth": float(entry["Depth"]),
                "R": float(entry["Yo"]) - float(entry["Depth"])}
    pts = entry["points"]
    return [{"X": float(x), "Y": float(y),
             "Movement": movement[i] if i < len(movement) else "Free"}
            for i, (x, y) in enumerate(pts)]


def _on_surface(slope_data, surface_def, circular):
    """A copy of the model whose specified surface is ``surface_def``."""
    sd = dict(slope_data)
    if circular:
        sd["circular"] = True
        sd["circles"] = [surface_def]
    else:
        sd["circular"] = False
        sd["non_circ"] = surface_def
    return sd


def yield_coefficient(slope_data, method, analysis="auto_search",
                      surface="circular", num_slices=40, seismic2=False,
                      composite=False, grid_seed=False, diagnostic=False,
                      cancel_check=None, fs_tol=None, tol=None, max_iter=None,
                      min_slip_depth=None, announce=True,
                      confirm_tol=CONFIRM_TOL, max_rounds=MAX_ROUNDS):
    """The seismic yield coefficient k_y: the seismic coefficient at which the
    factor of safety is exactly 1 (Duncan, Wright, and Brandon 2014, Sec. 10.1.3).

    Parameters
    ----------
    slope_data : dict
        The model. Its own seismic coefficient is not used; nothing in it is
        changed.
    method : str
        The solver's name for the method (``"spencer"``, ``"bishop"``, ...).
    analysis : {"auto_search", "single_surface"}
        ``"single_surface"`` finds k_y on the surface the model specifies (the
        first circle, or the non-circular surface). ``"auto_search"`` finds the
        surface with the minimum k_y: it searches for the critical surface with no
        seismic force, finds k_y on it, then searches again at that k. A search
        that finds a surface below F = 1 there moves the answer to that surface,
        and the rounds repeat until a search finds nothing more than
        ``confirm_tol`` below F = 1. A circular search starts every round from
        the model's circles. A non-circular search is a local one, so each
        round after the first searches twice, from the surface the last round
        found and from the model's entered surface, and keeps the lower: started
        only from the last surface, the rounds settled on RS2-68 Case 2 at
        k_y = 0.4366 where the entered surface leads to 0.4318, and started only
        from the entered surface they settle on VP104 at 0.1452 where the last
        surface leads to 0.1436 (both measured when the rounds stopped within
        0.002 of F = 1).
    surface : {"circular", "noncircular"}
    num_slices, seismic2, composite, grid_seed, diagnostic, cancel_check,
    fs_tol, tol, max_iter, min_slip_depth, announce
        As in :func:`xslope.search.run_lem_analysis`. ``seismic2`` finds k_y with
        the two-stage analysis: stage 1 has no seismic force, so each surface's
        undrained strengths are fixed and only stage 2 depends on k.
    confirm_tol : float
        How far below F = 1 a confirming search must find a surface for the
        answer to move to it.
    max_rounds : int
        The most confirming searches made.

    Returns
    -------
    dict
        The run bundle of :func:`xslope.search.run_lem_analysis`, solved on the
        k_y surface at k = k_y, so its factor of safety is 1 and its slices carry
        the seismic force at k_y. ``search`` is the search that found that
        surface (None on a single surface). ``options`` records ``ky=True``.
        ``yield`` holds the answer:

        - ``k_y``: the yield coefficient.
        - ``FS_k0``: F on the k_y surface with no seismic force (in a two-stage
          run, stage 2's F at k = 0, with the undrained strengths).
        - ``trace``: every (k, F) the root search passed through on that surface.
        - ``search``: True for the minimum-k_y search.
        - ``static``: the critical surface with no seismic force, its F, and its
          own k_y (``{"FS", "k_y", "surface"}``), the first of the book's
          candidate surfaces; None on a single surface.
        - ``rounds``: one entry per confirming search, ``{"k", "FS", "k_y"}``: the
          coefficient it searched at, the critical F it found, and k_y on that
          surface (None where it found nothing below F = 1).
        - ``converged``: False where ``max_rounds`` ran out first.
        - ``strength``: the two-stage strength choice, or None.

    Raises
    ------
    AnalysisError
        No k_y exists or none can be found: F is at or below 1 with no seismic
        force, F stays above 1 up to k = :data:`K_MAX`, the method finds no
        solution, or the surface cannot be built.
    """
    from .search import (AnalysisError, analyzed_surface_noun,
                         noncirc_left_to_right, run_lem_analysis)

    if analysis not in ("auto_search", "single_surface"):
        raise ValueError(f"Unknown analysis {analysis!r}. Use 'auto_search' or "
                         f"'single_surface'.")
    method = str(method or "").lower()
    circular = surface != "noncircular"
    search = analysis == "auto_search"
    noun = analyzed_surface_noun(search, circular)
    base = "circle" if circular else "surface"
    strength = seismic2_strength(seismic2)
    tag = (" (two-stage seismic, K_c = 1 envelope)" if strength == "kc1" else
           " (two-stage seismic)" if strength else "")
    made_under = {"analysis": analysis,
                  "surface": "circular" if circular else "noncircular",
                  "num_slices": num_slices, "rapid": False, "seismic2": seismic2,
                  "ky": True, "composite": composite, "grid_seed": grid_seed,
                  "diagnostic": diagnostic, "fs_tol": fs_tol, "tol": tol,
                  "max_iter": max_iter, "min_slip_depth": min_slip_depth}
    # A vertical seismic coefficient rides with k at the model's own ratio, the way
    # Eurocode 8 states it (k_v = ±0.5 k or ±0.33 k), so k_y is the horizontal
    # coefficient at which F = 1 with k_v = ratio·k_y. The ratio is the entered
    # k_v over the entered |k|; a k_v entered with no k sets no ratio.
    kv_ratio = 0.0
    k_v_in = float(slope_data.get("k_vertical") or 0.0)
    if k_v_in:
        k_h_in = abs(float(slope_data.get("k_seismic") or 0.0))
        if k_h_in == 0.0:
            raise AnalysisError(
                f"The model has a vertical seismic coefficient (k_v = {k_v_in:g}) "
                "and no horizontal one. The yield coefficient search holds k_v at "
                "the ratio k_v/k as it varies k, so enter the horizontal "
                "coefficient that ratio belongs to, or set k_v to 0.")
        kv_ratio = k_v_in / k_h_in
    made_under["kv_ratio"] = kv_ratio

    # The search returns its points left to right; pair each with its own flag.
    movement = [p.get("Movement") for p in noncirc_left_to_right(slope_data.get("non_circ"))]

    def say(text):
        if announce:
            print(text)

    def k_y_on(surface_def, where, hint=None):
        F = _OneSurface(slope_data, method,
                        surface_def if circular else None,
                        None if circular else surface_def,
                        num_slices, composite, seismic2, kv_ratio)
        k, f0 = _k_at_unity(F, where, hint=hint, cancel_check=cancel_check)
        if diagnostic:
            for kk, ff in F.trace:
                say(f"    k = {kk:.5f}: F = "
                    + ("no solution" if ff is None else f"{ff:.5f}"))
        return k, f0, F.trace

    def search_at(k):
        # A non-circular round searches from the last surface found and, once
        # that is no longer the entered one, from the entered surface as well,
        # and keeps the lower (see the docstring).
        starts = [None]
        if not circular:
            starts = [current_nc]
            if current_nc is not slope_data.get("non_circ"):
                starts.append(slope_data.get("non_circ"))
        b = None
        for start in starts:
            sd = dict(slope_data)
            sd["k_seismic"] = k
            sd["k_vertical"] = kv_ratio * k
            if not circular:
                sd["non_circ"] = start
            found = run_lem_analysis(sd, method, analysis="auto_search",
                                     surface=surface, num_slices=num_slices,
                                     seismic2=seismic2, composite=composite,
                                     grid_seed=grid_seed, diagnostic=diagnostic,
                                     cancel_check=cancel_check, fs_tol=fs_tol,
                                     tol=tol, max_iter=max_iter,
                                     min_slip_depth=min_slip_depth,
                                     announce=False)
            if b is None or found["results"]["FS"] < b["results"]["FS"]:
                b = found
        entry = b["search"]["fs_cache"][0]
        # The coefficient the search ran at, so a view that draws it can say so.
        b["search"]["k_seismic"] = k
        return (float(b["results"]["FS"]), _surface_of(entry, circular, movement),
                b["search"])

    static = None
    rounds = []
    converged = True
    if not search:
        if circular:
            if not (slope_data.get("circular") and slope_data.get("circles")):
                raise AnalysisError("A circular surface is required (no circles "
                                    "defined).")
            best_surface = slope_data["circles"][0]
        else:
            if not slope_data.get("non_circ"):
                raise AnalysisError("A non-circular surface is required (no "
                                    "non-circular points defined).")
            best_surface = slope_data["non_circ"]
        say(f"Finding the yield coefficient k_y with {method.upper()} on the "
            f"specified {noun}{tag}…")
        best_k, best_f0, best_trace = k_y_on(best_surface, f"the {noun}")
        best_search = None
    else:
        say(f"Finding the minimum yield coefficient k_y with {method.upper()} by "
            f"search ({'circular' if circular else 'non-circular'}){tag}…")
        current_nc = slope_data.get("non_circ")
        zero = zero_k_phrase({"strength": strength})
        f_static, s0, search0 = search_at(0.0)
        if f_static <= 1.0:
            raise AnalysisError(
                f"The {noun} {zero} has F = {f_static:.3f}. A slope that is not "
                f"stable before the earthquake has no yield coefficient.")
        k0, f0, trace0 = k_y_on(s0, f"the {noun} {zero}")
        static = {"FS": f_static, "k_y": k0, "surface": s0}
        best_surface, best_k, best_f0, best_trace, best_search = (
            s0, k0, f0, trace0, search0)
        current_nc = s0 if not circular else None
        converged = False
        for _ in range(max_rounds):
            k = best_k
            f_r, s_r, search_r = search_at(k)
            entry = {"k": k, "FS": f_r, "k_y": None}
            if f_r < 1.0:
                k_r, f0_r, trace_r = k_y_on(s_r, "the surface the search found",
                                            hint=k)
                entry["k_y"] = k_r
                if k_r < best_k:
                    best_surface, best_k, best_f0, best_trace, best_search = (
                        s_r, k_r, f0_r, trace_r, search_r)
                    if not circular:
                        current_nc = s_r
            rounds.append(entry)
            if diagnostic:
                say(f"  search at k_h = {k:.4f}: lowest F = {f_r:.4f}"
                    + ("" if entry["k_y"] is None
                       else f", k_y = {entry['k_y']:.4f} on that {base}"))
            if f_r >= 1.0 - confirm_tol:
                converged = True
                break

    # The answer as a run like any other: the k_y surface solved at k = k_y,
    # through the same entry point, so the slices, the result, and the views
    # that read them are exactly those of an ordinary run at that coefficient.
    sd = _on_surface(slope_data, best_surface, circular)
    sd["k_seismic"] = best_k
    sd["k_vertical"] = kv_ratio * best_k
    bundle = run_lem_analysis(sd, method, analysis="single_surface",
                              surface=surface, num_slices=num_slices,
                              seismic2=seismic2, composite=composite,
                              announce=False)
    if bundle.get("results") is None:
        raise AnalysisError(f"The method found no solution on the {noun} at "
                            f"k_y = {best_k:.4f}: {bundle.get('failure')}")
    bundle["search"] = best_search
    bundle["options"] = made_under
    bundle["yield"] = {"k_y": best_k, "FS_k0": best_f0, "trace": best_trace,
                       "search": search, "static": static, "rounds": rounds,
                       "converged": converged, "surface": best_surface,
                       "confirm_tol": confirm_tol, "strength": strength,
                       "kv_ratio": kv_ratio}
    say(yield_line(bundle["yield"], noun))
    return bundle


def zero_k_phrase(info):
    """How a yield coefficient result names the state at k = 0. In a two-stage
    run that is stage 2 with its undrained strengths and no seismic force, which
    is not the static (stage 1) factor of safety, so it says so."""
    if (info or {}).get("strength"):
        return "with its stage 2 strengths and no seismic force"
    return "with no seismic force"


def yield_line(info, surface_noun):
    """The lines a yield coefficient result is reported in: the Log's last
    lines, and the note on Studio's LEM Solution tab, which joins them. After a
    search, the first line is k_y on the static critical surface and the second
    is the lowest k_y the search found.

    ``info`` is a bundle's ``"yield"`` record; ``surface_noun`` names the surface
    (:func:`xslope.search.analyzed_surface_noun`).
    """
    zero = zero_k_phrase(info)
    base = surface_noun.split()[-1]
    k_y, f0 = info["k_y"], info["FS_k0"]
    static = info.get("static")
    if not (info.get("search") and static):
        return f"k_y = {k_y:.4f} on the {surface_noun} (F = {f0:.3f} {zero})."
    lines = [f"k_y = {static['k_y']:.4f} on the static critical {base} "
             f"(F = {static['FS']:.3f} {zero})."]
    if abs(static["k_y"] - k_y) < 1e-12:
        lines.append(f"k_y = {k_y:.4f} on the {base} with the lowest yield "
                     f"coefficient, which is the same {base}.")
    else:
        lines.append(f"k_y = {k_y:.4f} on the {base} with the lowest yield "
                     f"coefficient.")
    if not info.get("converged"):
        lines.append(f"The searches did not settle in "
                     f"{len(info.get('rounds') or [])} rounds, so a {base} with "
                     f"a lower k_y may exist.")
    return "\n".join(lines)


def model_at(slope_data, k):
    """The model at horizontal seismic coefficient ``k``, as a run made at a
    coefficient of its own (the yield coefficient, or a confirming search) was
    solved: a vertical coefficient rides with k at the model's own ratio k_v/|k|,
    as :func:`yield_coefficient` holds it. A model with no horizontal coefficient
    sets no ratio and keeps its k_v. A copy; ``slope_data`` is not changed."""
    sd = dict(slope_data)
    k_h = abs(float(slope_data.get("k_seismic") or 0.0))
    k_v = float(slope_data.get("k_vertical") or 0.0)
    sd["k_seismic"] = k
    if k_v and k_h:
        sd["k_vertical"] = k_v / k_h * abs(float(k))
    return sd


# ---------------------------------------------------------------------------
# The yield coefficient by the finite element method
# ---------------------------------------------------------------------------

#: The finite element yield coefficient's bisection stops when the bracket on k
#: is narrower than this.
FEM_K_TOL = 0.005

#: The first coefficient the finite element bisection tries as its failing edge
#: when none is given, doubled until a trial fails or K_MAX is passed.
FEM_K_START = 0.25


def _fem_direction(fem_data, direction):
    """The sign of the seismic body force: -1 drives the model toward -x (a
    left-facing slope), +1 toward +x. The finite element solver applies k
    exactly as signed (see the FEM overview page), so the direction is the
    caller's, or the sign of the model's own coefficient."""
    if direction in ("left", -1, -1.0):
        return -1.0
    if direction in ("right", 1, 1.0):
        return 1.0
    if direction is not None:
        raise ValueError(f"direction must be 'left' or 'right', not {direction!r}")
    k = float(fem_data.get("k_seismic") or 0.0)
    if k == 0.0:
        raise ValueError(
            "The finite element yield coefficient needs the direction of the "
            "seismic force: pass direction='left' for a slope that faces -x or "
            "'right' for one that faces +x, or give the model a signed k.")
    return float(np.sign(k))


def fem_yield_coefficient(fem_data, direction=None, k_stand=0.0, k_fail=None,
                          tol=FEM_K_TOL, k0=1.0, cancel_check=None,
                          progress_callback=None, **solve_kwargs):
    """The seismic yield coefficient k_y by the finite element method: the
    largest seismic coefficient under which the slope stands at full strength.

    A bisection on k. Each trial is one :func:`xslope.fem.solve_fem` at F = 1,
    with the seismic body force k·γ acting in ``direction`` and nothing reduced,
    so a trial asks whether the slope stands under that coefficient. A trial
    starts from the K0 in-situ field, as every strength-reduction trial does,
    and is read by the solver's own verdict (the hybrid criterion unless
    ``failure_criterion`` says otherwise): it stands when the solver reports it
    stable. It is the strength-reduction bisection run on the load instead of
    the strength, and each trial is the full-strength solve a ``fem_tilt``
    verification row locks.

    A trial that ends undecided (an exit in
    :data:`xslope.fem.SSRM_UNDECIDED_EXITS`) is not counted as a failure. As in
    the strength-reduction search it becomes the top of the bracket, and where
    one is still the top when the bracket closes, k_y is reported as a lower
    bound, the bracket's bottom.

    Parameters
    ----------
    fem_data : dict
        From :func:`xslope.fem.build_fem_data`. Its own seismic coefficient is
        replaced on every trial; nothing in it is changed.
    direction : {"left", "right"} or None
        The way the seismic force drives the slope: "left" for a slope whose
        face looks toward -x, "right" for +x. None takes the sign of the
        model's own k.
    k_stand, k_fail : float
        The starting bracket, as magnitudes. ``k_stand`` defaults to 0, the
        static case. With no ``k_fail`` the search tries FEM_K_START and doubles
        it until a trial fails. A ``k_stand`` above 0 is solved and must stand.
    tol : float
        The bisection stops when k_fail - k_stand is below this.
    k0 : float or None
        The at-rest ratio of the in-situ field each trial starts from.
    cancel_check : callable or None
        Called between trials; True cancels the run
        (:class:`xslope.search.AnalysisCancelled`).
    progress_callback : callable or None
        Called as ``progress_callback(trial_number, k, verdict)`` after each
        trial.
    **solve_kwargs
        Passed to every :func:`xslope.fem.solve_fem` trial (``max_iterations``,
        ``tension_srf``, ``failure_criterion``, ``fast_kernel``, ...).

    Returns
    -------
    dict
        - ``k_y``: the midpoint of the final bracket, or its bottom where the
          top trial was undecided.
        - ``k_stand``, ``k_fail``: the final bracket, as magnitudes.
        - ``lower_bound``: True where the top trial was undecided.
        - ``direction``: -1.0 or +1.0.
        - ``trials``: every trial, in order, as dicts with ``k``, ``verdict``,
          ``stable``, ``exit_reason``, ``iterations`` and ``max_displacement``.
        - ``solution``: the solve at ``k_stand`` (None when the bracket's
          bottom is the static case and was never solved), and
          ``failure_solution``: the solve at ``k_fail``.
    """
    from .fem import solve_fem, SSRM_UNDECIDED_EXITS
    from .search import AnalysisError, _check_cancel

    sign = _fem_direction(fem_data, direction)
    trials = []
    solved = {}

    def run(k):
        _check_cancel(cancel_check)
        fd = dict(fem_data)
        fd["k_seismic"] = sign * float(k)
        sol = solve_fem(fd, F=1.0, k0=k0, **solve_kwargs)
        stable = bool(sol.get("stable"))
        undecided = (not stable
                     and sol.get("exit_reason") in SSRM_UNDECIDED_EXITS)
        trials.append({
            "k": float(k),
            "verdict": sol.get("verdict"),
            "stable": stable,
            "undecided": undecided,
            "exit_reason": sol.get("exit_reason"),
            "iterations": int(sol.get("iterations") or 0),
            "max_displacement": float(sol.get("max_displacement") or 0.0),
        })
        solved[float(k)] = sol
        if progress_callback is not None:
            progress_callback(len(trials), float(k), sol.get("verdict"))
        return stable, undecided

    lo = abs(float(k_stand))
    if lo > 0.0:
        stands, _ = run(lo)
        if not stands:
            raise AnalysisError(
                f"The slope does not stand at k_stand = {lo:g}; start the "
                "bracket lower.")
    hi = None if k_fail is None else abs(float(k_fail))
    top_undecided = False
    if hi is None:
        k = max(FEM_K_START, 2.0 * lo)
        while True:
            stands, top_undecided = run(k)
            if not stands:
                hi = k
                break
            lo = k
            if k >= K_MAX:
                raise AnalysisError(
                    f"The slope still stands at k = {k:g}; its yield "
                    f"coefficient is above {K_MAX:g}.")
            k = min(2.0 * k, K_MAX)
    else:
        if not hi > lo:
            raise ValueError(f"k_fail = {hi:g} must be above k_stand = {lo:g}")
        stands, top_undecided = run(hi)
        if stands:
            raise AnalysisError(
                f"The slope stands at k_fail = {hi:g}; start the bracket "
                "higher.")

    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        stands, undecided = run(mid)
        if stands:
            lo = mid
        else:
            hi, top_undecided = mid, undecided

    if lo == 0.0 and 0.0 not in solved:
        stands, _ = run(0.0)
        if not stands:
            raise AnalysisError(
                "The slope does not stand at full strength with no seismic "
                "force, so it has no yield coefficient.")

    k_y = lo if top_undecided else 0.5 * (lo + hi)
    return {"k_y": k_y, "k_stand": lo, "k_fail": hi,
            "lower_bound": top_undecided, "direction": sign,
            "trials": trials, "solution": solved.get(lo),
            "failure_solution": solved.get(hi)}


# ---------------------------------------------------------------------------
# Newmark sliding-block displacement
# ---------------------------------------------------------------------------

#: The acceleration of gravity in each unit system's length per second squared.
GRAVITY = {"si": 9.81, "imperial": 32.174}

#: The unit a Newmark displacement is reported in, and how many of it make one
#: of the model's length units.
DISPLACEMENT_UNIT = {"si": ("cm", 100.0), "imperial": ("in", 12.0)}

#: A number as a record file writes it: a sign, digits with or without a point,
#: and a Fortran exponent (E or D). Matched rather than split, because some
#: fixed-width record files run negative values together with no space between.
_NUMBER = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?")

#: The sample-count line of a PEER .AT2 file: "NPTS=  5862, DT=   .0050 SEC"
#: (NGA-West2), or "5862   .0050   NPTS, DT" (the older PEER format).
_AT2_NEW = re.compile(r"NPTS\s*=\s*(\d+)\s*,?\s*DT\s*=\s*([-+.\dEe]+)", re.I)
_AT2_OLD = re.compile(r"^\s*(\d+)\s+([-+.\dEe]+)\s+NPTS\s*,\s*DT", re.I)


def _float(text):
    return float(text.replace("D", "E").replace("d", "e"))


def _record_out(path, fmt, header, time, accel, description=None):
    time = np.asarray(time, dtype=float)
    accel = np.asarray(accel, dtype=float)
    if len(time) < 2:
        raise ValueError(f"{os.path.basename(path)}: an acceleration record needs "
                         f"at least two samples.")
    steps = np.diff(time)
    if np.any(steps <= 0):
        i = int(np.argmax(steps <= 0)) + 1
        raise ValueError(f"{os.path.basename(path)}: the times must increase; "
                         f"sample {i + 1} is at {time[i]:g} s, after "
                         f"{time[i - 1]:g} s.")
    mean = float(steps.mean())
    uniform = bool(np.all(np.abs(steps - mean) <= 1e-6 * mean))
    header = [h for h in (header or []) if h]
    return {"name": os.path.basename(path), "path": os.path.abspath(path),
            "format": fmt, "header": header,
            "description": description or (header[0] if header else
                                            os.path.basename(path)),
            "time": time, "accel": accel, "npts": len(time),
            "dt": mean if uniform else None,
            "duration": float(time[-1] - time[0]),
            "pga": float(np.max(np.abs(accel)))}


def read_record(path):
    """Load a horizontal acceleration record, in g.

    Two formats are read:

    - a PEER ``.AT2`` file (the NGA databases): header lines, a line stating the
      number of points and the time step (``NPTS=  5862, DT=   .0050 SEC``, or
      the older ``5862  .0050  NPTS, DT``), then the accelerations in g, several
      to a line;
    - a text file of two columns, time in s and acceleration in g, separated by
      spaces, tabs, or commas. Lines before the first pair of numbers are its
      header (a leading ``#`` is dropped), and ``#`` lines after it are skipped.
      The times need not be evenly spaced.

    Returns a dict: ``name``, ``path``, ``format`` (``"AT2"`` or ``"text"``),
    ``header`` (the header lines), ``description`` (the line that names the
    record: an AT2 file's event and station line, else the first header line,
    else the file name), ``time`` and ``accel`` (arrays, s and g), ``npts``,
    ``dt`` (the time step, or None when the steps are uneven), ``duration``,
    and ``pga`` (the largest absolute acceleration, in g).
    """
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    name = os.path.basename(path)

    # --- PEER .AT2 ---
    for i, line in enumerate(lines[:12]):
        m = _AT2_NEW.search(line) or _AT2_OLD.search(line)
        if not m:
            continue
        npts, dt = int(m.group(1)), _float(m.group(2))
        header = [h.strip() for h in lines[:i]]
        units = " ".join(header).upper()
        if re.search(r"UNITS OF\s+(?:CM|MM|M|IN|FT)\b|CM/S|M/S", units):
            raise ValueError(f"{name}: the record is not in g ({units.strip()}); "
                             f"XSLOPE reads accelerations in g.")
        values = [_float(v) for ln in lines[i + 1:] for v in _NUMBER.findall(ln)]
        if len(values) < npts:
            raise ValueError(f"{name}: the header states {npts} points and the "
                             f"file holds {len(values)}.")
        if dt <= 0:
            raise ValueError(f"{name}: the header's time step is {dt:g} s.")
        # Both PEER layouts open with the database's title and name the event,
        # date, station, and component on the second line.
        desc = header[1] if len(header) > 1 else (header[0] if header else name)
        return _record_out(path, "AT2", header, np.arange(npts) * dt,
                           values[:npts], description=desc)

    # --- two columns: time (s), acceleration (g) ---
    header, time, accel = [], [], []
    for n, line in enumerate(lines, start=1):
        text = line.strip()
        if not text:
            continue
        if text.startswith("#"):
            if not time:
                header.append(text.lstrip("#").strip())
            continue
        parts = [p for p in re.split(r"[\s,;]+", text) if p]
        try:
            t, a = _float(parts[0]), _float(parts[1])
        except (IndexError, ValueError):
            if not time:
                header.append(text)
                continue
            raise ValueError(f"{name}, line {n}: expected a time and an "
                             f"acceleration, found {text!r}.")
        time.append(t)
        accel.append(a)
    if not time:
        raise ValueError(f"{name}: no samples found. A record is a PEER .AT2 "
                         f"file, or two columns of time (s) and acceleration (g).")
    return _record_out(path, "text", header, time, accel)


def scale_for_pga(record, target_pga):
    """The scale factor that brings ``record``'s peak acceleration to
    ``target_pga`` (both in g)."""
    pga = float(record["pga"])
    if pga <= 0:
        raise ValueError("The record's peak acceleration is zero.")
    return float(target_pga) / pga


def _slide(time, accel, ky, g, sign):
    """One direction of a rigid-block integration.

    The block slides downslope when ``sign * accel`` exceeds ``ky`` (both in g):
    its acceleration relative to the ground is then ``g * (sign * a - ky)``, and
    once moving it keeps moving until its relative velocity returns to zero. It
    never slides back. The record is taken as linear between samples, and each
    step is integrated exactly on that line, the moments sliding starts and stops
    included, so a record sampled at the corners of a piecewise-linear history
    gives that history's displacement whatever its time step.

    Returns the relative velocity and the displacement at every sample (model
    length per s, and model length), and the number of times the block started
    sliding.
    """
    n = len(time)
    v_out = np.zeros(n)
    d_out = np.zeros(n)
    v = d = 0.0
    sliding = False
    episodes = 0
    a_all = sign * np.asarray(accel, dtype=float)
    for i in range(1, n):
        h = time[i] - time[i - 1]
        a0, a1 = a_all[i - 1], a_all[i]
        if not sliding and a0 <= ky and a1 <= ky:
            v_out[i], d_out[i] = 0.0, d
            continue
        slope = (a1 - a0) / h            # g per s along the step
        s = 0.0                          # time into the step
        while s < h:
            if not sliding:
                a_s = a0 + slope * s
                if a_s > ky:
                    pass                 # sliding from here
                elif a1 > ky:
                    s = (ky - a0) / slope    # where the line crosses ky
                else:
                    break                # stuck for the rest of the step
                sliding = True
                episodes += 1
            # Relative acceleration on [s, h], linear: r(w) = r0 + q w, w = time
            # since s; velocity v + r0 w + q w²/2; the block stops at its first
            # positive root.
            r0 = g * (a0 + slope * s - ky)
            q = g * slope
            span = h - s
            w_stop = None
            A, B, C = 0.5 * q, r0, v
            if abs(A) < 1e-300:
                if B < 0:
                    w_stop = C / -B
            else:
                disc = B * B - 4.0 * A * C
                if disc >= 0.0:
                    root = disc ** 0.5
                    qq = -0.5 * (B + (root if B >= 0 else -root))
                    roots = [qq / A] + ([C / qq] if qq != 0.0 else [])
                    tiny = 1e-12 * h
                    pos = [w for w in roots if w > tiny]
                    if pos:
                        w_stop = min(pos)
            if w_stop is not None and w_stop <= span:
                w = w_stop
                d += v * w + 0.5 * r0 * w * w + q * w ** 3 / 6.0
                v = 0.0
                sliding = False
                s += w
                continue
            w = span
            d += v * w + 0.5 * r0 * w * w + q * w ** 3 / 6.0
            v = v + r0 * w + 0.5 * q * w * w
            s = h
        v_out[i], d_out[i] = v, d
    return v_out, d_out, episodes


def newmark_displacement(record, ky, scale=1.0, units="si"):
    """Permanent displacement of a rigid sliding block (Newmark 1965).

    The sliding mass is a rigid block on the ground. Whenever the ground
    acceleration exceeds the yield acceleration k_y·g it slides downslope, and it
    keeps sliding until its velocity relative to the ground returns to zero
    (Duncan, Wright, and Brandon 2014, Fig. 10.3). It never slides back up.

    Nothing in a record says which way its positive accelerations point on the
    slope being analyzed: its sign is the compass direction of the instrument's
    axis where it was recorded. So the record is run both ways — once with its
    positive accelerations pushing the block downslope, once with its negative
    ones — and the larger displacement is the answer.

    Parameters
    ----------
    record : dict
        From :func:`read_record`, or any dict carrying ``time`` (s) and
        ``accel`` (g).
    ky : float
        The yield coefficient, in g (:func:`yield_coefficient`).
    scale : float
        Multiplies every acceleration (:func:`scale_for_pga` gives the factor
        for a target peak).
    units : {"si", "imperial"}
        The model's unit system: it sets g (9.81 m/s² or 32.174 ft/s²) and the
        unit the displacement is reported in (cm or in).

    Returns
    -------
    dict
        ``D`` (the larger of the two directions) and ``D_unit``; ``governing``
        (+1 or −1, the sign of the accelerations that push the block in the
        governing direction); ``directions``, keyed +1 and −1, each with ``D``,
        ``episodes``, ``v_max``, and the histories ``v`` and ``d`` (in D_unit
        per s and D_unit); ``time`` and ``accel`` (the scaled record, g);
        ``ky``, ``scale``, ``pga`` (scaled, g), ``units``, ``g``; and
        ``record`` (the record's name, path, format, description, npts, dt,
        duration, and unscaled pga).
    """
    if units not in GRAVITY:
        raise ValueError(f"units must be one of {sorted(GRAVITY)}, not {units!r}.")
    ky = float(ky)
    if not ky > 0:
        raise ValueError(f"k_y must be positive; got {ky:g}.")
    scale = float(scale)
    if not scale > 0:
        raise ValueError(f"The scale factor must be positive; got {scale:g}.")
    time = np.asarray(record["time"], dtype=float)
    accel = scale * np.asarray(record["accel"], dtype=float)
    g = GRAVITY[units]
    unit, per = DISPLACEMENT_UNIT[units]
    directions = {}
    for sign in (+1, -1):
        v, d, episodes = _slide(time, accel, ky, g, sign)
        directions[sign] = {"D": float(d[-1] * per), "episodes": episodes,
                            "v_max": float(v.max() * per), "v": v * per,
                            "d": d * per}
    governing = +1 if directions[+1]["D"] >= directions[-1]["D"] else -1
    meta = {k: record.get(k) for k in ("name", "path", "format", "description",
                                       "header", "npts", "dt", "duration", "pga")}
    return {"D": directions[governing]["D"], "D_unit": unit,
            "governing": governing, "directions": directions,
            "time": time, "accel": accel, "ky": ky, "scale": scale,
            "pga": float(np.max(np.abs(accel))), "units": units, "g": g,
            "record": meta}


def direction_phrase(sign):
    """How a result names a direction of the record."""
    return ("its positive accelerations pushing the block downslope" if sign > 0
            else "its negative accelerations pushing the block downslope")


def newmark_line(result):
    """The sentence a Newmark result is reported in: the Log, the results tab,
    and the report's summary."""
    u = result["D_unit"]
    gov = result["governing"]
    other = result["directions"][-gov]
    main = result["directions"][gov]
    return (f"Newmark displacement D = {result['D']:.2f} {u} at k_y = "
            f"{result['ky']:.4f} under {result['record']['name']} "
            f"(peak {result['pga']:.3f} g), with {direction_phrase(gov)} "
            f"({main['episodes']} slip episodes); the other way it is "
            f"{other['D']:.2f} {u}.")
