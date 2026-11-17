"""Seismic slope stability beyond the pseudostatic factor of safety.

The pseudostatic force itself lives in every limit equilibrium method (a
horizontal force kW through each slice's center of gravity, directed by the
failure surface) and the two-stage analysis in :mod:`xslope.advanced`. This module
holds what is computed FROM them:

- :func:`yield_coefficient`: the seismic yield coefficient k_y, the seismic
  coefficient at which the factor of safety is exactly 1 (Duncan, Wright, and
  Brandon 2014, Sec. 10.1.3), on a given surface or by search.

Reference: Duncan, J.M., Wright, S.G., and Brandon, T.L. (2014), *Soil Strength
and Slope Stability*, 2nd ed., Wiley, Chapter 10.
"""
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
CONFIRM_TOL = 0.002

#: The most confirming searches the minimum-k_y search makes.
MAX_ROUNDS = 6


class _OneSurface:
    """F(k) on one fixed surface.

    The slices are cut once. The slicer computes the seismic force of each slice
    as |k| times its soil weight (:mod:`xslope.slice`), so F at any k is a solve
    of the same slices with that column scaled: no re-slicing, and the pore
    pressures, loads, and strengths stay exactly as cut. Every solve starts from
    a fresh copy, because the solvers and the two-stage analysis write into the
    table they are given. ``trace`` keeps every (k, F) the root search passed
    through, which is what the Log and the report show.
    """

    def __init__(self, slope_data, method, circle, non_circ, num_slices,
                 composite, seismic2):
        from .search import AnalysisError
        sd = dict(slope_data)
        sd["k_seismic"] = 1.0
        ok, out = generate_slices(sd, circle=circle, non_circ=non_circ,
                                  num_slices=num_slices, composite=composite,
                                  debug=False)
        if not ok:
            raise AnalysisError(str(out))
        self.df, self.failure_surface = out
        self.kw1 = self.df["kw"].to_numpy(copy=True)
        self.method = method
        self.strength = seismic2_strength(seismic2)
        self.func = getattr(solve, method)
        self.trace = []

    def __call__(self, k):
        """(F, None) at coefficient k, or (None, the solver's reason)."""
        df = self.df.copy()
        df["kw"] = k * self.kw1
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
        ``confirm_tol`` below F = 1.
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
    strength = seismic2_strength(seismic2)
    tag = (" (two-stage seismic, K_c = 1 envelope)" if strength == "kc1" else
           " (two-stage seismic)" if strength else "")
    made_under = {"analysis": analysis,
                  "surface": "circular" if circular else "noncircular",
                  "num_slices": num_slices, "rapid": False, "seismic2": seismic2,
                  "ky": True, "composite": composite, "grid_seed": grid_seed,
                  "diagnostic": diagnostic, "fs_tol": fs_tol, "tol": tol,
                  "max_iter": max_iter, "min_slip_depth": min_slip_depth}
    # The search returns its points left to right; pair each with its own flag.
    movement = [p.get("Movement") for p in noncirc_left_to_right(slope_data.get("non_circ"))]

    def say(text):
        if announce:
            print(text)

    def k_y_on(surface_def, where, hint=None):
        F = _OneSurface(slope_data, method,
                        surface_def if circular else None,
                        None if circular else surface_def,
                        num_slices, composite, seismic2)
        k, f0 = _k_at_unity(F, where, hint=hint, cancel_check=cancel_check)
        if diagnostic:
            for kk, ff in F.trace:
                say(f"    k = {kk:.5f}: F = "
                    + ("no solution" if ff is None else f"{ff:.5f}"))
        return k, f0, F.trace

    def search_at(k):
        sd = dict(slope_data)
        sd["k_seismic"] = k
        if not circular:
            sd["non_circ"] = current_nc
        b = run_lem_analysis(sd, method, analysis="auto_search",
                             surface=surface, num_slices=num_slices,
                             seismic2=seismic2, composite=composite,
                             grid_seed=grid_seed, diagnostic=diagnostic,
                             cancel_check=cancel_check, fs_tol=fs_tol, tol=tol,
                             max_iter=max_iter, min_slip_depth=min_slip_depth,
                             announce=False)
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
        say(f"The {noun} {zero} has F = {f_static:.3f}; on it, "
            f"F = 1 at k = {k0:.4f}.")
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
            if f_r >= 1.0 - confirm_tol:
                say(f"The search at k = {k:.4f} finds F = {f_r:.4f}, within "
                    f"{confirm_tol:g} of 1: k_y has settled.")
                converged = True
                break
            say(f"The search at k = {k:.4f} finds F = {f_r:.4f}; on that "
                f"surface, F = 1 at k = {entry['k_y']:.4f}.")
        if not converged:
            say(f"k_y did not settle in {max_rounds} confirming searches; the "
                f"lowest found is reported.")

    # The answer as a run like any other: the k_y surface solved at k = k_y,
    # through the same entry point, so the slices, the result, and the views
    # that read them are exactly those of an ordinary run at that coefficient.
    sd = _on_surface(slope_data, best_surface, circular)
    sd["k_seismic"] = best_k
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
                       "confirm_tol": confirm_tol, "strength": strength}
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
    """The sentence a yield coefficient result is reported in: the Log's last
    line, the note on Studio's LEM Solution tab, and the report's summary.

    ``info`` is a bundle's ``"yield"`` record; ``surface_noun`` names the surface
    (:func:`xslope.search.analyzed_surface_noun`).
    """
    zero = zero_k_phrase(info)
    text = (f"k_y = {info['k_y']:.3f} on the {surface_noun}, which has "
            f"F = {info['FS_k0']:.3f} {zero}.")
    static = info.get("static")
    if static:
        text += (f" The {surface_noun} {zero} has "
                 f"F = {static['FS']:.3f} and k_y = {static['k_y']:.3f}.")
    if info.get("search") and not info.get("converged"):
        text += (f" The confirming searches did not settle in "
                 f"{len(info.get('rounds') or [])} rounds; this is the lowest "
                 f"k_y they found.")
    return text
