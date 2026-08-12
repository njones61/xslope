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

import time

import numpy as np
from shapely.geometry import LineString, Point

from . import solve
from .advanced import rapid_drawdown, validate_rapid_drawdown
# Re-exported beside the search it feeds: a starting set for a model that carries
# no circles is a search input, and a caller reaching for one is already here.
# The implementation lives in generators.py, which stays free of search's imports
# so the preflight remedy can build a starting set without pulling in the solvers.
from .generators import (generate_noncircular_surface, generate_starting_circles,
                         rank_weak_zones, slope_geometry)
from .preflight import preflight
from .slice import generate_slices, get_y_from_intersection

# A valid limit-equilibrium failure surface must have a meaningful net gravitational
# driving force. A surface under flat ground (e.g. a circle in a reservoir bottom)
# has |sum W sin(alpha)| ~ 0 — no failure mechanism — yet with high pore pressure the
# factor of safety can blow up to a huge negative/near-zero value and win the search
# minimum. Reject surfaces whose |sum W sin(alpha)| is below this fraction of the
# total weight. Valid criticals run 0.2-0.6; degenerate flat circles run ~5e-5.
MIN_DRIVING_FRAC = 0.01

# A surface with more than this fraction of slices in base tension (negative
# effective normal) is rejected as a non-failure surface. Valid criticals run
# 0-10%; the degenerate high-pore-pressure circles run 40-80%.
MAX_BASE_TENSION_FRAC = 0.25

#: Short names for the disclosure sentence, so it reads as a method rather than
#: as a solver key. Anything not named here is capitalized.
_METHOD_LABEL = {
    "oms": "The Ordinary Method of Slices",
    "bishop": "Bishop's method",
    "janbu": "The Janbu method",
    "corps": "The Corps of Engineers method",
    "lowe": "The Lowe & Karafiath method",
    "spencer": "Spencer",
    "mprice": "Morgenstern-Price",
}


def _method_label(method_name):
    return _METHOD_LABEL.get(str(method_name).lower(), str(method_name).capitalize())


class UnsolvedTrials:
    """The trial surfaces a search's method could not solve, kept and reported.

    A search scores an inadmissible trial and an UNSOLVABLE one the same way —
    ``fs_fail``, dropped, gone — and the two are not the same fact. A trial
    rejected on geometry was never a candidate. A trial the solver could not
    answer on WAS a candidate: the search looked at it, got no number, and moved
    on, and if it was the lowest surface in the model the reported minimum is the
    minimum of what the method could solve, not of the model. That is the
    honest-search hole this closes: the counts are collected here and disclosed
    on the search's own output, so the reported minimum is never quoted as global
    while surfaces that outrank it went unanswered.

    The ranking is by the MOMENT factor of safety, both sides — the unsolved
    trials against the reported minimum's own moment answer, never against the
    reported factor of safety, which is a different method's number and cannot
    order surfaces beside Bishop's. It is free where the method already computed
    it (Spencer's Bishop-seeded restart records it, and at phi = 0 its existence
    test knows it in closed form); otherwise it is one Bishop solve per unsolved
    trial, and those are a minority of the sweep.

    Counting is per unique trial SURFACE. The depth optimizer re-evaluates its
    own best depth at every step, so counting solver calls would report a surface
    as several.
    """

    def __init__(self, method_name):
        self.method = str(method_name).lower()
        self.attempted = 0          # unique surfaces handed to the solver
        self.unsolved = 0           # of those, the ones it returned no answer for
        self.no_solution = 0        # of those, the ones that admit no solution
        self.not_converged = 0      # of those, the ones the iteration failed on
        self.inadmissible = 0       # of those, solved but refused on base tension
        self.moment_fs = []         # moment factor of safety of each unsolved trial
        self.reported_moment_fs = None
        self.lower_by_moment = 0
        self._seen = set()

    @staticmethod
    def _key(x, y, d):
        return (round(float(x), 9), round(float(y), 9), round(float(d), 9))

    def record(self, x, y, d, solved, df_slices, message=None):
        """Log one trial surface the solver was given. Repeats are ignored."""
        key = self._key(x, y, d)
        if key in self._seen:
            return
        self._seen.add(key)
        self.attempted += 1
        if solved:
            return
        self.unsolved += 1
        kind = solve.spencer_failure_kind(message)
        setattr(self, kind, getattr(self, kind) + 1)
        fs = self._moment_fs(df_slices)
        if fs is not None:
            self.moment_fs.append(fs)

    @staticmethod
    def _moment_fs(df_slices):
        """This surface's moment factor of safety, reused where it was recorded.

        ``solve.spencer`` writes ``moment_fs`` onto the slice table whenever it
        has computed one — its Bishop-seeded restart, and the phi = 0 existence
        test, which knows it in closed form. Absent that, Bishop is solved here,
        ON A COPY: Bishop writes its own effective normal into ``n_eff``, and the
        table this may be handed is the REPORTED one, whose ``n_eff`` belongs to
        the method that solved it and is printed and plotted as that method's.
        """
        if df_slices is None:
            return None
        fs = df_slices.attrs.get('moment_fs')
        if fs is not None and np.isfinite(fs) and fs > 0:
            return float(fs)
        try:
            ok, res = solve.bishop(df_slices.copy())
        except Exception:
            return None
        if ok and np.isfinite(res.get('FS', np.nan)) and res['FS'] > 0:
            return float(res['FS'])
        return None

    def finish(self, critical_df):
        """Rank the unsolved trials against the reported minimum. Idempotent."""
        if not self.unsolved or not self.moment_fs:
            return
        self.reported_moment_fs = self._moment_fs(critical_df)
        if self.reported_moment_fs is None:
            self.lower_by_moment = 0
            return
        self.lower_by_moment = sum(1 for f in self.moment_fs
                                   if f < self.reported_moment_fs)

    def sentence(self):
        """The disclosure line, or '' when every trial surface was solved."""
        if not self.unsolved:
            return ""
        what = _method_label(self.method)
        parts = [(self.no_solution, "admit no solution"),
                 (self.not_converged, "failed to converge"),
                 (self.inadmissible, "solved only with anomalous base tension")]
        breakdown = ", ".join(f"{n} {t}" for n, t in parts if n)
        line = (f"{what} could not solve {self.unsolved} of {self.attempted} "
                f"trial surfaces ({breakdown})")
        if self.reported_moment_fs is None:
            return line + "; their ranking could not be measured."
        return (f"{line}; {self.lower_by_moment} of them rank lower than the "
                f"reported minimum by the moment measure.")

    def as_dict(self):
        """The counts, for a caller that prints them somewhere else."""
        return {"method": self.method,
                "attempted": self.attempted,
                "unsolved": self.unsolved,
                "no_solution": self.no_solution,
                "not_converged": self.not_converged,
                "inadmissible": self.inadmissible,
                "lower_by_moment": self.lower_by_moment,
                "reported_moment_fs": self.reported_moment_fs,
                "unsolved_moment_fs": list(self.moment_fs),
                "sentence": self.sentence()}


def _with_water_loads_once(slope_data):
    """The model to search, with any automatic water load already derived.

    In automatic mode (``main!D23``) the weight of standing water is not typed by
    the user: the slicer derives it, from the water definition and the ground
    surface, every time it is called. Under a search that is once per TRIAL
    SURFACE -- several thousand times for one answer -- and every one of them
    returns the same blocks, because neither the water line nor the ground
    surface is a thing a search varies. It moves circles.

    So it is derived here instead, once, at the entry to the search, and the
    derived model is what every trial is sliced against.
    :func:`xslope.water.with_water_loads` is idempotent -- a model that already
    carries the derivation is returned by identity -- so the slicer's own
    defensive call becomes a dictionary lookup and the loads reaching the slices
    are the same objects they always were.

    The derivation itself is sub-millisecond, and it is not what this saves. In
    automatic mode ``with_water_loads`` hands back a shallow COPY, so under the
    old arrangement every trial was sliced against a model dict that had existed
    for a few milliseconds -- and anything the slicer memoizes INTO the model it
    is given went into the throwaway with it. The expensive one is
    ``_water_table_profile``, the u = 0 contour traced through a seepage mesh for
    the gamma/gamma_sat weight split: 201 stations, each bisected against the
    field. On a dam with saturated unit weights and a seepage solution that is
    ~6,000 mesh point-locations, rebuilt per trial. Deriving once means the model
    survives the search, and the profile is traced once with it: measured on the
    Rocscience VP72 dam, a full search fell from 349 s to 20 s with the reported
    factor of safety and critical circle bit-identical.

    Deliberately NOT a cache. Nothing is memoized and nothing survives the call:
    the reuse is scoped to one search, over a model that cannot change while the
    search is running, which is the only window in which reuse is provably safe.
    Studio edits the water definition between runs -- a piezometric line dragged,
    a reservoir level retyped -- and a cache keyed on the model would answer the
    next run with the last run's reservoir.

    In manual mode this returns the caller's own model unchanged, by identity.
    """
    from .water import with_water_loads
    return with_water_loads(slope_data)


def _net_driving_too_small(df_slices):
    """True if the surface has negligible net gravitational driving force (a flat,
    non-failure surface) and should be rejected as a search candidate."""
    # Degenerate circles (e.g. grazing a stepped wall face) can yield an empty
    # slice table with no columns: reject those outright.
    if df_slices is None or len(df_slices) == 0 or 'w' not in df_slices.columns:
        return True
    W = df_slices['w'].values
    total = float(np.abs(W).sum())
    if total <= 0:
        return True
    driving = abs(float((W * np.sin(np.radians(df_slices['alpha'].values))).sum()))
    return driving < MIN_DRIVING_FRAC * total


MIN_THICKNESS_FRAC = 0.05


def _too_thin(df_slices, h_ground, frac):
    """True if the sliding mass is thinner than ``frac`` of the slope height.

    GRID MODE ONLY (frac=0 disables). A global sweep of a benched geometry finds
    surficial skins hugging any steep face — on a c~0 face they undercut every
    real mechanism. They are usually a raveling mode, not the slope failure the
    search is after, and Slide2 suppresses them the same way (its 'minimum
    depth' filter). But thin criticals can also be THE answer (a submerged c=0
    slope fails as an infinite-slope skin at zero depth), which is why this must
    never apply to the default seeded search: there the user's circles define
    the mechanism of interest, thin or not."""
    if frac <= 0:
        return False
    if df_slices is None or len(df_slices) == 0 or 'y_ct' not in df_slices.columns:
        return True
    thick = float((df_slices['y_ct'] - df_slices['y_cb']).max())
    return thick < frac * h_ground


def _too_shallow(df_slices, min_slip_depth):
    """True if the deepest point of the slip surface lies less than
    ``min_slip_depth`` below the ground surface — an ABSOLUTE minimum-slip-depth
    filter (Slide2's 'minimum depth'), and the search-side twin of solve_fem's
    ``min_slip_depth`` in SSRM. ``None`` or a non-positive value disables it.

    Unlike ``_too_thin`` (a fraction of slope height, grid-mode only), this
    applies to EVERY trial when set, so a caller can suppress surficial
    cohesionless skins (FS = tan phi / tan beta) on demand while leaving the
    default seeded search untouched."""
    if min_slip_depth is None or min_slip_depth <= 0:
        return False
    if df_slices is None or len(df_slices) == 0 or 'y_ct' not in df_slices.columns:
        return True
    max_depth = float((df_slices['y_ct'] - df_slices['y_cb']).max())
    return max_depth < min_slip_depth


def _base_tension_too_extensive(df_slices):
    """True if a large fraction of slices carry a negative effective base normal
    (base in tension) — a non-physical surface, e.g. a circle through a
    high-pore-pressure zone where the moment methods then return a negative/near-zero
    factor of safety. Valid criticals run ~0%; the degenerate cases run 60-80%. A
    few negative base normals are normal and are NOT rejected (extent only)."""
    if 'n_eff' not in df_slices.columns:
        return False
    N = df_slices['n_eff'].values
    return len(N) > 0 and float((N < 0).mean()) > MAX_BASE_TENSION_FRAC


# === Optional search-window constraints (Slide2 "search limits" semantics) ===
# All three default to None, which is exactly the historical unconstrained search
# (the helpers below short-circuit to True and add no behaviour). When set, they
# confine the adaptive circular search to a region of the geometry so it settles on
# a chosen local minimum instead of the global one — the LEM analog of RS2's SSR
# Polygon Search Area / Slide2's entry-and-exit and slip-center limits.


def _center_in_box(x, y, center_box):
    """True if the trial circle CENTER ``(x, y)`` lies within ``center_box`` —
    an ``(x1, y1, x2, y2)`` rectangle whose corners may be given in any order.
    ``None`` disables the box (every center passes)."""
    if center_box is None:
        return True
    x1, y1, x2, y2 = center_box
    return (min(x1, x2) <= x <= max(x1, x2)) and (min(y1, y2) <= y <= max(y1, y2))


def _tangent_depth_ok(depth, tangent_depth):
    """True if the circle's lowest point — elevation ``depth`` (= Yo - R, the
    tangent line of the arc) — lies within the ``tangent_depth`` = ``(ymin, ymax)``
    elevation band. ``None`` disables the band."""
    if tangent_depth is None:
        return True
    ymin, ymax = tangent_depth
    return min(ymin, ymax) <= depth <= max(ymin, ymax)


def _endpoints_in_ranges(failure_surface, entry_range, exit_range):
    """True if the failure surface's two trace endpoints satisfy the entry/exit
    x-range limits. Following Slide2's convention the ENTRY endpoint is the
    crest-side (higher-ground) trace point and the EXIT endpoint is the toe-side
    (lower-ground) trace point, so the limits are independent of slope facing.
    Either range may be ``None`` (that end is unconstrained); both ``None`` passes
    every surface."""
    if entry_range is None and exit_range is None:
        return True
    if failure_surface is None:
        return False
    coords = list(failure_surface.coords)
    if len(coords) < 2:
        return False
    (xl, yl), (xr, yr) = coords[0], coords[-1]
    # crest side = the endpoint sitting on higher ground; toe side = the lower one.
    if yl >= yr:
        entry_x, exit_x = xl, xr
    else:
        entry_x, exit_x = xr, xl
    if entry_range is not None:
        a, b = entry_range
        if not (min(a, b) <= entry_x <= max(a, b)):
            return False
    if exit_range is not None:
        a, b = exit_range
        if not (min(a, b) <= exit_x <= max(a, b)):
            return False
    return True


def _validate_search_bounds(center_box, entry_range, exit_range, tangent_depth):
    """Shape-check the optional search-window tuples, raising a clear ValueError on
    a malformed bound. ``None`` bounds are left untouched."""
    if center_box is not None and len(tuple(center_box)) != 4:
        raise ValueError("center_box must be a 4-tuple (x1, y1, x2, y2)")
    for name, rng in (("entry_range", entry_range), ("exit_range", exit_range),
                      ("tangent_depth", tangent_depth)):
        if rng is not None and len(tuple(rng)) != 2:
            raise ValueError(f"{name} must be a 2-tuple (lo, hi)")


class AnalysisCancelled(Exception):
    """Raised by a search/reliability run when its ``cancel_check`` asks it to stop.

    Cooperative cancellation: long loops call ``cancel_check()`` at iteration
    boundaries and raise this if it returns True, so a caller (e.g. the GUI's
    worker thread) can abort cleanly without killing the thread mid-computation."""


class AnalysisError(Exception):
    """Raised when an analysis cannot be run at all, with the reason as its message.

    The reason is written for the person who asked for the run — "A circular
    surface is required (no circles defined)." — because that is where it is
    shown: Studio puts it in the Run LEM box, and the report prints it where the
    method would have been documented. Not a failure to converge, which is an
    answer about the model; this is the run never having happened.
    """


def _check_cancel(cancel_check):
    if cancel_check is not None and cancel_check():
        raise AnalysisCancelled()


def _grid_seed_circles(slope_data, method_name, num_slices=20, fs_fail=9999,
                       rapid=False, composite=False, nx=10, ny=5, n_tangents=6,
                       keep=4, diagnostic=False, cancel_check=None, circle_cache=None,
                       min_slip_depth=None, center_box=None, entry_range=None,
                       exit_range=None, tangent_depth=None, tally=None):
    """Coarse grid-and-tangent sweep that SEEDS the adaptive circular search.

    The adaptive 9-point search is a local optimizer: it refines whatever
    neighborhood its starting circles put it in, and a single seed in the wrong
    family converges to a local minimum with no warning (a shallow circle in a
    fill can read 20%+ high while the true critical rides the base of the
    foundation). This sweep is the global stage that protects against that.

    It is Slide2's grid-and-tangent model. Centers are laid out on an nx-by-ny
    grid over a box derived from the slope geometry — spanning the slope
    horizontally (extended 0.75H past the toe and crest) and sitting 0.25H to
    2H above the crest, which brackets where the center of a circle that cuts
    the slope can be. For each center, one trial circle is made TANGENT to each
    of ``n_tangents`` elevations between the bottom of the model and mid-slope.
    Sweeping tangent elevations, not depths-below-center, is what finds
    competing families: the shallow face circle and the deep base-tangent
    circle live at similar centers but different tangent lines, and a local
    depth optimizer started in one family never sees the other.

    Every trial is solved (at a reduced slice count) with the same degenerate-
    surface guards as the search proper. The best circle from each tangent
    family within 25% of the global best FS is returned as a seed, deepest
    families first among ties, capped at ``keep``. The caller merges these with
    any user-entered circles and runs the normal adaptive search from all of
    them.

    ``tally`` is the caller's :class:`UnsolvedTrials`, if it keeps one: the sweep
    hands the solver thousands of surfaces and its failures belong in the same
    disclosure as the refinement's.
    """
    solver = getattr(solve, method_name)
    coords = list(slope_data['ground_surface'].coords)
    ys_g = [y for _, y in coords]
    y_hi, y_lo = max(ys_g), min(ys_g)
    H = y_hi - y_lo
    if H <= 0:
        return []
    right_facing = coords[0][1] > coords[-1][1]
    tol_y = 1e-6 * max(1.0, H)
    his = [x for x, y in coords if y > y_hi - tol_y]
    los = [x for x, y in coords if y < y_lo + tol_y]
    x_crest = max(his) if right_facing else min(his)
    x_toe = min(los) if right_facing else max(los)

    # The box must scale with the FULL model depth, not just the slope height: a
    # circle tangent to the base of a deep foundation spans the whole model, and
    # its center sits several slope-heights above the crest (VP75's critical
    # centers run to ~2.5x the crest-to-floor depth above the toe).
    depth_floor = slope_data['domain_polygon'].bounds[1]
    Hm = max(H, y_hi - depth_floor)
    xs = np.linspace(min(x_toe, x_crest) - 0.75 * Hm, max(x_toe, x_crest) + 0.75 * Hm, nx)
    ys = np.linspace(y_hi + 0.25 * Hm, y_hi + 2.5 * Hm, ny)
    y_mid = 0.5 * (y_lo + y_hi)
    tangents = list(np.linspace(depth_floor, y_mid, n_tangents))
    if composite:
        # one below-floor family: its circles truncate at the base and run along it
        tangents.insert(0, depth_floor - 0.5 * H)

    if rapid:
        from .advanced import rapid_drawdown

    sweep_slices = min(num_slices, 20)
    best_by_tan = {}
    n_valid = 0
    for yt in tangents:
        for yc in ys:
            _check_cancel(cancel_check)
            for xc in xs:
                R = yc - yt
                if R <= 0:
                    continue
                # Honour the optional search window: centers outside center_box and
                # tangent lines outside tangent_depth are not candidates (skipped,
                # never clamped). None bounds pass everything -> unchanged sweep.
                if not _center_in_box(xc, yc, center_box):
                    continue
                if not _tangent_depth_ok(yt, tangent_depth):
                    continue
                circle = {'Xo': float(xc), 'Yo': float(yc), 'Depth': float(yt), 'R': float(R)}
                success, result = generate_slices(slope_data, circle=circle,
                                                  num_slices=sweep_slices,
                                                  composite=composite,
                                                  check_inputs=False)
                if not success:
                    continue
                df_slices, failure_surface = result
                if (_net_driving_too_small(df_slices) or _too_thin(df_slices, H, MIN_THICKNESS_FRAC)
                        or _too_shallow(df_slices, min_slip_depth)
                        or not _endpoints_in_ranges(failure_surface, entry_range, exit_range)):
                    continue
                if rapid:
                    ok, solver_result = rapid_drawdown(df_slices, method_name, debug_level=0)
                else:
                    ok, solver_result = solver(df_slices)
                if tally is not None:
                    tally.record(circle['Xo'], circle['Yo'], circle['Depth'],
                                 ok, df_slices, None if ok else solver_result)
                if not ok:
                    continue
                FS = solver_result['FS']
                if not (0 < FS < fs_fail) or _base_tension_too_extensive(df_slices):
                    continue
                n_valid += 1
                if circle_cache is not None:
                    circle_cache.append({"Xo": circle['Xo'], "Yo": circle['Yo'],
                                         "Depth": circle['Depth'], "R": circle['R'],
                                         "FS": FS, "failure_surface": failure_surface})
                if yt not in best_by_tan or FS < best_by_tan[yt][0]:
                    best_by_tan[yt] = (FS, circle)

    if not best_by_tan:
        return []
    fs_best = min(fs for fs, _ in best_by_tan.values())
    ranked = sorted(best_by_tan.items(), key=lambda kv: kv[1][0])
    seeds = [circ for yt, (fs, circ) in ranked if fs <= 1.25 * fs_best][:keep]
    print(f"[grid seed] swept {len(tangents) * nx * ny} circles ({n_valid} valid), "
          f"best FS={fs_best:.4f}; seeding {len(seeds)} tangent "
          f"famil{'y' if len(seeds) == 1 else 'ies'}")
    if diagnostic:
        for fs, c in [v for _, v in ranked]:
            tag = 'seed' if c in seeds else '    '
            print(f"  [{tag}] tangent y={c['Depth']:8.2f}  FS={fs:7.4f}  "
                  f"center=({c['Xo']:.1f}, {c['Yo']:.1f})")
    return seeds


# ---------------------------------------------------------------------------
# The model's own search window
# ---------------------------------------------------------------------------

# The circles-sheet window keys (v19) that this module turns into search kwargs.
# Named here so the reading below and any caller inspecting a model talk about
# the same set.
SEARCH_WINDOW_KEYS = ('entry_x_min', 'entry_x_max', 'exit_x_min', 'exit_x_max',
                      'center_box_x_min', 'center_box_y_min',
                      'center_box_x_max', 'center_box_y_max',
                      'max_tangent_depth', 'min_slip_depth')


def file_search_window(slope_data, already=()):
    """The model's own search window as :func:`circular_search` keyword arguments.

    A model may declare a search window on its circles sheet (v19,
    ``slope_data['search_window']``) — entry and exit ranges, a center box, a
    maximum tangent depth, a minimum slip depth. It is how an engineer says which
    mechanism is under study, and every path that searches on the model's behalf
    reads it here so that a windowed model gets the same surface family from a
    sweep, a reliability run and Studio's Run LEM button as it does from a single
    scripted search.

    A limit is forwarded only when the file supplies BOTH ends of it — a
    half-declared range is not a window, and inventing the missing end would
    constrain the search in a direction nobody asked for. Anything the caller has
    already set (``already``, any container supporting ``in``: the kwargs dict
    being built, or a sequence of names) wins and is left alone.

    "Max tangent depth" is the one limit that is not stored as a pair: it is a
    single ELEVATION, the deepest a circle's tangent point may reach.
    ``circular_search`` takes a band, so it is paired with the top of the ground
    surface — the only upper end that adds no constraint of its own — and dropped
    when the model has no ground surface to read or the declared floor is above
    it.

    Returns a dict of ``circular_search`` kwargs, empty when the file declares no
    usable window. Pass it through :func:`noncircular_search_opts` before handing
    it to a non-circular search, which understands only one of these limits.
    """
    sw = (slope_data.get('search_window') or {}) if isinstance(slope_data, dict) else {}
    if not sw:
        return {}
    kw = {}

    def pair(name, lo, hi):
        if name in already:
            return
        if sw.get(lo) is not None and sw.get(hi) is not None:
            kw[name] = (float(sw[lo]), float(sw[hi]))

    pair('entry_range', 'entry_x_min', 'entry_x_max')
    pair('exit_range', 'exit_x_min', 'exit_x_max')
    box = ('center_box_x_min', 'center_box_y_min',
           'center_box_x_max', 'center_box_y_max')
    if 'center_box' not in already and all(sw.get(k) is not None for k in box):
        kw['center_box'] = tuple(float(sw[k]) for k in box)
    if 'tangent_depth' not in already and sw.get('max_tangent_depth') is not None:
        gs = slope_data.get('ground_surface')
        try:
            y_top = max(y for _x, y in gs.coords)
        except Exception:                                  # noqa: BLE001
            y_top = None
        if y_top is not None and y_top > float(sw['max_tangent_depth']):
            kw['tangent_depth'] = (float(sw['max_tangent_depth']), float(y_top))
    if 'min_slip_depth' not in already and sw.get('min_slip_depth') is not None:
        kw['min_slip_depth'] = float(sw['min_slip_depth'])
    return kw


# The only window limit a non-circular search has a notion of. Entry, exit,
# center and tangent limits are all shapes of a CIRCLE, and there is no circle to
# constrain on that branch.
NONCIRCULAR_WINDOW_KEYS = ('min_slip_depth',)


def noncircular_search_opts(search_opts):
    """The subset of a search window that :func:`noncircular_search` accepts.

    A non-circular search takes ``min_slip_depth`` and has no notion of the
    circle-shaped limits, so those are DROPPED rather than raised on — the same
    line Studio draws between the two branches. Returns a new dict; the caller's
    is untouched.
    """
    return {k: v for k, v in dict(search_opts or {}).items()
            if k in NONCIRCULAR_WINDOW_KEYS}


def circular_search(slope_data, method_name, rapid=False, tol=1e-2, fs_tol=5e-4, max_iter=50,
                    shrink_factor=0.5, fs_fail=9999, min_grid_frac=0.03, depth_tol_frac=0.03,
                    diagnostic=False, num_slices=40, cancel_check=None, composite=False,
                    seed='circles', min_slip_depth=None, center_box=None, entry_range=None,
                    exit_range=None, tangent_depth=None, unsolved_out=None):
    """
    Global 9-point circular search with adaptive grid refinement.

    Convergence is driven by the factor of safety, not the grid spacing: the
    center grid keeps refining until the best FS stops improving — specifically
    until two successive refinement levels each gain less than ``fs_tol`` (so the
    reported FS is stable to ~3 decimal places). FS convergence is only accepted
    once the center grid has refined below ``min_grid_frac`` of the domain height,
    so a coarse-grid FS plateau cannot be mistaken for the true minimum. This grid
    gate is deliberately loose (3% of the domain height): near the critical circle
    the FS surface is flat, so once the FS has plateaued the center is already
    located well enough — tightening the gate only spends extra refinement levels
    that don't move the reported FS. ``tol`` is a geometric backstop on the grid
    spacing and ``max_iter`` caps the count. Keying the stop on FS (not an absolute
    length like ``grid < 0.01``, which means different things on a 20 ft vs a
    500 ft slope) makes it scale-invariant.

    ``seed`` selects the starting circles. ``'circles'`` (default) refines from
    the circles sheet, exactly as before. ``'grid'`` first runs a coarse
    grid-and-tangent sweep over a box derived from the slope geometry and seeds
    the refinement with the best circle from each competing tangent family (plus
    any user circles) — see _grid_seed_circles. Use it when a single starting
    circle risks trapping the search in a local minimum, e.g. an embankment on a
    layered foundation where a shallow fill circle and a deep foundation circle
    compete.

    ``composite`` (default False) lets trial circles dip below the bottom of the
    model, where they are truncated into a composite surface that runs along the
    base (see slice.CompositeSurface). Turn it on when the base is a real
    impenetrable boundary with a weak seam or soft layer on it — the critical
    mechanism there rides the base, and a search clamped to the floor can only ever
    reach a circle tangent to it. Leave it off when the floor is just ``max_depth``,
    a search bound rather than bedrock.

    ``min_slip_depth`` (default None = off) is an absolute minimum-slip-depth filter:
    trial surfaces whose deepest point is less than this far below the ground surface
    are rejected, so a shallow surficial skin (FS = tan phi / tan beta on a
    cohesionless face) cannot win. This is Slide2's "minimum depth" filter and the
    search-side twin of solve_fem's ``min_slip_depth`` in SSRM — set the same value in
    both to keep an LEM/SSRM comparison on the same mechanism. Off by default, so the
    reported surface is the true minimum unless the caller opts in.

    ``center_box``, ``entry_range``, ``exit_range`` and ``tangent_depth`` are optional
    SEARCH-WINDOW constraints, all default None = today's unconstrained search
    (byte-for-byte identical). They confine the adaptive refinement to a region of the
    geometry so it settles on a chosen LOCAL minimum rather than the global one — the
    LEM analog of RS2's SSR Polygon Search Area and Slide2's slip-center / entry-and-
    exit limits (a benched slope has several competing minima; see RS2-61):

      * ``center_box=(x1, y1, x2, y2)`` confines candidate circle CENTERS to a
        rectangle (corners in any order). The refined grid stays inside the box —
        grid points that fall outside are dropped, so the search cannot walk out of
        it. A starting circle whose center lies outside the box is clamped to the box
        for its launch only (a seed, not a reported candidate).
      * ``entry_range=(xa, xb)`` / ``exit_range=(xc, xd)`` confine the surface trace
        endpoints, in x, to the given ranges. ENTRY is the crest-side (higher-ground)
        endpoint and EXIT the toe-side (lower-ground) one, independent of facing. A
        trial whose endpoints fall outside is REJECTED (scored fs_fail), never
        clamped, so the reported minimum genuinely honours the window.
      * ``tangent_depth=(ymin, ymax)`` confines the circle's lowest point (its tangent
        elevation) to an elevation band; out-of-band trials are rejected.

    All four survive the grid refinement and the grid-seed sweep. A malformed bound
    raises ValueError.

    ``unsolved_out``, when a dict is passed, receives the counts of the trial
    surfaces THE METHOD COULD NOT SOLVE (see :class:`UnsolvedTrials`), which the
    search also prints as one line when there are any. A search whose method
    answered on every admissible trial prints nothing extra and fills the dict
    with zeros, so clean searches read exactly as they always have.

    Returns:
        list of dict: sorted fs_cache by FS
        bool: convergence flag
        list of dict: search path
        list of dict: circle_cache - all circles tested during search
    """

    _validate_search_bounds(center_box, entry_range, exit_range, tangent_depth)

    if seed != 'grid' and not slope_data.get('circles'):
        raise ValueError(
            "Circular search requires at least one circle defined in the input. "
            "Add circle data to the 'circles' sheet in the input template, "
            "or run with seed='grid' to generate starting circles automatically.")

    # Preflight once, here, rather than once per trial surface: the model does not
    # change between trials, so a per-trial check would be the same answer several
    # thousand times over -- and a refusal raised inside the sweep would be reported
    # as "no viable surface", blaming geometry for a missing input.
    # A rapid-drawdown search is checked as a rapid run: it inherits every
    # limit-equilibrium rule and adds the stage-2 water, load and strength
    # requirements, which is exactly what validate_rapid_drawdown below only
    # PRINTS -- and only on this path.
    preflight(slope_data, 'rapid' if rapid else 'lem',
              {'surface': 'circular', 'search': True,
               'method': method_name}).raise_for_errors()

    if rapid:
        validate_rapid_drawdown(slope_data)

    slope_data = _with_water_loads_once(slope_data)

    solver = getattr(solve, method_name)
    circle_cache = []  # Store ALL circles tested for plotting
    tally = UnsolvedTrials(method_name)

    def _report_unsolved(critical_df):
        """Close the tally, disclose it, and hand the counts to the caller."""
        tally.finish(critical_df)
        if tally.unsolved:
            print(f"[⚠️ unsolved trials] {tally.sentence()}")
        if unsolved_out is not None:
            unsolved_out.update(tally.as_dict())

    start_time = time.time()  # Start timing

    ground_surface = slope_data['ground_surface']
    ground_surface = LineString([(x, y) for x, y in ground_surface.coords])
    y_max = max(y for _, y in ground_surface.coords)
    h_ground = y_max - min(y for _, y in ground_surface.coords)   # slope height, for the skin-slide guard
    # The deepest allowable failure-surface elevation is the bottom of the domain
    # polygon (not max_depth, which is only a profile-sheet input). Circles below an
    # irregular bottom are rejected by the containment check in generate_slices.
    depth_floor = slope_data['domain_polygon'].bounds[1]
    y_min = depth_floor
    delta_y = y_max - y_min

    # composite=True lets trial circles dip BELOW the floor, where generate_slices
    # truncates them into a composite surface that runs along the floor (see
    # slice.CompositeSurface). This is what a weak seam or a soft layer sitting on
    # bedrock actually fails along, and a search clamped to the floor can never find
    # it: the best it can do is a circle tangent to the base. Off by default,
    # because on a profile-line model the floor is `max_depth` — a search bound the
    # user picked, not a real impenetrable boundary, and truncating there would be
    # meaningless. The lower bound is generous; a circle that goes too deep flattens
    # out along the floor and its FS climbs again, so the optimizer turns back.
    search_floor = depth_floor - delta_y if composite else depth_floor
    min_grid = delta_y * min_grid_frac   # required center-grid resolution before FS convergence

    # The skin-slide filter runs only in grid mode: a global sweep must not report
    # a raveling sliver as the critical, but the seeded search must stay able to
    # follow a user's circles to a genuinely thin critical (e.g. a submerged c=0
    # slope, whose infinite-slope skin IS the answer). See _too_thin.
    thin_frac = MIN_THICKNESS_FRAC if seed == 'grid' else 0.0

    if seed == 'grid':
        # Global stage: sweep a geometry-derived center grid against tangent
        # elevations and seed the adaptive search with the best circle from each
        # competing family (see _grid_seed_circles). User circles, if any, stay in
        # as extra seeds. This is the protection against the local-minimum trap —
        # a single user seed in the wrong family otherwise converges 20%+ high
        # with no warning.
        _check_cancel(cancel_check)
        circles = _grid_seed_circles(
            slope_data, method_name, num_slices=num_slices, fs_fail=fs_fail,
            rapid=rapid, composite=composite, diagnostic=diagnostic,
            cancel_check=cancel_check, circle_cache=circle_cache,
            min_slip_depth=min_slip_depth, center_box=center_box,
            entry_range=entry_range, exit_range=exit_range, tangent_depth=tangent_depth,
            tally=tally)
        circles = circles + list(slope_data.get('circles') or [])
        if not circles:
            _report_unsolved(None)
            return [], False, [], circle_cache
    else:
        circles = slope_data['circles']

    def optimize_depth(x, y, depth_guess, depth_step_init, depth_shrink_factor, tol_frac, fs_fail, circle_cache, diagnostic=False):
        depth_step = min(10.0, depth_step_init)
        # The floor bounds the SEARCH RANGE, not a geometrically valid seed's own
        # depth. `Depth` is the elevation of the circle's LOWEST POINT — its nadir at
        # x = Xo — and for a SKIMMING circle (a very large radius whose arc runs nearly
        # parallel to the face) that nadir sits far outside the model, hundreds of
        # metres below the domain floor, while the failure SURFACE it cuts stays
        # entirely inside. What decides whether a trial is admissible is
        # generate_slices' containment test on that surface, never the nadir — so
        # clamping the seed to the floor destroyed every skimming seed handed to a
        # search. On Talbingo the generated skim solves DIRECTLY to FS = 1.6696
        # against the face-parallel closed form tan(phi)/tan(beta) = 1.6687, but the
        # search seeded with it returned 1.8135, because the very first refinement
        # pulled Depth from -54.1 up to the floor at 0.0 and turned the skin into an
        # ordinary toe circle it could never climb back out of.
        #
        # So a seed that actually slices keeps its own depth, and the probe range
        # simply extends down to it. The extra slice is only ever built for a seed
        # BELOW the floor, which on an ordinary model never happens.
        floor = search_floor
        if depth_guess < search_floor:
            seed_ok, _seed_res = generate_slices(
                slope_data,
                circle={'Xo': x, 'Yo': y, 'Depth': depth_guess, 'R': y - depth_guess},
                num_slices=num_slices, composite=composite, check_inputs=False)
            if seed_ok:
                floor = depth_guess
        best_depth = max(depth_guess, floor)
        best_fs = fs_fail
        # These three are only rebound inside the loop below, so they must be seeded
        # here: a degenerate depth_step (<= depth_tol on entry) skips the loop
        # entirely, and the return would otherwise raise UnboundLocalError instead of
        # reporting an honest fs_fail.
        best_df = None
        best_surface = None
        best_solver_result = None
        depth_tol = depth_step * tol_frac
        iterations = 0

        while depth_step > depth_tol:
            depths = [
                max(best_depth - depth_step, floor),
                best_depth,
                best_depth + depth_step
            ]
            fs_results = []
            for d in depths:
                test_circle = {'Xo': x, 'Yo': y, 'Depth': d, 'R': y - d}
                success, result = generate_slices(slope_data, circle=test_circle,
                                                  num_slices=num_slices, composite=composite,
                                                  check_inputs=False)
                if not success:
                    FS = fs_fail
                    df_slices = None
                    failure_surface = None
                    solver_result = None
                else:
                    df_slices, failure_surface = result
                    if (_net_driving_too_small(df_slices) or _too_thin(df_slices, h_ground, thin_frac)
                            or _too_shallow(df_slices, min_slip_depth)
                            or not _tangent_depth_ok(d, tangent_depth)
                            or not _endpoints_in_ranges(failure_surface, entry_range, exit_range)):
                        # degenerate (near-zero driving; grid mode also drops skins),
                        # or the trial falls outside the optional search window
                        FS = fs_fail
                        solver_result = None
                    else:
                        if rapid:
                            solver_success, solver_result = rapid_drawdown(df_slices, method_name, debug_level=0)
                        else:
                            solver_success, solver_result = solver(df_slices)
                        # The two ways a trial scores fs_fail are different facts,
                        # and only this one is about the METHOD: the surface was
                        # admissible and the solver returned no answer on it. It
                        # still scores fs_fail — the search is unchanged — but it
                        # is counted, and disclosed at the end.
                        tally.record(x, y, d, solver_success, df_slices,
                                     None if solver_success else solver_result)
                        FS = solver_result['FS'] if solver_success else fs_fail
                        if solver_success and _base_tension_too_extensive(df_slices):
                            FS = fs_fail  # degenerate surface (base mostly in tension)
                        if not (FS > 0):  # reject non-positive / NaN factor of safety
                            FS = fs_fail
                fs_results.append((FS, d, df_slices, failure_surface, solver_result))

                # Add to circle_cache for plotting all tested circles
                if FS != fs_fail:
                    circle_cache.append({
                        "Xo": x,
                        "Yo": y,
                        "Depth": d,
                        "R": y - d,
                        "FS": FS,
                        "failure_surface": failure_surface
                    })

            fs_results.sort(key=lambda t: t[0])
            best_fs, best_depth, best_df, best_surface, best_solver_result = fs_results[0]

            if all(FS == fs_fail for FS, *_ in fs_results):
                if diagnostic:
                    print(f"[❌ all fail] x={x:.2f}, y={y:.2f}")
                return best_depth, fs_fail, None, None, None

            if diagnostic:
                print(f"[✓ depth-opt] x={x:.2f}, y={y:.2f}, depth={best_depth:.2f}, FS={best_fs:.4f}, step={depth_step:.2f}")

            depth_step *= depth_shrink_factor
            iterations += 1
            if iterations > 50:
                if diagnostic:
                    print(f"[⚠️ warning] depth iterations exceeded at (x={x:.2f}, y={y:.2f})")
                break

        return best_depth, best_fs, best_df, best_surface, best_solver_result

    def evaluate_grid(x0, y0, grid_size, depth_guess, slope_data, diagnostic=False, fs_cache=None, circle_cache=None):
        if fs_cache is None:
            fs_cache = {}

        Xs = [x0 - grid_size, x0, x0 + grid_size]
        Ys = [y0 - grid_size, y0, y0 + grid_size]
        points = [(x, y) for y in Ys for x in Xs]

        for i, (x, y) in enumerate(points):
            # center_box confines candidate centers: grid points outside the box are
            # dropped, so the refined grid can never walk out of it. None = no box.
            if not _center_in_box(x, y, center_box):
                if diagnostic:
                    print(f"[out of box] grid pt {i + 1}/9 at (x={x:.2f}, y={y:.2f}) — skipped")
                continue
            if (x, y) in fs_cache:
                result = fs_cache[(x, y)]
                if diagnostic:
                    print(f"[cache hit] grid pt {i + 1}/9 at (x={x:.2f}, y={y:.2f}) → FS={result['FS']:.4f}")
                continue

            depth_step_init = grid_size * 0.75
            d, FS, df_slices, failure_surface, solver_result = optimize_depth(
                x, y, depth_guess, depth_step_init, depth_shrink_factor=0.25, tol_frac=0.01, fs_fail=fs_fail,
                circle_cache=circle_cache, diagnostic=diagnostic
            )

            fs_cache[(x, y)] = {
                "Xo": x,
                "Yo": y,
                "Depth": d,
                "FS": FS,
                "slices": df_slices,
                "failure_surface": failure_surface,
                "solver_result": solver_result
            }

            if diagnostic:
                print(f"[grid pt {i + 1}/9] x={x:.2f}, y={y:.2f} → FS={FS:.4f} at d={d:.2f}")

        sorted_fs = sorted(fs_cache.items(), key=lambda item: item[1]['FS'])
        best_point = sorted_fs[0][1]
        best_index = list(fs_cache.keys()).index((best_point['Xo'], best_point['Yo']))

        if diagnostic:
            print(f"[★ grid best {best_index + 1}/9] FS={best_point['FS']:.4f} at (x={best_point['Xo']:.2f}, y={best_point['Yo']:.2f})")

        return fs_cache, best_point

    # === Step 1: Evaluate starting circles ===
    all_starts = []
    fs_cache = {}  # Shared cache for all starting circles
    for i, start_circle in enumerate(circles):
        _check_cancel(cancel_check)
        x0 = start_circle['Xo']
        y0 = start_circle['Yo']
        depth_guess = start_circle['Depth']
        # With a search window, a seed circle from the input may sit outside it. The
        # SEED is only a launch point, not a reported candidate, so clamp its center
        # into center_box (and its tangent into the tangent_depth band) so the first
        # grid begins inside the window. Candidates are still rejected, never clamped.
        if center_box is not None:
            x1, y1, x2, y2 = center_box
            x0 = min(max(x0, min(x1, x2)), max(x1, x2))
            y0 = min(max(y0, min(y1, y2)), max(y1, y2))
        if tangent_depth is not None:
            tlo, thi = min(tangent_depth), max(tangent_depth)
            depth_guess = min(max(depth_guess, tlo), thi)
        r0 = y0 - depth_guess
        if diagnostic:
            print(f"\n[⏱ starting circle {i+1}] x={x0:.2f}, y={y0:.2f}, r={r0:.2f}")
        grid_size = r0 * 0.15
        fs_cache, best_point = evaluate_grid(x0, y0, grid_size, depth_guess, slope_data, diagnostic=diagnostic, fs_cache=fs_cache, circle_cache=circle_cache)
        all_starts.append((start_circle, best_point))

    all_starts.sort(key=lambda t: t[1]['FS'])

    def refine(start_circle, best_start, launch_label=""):
        """Adaptive 9-point refinement from one start. Returns (fs, converged, path)."""
        nonlocal fs_cache
        x0 = best_start['Xo']
        y0 = best_start['Yo']
        depth_guess = best_start['Depth']
        grid_size = (y0 - depth_guess) * 0.15
        best_fs = best_start['FS']

        # Include initial jump from the seed circle to the best point on its grid
        path = [
            {"x": start_circle['Xo'], "y": start_circle['Yo'], "FS": None},
            {"x": x0, "y": y0, "FS": best_fs}
        ]
        converged = False

        if diagnostic:
            print(f"\n[✅ launch grid{launch_label}] Starting refinement from FS={best_fs:.4f} at ({x0:.2f}, {y0:.2f})")

        fs_level_start = best_fs   # best FS at the start of the current grid resolution
        small_levels = 0           # consecutive refinements that barely improved FS

        for iteration in range(max_iter):
            _check_cancel(cancel_check)
            print(f"[🔁 iteration {iteration+1}{launch_label}] center=({x0:.2f}, {y0:.2f}), FS={best_fs:.4f}, grid={grid_size:.4f}")
            fs_cache, best_point = evaluate_grid(x0, y0, grid_size, depth_guess, slope_data, diagnostic=diagnostic, fs_cache=fs_cache, circle_cache=circle_cache)

            if best_point['FS'] < best_fs:
                # Found a better center at this resolution — move there and keep exploring.
                best_fs = best_point['FS']
                x0 = best_point['Xo']
                y0 = best_point['Yo']
                depth_guess = best_point['Depth']
                path.append({"x": x0, "y": y0, "FS": best_fs})
                continue

            # No improvement at this resolution. Refine the grid, and converge once the
            # FS gain over a refinement level has been negligible twice in a row (so the
            # third decimal of FS is stable). The grid floor `tol` is only a backstop.
            level_gain = fs_level_start - best_fs
            if level_gain < fs_tol and grid_size < min_grid:
                small_levels += 1
                if small_levels >= 2:
                    converged = True
                    elapsed = time.time() - start_time
                    print(f"[✅ converged] Iter={iteration+1}, FS={best_fs:.4f} (ΔFS<{fs_tol}) at (x={x0:.2f}, y={y0:.2f}, depth={depth_guess:.2f}), elapsed time={elapsed:.2f} seconds")
                    break
            else:
                small_levels = 0
            fs_level_start = best_fs
            grid_size *= shrink_factor

            if grid_size < tol:
                converged = True
                elapsed = time.time() - start_time
                print(f"[✅ converged: grid floor] Iter={iteration+1}, FS={best_fs:.4f} at (x={x0:.2f}, y={y0:.2f}, depth={depth_guess:.2f}), elapsed time={elapsed:.2f} seconds")
                break

        if not converged and diagnostic:
            print(f"\n[❌ max iterations reached] FS={best_fs:.4f} at (x={x0:.2f}, y={y0:.2f})")
        return best_fs, converged, path

    # Grid seeding refines EVERY surviving family, not just the best start —
    # competing families can differ by only a few percent at the coarse stage yet
    # refine to very different minima, and launching only the best start is
    # exactly the local-minimum trap the grid exists to remove. Seeded mode
    # ('circles') keeps the original single-launch behavior unchanged.
    if seed == 'grid':
        fs_cut = all_starts[0][1]['FS'] * 1.15
        launches = [t for t in all_starts if t[1]['FS'] <= fs_cut][:4]
    else:
        launches = all_starts[:1]

    search_path = []
    best_fs, converged = None, False
    for i, (sc, bs) in enumerate(launches):
        label = f" (family {i+1}/{len(launches)})" if len(launches) > 1 else ""
        fs_i, conv_i, path_i = refine(sc, bs, launch_label=label)
        search_path.extend(path_i)
        if best_fs is None or fs_i < best_fs:
            best_fs, converged = fs_i, conv_i

    sorted_fs_cache = sorted(fs_cache.values(), key=lambda d: d['FS'])
    _report_unsolved(sorted_fs_cache[0].get('slices') if sorted_fs_cache else None)
    return sorted_fs_cache, converged, search_path, circle_cache

def noncircular_search(slope_data, method_name, rapid=False, diagnostic=True, movement_distance=4.0, shrink_factor=0.8, fs_tol=0.001, max_iter=100, move_tol=0.1, num_slices=30, max_base_angle=65.0, cancel_check=None, min_slip_depth=None):
    """
    Non-circular search using the specified solver.

    A geometric admissibility guard (`max_base_angle`) rejects trial surfaces with
    an over-steep base segment. Without it the coordinate-descent search drives the
    points into a near-vertical base running up to the toe, which the rigorous
    methods (Spencer especially) score as a spurious low minimum. Capping the base
    inclination keeps the search on physically realistic surfaces. (The convergence
    criterion itself is kept absolute, but the guard removes the degeneracy that
    blocked refining it.)

    Parameters:
    -----------
    data : dict
        Input data dictionary containing all necessary parameters
    method_name : str
        The method name to use (e.g., 'lowe', 'spencer')
    diagnostic : bool
        If True, print diagnostic information during search
    movement_distance : float
        Initial distance to move points in each iteration
    shrink_factor : float
        Factor to reduce movement_distance by when no improvement is found
    fs_tol : float
        Factor of safety convergence tolerance
    max_iter : int
        Maximum number of iterations
    move_tol : float
        Minimum movement distance for convergence (AND logic with fs_tol)
    max_base_angle : float
        Maximum allowed base inclination (degrees) for any slice. Trial surfaces
        with a steeper base are rejected as inadmissible. Default 65 deg, the
        active-wedge angle (45 + phi/2) for phi ~ 40 deg, near the steep end of
        realistic soils; steeper bases are geometrically unrealistic slip surfaces.

    Returns:
    --------
    tuple : (fs_cache, converged, search_path)
        fs_cache : dict of all evaluated surfaces and their FS values
        converged : bool indicating if search converged
        search_path : list of surfaces evaluated during search
    """
    if not slope_data.get('non_circ'):
        raise ValueError(
            "Non-circular search requires a non-circular surface defined in the "
            "input. Add surface point data to the 'non-circ' sheet in the input "
            "template.")

    # Gate once at the entry, not per trial -- see circular_search. This is also
    # what turns "the starting surface is not viable" into the real reason when the
    # method cannot use a non-circular surface at all.
    preflight(slope_data, 'rapid' if rapid else 'lem',
              {'surface': 'noncircular', 'search': True,
               'method': method_name}).raise_for_errors()

    if rapid:
        validate_rapid_drawdown(slope_data)

    slope_data = _with_water_loads_once(slope_data)

    # Get the solver function from solve module
    solver = getattr(solve, method_name)
    def move_point(points, i, dx, dy, movement_type, ground_surface, depth_floor):
        """Move a point while respecting constraints"""
        # Get current point
        point = points[i]
        
        # Calculate new position
        new_x = point[0] + dx
        new_y = point[1] + dy
        
        # For endpoints, ensure they stay on ground surface
        if i == 0 or i == len(points)-1:
            # Create vertical line at new_x
            vertical_line = LineString([(new_x, 0), (new_x, 1000)])  # Arbitrary high y value
            intersection = ground_surface.intersection(vertical_line)
            y = get_y_from_intersection(intersection)
            if y is None:
                return False
            new_y = y
        else:
            # For middle points, ensure they stay below ground surface but above the
            # domain floor (surfaces leaving an irregular bottom are rejected later
            # by the containment check in generate_slices).
            if new_y > ground_surface.interpolate(ground_surface.project(Point(new_x, new_y))).y:
                return False
            if new_y < depth_floor:
                return False
        
        # Check x-ordering constraints
        if i > 0 and new_x <= points[i-1][0]:  # Don't move past left neighbor
            return False
        if i < len(points)-1 and new_x >= points[i+1][0]:  # Don't move past right neighbor
            return False
        
        # Update point
        points[i] = [new_x, new_y]
        return True

    def evaluate_surface(points, distance, fs_cache=None):
        """Evaluate factor of safety for current surface configuration"""
        if fs_cache is None:
            fs_cache = {}
            
        # Create non_circ format from points
        non_circ = [{'X': x, 'Y': y, 'Movement': movements[i]} for i, (x, y) in enumerate(points)]

        # Geometric admissibility on the surface itself: slicing breaks at
        # every vertex, so a steep segment with real width gets judged by the
        # per-slice alpha guard below - but generate_slices merges slices
        # thinner than its min_width (1e-2), so a near-vertical segment
        # narrower than that is absorbed into its neighbour and its angle
        # diluted out of the alpha table. Reject exactly those hidden risers.
        seg = np.diff(np.array([[p_['X'], p_['Y']] for p_ in non_circ]), axis=0)
        seg_ang = np.degrees(np.arctan2(np.abs(seg[:, 1]),
                                        np.maximum(np.abs(seg[:, 0]), 1e-300)))
        hidden_riser = (np.abs(seg[:, 0]) < 1e-2) & (seg_ang > max_base_angle)
        if hidden_riser.any():
            return float('inf'), None, None, None, fs_cache

        # Generate slices and compute FS
        success, result = generate_slices(slope_data, non_circ=non_circ,
                                          num_slices=num_slices, check_inputs=False)
        if not success:
            return float('inf'), None, None, None, fs_cache
            
        df_slices, failure_surface = result

        # Admissibility guard: reject surfaces with an over-steep base segment. A
        # near-vertical base (typically running up to the toe) is geometrically
        # unrealistic and is the spurious low-FS minimum the coordinate-descent
        # search would otherwise slide into (worst for Spencer/Lowe).
        if np.abs(df_slices['alpha'].values).max() > max_base_angle:
            return float('inf'), None, None, None, fs_cache

        # Reject surfaces with negligible net driving force (flat / non-failure).
        if _net_driving_too_small(df_slices):
            return float('inf'), None, None, None, fs_cache

        # Optional minimum-slip-depth filter: drop surficial skins on demand.
        if _too_shallow(df_slices, min_slip_depth):
            return float('inf'), None, None, None, fs_cache

        if rapid:
            solver_success, solver_result = rapid_drawdown(df_slices, method_name, debug_level=0)
        else:
            solver_success, solver_result = solver(df_slices)
        FS = solver_result['FS'] if solver_success else float('inf')
        if solver_success and _base_tension_too_extensive(df_slices):
            FS = float('inf')  # degenerate surface (base mostly in tension)
        if not (FS > 0):  # reject non-positive / NaN factor of safety
            FS = float('inf')
        
        # Cache result
        key = tuple(map(tuple, points))
        fs_cache[key] = {
            'points': points.copy(),
            'FS': FS,
            'slices': df_slices,
            'failure_surface': failure_surface,
            'solver_result': solver_result
        }
        
        return FS, df_slices, failure_surface, solver_result, fs_cache

    # Get initial surface from non_circ data
    non_circ = slope_data['non_circ']
    points = np.array([[p['X'], p['Y']] for p in non_circ])
    movements = [p['Movement'] for p in non_circ]
    ground_surface = slope_data['ground_surface']
    
    # Initialize cache and search path
    fs_cache = {}
    search_path = []
    
    # Evaluate initial surface
    FS, df_slices, failure_surface, solver_result, fs_cache = evaluate_surface(
        points, movement_distance, fs_cache)
    
    # A dead starting surface (slice generation or the solver failed on it)
    # leaves failure_surface=None; the search loop can't recover from that, so
    # report "no valid surface" instead of crashing on best_surface.coords.
    if failure_surface is None or not np.isfinite(FS):
        print("[❌ noncircular_search] the starting surface is not viable "
              "(slice generation or the solver failed on it) — adjust the "
              "non_circ starting points.")
        return [], False, []

    # Initialize best surface with initial evaluation
    best_points = points.copy()
    best_fs = FS
    best_df = df_slices
    best_surface = failure_surface
    best_solver_result = solver_result
    
    # Track convergence
    converged = False
    start_time = time.time()
    prev_fs = best_fs
    
    if diagnostic:
        print(f"\n[✅ starting search] Initial FS={best_fs:.4f}\n")
        print("Initial failure surface:")
        for i, point in enumerate(points):
            print(f"Point {i}: ({point[0]:.2f}, {point[1]:.2f})")
        print("\nGround surface:")
        for i, point in enumerate(ground_surface.coords):
            print(f"Point {i}: ({point[0]:.2f}, {point[1]:.2f})")
    
    # Main search loop
    for iteration in range(max_iter):
        _check_cancel(cancel_check)
        improved = False
        
        if diagnostic:
            print(f"\nIteration {iteration + 1}")
            print("Current surface points:")
            for i, point in enumerate(best_surface.coords):
                print(f"Point {i}: ({point[0]:.2f}, {point[1]:.2f})")
        
        # Try moving each point
        for i in range(len(points)):
            # Try both positive and negative directions
            for direction in [-1, 1]:
                test_points = points.copy()
                
                # Get movement direction based on point type
                if i == 0 or i == len(points)-1:  # End points
                    dx = direction * movement_distance
                    dy = 0  # y will be determined by ground surface
                elif movements[i] == 'Horiz':
                    dx = direction * movement_distance
                    dy = 0
                elif movements[i] == 'Free':
                    # For free points, move perpendicular to tangent
                    dx_tangent = points[i+1][0] - points[i-1][0] if i > 0 and i < len(points)-1 else 1
                    dy_tangent = points[i+1][1] - points[i-1][1] if i > 0 and i < len(points)-1 else 0
                    length = np.sqrt(dx_tangent**2 + dy_tangent**2)
                    if length > 0:
                        dx = -dy_tangent/length * direction * movement_distance
                        dy = dx_tangent/length * direction * movement_distance
                    else:
                        dx = direction * movement_distance
                        dy = 0
                else:  # Fixed
                    continue
                
                # Try to move the point
                if move_point(test_points, i, dx, dy, movements[i], ground_surface, slope_data['domain_polygon'].bounds[1]):
                    # Evaluate new surface
                    FS, df_slices, failure_surface, solver_result, fs_cache = evaluate_surface(
                        test_points, movement_distance, fs_cache)
                    
                    if FS < best_fs:
                        best_fs = FS
                        best_points = test_points.copy()
                        best_df = df_slices
                        best_surface = failure_surface
                        best_solver_result = solver_result
                        improved = True
                        if diagnostic:
                            print(f"[✓ improved] iter={iteration}, point={i}, FS={FS:.4f}")
        
        # print iteration results
        print(f"iteration {iteration+1} FS={best_fs:.4f}")
        
        # Check convergence based on FS change and movement distance (AND logic)
        fs_change = abs(best_fs - prev_fs)
        if fs_change < fs_tol and movement_distance < move_tol:
            converged = True
            if diagnostic:
                print(f"[✓ converged] FS change {fs_change:.6f} < tolerance {fs_tol} and movement_distance {movement_distance:.4f} < move_tol {move_tol}")
            break
        prev_fs = best_fs
        
        if not improved or fs_change < fs_tol:
            movement_distance *= shrink_factor
            if True:
                print(f"[↘️ shrinking] movement_distance={movement_distance:.4f}")
        
        points = best_points.copy()
        
    end_time = time.time()
    elapsed = end_time - start_time
    
    if converged:
        print(f"\n[✅ converged] Iter={iteration+1}, FS={best_fs:.4f}, elapsed time={elapsed:.2f} seconds")
    else:
        print(f"\n[❌ max iterations reached] FS={best_fs:.4f}, elapsed time={elapsed:.2f} seconds")

    sorted_fs_cache = sorted(fs_cache.values(), key=lambda d: d['FS'])
    return sorted_fs_cache, converged, search_path


# ---------------------------------------------------------------------------
# One analysis, run once
#
# A limit equilibrium analysis is a method, a surface family, and the choice
# between searching for the critical surface and solving the one the input
# specifies. Studio's Run LEM button and the report both need that done, and the
# report needs it for a method the run never solved: the owner's rule is that
# every method the report is asked for is RUN, each searching for its own
# critical surface, rather than borrowed from another method's answer. Both come
# through here, so a method run from the report is run exactly as the same method
# ticked in the Run dialog would have been — the file's own search window
# included.
# ---------------------------------------------------------------------------

def analysis_search_kwargs(slope_data, circular=True, fs_tol=None, tol=None,
                           max_iter=None, min_slip_depth=None, announce=True):
    """The search keyword arguments one analysis is run under.

    The caller's tolerances where it states them, and the model's own window
    (:func:`file_search_window`) for everything it does not: a windowed model is
    searched inside its window whoever started the search. ``tol`` is circular
    only — a non-circular search has no geometric tolerance — and the window is
    reduced to the limits that branch understands
    (:func:`noncircular_search_opts`).
    """
    kw = {}
    if fs_tol is not None:
        kw["fs_tol"] = fs_tol
    if max_iter is not None:
        kw["max_iter"] = max_iter
    if min_slip_depth is not None:
        kw["min_slip_depth"] = min_slip_depth
    if circular and tol is not None:
        kw["tol"] = tol
    win = file_search_window(slope_data, already=kw)
    if not circular:
        win = noncircular_search_opts(win)
    if win and announce:
        print(f"Applying the file's search window: {', '.join(sorted(win))}.")
    kw.update(win)
    return kw


def run_lem_analysis(slope_data, method, analysis="auto_search", surface="circular",
                     num_slices=40, rapid=False, composite=False, grid_seed=False,
                     diagnostic=False, cancel_check=None, fs_tol=None, tol=None,
                     max_iter=None, min_slip_depth=None, announce=True):
    """Run ONE method on this model and return the bundle every consumer reads.

    Parameters
    ----------
    slope_data : dict
        The model, as :func:`xslope.fileio.load_slope_data` returns it.
    method : str
        The solver's own name for the method (``"spencer"``, ``"bishop"``, ...).
    analysis : {"auto_search", "single_surface"}
        Whether to search for the critical surface or solve the surface the input
        specifies. A search finds the critical surface FOR THIS METHOD: run the
        same model under two methods and the two searches settle on two surfaces.
    surface : {"circular", "noncircular"}
        Which family to run over, for a model that defines both.
    num_slices, rapid, composite, grid_seed, diagnostic
        The run options, under the names Studio's Run LEM dialog gives them.
    cancel_check : callable, optional
        Polled at iteration boundaries; raises :class:`AnalysisCancelled`.
    fs_tol, tol, max_iter, min_slip_depth
        Search tolerances and limits; anything left None is taken from the model's
        own search window where it declares one (:func:`analysis_search_kwargs`).
    announce : bool
        Whether to print what is being run, which is how Studio's Log pane and a
        script's console follow a solve that takes minutes.

    Returns
    -------
    dict
        ``{"slice_df", "failure_surface", "results", "search", "method",
        "options"}`` — the bundle the report, the plots and Studio's result views
        all read. ``search`` is None for a single-surface run and a dict
        describing the search otherwise — including ``unsolved``, the count of
        trial surfaces the method could not answer on and how they rank against
        the reported minimum (:class:`UnsolvedTrials`; None on the non-circular
        branch, which does not yet distinguish them). ``results`` is None where the solver
        returned no solution on a surface that was otherwise built, with the
        solver's reason on ``failure``: the method ran and did not converge,
        which is an answer about the model rather than a run that never happened.

        ``options`` is what this run was made under, recorded here because here
        is the only place that knows it. A second method run beside this one is
        run from that record rather than from anything read back out of the
        answer — a slice count in particular, which cannot be recovered from the
        slice table: the slicer splits at material and water boundaries, so a run
        of 20 comes back as 21 rows, and a neighbour re-run at 21 is a different
        discretization on the same surface.

    Raises
    ------
    AnalysisError
        The run could not be made at all — no surface of the requested family, a
        surface that cannot be sliced, a search that found nothing valid.
    AnalysisCancelled
        ``cancel_check`` asked it to stop.
    """
    method = str(method or "").lower()
    circular = surface != "noncircular"
    tag = " (rapid drawdown)" if rapid else ""
    # Recorded before anything is solved, so what comes back says how it was made
    # whichever branch made it.
    made_under = {"analysis": analysis,
                  "surface": "circular" if circular else "noncircular",
                  "num_slices": num_slices, "rapid": rapid,
                  "composite": composite, "grid_seed": grid_seed,
                  "diagnostic": diagnostic, "fs_tol": fs_tol, "tol": tol,
                  "max_iter": max_iter, "min_slip_depth": min_slip_depth}
    if analysis == "auto_search":
        kw = analysis_search_kwargs(slope_data, circular=circular, fs_tol=fs_tol,
                                    tol=tol, max_iter=max_iter,
                                    min_slip_depth=min_slip_depth,
                                    announce=announce)
        seed = "grid" if grid_seed else "circles"
        if circular:
            if seed != "grid" and not (slope_data.get("circular")
                                       and slope_data.get("circles")):
                raise AnalysisError("Circular search needs at least one starting "
                                    "circle (or turn on grid seeding).")
            if announce:
                print(f"Searching for the critical circular surface with "
                      f"{method.upper()}{tag}…")
            unsolved = {}
            fs_cache, converged, search_path, circle_cache = circular_search(
                slope_data, method, rapid=rapid, num_slices=num_slices,
                diagnostic=diagnostic, cancel_check=cancel_check,
                composite=composite, seed=seed, unsolved_out=unsolved, **kw)
            search = {"kind": "circular", "fs_cache": fs_cache,
                      "search_path": search_path, "circle_cache": circle_cache,
                      "unsolved": unsolved}
        else:
            if not slope_data.get("non_circ"):
                raise AnalysisError("Non-circular search needs a starting "
                                    "non-circular surface.")
            if announce:
                print(f"Searching for the critical non-circular surface with "
                      f"{method.upper()}{tag}…")
            fs_cache, converged, search_path = noncircular_search(
                slope_data, method, rapid=rapid, num_slices=num_slices,
                diagnostic=diagnostic, cancel_check=cancel_check, **kw)
            # The unsolved-trial disclosure is the circular search's; the
            # non-circular branch scores a solver failure the same way it scores
            # a geometric rejection and cannot yet tell them apart.
            search = {"kind": "noncircular", "fs_cache": fs_cache,
                      "search_path": search_path, "circle_cache": None,
                      "unsolved": None}
        if not fs_cache:
            raise AnalysisError("Search produced no valid surfaces.")
        critical = fs_cache[0]
        results = critical.get("solver_result")
        if not isinstance(results, dict):
            raise AnalysisError("Search found no surface with a valid solution.")
        if announce:
            tail = "" if converged else "  (search did not fully converge)"
            print(f"Critical FS = {results.get('FS'):.3f}{tail}")
        return {"slice_df": critical.get("slices"),
                "failure_surface": critical.get("failure_surface"),
                "results": results, "search": search, "method": method,
                "options": made_under}

    # --- the surface the input specifies ---
    if circular:
        circle = (slope_data["circles"][0]
                  if slope_data.get("circular") and slope_data.get("circles")
                  else None)
        if circle is None:
            raise AnalysisError("A circular surface is required (no circles "
                                "defined).")
        non_circ = None
        if announce:
            print(f"Running {method.upper()} — single circular surface "
                  f"(Xo={circle.get('Xo')}, Yo={circle.get('Yo')}, "
                  f"R={circle.get('R'):.3g}), {num_slices} slices{tag}…")
    else:
        non_circ = slope_data.get("non_circ") or None
        if not non_circ:
            raise AnalysisError("A non-circular surface is required (no "
                                "non-circular points defined).")
        circle = None
        if announce:
            print(f"Running {method.upper()} — single non-circular surface, "
                  f"{num_slices} slices{tag}…")

    ok, out = generate_slices(slope_data, circle=circle, non_circ=non_circ,
                              num_slices=num_slices, composite=composite)
    if not ok:
        raise AnalysisError(str(out))
    slice_df, failure_surface = out
    if announce:
        print(f"Generated {len(slice_df)} slices; solving…")
    results = solve.solve_selected(method, slice_df, rapid=rapid)
    bundle = {"slice_df": slice_df, "failure_surface": failure_surface,
              "results": results if isinstance(results, dict) else None,
              "search": None, "method": method, "options": made_under}
    if bundle["results"] is None:
        # The surface was built and the method was given it: it ran, and it did
        # not converge. The bundle keeps the surface so that answer can be
        # reported as the answer it is, with the solver's own reason beside it.
        bundle["failure"] = str(results)
    return bundle
