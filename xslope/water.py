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

"""Water: where the pool stands, and what it weighs.

Ponded water is stored implicitly by most commercial slope-stability programs —
SLOPE/W and Slide2 both derive the weight of water standing on the slope from the
water surface, never as an object of its own — and from template v22 so does
xslope: the engine measures the reservoir at solve time from the model's own water
definition. (RS2 is the exception in both directions: it stores an explicit ponded
water load object, and its whole-domain piezometric surface would not support the
derivation below, so rs2.py imports its models in manual mode.)

This module holds the geometry that measures it:

  - material_above_ground_dload(ground, upper_line, unit_weight): the weight of
    whatever fills the gap between a line and the ground, as dload blocks. Both
    ponded water (upper_line = the water surface, unit_weight = gamma_w) and a
    GeoStudio-style surcharge (upper_line = the fill's top, unit_weight = the
    fill's unit weight) are this same object.
  - ponded_water_dload(ground, piezo_line, gamma_water): the ponded-water name
    for the above, for readers where the intent is only ever water.
  - _y_on(line, x): elevation of a left-to-right polyline at x.

The dload convention is xslope's own: each block is a list of ``{'X','Y','Normal'}``
points along the loaded surface, ``Normal`` the pressure (perpendicular
intensity) at that point, tapering to zero at a water/ground crossing so the load
ends where the water does. It was validated exact against the dry-buoyant oracle
(a fully submerged still-water slope reads identically to the same slope run dry
with gamma' = gamma_sat - gamma_w).

Where the pool stands
---------------------

The blocks above need a water *surface* to measure down from, and a model can
state where the water stands in two ways: with a piezometric line, or with the
seepage boundary conditions. The second half of this module reads either one:

  - :func:`seep_bc_water_line` turns a reservoir or head boundary block drawn on
    the ground surface into the water line it implies. This is the one piece of
    geometry the package did not already have.
  - :func:`water_line_for_stage` applies the precedence between the two sources —
    seepage boundary conditions where a seepage analysis is defined, otherwise the
    stage's piezometric line — and returns the line with the sheet it came from.
  - :func:`derive_water_loads` measures the load: the water line against the ground
    surface, through :func:`material_above_ground_dload`, with the peak pressure
    and reach reported alongside the blocks so a caller can say what it is about
    to add before it adds it.

All three are pure: they read a model and return a description, and nothing here
writes to a model or to a file. Three callers consume them. The preflight remedy
(:mod:`xslope.remedies`) offers the derived blocks for a manual-mode model that
is missing its reservoir; the automatic water-load mode derives the same blocks at
solve time; and the vendor importers (:mod:`xslope.geostudio`,
:mod:`xslope.slide2`) call :func:`derive_water_loads` to *report* what the engine
will apply to the model they just built, rather than writing a load of their own.
One derivation, so a load offered by a button, a load described by an import and a
load applied by the engine can never be different loads.

Where a boundary condition's level is time-varying, its value comes from the
transient march's own series evaluator, so a derived load and the seepage field
it accompanies can never disagree about where the pool was at that moment.

The engine's side of it
-----------------------

``main!D23`` (**Water loads: auto | manual**) decides who supplies the load.

  - :func:`water_loads_mode` reads it. A model that says nothing — every template
    version through v21 — means ``manual``, because those files carry the water
    load already typed in and deriving a second one would double the reservoir.
  - :func:`with_water_loads` is what an engine calls. In manual mode it hands
    back the caller's own model, unchanged and by identity; in auto mode a copy
    carrying the derived blocks in keys of their own. Both the LEM slice forces
    and the FEM tractions go through it, so the two engines cannot end up
    carrying different water.
  - :func:`match_water_blocks` answers when a transcribed block *is* the derived
    reservoir. The auto-mode double-count warning, the switch-to-auto remedy and
    the corpus migration all ask that question, and all three ask it here.
"""

from shapely.geometry import LineString, Point


def _y_on(line, x):
    """Elevation of a left-to-right polyline at x, or None if x is off its ends."""
    if line is None or line.is_empty:
        return None
    coords = list(line.coords)
    if not coords or x < coords[0][0] or x > coords[-1][0]:
        return None
    for (x1, y1), (x2, y2) in zip(coords, coords[1:]):
        if x1 <= x <= x2:
            if x2 == x1:
                return max(y1, y2)
            return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
    return None


def _ground_walk(ground_surface, upper_line):
    """The ground's own vertices in order, with the upper line's vertices spliced in.

    Sampling the ground by its *x* values -- the obvious thing, and what this
    function replaces -- loses a **vertical face**: two vertices share one x, a set
    of x values keeps one of them, and the load on the face disappears into a
    fictitious diagonal drawn to the next vertex. A stepped wall or a bench is
    exactly that shape, so the walk is by vertex and not by coordinate.

    Splicing the upper line's vertices into the segment that contains them keeps
    the other half of the old behaviour: a bend in the water surface between two
    ground vertices still becomes a sample, so the pressure diagram bends where the
    water does. And because every sample now lies on the ground's own polyline, two
    consecutive samples always share one ground segment -- which is what lets the
    caller interpolate a crossing in x *and* y and land exactly on the ground.
    """
    g = list(ground_surface.coords)
    extra = sorted({float(x) for x, _ in upper_line.coords})
    out = []
    for k, (x1, y1) in enumerate(g):
        out.append((float(x1), float(y1)))
        if k + 1 >= len(g):
            break
        x2, y2 = g[k + 1]
        if x2 <= x1:                                  # vertical or repeated: nothing between
            continue
        for x in extra:
            if x1 < x < x2:
                out.append((x, y1 + (y2 - y1) * (x - x1) / (x2 - x1)))
    return out


def material_above_ground_dload(ground_surface, upper_line, unit_weight):
    """The weight of whatever fills the gap between a line and the ground, as dloads.

    Two GeoStudio things are this same object, and both would otherwise be lost:

    * **Ponded water.** GeoStudio has no ponded-water object -- where the piezometric
      line rises above the ground, SLOPE/W simply carries the water's weight. Pass the
      piezo line and gamma_w.
    * **A surcharge.** GeoStudio's ``<Surcharge>`` is a body of fill between a drawn
      line and the ground, and its ``<Pressure>`` is the fill's UNIT WEIGHT, not a
      pressure -- so the load varies with how deep the fill is. (Verified against
      SLOPE/W's own per-slice surcharge forces: a 20 kN/m3 fill under a line running
      1 m to 2 m above flat ground gives 300 kN/m, not 20 x 10 = 200.)

    Returns a list of dload blocks -- one per stretch where the line is above ground --
    each a list of ``{'X','Y','Normal'}``, with ``Normal`` the pressure
    ``unit_weight * depth``, zero at the crossing so the load tapers out correctly.
    """
    if not upper_line or ground_surface is None or ground_surface.is_empty:
        return []
    piezo = LineString(upper_line)
    gamma_water = unit_weight

    samples = []
    for x, yg in _ground_walk(ground_surface, piezo):
        yp = _y_on(piezo, x)
        if yp is not None:
            samples.append((x, yg, yp - yg))          # depth of water over the ground
    if not samples:
        return []

    blocks, run = [], []
    for i, (x, yg, d) in enumerate(samples):
        if d > 1e-9:
            if not run and i > 0:                     # entering the water: add the waterline
                px, pyg, pd = samples[i - 1]
                t = -pd / (d - pd) if (d - pd) else 0.0
                run.append({"X": px + t * (x - px), "Y": pyg + t * (yg - pyg),
                            "Normal": 0.0})
            run.append({"X": x, "Y": yg, "Normal": gamma_water * d})
        elif run:                                     # leaving the water: close at the waterline
            px, pyg, pd = samples[i - 1]
            t = pd / (pd - d) if (pd - d) else 0.0
            run.append({"X": px + t * (x - px), "Y": pyg + t * (yg - pyg),
                        "Normal": 0.0})
            if len(run) >= 2:
                blocks.append(_drop_redundant(run))
            run = []
    if len(run) >= 2:
        blocks.append(_drop_redundant(run))
    return blocks


def _drop_redundant(block):
    """The block's shape, without the vertices that only record where it was sampled.

    The derivation samples wherever *either* line has a vertex, so a bend in the
    water surface over a straight reach of ground -- or the reverse -- leaves a
    point that carries no information: it sits on the straight line between its
    neighbours in position and in pressure alike. Harmless to the statics, and not
    harmless to everything else. A dload vertex places a slice boundary and seeds a
    mesh node, so a derived block that is sampled more finely than the load it
    describes discretises the model differently from the transcribed block it
    replaces. Dropping them is what lets a conversion be a change of expression
    rather than of model.

    Only exact redundancy goes, judged against the block's own scale so the test
    means the same thing in feet and in metres.
    """
    pts = list(block)
    if len(pts) < 3:
        return pts
    ys = [abs(float(p["Y"])) for p in pts]
    ps = [abs(float(p["Normal"])) for p in pts]
    ytol = 1e-9 * max(max(ys), 1.0)
    ptol = 1e-9 * max(max(ps), 1.0)
    out = [pts[0]]
    for prev, cur, nxt in zip(pts, pts[1:], pts[2:]):
        x1, x3 = float(prev["X"]), float(nxt["X"])
        if x3 == x1:                                  # a vertical run: never redundant
            out.append(cur)
            continue
        t = (float(cur["X"]) - x1) / (x3 - x1)
        y = float(prev["Y"]) + t * (float(nxt["Y"]) - float(prev["Y"]))
        p = float(prev["Normal"]) + t * (float(nxt["Normal"]) - float(prev["Normal"]))
        if abs(y - float(cur["Y"])) > ytol or abs(p - float(cur["Normal"])) > ptol:
            out.append(cur)
    out.append(pts[-1])
    return out


def ponded_water_dload(ground_surface, piezo_line, gamma_water):
    """Water standing above the ground surface, as a distributed load.

    A thin name for :func:`material_above_ground_dload` -- ponded water is just the
    material above the ground being water.
    """
    return material_above_ground_dload(ground_surface, piezo_line, gamma_water)


# ===========================================================================
# Where the pool stands: reading a water surface out of a model
# ===========================================================================

#: A boundary vertex counts as drawn ON the ground surface when it sits this
#: close to it, as a fraction of the model's own width. Boundary polylines are
#: normally traced along the ground's own vertices, so this only has to absorb
#: transcription round-off -- and it must stay tight, because the test is what
#: keeps a deep aquifer head boundary from being read as a reservoir.
_ON_GROUND_FRAC = 1e-6


def _as_coords(obj):
    """``[(x, y), ...]`` from a shapely geometry, a list of pairs, or ``None``."""
    if obj is None:
        return []
    if hasattr(obj, "is_empty") and obj.is_empty:
        return []
    if hasattr(obj, "coords"):
        return [(float(x), float(y)) for x, y in obj.coords]
    out = []
    for p in obj:
        try:
            if isinstance(p, dict):
                out.append((float(p["X"]), float(p["Y"])))
            else:
                out.append((float(p[0]), float(p[1])))
        except (TypeError, ValueError, KeyError, IndexError):
            continue
    return out


def _ascending(pts):
    """True when a polyline's x values never decrease.

    The ordering precondition every derivation here depends on: :func:`_y_on`
    walks a polyline's segments in order and returns ``None`` the moment x runs
    backwards, so a line entered right-to-left would silently produce no water
    rather than the reservoir it describes. Checked and reported, never repaired
    in passing -- reversing a line is its own explicit act.
    """
    return all(b[0] >= a[0] for a, b in zip(pts, pts[1:]))


def series_value(slope_data, name, time=0.0):
    """A ``tseep`` time series evaluated at ``time``, or ``None`` if there is none.

    The evaluation goes through the transient march's own
    :func:`xslope.seep._eval_series`, deliberately: a derived water load and the
    seepage field it accompanies must not be able to disagree about where the pool
    stood at that moment, and the only way to guarantee that is one implementation
    of the interpolation.
    """
    ts = (slope_data or {}).get("tseep") or {}
    vals = (ts.get("series") or {}).get(name)
    times = ts.get("times")
    if not vals or times is None:
        return None
    ta, va = [], []
    for t, v in zip(times, vals):
        if v is None:
            continue
        try:
            ta.append(float(t))
            va.append(float(v))
        except (TypeError, ValueError):
            continue
    if not ta:
        return None
    from .seep import _eval_series             # local: seep.py is a heavy import
    return float(_eval_series((ta, va), 0.0 if time is None else float(time)))


def bc_pool_levels(bc, ground_surface, slope_data=None, time=0.0):
    """The pool levels a boundary-condition set describes, with where each applies.

    Returns ``[(level, [x, ...], label), ...]`` -- one entry per specified-head
    block that describes standing water: its level, the x values of the vertices it
    is anchored at, and a label naming it in the sheet's vocabulary.

    Two filters decide what counts, and the first is the important one:

    * **The block must be drawn on the ground surface.** A head boundary along a
      deep aquifer, a cutoff, or a model side is a groundwater condition, not a
      reservoir; reading its level as a water surface would flood the section with
      a pool nobody drew. A block qualifies only where at least one of its vertices
      lies on the ground surface.
    * **Its level must resolve to a number.** A block whose VALUE cell names a
      ``tseep`` series is evaluated at ``time`` (stage 1 reads t = 0); a series that
      does not exist yields nothing rather than a guess.

    Flux boundaries and exit faces never describe standing water and are not read.
    """
    ground = _as_coords(ground_surface)
    if len(ground) < 2:
        return []
    line = LineString(ground)
    xs = [p[0] for p in ground]
    tol = max(1e-12, _ON_GROUND_FRAC * (max(xs) - min(xs)))
    out = []
    for n, blk in enumerate(((bc or {}).get("specified_heads") or [])):
        raw = blk.get("head")
        if isinstance(raw, str):
            level = series_value(slope_data, raw.strip(), time)
            named = f"the '{raw.strip()}' series at t = {0.0 if time is None else time:g}"
        else:
            try:
                level = float(raw)
            except (TypeError, ValueError):
                level = None
            named = None
        if level is None:
            continue
        anchors = [c[0] for c in (blk.get("coords") or [])
                   if line.distance(Point(float(c[0]), float(c[1]))) <= tol]
        if not anchors:
            continue
        kind = str(blk.get("kind") or "head").strip().lower()
        label = (f"{'reservoir' if kind == 'reservoir' else 'head'} boundary #{n + 1}"
                 f" at elevation {level:g}")
        if named:
            label += f" (from {named})"
        out.append((float(level), sorted(anchors), label))
    return out


def seep_bc_water_line(ground_surface, bc, slope_data=None, time=0.0):
    """The water surface a seepage boundary set implies, as ``[(x, y), ...]``.

    The boundary conditions are the authority on where water stands: a reservoir or
    head block traced along the upstream face states the pool elevation directly,
    and no seepage solution is needed to read it. This turns that statement into a
    line the load derivation can measure against -- flat at the boundary's level
    over the reach of ground the pool actually covers, and lying on the ground
    everywhere else.

    The covered reach is found by walking outward from the boundary's own vertices
    for as long as the ground stays below the level, so the pool ends where the
    ground rises out of it rather than at the boundary's last vertex. Two blocks
    that overlap take the higher level.

    **The shoreline is a sample of its own.** The line is sampled at every ground
    vertex, at the boundary's own vertices, and at every point where the ground
    crosses a pool's level -- the last of those because otherwise the pool's edge
    lands on the next ground vertex instead of the waterline, and the load runs up
    the face above the water. Where a modeller happens to trace the boundary block
    to the shoreline the two agree; where the level is time-varying, or the block is
    drawn up the whole face as RS2 draws it, they do not, and the pool was
    over-measured by the height of the face above the water.

    Returns ``[]`` when no block describes standing water above the ground.
    """
    ground = _as_coords(ground_surface)
    if len(ground) < 2 or not _ascending(ground):
        return []
    pools = bc_pool_levels(bc, ground_surface, slope_data, time)
    if not pools:
        return []

    shore = set()
    for level, _anchors, _label in pools:
        for (x1, y1), (x2, y2) in zip(ground, ground[1:]):
            if (y1 - level) * (y2 - level) < 0 and y2 != y1:
                shore.add(x1 + (level - y1) / (y2 - y1) * (x2 - x1))
    xs = sorted({p[0] for p in ground} | shore
                | {x for _, anchors, _ in pools for x in anchors})
    line = LineString(ground)
    yg = [_y_on(line, x) for x in xs]
    if any(y is None for y in yg):
        xs = [x for x, y in zip(xs, yg) if y is not None]
        yg = [y for y in yg if y is not None]
    if len(xs) < 2:
        return []
    span = max(xs) - min(xs)
    tol = max(1e-12, _ON_GROUND_FRAC * span)

    level_at = [None] * len(xs)
    for level, anchors, _label in pools:
        covered = set()
        for xa in anchors:
            # nearest sample to the anchor, then walk both ways while submerged
            i = min(range(len(xs)), key=lambda k: abs(xs[k] - xa))
            if yg[i] > level + tol:
                continue
            covered.add(i)
            for step in (-1, 1):
                j = i + step
                while 0 <= j < len(xs) and yg[j] <= level + tol:
                    covered.add(j)
                    j += step
        for i in covered:
            if level_at[i] is None or level > level_at[i]:
                level_at[i] = level

    if not any(l is not None and l - g > tol for l, g in zip(level_at, yg)):
        return []
    return [(x, g if l is None else max(l, g)) for x, g, l in zip(xs, yg, level_at)]


def water_line_for_stage(slope_data, stage=1, time=None):
    """Where the water stands for one stage, and which sheet says so.

    Applies the load-source precedence: **seepage boundary conditions wherever a
    seepage analysis is defined** (``seep bc`` for stage 1, ``seep bc (2)`` for
    stage 2), **otherwise the stage's piezometric line** (Line 1, or Line 2 for
    stage 2). The boundary conditions win because they are where the modeller
    stated the pool for the analysis that will actually be run; the piezometric
    line is the answer for a model that has no seepage analysis at all.

    ``time`` selects the moment for a time-varying boundary level. It defaults to
    t = 0 for stage 1 and to the ``tseep`` sheet's ``stage_2`` time for stage 2,
    which are the two moments the staged analysis solves at.

    Returns a dict with ``points`` (the water line, possibly empty), ``source``
    (where it came from, named the way the interface labels it) and ``reason`` (why
    there is no line, when there is none).
    """
    sd = slope_data or {}
    ground = _as_coords(sd.get("ground_surface"))
    if len(ground) < 2:
        return {"points": [], "source": "", "reason": "the model defines no ground surface"}
    if not _ascending(ground):
        return {"points": [], "source": "",
                "reason": "the ground surface's X values are not in increasing order"}

    bc = sd.get("seepage_bc" if stage == 1 else "seepage_bc2") or {}
    if time is None:
        ts = sd.get("tseep") or {}
        time = 0.0 if stage == 1 else ts.get("stage_2", 0.0)
    which = "the seepage head boundaries" if stage == 1 else \
        "the stage-2 seepage head boundaries"
    if bc.get("specified_heads"):
        pts = seep_bc_water_line(sd.get("ground_surface"), bc, sd, time)
        levels = bc_pool_levels(bc, sd.get("ground_surface"), sd, time)
        if pts:
            names = ", ".join(label for _l, _a, label in levels)
            return {"points": pts, "source": f"{which} ({names})", "reason": ""}
        return {"points": [], "source": which,
                "reason": f"{which} do not put water above the ground surface"}

    key = "piezo_line" if stage == 1 else "piezo_line2"
    name = "Piezometric Line 1" if stage == 1 else "Piezometric Line 2"
    pz = _as_coords(sd.get(key))
    if len(pz) < 2:
        return {"points": [], "source": "",
                "reason": (f"the model defines neither seepage head boundaries nor "
                           f"{name}, so nothing says where the water stands")}
    if not _ascending(pz):
        return {"points": [], "source": name,
                "reason": (f"{name} does not run left to right, and a water surface "
                           f"can only be measured along a line whose X values "
                           f"increase (Piezometric lines; piezo sheet)")}
    return {"points": pz, "source": name, "reason": ""}


#: How deep a pool must be before it is a pool, as a fraction of the section's own
#: vertical extent. Coordinates are transcribed, digitised off figures, or produced
#: by a parametric frame, and a water line meant to *meet* the ground -- a phreatic
#: surface exiting at a toe, a piezometric line drawn along a flat foreshore, a
#: piezo tail carried past the model's edge -- lands a few millimetres to a few
#: centimetres off it. Those residuals are the coordinates' own precision, not
#: standing water. The Wave 4b corpus audit measured both populations: the deepest
#: residual is 8.0e-4 of its section's height (vp072b, 72 mm at a downstream toe on
#: a 90 ft dam) and the shallowest real pool is 3.96e-2 (vp034, 4.2 ft of tailwater
#: on a 106 ft section) -- a factor of 49 between them, with nothing in between.
#: The fence sits at 1e-3: 1.25x above the deepest residual, 40x below the
#: shallowest pool.
POOL_MIN_DEPTH_FRAC = 1e-3


def _vertical_face_ambiguity(ground, blocks, levels):
    """Vertical ground faces the derivation loaded while two pools were in play.

    A water surface is a function of x, and a **vertical face is the one place two
    pools can meet at the same x**: a dam's stepped apron with the reservoir behind
    and the tailrace in front, a dewatered trench cut into a seabed. On either side
    of such a face the water stands at a different elevation, the single-valued
    water line can only carry one of them, and the derivation loads the face to the
    higher -- flooding the dry side.

    Reported rather than guessed at. Where a section has one pool the question does
    not arise and this is empty; where it has two and the ground has a vertical face
    between them, the derived load on that face is not trustworthy and the model
    should keep its water loads typed in.
    """
    if len({round(float(v), 9) for v in levels}) < 2:
        return []
    loaded = {(round(float(p["X"]), 9), round(float(p["Y"]), 9))
              for blk in blocks for p in blk}
    out = []
    for (x1, y1), (x2, y2) in zip(ground, ground[1:]):
        if x2 != x1 or y2 == y1:
            continue
        if (round(x1, 9), round(y1, 9)) in loaded or (round(x2, 9), round(y2, 9)) in loaded:
            out.append(
                f"the vertical face at x = {x1:.6g}, from elevation {min(y1, y2):.6g} "
                f"to {max(y1, y2):.6g}, was loaded from a water surface that stands "
                f"at more than one elevation on this section -- the water line cannot "
                f"say which pool wets the face, so it used the one above it")
    return out


def derive_water_loads(slope_data, stage=1, time=None):
    """The distributed load the standing water applies, derived from the model.

    The whole derivation in one call: find where the water stands
    (:func:`water_line_for_stage`), measure it against the ground surface
    (:func:`material_above_ground_dload`), and report the result with enough
    detail to describe it before anything is changed -- which is what lets a
    remedy state *"1 block, 262 peak, over x = 346.7 to 594"* while the user is
    still deciding.

    Blocks shallower than :data:`POOL_MIN_DEPTH_FRAC` of the section's vertical
    extent are dropped as coordinate round-off rather than water, and what was
    dropped is reported in ``negligible`` so the decision is never silent.

    Returns a dict:

    ``blocks``      dload blocks in xslope's own form (``X``/``Y``/``Normal``)
    ``dirs``        the matching direction words -- water always acts ``normal``
    ``source``      which sheet the water line came from
    ``peak``        the largest pressure in the derived blocks
    ``reach``       ``(x_first, x_last)`` of the loaded stretch
    ``resultant``   the total load per unit width, for comparing against a
                    transcribed block
    ``negligible``  one description per block dropped as too shallow to be a pool
    ``ambiguous``   one description per vertical ground face loaded while the
                    section carried water at more than one elevation
    ``reason``      why nothing was derived, when nothing was
    """
    sd = slope_data or {}
    empty = {"blocks": [], "dirs": [], "source": "", "peak": 0.0, "reach": None,
             "resultant": 0.0, "negligible": [], "ambiguous": [], "reason": ""}
    try:
        gamma_w = float(sd.get("gamma_water"))
    except (TypeError, ValueError):
        gamma_w = None
    if not gamma_w or gamma_w <= 0:
        return dict(empty, reason="Unit weight of water is blank or not positive, "
                                  "so the water has no weight to apply (Global "
                                  "parameters; main D10)")

    found = water_line_for_stage(sd, stage=stage, time=time)
    if not found["points"]:
        return dict(empty, source=found["source"], reason=found["reason"])

    ground = _as_coords(sd.get("ground_surface"))
    blocks = material_above_ground_dload(LineString(ground), found["points"], gamma_w)
    if not blocks:
        return dict(empty, source=found["source"],
                    reason=(f"{found['source']} does not stand above the ground "
                            f"surface anywhere, so there is no water load to add"))

    bc = sd.get("seepage_bc" if stage == 1 else "seepage_bc2") or {}
    ambiguous = _vertical_face_ambiguity(
        ground, blocks,
        [lvl for lvl, _a, _lab in bc_pool_levels(
            bc, sd.get("ground_surface"), sd,
            time if time is not None else (0.0 if stage == 1
                                           else (sd.get("tseep") or {}).get("stage_2", 0.0)))])
    ys = [y for _x, y in ground] + [y for _x, y in found["points"]]
    extent = max(ys) - min(ys)
    fence = gamma_w * POOL_MIN_DEPTH_FRAC * extent
    kept, negligible = [], []
    for blk in blocks:
        blk_peak = max(float(p["Normal"]) for p in blk)
        if blk_peak > fence:
            kept.append(blk)
        else:
            bx = [float(p["X"]) for p in blk]
            negligible.append(
                f"{blk_peak / gamma_w:.4g} deep over x = {min(bx):.6g} to {max(bx):.6g} "
                f"-- below {POOL_MIN_DEPTH_FRAC:g} of the section's {extent:.6g} "
                f"height, so it is where the water line meets the ground rather than "
                f"a pool standing on it")
    if not kept:
        return dict(empty, source=found["source"], negligible=negligible,
                    ambiguous=ambiguous,
                    reason=(f"{found['source']} stands above the ground surface only "
                            f"by less than {POOL_MIN_DEPTH_FRAC:g} of the section's "
                            f"height, which is coordinate round-off rather than a pool"))

    peak = max(float(p["Normal"]) for blk in kept for p in blk)
    xs = [float(p["X"]) for blk in kept for p in blk]
    resultant = sum(block_resultant(blk) for blk in kept)
    return {"blocks": kept, "dirs": ["normal"] * len(kept),
            "source": found["source"], "peak": peak, "reach": (min(xs), max(xs)),
            "resultant": resultant, "negligible": negligible,
            "ambiguous": ambiguous, "reason": ""}


# ===========================================================================
# The mode, and the loads the engine derives under it
# ===========================================================================

#: The two values of ``main!D23`` (Water loads), and of ``slope_data['water_loads']``.
WATER_LOAD_MODES = ("auto", "manual")

#: Where the derivation puts its blocks, per stage. They are kept in their own
#: keys rather than appended to ``dloads`` for one reason: a derived load and a
#: load the user typed must stay tellable apart -- by the plots, which draw them
#: differently; by the double-count rule, which compares one against the other;
#: and by the writer, which saves the user's sheet and never the derivation.
DERIVED_KEYS = {1: "dloads_derived", 2: "dloads2_derived"}

#: The marker a derived model carries, so deriving twice is a no-op rather than a
#: second reservoir. It also carries the per-stage description of what was
#: derived, which is what lets a report say where the load came from.
DERIVED_META_KEY = "water_derived"


def water_loads_mode(slope_data):
    """``"auto"`` or ``"manual"`` for a model -- who supplies the water load.

    ``auto`` means the engine derives the ponded-water surface load from the
    model's own water definition at solve time, and the dloads sheets carry
    non-water loads only. ``manual`` means the load is whatever the user typed,
    which is what every template version through v21 meant and therefore what a
    model that says nothing means.
    """
    mode = (slope_data or {}).get("water_loads")
    if mode is None:
        return "manual"
    mode = str(mode).strip().lower()
    if mode not in WATER_LOAD_MODES:
        raise ValueError(
            f"Water loads is {mode!r}, which is not a mode. Expected one of: "
            f"{', '.join(WATER_LOAD_MODES)} (main D23).")
    return mode


def with_water_loads(slope_data, time=None):
    """The model an engine should actually solve, water loads included.

    In **manual** mode this is the caller's own model, returned unchanged and by
    identity: nothing is derived, nothing is copied, and the analysis is
    bit-identical to the one that ran before automatic water loads existed.

    In **auto** mode it is a shallow copy carrying two extra keys --
    ``dloads_derived`` and ``dloads2_derived`` (:data:`DERIVED_KEYS`) -- holding
    the ponded-water blocks derived from the water definition for stage 1 and
    stage 2, plus :data:`DERIVED_META_KEY` describing where each came from. The
    caller's model is never mutated, and the user's own ``dloads`` are passed
    through untouched: under the automatic semantics those are non-water loads,
    and a load the user typed is never removed by an engine.

    ``time`` fixes the moment a time-varying boundary level is read at. Left
    None, stage 1 is read at the moment the run is solving for -- the transient
    frame's own time when one has been selected (``slope_data['seep_u_time']``,
    set by :func:`xslope.seep.select_transient_frame_u`), otherwise t = 0 -- and
    stage 2 at the ``tseep`` sheet's ``stage_2`` time. One derivation per moment,
    so the water pressing on the slope and the pore-pressure field inside it can
    never disagree about where the pool stood.

    Calling this on a model that already carries the derivation returns it
    unchanged, so an engine may call it defensively without deriving twice.
    """
    sd = slope_data or {}
    if DERIVED_META_KEY in sd:
        return slope_data
    if water_loads_mode(sd) != "auto":
        return slope_data
    if time is None:
        time = sd.get("seep_u_time")
    out = dict(sd)
    meta = {"mode": "auto"}
    for stage in (1, 2):
        found = derive_water_loads(sd, stage=stage,
                                   time=time if stage == 1 else None)
        out[DERIVED_KEYS[stage]] = found["blocks"]
        meta[stage] = {k: found[k] for k in
                       ("source", "peak", "reach", "resultant", "negligible",
                        "ambiguous", "reason")}
    out[DERIVED_META_KEY] = meta
    return out


def derived_blocks(slope_data, stage=1):
    """The derived water blocks a model carries for one stage (``[]`` in manual mode).

    Reads what :func:`with_water_loads` put there rather than deriving again, so a
    plot and the solve it illustrates show the same load.
    """
    return list((slope_data or {}).get(DERIVED_KEYS.get(stage, "")) or [])


def has_derived_loads(slope_data):
    """True when a model carries any derived water load, either stage."""
    return bool(derived_blocks(slope_data, 1) or derived_blocks(slope_data, 2))


# ===========================================================================
# Matching a transcribed block against the derivation
# ===========================================================================

#: How close a transcribed block must be to the derivation to count as the same
#: load. Every measure below is relative -- to the derived peak, the derived
#: resultant, and the derived reach -- so the number travels between unit systems
#: and model sizes.
#:
#: Read off the Wave 4b corpus audit rather than chosen. Across every input file
#: in the repository that carries both a transcribed water load and a derivable
#: pool -- 52 file-stages spanning piezometric lines, seepage boundary levels,
#: reservoir series and both rapid-drawdown stages -- the worst disagreement is
#: **2.4e-4**, and it is vp069, whose transcribed shoreline vertex is rounded to
#: a tenth of a foot. The median is 0: a faithful transcription and the
#: derivation agree bit for bit, because both measure the same two polylines.
#: 2% is therefore a fence at 82x the worst genuine residual, not a fitted
#: number, and nothing in the corpus sits between the two.
WATER_MATCH_TOL = 0.02


def block_resultant(block):
    """The total force per unit width of one dload block (the area under it).

    Integrated along the loaded surface, not along x, because that is where the
    pressure acts -- the two differ by sec(beta) on a sloping face.
    """
    total = 0.0
    pts = list(block or [])
    for a, b in zip(pts, pts[1:]):
        try:
            length = ((float(b["X"]) - float(a["X"])) ** 2
                      + (float(b["Y"]) - float(a["Y"])) ** 2) ** 0.5
            total += 0.5 * (float(a["Normal"]) + float(b["Normal"])) * length
        except (TypeError, ValueError, KeyError):
            continue
    return total


def _block_xs(block):
    return [float(p["X"]) for p in (block or []) if "X" in p]


def _p_at(block, x):
    """Pressure of a dload block at x, linearly interpolated, zero off its ends."""
    pts = sorted(((float(p["X"]), float(p["Normal"])) for p in (block or [])),
                 key=lambda t: t[0])
    if len(pts) < 2 or x < pts[0][0] or x > pts[-1][0]:
        return 0.0
    for (x1, p1), (x2, p2) in zip(pts, pts[1:]):
        if x1 <= x <= x2:
            if x2 == x1:
                return max(p1, p2)
            return p1 + (p2 - p1) * (x - x1) / (x2 - x1)
    return 0.0


def compare_water_block(block, derived_block):
    """How far one transcribed block is from one derived block, in three measures.

    Returns ``{"pressure", "resultant", "reach", "worst"}`` -- each a relative
    difference, and ``worst`` the largest of them. All three are checked because
    each can be right while another is wrong: two blocks can carry the same total
    force over different reaches, or the same reach at different depths.

    ``pressure`` is the largest pointwise gap between the two pressure diagrams,
    sampled at the union of both blocks' vertices (which is every point at which
    either can bend), relative to the derived peak.
    """
    du = block_resultant(derived_block)
    peak = max([abs(float(p["Normal"])) for p in (derived_block or [])] or [0.0])
    xs_d, xs_u = _block_xs(derived_block), _block_xs(block)
    if not xs_d or not xs_u or peak <= 0.0:
        return {"pressure": float("inf"), "resultant": float("inf"),
                "reach": float("inf"), "worst": float("inf")}
    span = max(xs_d) - min(xs_d)
    span = span if span > 0 else 1.0
    reach = max(abs(min(xs_u) - min(xs_d)), abs(max(xs_u) - max(xs_d))) / span
    pressure = max(abs(_p_at(block, x) - _p_at(derived_block, x))
                   for x in sorted(set(xs_d) | set(xs_u))) / peak
    resultant = abs(block_resultant(block) - du) / abs(du) if du else float("inf")
    return {"pressure": pressure, "resultant": resultant, "reach": reach,
            "worst": max(pressure, resultant, reach)}


def match_water_blocks(blocks, derived, tol=WATER_MATCH_TOL):
    """Which transcribed blocks reproduce the derivation, and how closely.

    The one piece of logic behind three consumers, so they can never disagree
    about what "this block is the reservoir" means: the auto-mode double-count
    warning (a user block that matches the derivation is the reservoir entered
    twice), the switch-to-auto remedy (the blocks it may remove are exactly the
    ones that match), and the corpus migration (same rule, applied per file).

    Each derived block is paired with its closest transcribed block, best first,
    and the pair is a match when every measure of :func:`compare_water_block` is
    within ``tol``. Nothing is matched twice.

    Returns ``{"pairs": [(i, j, measures)], "unmatched": [i, ...],
    "missing": [j, ...], "worst": float}`` -- ``i`` indexing ``blocks``, ``j``
    indexing ``derived``, ``unmatched`` the transcribed blocks that are not the
    reservoir (a surcharge, a footing: user data, never removed), and ``missing``
    the derived blocks nothing transcribed.
    """
    blocks = list(blocks or [])
    derived = list(derived or [])
    scored = []
    for j, d in enumerate(derived):
        for i, b in enumerate(blocks):
            scored.append((compare_water_block(b, d)["worst"], i, j))
    scored.sort(key=lambda t: t[0])
    pairs, used_b, used_d = [], set(), set()
    for score, i, j in scored:
        if i in used_b or j in used_d or score > tol:
            continue
        used_b.add(i)
        used_d.add(j)
        pairs.append((i, j, compare_water_block(blocks[i], derived[j])))
    pairs.sort(key=lambda t: t[0])
    return {"pairs": pairs,
            "unmatched": [i for i in range(len(blocks)) if i not in used_b],
            "missing": [j for j in range(len(derived)) if j not in used_d],
            "worst": max([p[2]["worst"] for p in pairs], default=0.0)}


# ---------------------------------------------------------------------------
# The global water table
#
# One phreatic surface per problem, independent of any material's pore-pressure
# option (the "sidecar" model). It governs unit weight — gamma_sat below it,
# moist gamma above — and the effective overburden an overburden-dependent
# pullout law reads. Base pore pressure stays with the per-material u option.
#
# Source precedence: a seepage solution's u = 0 contour beats a hand-drawn
# piezometric line. Rapid drawdown deliberately keys weight to the PRE-drawdown
# (stage-1) surfaces, never the staged ones — the premise of rapid drawdown is
# that the soil stays saturated while the pore pressures fall — so only
# 'seep_u' and 'piezo_line' are read here, never their stage-2 twins.
# ---------------------------------------------------------------------------

def seep_water_table_profile(slope_data):
    """(xs, ys) sample of the u = 0 contour of the model's seepage solution, or
    None when the model carries no seepage field.

    The contour is root-found on the UNCLAMPED signed field by bisection, on a
    201-point grid spanning the mesh. A column that is wet at its top reads the
    mesh top; one that is dry at its bottom reads the mesh bottom, so a fully
    saturated or fully unsaturated column falls out of the same code.

    Cached on ``slope_data['_water_table_profile']`` — the slice generator and
    the reinforcement pullout law ask for the same surface, and bisecting a
    finite element field 201 times is not something to repeat per trial surface.
    """
    import numpy as np
    from .mesh import interpolate_at_point

    cache = slope_data.get('_water_table_profile')
    if cache is not None:
        return cache
    mesh = slope_data.get('mesh')
    seep_u = slope_data.get('seep_u')
    if mesh is None or seep_u is None:
        return None

    seep_u = np.asarray(seep_u, dtype=float)
    nodes = np.asarray(mesh['nodes'], dtype=float)
    xs_grid = np.linspace(nodes[:, 0].min(), nodes[:, 0].max(), 201)
    y_lo_all, y_hi_all = nodes[:, 1].min(), nodes[:, 1].max()
    wt = np.full(xs_grid.shape, np.nan)
    for k, xg in enumerate(xs_grid):
        def _u(yy):
            val, found = interpolate_at_point(
                mesh['nodes'], mesh['elements'], mesh['element_types'],
                seep_u, (xg, yy), return_found=True)
            return (val, found)
        u_lo, f_lo = _u(y_lo_all + 1e-6)
        u_hi, f_hi = _u(y_hi_all - 1e-6)
        if not (f_lo or f_hi):
            continue
        if f_hi and u_hi >= 0:
            wt[k] = y_hi_all          # fully saturated column
            continue
        if f_lo and u_lo <= 0:
            wt[k] = y_lo_all          # fully unsaturated column
            continue
        lo, hi = y_lo_all, y_hi_all
        for _ in range(30):           # bisect u(y) = 0
            mid = 0.5 * (lo + hi)
            um, fm = _u(mid)
            if not fm:
                hi = mid
                continue
            if um > 0:
                lo = mid
            else:
                hi = mid
        wt[k] = 0.5 * (lo + hi)
    cache = (xs_grid, wt)
    slope_data['_water_table_profile'] = cache
    return cache


def water_table_y(slope_data, x):
    """Water-table elevation at x (scalar or array), NaN where the model does
    not define one there.

    Follows the precedence above: the seepage solution's u = 0 contour first,
    the piezometric line second. Off the ends of a piezometric line the answer
    is NaN rather than the nearest endpoint, so a column outside the line's own
    extent is treated as undefined instead of silently flooded or drained.
    """
    import numpy as np
    x = np.asarray(x, dtype=float)
    prof = seep_water_table_profile(slope_data)
    if prof is not None:
        xs_grid, wt = prof
        return np.interp(x, xs_grid, wt)
    piezo = slope_data.get('piezo_line') or []
    if len(piezo) >= 2:
        pts = np.asarray(piezo, dtype=float)
        order = np.argsort(pts[:, 0])
        return np.interp(x, pts[order, 0], pts[order, 1],
                         left=np.nan, right=np.nan)
    return np.full(x.shape, np.nan)
