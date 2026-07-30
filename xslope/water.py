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

Ponded water is the one thing every commercial slope-stability program stores
implicitly — SLOPE/W, Slide2 and RS2 all derive the weight of water standing on
the slope from the water surface, never as an object of its own. xslope needs it
as an explicit distributed load, or the reservoir is simply lost from the
statics. This module holds the geometry that recovers it, so
geostudio.py, slide2.py and rs2.py all convert water the same way:

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
writes to a model or to a file. Two callers consume them. The preflight remedy
(:mod:`xslope.remedies`) offers the derived blocks for a manual-mode model that
is missing its reservoir, and the automatic water-load mode derives the same
blocks at solve time — one derivation, so a load offered by a button and a load
applied by the engine can never be different loads.

Where a boundary condition's level is time-varying, its value comes from the
transient march's own series evaluator, so a derived load and the seepage field
it accompanies can never disagree about where the pool was at that moment.
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

    # Sample where either line has a vertex, so no break in either is missed.
    xs = sorted({x for x, _ in ground_surface.coords} | {x for x, _ in piezo.coords})
    samples = []
    for x in xs:
        yg, yp = _y_on(ground_surface, x), _y_on(piezo, x)
        if yg is not None and yp is not None:
            samples.append((x, yg, yp - yg))          # depth of water over the ground
    if not samples:
        return []

    blocks, run = [], []
    for i, (x, yg, d) in enumerate(samples):
        if d > 1e-9:
            if not run and i > 0:                     # entering the water: add the waterline
                px, pyg, pd = samples[i - 1]
                t = -pd / (d - pd) if (d - pd) else 0.0
                xw = px + t * (x - px)
                run.append({"X": xw, "Y": _y_on(ground_surface, xw) or pyg, "Normal": 0.0})
            run.append({"X": x, "Y": yg, "Normal": gamma_water * d})
        elif run:                                     # leaving the water: close at the waterline
            px, pyg, pd = samples[i - 1]
            t = pd / (pd - d) if (pd - d) else 0.0
            xw = px + t * (x - px)
            run.append({"X": xw, "Y": _y_on(ground_surface, xw) or yg, "Normal": 0.0})
            if len(run) >= 2:
                blocks.append(run)
            run = []
    if len(run) >= 2:
        blocks.append(run)
    return blocks


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
    that overlap take the higher level. The returned line is sampled at every
    ground vertex, which makes the crossing that
    :func:`material_above_ground_dload` interpolates exact rather than approximate.

    Returns ``[]`` when no block describes standing water above the ground.
    """
    ground = _as_coords(ground_surface)
    if len(ground) < 2 or not _ascending(ground):
        return []
    pools = bc_pool_levels(bc, ground_surface, slope_data, time)
    if not pools:
        return []

    xs = sorted({p[0] for p in ground} | {x for _, anchors, _ in pools for x in anchors})
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

    Returns a dict with ``points`` (the water line, possibly empty), ``source`` (a
    sheet-vocabulary description of where it came from) and ``reason`` (why there
    is no line, when there is none).
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
    sheet = "seep bc" if stage == 1 else "seep bc (2)"
    if bc.get("specified_heads"):
        pts = seep_bc_water_line(sd.get("ground_surface"), bc, sd, time)
        levels = bc_pool_levels(bc, sd.get("ground_surface"), sd, time)
        if pts:
            names = ", ".join(label for _l, _a, label in levels)
            return {"points": pts, "source": f"{sheet} sheet, {names}", "reason": ""}
        return {"points": [], "source": f"{sheet} sheet",
                "reason": (f"the {sheet} sheet's head boundaries do not put water above "
                           f"the ground surface")}

    key = "piezo_line" if stage == 1 else "piezo_line2"
    name = "Piezometric Line 1" if stage == 1 else "Piezometric Line 2"
    pz = _as_coords(sd.get(key))
    if len(pz) < 2:
        return {"points": [], "source": "",
                "reason": (f"the model defines neither seepage head boundaries nor "
                           f"{name}, so nothing says where the water stands")}
    if not _ascending(pz):
        return {"points": [], "source": f"piezo sheet, {name}",
                "reason": (f"piezo sheet, {name} does not run left to right, and a "
                           f"water surface can only be measured along a line whose X "
                           f"values increase")}
    return {"points": pz, "source": f"piezo sheet, {name}", "reason": ""}


def derive_water_loads(slope_data, stage=1, time=None):
    """The distributed load the standing water applies, derived from the model.

    The whole derivation in one call: find where the water stands
    (:func:`water_line_for_stage`), measure it against the ground surface
    (:func:`material_above_ground_dload`), and report the result with enough
    detail to describe it before anything is changed -- which is what lets a
    remedy state *"1 block, 262 peak, over x = 346.7 to 594"* while the user is
    still deciding.

    Returns a dict:

    ``blocks``      dload blocks in xslope's own form (``X``/``Y``/``Normal``)
    ``dirs``        the matching direction words -- water always acts ``normal``
    ``source``      which sheet the water line came from
    ``peak``        the largest pressure in the derived blocks
    ``reach``       ``(x_first, x_last)`` of the loaded stretch
    ``resultant``   the total load per unit width, for comparing against a
                    transcribed block
    ``reason``      why nothing was derived, when nothing was
    """
    sd = slope_data or {}
    empty = {"blocks": [], "dirs": [], "source": "", "peak": 0.0,
             "reach": None, "resultant": 0.0, "reason": ""}
    try:
        gamma_w = float(sd.get("gamma_water"))
    except (TypeError, ValueError):
        gamma_w = None
    if not gamma_w or gamma_w <= 0:
        return dict(empty, reason="main sheet D10 (Unit weight of water) is blank or "
                                  "not positive, so the water has no weight to apply")

    found = water_line_for_stage(sd, stage=stage, time=time)
    if not found["points"]:
        return dict(empty, source=found["source"], reason=found["reason"])

    ground = _as_coords(sd.get("ground_surface"))
    blocks = material_above_ground_dload(LineString(ground), found["points"], gamma_w)
    if not blocks:
        return dict(empty, source=found["source"],
                    reason=(f"{found['source']} does not stand above the ground "
                            f"surface anywhere, so there is no water load to add"))

    peak = max(float(p["Normal"]) for blk in blocks for p in blk)
    xs = [float(p["X"]) for blk in blocks for p in blk]
    resultant = 0.0
    for blk in blocks:
        for a, b in zip(blk, blk[1:]):
            resultant += 0.5 * (a["Normal"] + b["Normal"]) * (
                (b["X"] - a["X"]) ** 2 + (b["Y"] - a["Y"]) ** 2) ** 0.5
    return {"blocks": blocks, "dirs": ["normal"] * len(blocks),
            "source": found["source"], "peak": peak,
            "reach": (min(xs), max(xs)), "resultant": resultant, "reason": ""}
