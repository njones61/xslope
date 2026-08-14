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

"""Generators -- starting surfaces derived from the model's own geometry.

A limit-equilibrium search has to start somewhere, and where it starts decides
what it finds: the adaptive search is a local optimiser, so a starting circle in
the wrong family converges to a local minimum and reports it without complaint.
Choosing that starting set is a small piece of standard practice -- a circle
through the toe, one tangent to the base of each layer, a flat one skimming a
cohesionless face -- and it was written down only as prose, to be re-derived by
hand for every model. This module is that practice as one implementation.

Three callers share it: the preflight remedy for a model with no failure surface
(:mod:`xslope.remedies`), a **Generate** button in Studio, and the analysis skill,
which cites this function instead of restating the arithmetic.

What it needs, and what it builds
---------------------------------

:func:`slope_geometry` reads the ground surface and reports the slope faces it
contains -- toe, crest, height, facing direction, and the face's own segments.
Nothing in the package derived a toe and a crest before this; the search's grid
seeding took the extreme high and low points, which is right for a single face
and wrong for a dam, whose two faces share one crest and have a toe each.

:func:`generate_starting_circles` then applies the rules:

* **Center.** ``Xo`` halfway between toe and crest, ``Yo`` at the toe elevation
  plus twice the slope height.
* **A toe circle** -- one passing *through* the toe, ``R = dist(center, toe)``.
  This is deliberately not a circle whose bottom sits at the toe *elevation*; the
  two are different surfaces and the distinction is easy to lose.
* **A circle tangent to the base of each distinct material layer**, with ``Depth``
  set to that layer's base elevation.
* **A skimming circle on a cohesionless face.** Where the material exposed on the
  face has ``c = 0``, the critical surface is an arbitrarily shallow face-parallel
  slide, and a toe-and-base seeded search converges to a deep local minimum with a
  non-conservatively high answer. A large-radius circle approximates the plane;
  its center lands far outside the model, which is expected and correct.
* **Daylighting.** Every generated circle must cut the ground surface twice,
  inside the model rather than at a vertical edge, and must not dip below the
  bottom of the domain -- so the set can be handed straight to a search.

``Depth`` is an **elevation**, not a depth below ground: it is the elevation of
the circle's lowest point, and the radius is ``Yo - Depth``. That is the circles
sheet's own convention and a persistent source of confusion.

The generated circles are a *starting set*, not an answer. They exist to give the
search a seed in every family that could win.

The non-circular surface
------------------------

Some slopes fail along a weak layer rather than along their own geometry, and no
circle passes through that mechanism: the surface runs flat inside the seam for
most of its length and turns up sharply at each end. :func:`rank_weak_zones` and
:func:`generate_noncircular_surface` build that shape.

Which layer is the weak one is answered by ranking every zone on the shear
strength it can mobilise **at the stress it actually carries** -- not on cohesion
and not on friction angle, neither of which is comparable between materials. That
one quantity spans every strength model the ``mat`` sheet offers, because each
reduces to a strength at a normal stress: undrained ``cp`` is already one,
Hoek-Brown is linearised at that stress by :func:`xslope.hoekbrown.hb_tangent`, a
power envelope is evaluated on it, and an ``elastic`` material cannot fail and is
left out.

When one zone is clearly the weakest, the generator seeds on it and says why. When
two are comparable, it returns the ranked candidates rather than guessing -- and
that is not a shortfall of the ranking, it is a property of the sections. Guo &
Griffiths' embankment-over-foundation pair fails deep on one set of numbers and
shallow on another with the *same* two zones, and the corpus contains both.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

#: A face qualifies as its own slope when its relief reaches this fraction of the
#: tallest face's. A dam has two faces worth seeding; a 0.3 m step at the far end
#: of the section is not a slope, and seeding it would spend the search's budget
#: on a mechanism nobody is asking about.
_FACE_RELIEF_FRAC = 0.25

#: Radius-to-chord ratio for a skimming circle. The skill's measured band is 15-20:
#: below about 10 the arc is too curved to find a face-parallel mechanism, and
#: above about 25 ``generate_slices`` rejects the result as a flat arc.
_SKIM_K = 15.0

#: How far a skimming circle is sunk below the face segment it skims, as a
#: fraction of the segment's length. The construction passes exactly through the
#: segment's two endpoints, and a surface that meets the ground exactly at a
#: vertex is a grazing intersection: the slicer finds one crossing where there
#: should be two and refuses the surface. Sinking it a thousandth of the segment
#: makes both crossings ordinary ones, and moves the arc by far less than a slice
#: width.
_SKIM_SINK = 1e-3

#: Center elevations to try, as multiples of the slope height above the toe. The
#: first is the rule; the rest are a fallback for a section transcribed too narrow
#: for it, where a standard circle would daylight past the model's own edge. A
#: lowered center is reported, because a section that needs one is cropped -- the
#: real repair is to widen the geometry.
_CENTER_LADDER = (2.0, 1.5, 1.0, 0.75)

#: Cohesion at or below this fraction of the model's own stress scale
#: (gamma * H) counts as cohesionless. Relative, so it travels between unit
#: systems rather than being a number that means one thing in kPa and another
#: in psf.
_COHESION_FRAC = 1e-4


@dataclass(frozen=True)
class SlopeFace:
    """One face of a slope: where it starts, where it ends, and how steep it is."""

    toe: tuple
    crest: tuple
    height: float
    right_facing: bool
    segments: tuple

    @property
    def steepest(self):
        """The segment with the largest inclination -- the one that governs.

        On a benched face this is a single bench, deliberately: chording crest to
        toe averages the benches away and misses the mechanism entirely.
        """
        return max(self.segments, key=lambda s: _inclination(s))

    @property
    def steepest_angle(self):
        """Inclination of :attr:`steepest`, in degrees from horizontal."""
        return math.degrees(math.atan(_inclination(self.steepest)))

    def __str__(self):
        return (f"{'right' if self.right_facing else 'left'}-facing face, toe "
                f"({self.toe[0]:g}, {self.toe[1]:g}) to crest ({self.crest[0]:g}, "
                f"{self.crest[1]:g}), height {self.height:g}")


@dataclass(frozen=True)
class SlopeGeometry:
    """The slope facts a generator, a plot frame or a search window needs."""

    ground: tuple
    faces: tuple
    floor: float

    @property
    def primary(self):
        """The tallest face -- the one a single-face rule means by "the slope"."""
        return max(self.faces, key=lambda f: f.height)

    @property
    def toe(self):
        return self.primary.toe

    @property
    def crest(self):
        return self.primary.crest

    @property
    def height(self):
        return self.primary.height

    @property
    def right_facing(self):
        return self.primary.right_facing

    @property
    def extent(self):
        """``(x_min, x_max)`` of the ground surface -- the model's own width."""
        xs = [p[0] for p in self.ground]
        return (min(xs), max(xs))


def _inclination(segment):
    (x1, y1), (x2, y2) = segment
    dx = abs(x2 - x1)
    return abs(y2 - y1) / dx if dx > 0 else float("inf")


def _ground_coords(slope_data):
    gs = (slope_data or {}).get("ground_surface")
    if gs is None or (hasattr(gs, "is_empty") and gs.is_empty):
        return []
    if hasattr(gs, "coords"):
        return [(float(x), float(y)) for x, y in gs.coords]
    return [(float(p[0]), float(p[1])) for p in gs]


def slope_geometry(slope_data):
    """The slope faces of a model's ground surface.

    A face is a run of the ground surface that gains or loses elevation in one
    direction, with flat benches inside it treated as part of the face rather than
    as its end -- which is what lets a benched face read as one slope. A dam
    therefore reports two faces sharing a crest, and a single embankment reports
    one.

    Faces whose relief is a small fraction of the tallest one's are dropped: they
    are steps and shoulders, not slopes.

    Returns
    -------
    SlopeGeometry

    Raises
    ------
    ValueError
        When the ground surface is missing, too short, or entirely flat -- there is
        no toe and no crest to derive, and inventing one would be worse than
        saying so.
    """
    ground = _ground_coords(slope_data)
    if len(ground) < 2:
        return _fail("this model defines no ground surface, so there is no toe or "
                     "crest to derive a starting surface from")
    if any(b[0] < a[0] for a, b in zip(ground, ground[1:])):
        ground = sorted(ground, key=lambda p: p[0])
    relief = max(p[1] for p in ground) - min(p[1] for p in ground)
    if relief <= 0:
        return _fail("the ground surface is flat, so the model has no slope to build "
                     "a starting surface for")

    tol = 1e-9 * max(1.0, relief)
    runs, run, direction = [], [], 0
    for a, b in zip(ground, ground[1:]):
        dy = b[1] - a[1]
        step = 0 if abs(dy) <= tol else (1 if dy > 0 else -1)
        if step == 0:
            if run:
                run.append((a, b))          # a bench inside a face, kept
            continue
        if direction == 0 or step == direction:
            run.append((a, b))
            direction = step
        else:
            runs.append((direction, _trim_flats(run)))
            run, direction = [(a, b)], step
    if run:
        runs.append((direction, _trim_flats(run)))

    faces = []
    for step, segs in runs:
        segs = tuple(s for s in segs)
        if not segs:
            continue
        start, end = segs[0][0], segs[-1][1]
        height = abs(end[1] - start[1])
        if height <= tol:
            continue
        toe, crest = (start, end) if end[1] > start[1] else (end, start)
        faces.append(SlopeFace(toe=toe, crest=crest, height=height,
                               right_facing=(step < 0),
                               segments=tuple(s for s in segs
                                              if abs(s[1][1] - s[0][1]) > tol)))
    if not faces:
        return _fail("no slope face could be found on the ground surface")
    tallest = max(f.height for f in faces)
    faces = [f for f in faces if f.height >= _FACE_RELIEF_FRAC * tallest]

    floor = _model_floor(slope_data, ground)
    return SlopeGeometry(ground=tuple(ground), faces=tuple(faces), floor=floor)


def _trim_flats(run):
    """Drop trailing flat segments, so a face ends where it stops climbing."""
    out = list(run)
    while out and abs(out[-1][1][1] - out[-1][0][1]) <= 0:
        out.pop()
    return out


def _fail(message):
    raise ValueError(message)


def _model_floor(slope_data, ground):
    """The bottom of the model: the domain polygon's floor, or the deepest zone."""
    dom = (slope_data or {}).get("domain_polygon")
    if dom is not None and hasattr(dom, "bounds"):
        return float(dom.bounds[1])
    lows = []
    for poly in (slope_data or {}).get("polygons") or []:
        shape = poly.get("polygon") if isinstance(poly, dict) else poly
        if shape is not None and hasattr(shape, "bounds"):
            lows.append(float(shape.bounds[1]))
    if lows:
        return min(lows)
    md = (slope_data or {}).get("max_depth")
    try:
        return float(md)
    except (TypeError, ValueError):
        return min(p[1] for p in ground)


# ---------------------------------------------------------------------------
# Circle geometry
# ---------------------------------------------------------------------------

def _crossings(xo, yo, r, ground):
    """Where a circle daylights on the ground surface, as x values, left to right.

    Solved segment by segment rather than by intersecting a polygonal buffer: a
    skimming circle's radius runs to twenty times the model's width, and a
    buffered approximation of it would miss the crossing by more than the depth it
    is meant to resolve.

    Only crossings **below the circle's center** count, which is the slicer's own
    rule (``slice.circle_polyline_intersections``): the failure surface is the
    bottom semicircle, so a crossing above the equator is not a point the analysed
    surface reaches. Asking the question the same way is what makes a circle this
    generator keeps a circle ``generate_slices`` accepts.
    """
    xs = []
    for (x1, y1), (x2, y2) in zip(ground, ground[1:]):
        dx, dy = x2 - x1, y2 - y1
        fx, fy = x1 - xo, y1 - yo
        a = dx * dx + dy * dy
        if a <= 0:
            continue
        b = 2 * (fx * dx + fy * dy)
        c = fx * fx + fy * fy - r * r
        disc = b * b - 4 * a * c
        if disc < 0:
            continue
        root = math.sqrt(disc)
        for t in ((-b - root) / (2 * a), (-b + root) / (2 * a)):
            if -1e-12 <= t <= 1 + 1e-12 and (y1 + t * dy) < yo:
                xs.append((x1 + t * dx, y1 + t * dy))
    xs.sort()
    merged = []
    for p in xs:
        if not merged or abs(p[0] - merged[-1][0]) > 1e-9:
            merged.append(p)
    return merged


def _analysed_span(xo, yo, r, ground):
    """The two crossings the slicer would actually analyse, as ``(x_left, x_right)``.

    A circle that dips in and out of a benched face crosses the ground more than
    twice, and ``slice.get_sorted_intersections`` keeps the first pair when the
    left end is higher and the last pair otherwise -- the pair on the crest side.
    Mirrored here rather than approximated, because the extreme pair spans ground
    the analysed surface never touches, and testing containment against *that*
    would reject exactly the skimming circles this generator exists to produce.
    """
    pts = _crossings(xo, yo, r, ground)
    if len(pts) < 2:
        return None
    if len(pts) == 2:
        pair = pts
    elif pts[0][1] > pts[-1][1]:
        pair = pts[:2]
    else:
        pair = pts[-2:]
    return pair[0][0], pair[1][0]


def _arc_min_y(xo, yo, r, x_left, x_right):
    """The lowest elevation the circle reaches between two x values."""
    if x_left <= xo <= x_right:
        return yo - r
    x = x_left if abs(x_left - xo) < abs(x_right - xo) else x_right
    d = min(abs(x - xo), r)
    return yo - math.sqrt(max(0.0, r * r - d * d))


def _arc_points(xo, yo, r, x_left, x_right, n=64):
    """The failure arc between two crossings, sampled evenly in angle.

    Sampled in angle rather than in x because the arc is near-vertical where it
    meets the ground, and an x-sampled polyline would cut the corner exactly where
    the containment test needs to be tight.
    """
    thetas = []
    for x in (x_left, x_right):
        dx = max(-r, min(r, x - xo))
        thetas.append(-math.acos(max(-1.0, min(1.0, dx / r))))   # lower half
    t0, t1 = sorted(thetas)
    return [(xo + r * math.cos(t0 + (t1 - t0) * i / n),
             yo + r * math.sin(t0 + (t1 - t0) * i / n)) for i in range(n + 1)]


def _daylights(circle, ground, floor, tol, domain=None):
    """True when a circle cuts the ground twice, inside the model and inside the domain.

    Three ways a trial circle is not a usable starting surface, and all three are
    rejected here rather than handed to a search that would have to score them:
    it never reaches the ground surface; it leaves through a vertical model edge
    rather than daylighting on the ground; or its arc runs outside the domain
    polygon, which is the refusal ``generate_slices`` raises. Testing against the
    domain rather than against the flat bottom of its bounding box is what makes
    this right on a model whose floor slopes.
    """
    xo, yo, r = circle["Xo"], circle["Yo"], circle["R"]
    if r <= 0:
        return False
    span = _analysed_span(xo, yo, r, ground)
    if span is None:
        return False
    x0, x1 = span
    gx0, gx1 = ground[0][0], ground[-1][0]
    if x0 <= gx0 + tol or x1 >= gx1 - tol:
        return False
    if (x1 - x0) <= tol:
        return False
    if domain is not None and hasattr(domain, "covers"):
        from shapely.geometry import LineString
        arc = LineString(_arc_points(xo, yo, r, x0, x1))
        try:
            return bool(domain.buffer(tol).covers(arc))
        except Exception:                       # pragma: no cover - defensive
            pass
    return _arc_min_y(xo, yo, r, x0, x1) >= floor - tol


def _skimming_circle(a, b, k=_SKIM_K):
    """A large-radius circle whose arc skims just under the face segment ``a``-``b``.

    The center sits on the outward side of the face, so the arc sags into the
    slope. ``k = R / L``; the center lands far outside the model, which is what
    makes the arc nearly planar.
    """
    (x1, y1), (x2, y2) = a, b
    mx, my = 0.5 * (x1 + x2), 0.5 * (y1 + y2)
    dx, dy = x2 - x1, y2 - y1
    length = math.hypot(dx, dy)
    if length <= 0:
        return None
    nx, ny = -dy / length, dx / length
    if ny < 0:
        nx, ny = -nx, -ny                       # the normal must point out of the slope
    r = k * length
    h = math.sqrt(max(0.0, r * r - (length / 2.0) ** 2)) + _SKIM_SINK * length
    cx, cy = mx + h * nx, my + h * ny
    return {"Xo": float(cx), "Yo": float(cy), "R": float(r), "Depth": float(cy - r)}


# ---------------------------------------------------------------------------
# Material facts the rules need
# ---------------------------------------------------------------------------

def _layer_base_elevations(slope_data, tol):
    """The distinct base elevations of the model's material zones, deepest first.

    One circle is seeded tangent to each: the base of a layer is where a
    mechanism that cannot cut the layer below it has to run, and those are the
    competing families a depth optimiser started in one of them never sees.
    """
    lows = []
    for poly in (slope_data or {}).get("polygons") or []:
        shape = poly.get("polygon") if isinstance(poly, dict) else poly
        if shape is None or not hasattr(shape, "bounds"):
            continue
        lows.append(float(shape.bounds[1]))
    out = []
    for y in sorted(lows):
        if not out or abs(y - out[-1]) > tol:
            out.append(y)
    return out


def _material_at(slope_data, point):
    """The material whose zone contains a point, or ``None``."""
    from shapely.geometry import Point
    materials = (slope_data or {}).get("materials") or []
    p = Point(point[0], point[1])
    for poly in (slope_data or {}).get("polygons") or []:
        shape = poly.get("polygon") if isinstance(poly, dict) else poly
        mat_id = poly.get("mat_id") if isinstance(poly, dict) else None
        if shape is None or mat_id is None:
            continue
        try:
            if shape.contains(p):
                return materials[int(mat_id)]
        except (IndexError, TypeError, ValueError):
            return None
    return None


def _is_cohesionless(material, scale):
    """True when a material carries no cohesion on the Mohr-Coulomb envelope.

    Only the linear envelope is read: a curved one (``hb``, ``pow``) has no ``c``
    to be zero, and an undrained ``cp`` material is cohesion and nothing else.
    """
    if not material:
        return False
    option = str(material.get("option") or "mc").strip().lower()
    if option not in ("mc", ""):
        return False
    try:
        c = float(material.get("c") or 0.0)
        phi = float(material.get("phi") or 0.0)
    except (TypeError, ValueError):
        return False
    return abs(c) <= _COHESION_FRAC * scale and phi > 0


def _face_material(slope_data, face, scale):
    """The material exposed on a face, sampled just inside its steepest segment."""
    (x1, y1), (x2, y2) = face.steepest
    mx, my = 0.5 * (x1 + x2), 0.5 * (y1 + y2)
    dx, dy = x2 - x1, y2 - y1
    length = math.hypot(dx, dy)
    if length <= 0:
        return None
    nx, ny = -dy / length, dx / length
    if ny > 0:                                  # step INTO the slope, not out of it
        nx, ny = -nx, -ny
    step = max(1e-9, 1e-3 * face.height)
    return _material_at(slope_data, (mx + step * nx, my + step * ny))


# ---------------------------------------------------------------------------
# The generator
# ---------------------------------------------------------------------------

def generate_starting_circles(slope_data, faces="significant", skim=True, report=False):
    """A starting set of trial circles, derived from the model's own geometry.

    Parameters
    ----------
    slope_data : dict
        A loaded model. Nothing is mutated; only the geometry and the material
        table are read.
    faces : {"significant", "primary"}
        Which faces to seed. ``"significant"`` (the default) seeds every face
        whose relief is a reasonable fraction of the tallest -- a dam gets a set on
        each of its two faces, because either can be the critical one.
        ``"primary"`` seeds the tallest face alone.
    skim : bool
        Add the large-radius skimming circle where the face material is
        cohesionless. On by default: without it a search on a ``c = 0`` face
        converges to a deep local minimum and reports a non-conservatively high
        factor of safety.
    report : bool
        Return a dict carrying the circles, a one-line summary and -- when nothing
        could be generated -- the reason, instead of the bare list. This is what a
        remedy or a Generate button needs, since it has to say what it will do
        before it does it, and dim itself with a reason when it cannot.

    Returns
    -------
    list of dict
        Circles in the sheet's own form: ``Xo``, ``Yo``, ``Depth`` (the elevation
        of the lowest point) and ``R``. Or, with ``report=True``, a dict with
        ``circles``, ``summary`` and ``reason``.

    Examples
    --------
    >>> circles = generate_starting_circles(slope_data)          # doctest: +SKIP
    >>> slope_data["circles"] = circles                          # doctest: +SKIP
    """
    try:
        geom = slope_geometry(slope_data)
    except ValueError as exc:
        if report:
            return {"circles": [], "summary": "", "reason": str(exc)}
        return []

    ground = list(geom.ground)
    span = geom.extent[1] - geom.extent[0]
    tol = 1e-6 * max(span, 1.0)
    chosen = [geom.primary] if faces == "primary" else list(geom.faces)
    scale = _stress_scale(slope_data, geom)
    bases = _layer_base_elevations(slope_data, tol=1e-6 * max(geom.height, 1.0))

    domain = (slope_data or {}).get("domain_polygon")
    circles, notes, rejected = [], [], 0
    for face in chosen:
        cohesionless = skim and _is_cohesionless(_face_material(slope_data, face, scale),
                                                 scale)
        for mult in _CENTER_LADDER:
            xo = 0.5 * (face.toe[0] + face.crest[0])
            yo = face.toe[1] + mult * face.height

            candidates = []
            r_toe = math.hypot(xo - face.toe[0], yo - face.toe[1])
            candidates.append(("toe", {"Xo": xo, "Yo": yo, "R": r_toe,
                                       "Depth": yo - r_toe}))
            for base in bases:
                r = yo - base
                if r <= 0 or base >= face.crest[1]:
                    continue
                candidates.append(("base", {"Xo": xo, "Yo": yo, "R": r, "Depth": base}))
            if cohesionless:
                sk = _skimming_circle(*face.steepest)
                if sk is not None:
                    candidates.append(("skim", sk))

            kept, skimmed, dropped = [], False, 0
            for kind, circle in candidates:
                if not _daylights(circle, ground, geom.floor, tol, domain):
                    dropped += 1
                    continue
                if any(abs(c["Xo"] - circle["Xo"]) <= tol
                       and abs(c["Yo"] - circle["Yo"]) <= tol
                       and abs(c["Depth"] - circle["Depth"]) <= tol
                       for c in circles + kept):
                    continue
                kept.append({k: float(v) for k, v in circle.items()})
                skimmed = skimmed or kind == "skim"
            if not kept:
                continue
            circles.extend(kept)
            rejected += dropped
            note = (f"{len(kept)} on the {'right' if face.right_facing else 'left'}"
                    f"-facing face (toe at x = {face.toe[0]:g}, height "
                    f"{face.height:g})")
            if skimmed:
                note += (f", one of them skimming its {face.steepest_angle:.0f} degree "
                         f"cohesionless face")
            if mult != _CENTER_LADDER[0]:
                note += (f", with the center at {mult:g} x the slope height rather than "
                         f"the usual {_CENTER_LADDER[0]:g} -- the section is not wide "
                         f"enough for a standard circle to daylight on the ground "
                         f"inside it")
            notes.append(note)
            break

    if not report:
        return circles
    if not circles:
        return {"circles": [], "summary": "",
                "reason": ("no trial circle derived from this slope daylights on the "
                           "ground surface inside the model: every candidate leaves "
                           "through a vertical edge of the section or runs outside the "
                           "domain. The section needs flat ground beyond the toe and "
                           "the crest -- about twice the slope height at each end -- "
                           "for a failure surface to exit on the ground")}
    summary = "; ".join(notes)
    if rejected:
        summary += (f" ({rejected} candidate{'' if rejected == 1 else 's'} dropped for "
                    f"not daylighting inside the model)")
    return {"circles": circles, "summary": summary, "reason": ""}


# ---------------------------------------------------------------------------
# Non-circular: finding the weak zone, and the surface that tracks it
# ---------------------------------------------------------------------------

#: How much weaker the weakest zone must be than the next before the generator
#: picks it without asking: ``tau_weakest / tau_next <= _SEPARATION``.
#:
#: Tuned against the corpus's own weak-seam rows rather than chosen, by seeding a
#: search on every ranked zone of every corpus row that carries a non-circular
#: answer and comparing what it converges to against that row's standing value.
#: The two sides of the bracket:
#:
#: * **Below.** Every ranking at or under 0.56 picks a zone whose seeded search
#:   reaches or beats the row's own answer -- Griffiths & Lane's dipping band at
#:   0.20, ACADS problem 4's weak layer at 0.24, problem 3(b)'s at 0.33, the
#:   UTEXAS soft-clay sample at 0.35, ``vp050``'s soft upper layer at 0.56.
#: * **Above.** The first ranking that picks the WRONG zone is 0.63, on Guo &
#:   Griffiths' embankment-over-foundation pair (``vp103c``/``vp103d``), where the
#:   ranking prefers the embankment and three of the four cases fail deep through
#:   the foundation. ``vp103b`` at 0.67, ``vp103a`` at 0.71 and Clarence Cannon Dam
#:   (``vp034``) at 0.74 fail the same way.
#:
#: 0.60 is the middle of that gap. The ranking is not a perfect discriminator and
#: the threshold cannot make it one: ``vp020``'s seam ranks at 0.75 and is right,
#: while ``vp034`` ranks at 0.74 and is wrong, so no cut separates them. That is
#: the argument for asking rather than guessing, not against it -- and it is why
#: this is not tuned until nothing is ambiguous. Five corpus rows have two
#: genuinely comparable candidates, and on three of them the SECOND-ranked zone is
#: the one carrying the mechanism. Those rows are what the picker is for.
_SEPARATION = 0.60

#: Where the track runs inside the weak zone, as a fraction of the zone's local
#: thickness above its base. Not zero: a surface lying exactly on a material
#: boundary is a slicing hazard, and the base is where the strength contrast the
#: mechanism exploits actually is. Measured across the corpus's weak-seam rows,
#: the converged minimum is flat between about 0.05 and 0.25 of the thickness and
#: rises outside it -- 0.10 is the middle of that plateau and is also the offset
#: the hand-built corpus surfaces use.
_OFFSET_FRAC = 0.10

#: Ceiling on the track's offset above the zone base, as a fraction of the slope
#: height. ``_OFFSET_FRAC`` of the local thickness is the rule, and it is the right
#: rule for a seam; on a zone tens of units thick the same fraction would lift the
#: track clear of the contact whose weakness the mechanism is exploiting. The
#: ceiling keeps "just above the base" meaning just above the base at any zone
#: thickness, which is what the search needs to be started on.
_OFFSET_CAP = 0.02

#: Evenly spaced stations along the track, on top of the zone's own vertices.
#: Deliberately few. The search is a coordinate descent that moves one point at a
#: time, so a densely sampled track is not a better starting surface -- it is a
#: greedier one, and it stalls: measured across the weak-seam rows, going from four
#: stations to six moves ``vp009`` from 0.615 to 0.700 and ``vp103a`` from 1.264 to
#: 1.279, both the wrong way. The corpus's own hand-drawn surfaces carry two to
#: four interior points, and there is a reason for that.
#:
#: A count of stations LAID, not of points kept: one that lands on a level stretch
#: of the track is dropped again, because there it is a move the search cannot
#: make. See :func:`_drop_inert_points`.
_TRACK_POINTS = 4

#: How far an interior point may lie off the straight line through its neighbours
#: and still count as sitting on it, as a fraction of the slope height. Set at
#: float-noise scale on purpose: a point this close to the chord was put there by
#: the even-station sampling, never by a bend in the zone's own base. Measured
#: across the corpus, the point that comes nearest to being dropped and is not --
#: the flattest surviving bend on any model the button reaches -- clears the
#: tolerance by a factor of 4,000, so nothing real is anywhere near it.
_COLLINEAR_TOL = 1e-6

#: Ramp inclinations to try, as offsets applied to the textbook wedge angles
#: (``45 - phi/2`` leaving at the toe, ``45 + phi/2`` entering at the crest).
#: The first entry is the rule; the rest flatten it, for an end where the textbook
#: ramp would leave the model or cut back over the track.
_RAMP_LADDER = (1.0, 0.75, 0.5, 0.3)

#: Hard bounds on a ramp inclination, in degrees. The upper bound is below
#: ``noncircular_search``'s own ``max_base_angle`` (65 deg), which rejects any
#: trial surface with a steeper base: a starting surface the search refuses is
#: worse than no starting surface at all.
_RAMP_MIN_DEG, _RAMP_MAX_DEG = 5.0, 60.0

#: Stations used to sample a zone when ranking it and when hunting its outcrop.
_ZONE_SAMPLES = 25


@dataclass(frozen=True)
class WeakZone:
    """One material zone, with the comparable strength that ranks it.

    ``tau`` is the mobilisable shear strength at the stress the zone actually
    carries -- not a cohesion and not a friction angle, neither of which can be
    compared across materials. It is what makes a ``c = 0, phi = 35`` sand and a
    ``c = 50, phi = 0`` clay rankable against each other at all.
    """

    index: int
    mat_id: int
    name: str
    option: str
    tau: float
    polygon: object
    x_range: tuple
    y_range: tuple

    def __str__(self):
        return (f"{self.name} ({self.tau:.3g} at x = {self.x_range[0]:g} to "
                f"{self.x_range[1]:g}, y = {self.y_range[0]:g} to {self.y_range[1]:g})")


def _vertical_span(polygon, x):
    """``(y_bottom, y_top)`` of a polygon on the vertical line at ``x``."""
    from shapely.geometry import LineString
    x0, y0, x1, y1 = polygon.bounds
    if not (x0 - 1e-12 <= x <= x1 + 1e-12):
        return None
    try:
        cut = polygon.intersection(LineString([(x, y0 - 1.0), (x, y1 + 1.0)]))
    except Exception:                                   # pragma: no cover
        return None
    if cut.is_empty:
        return None
    ys = []
    for part in (cut.geoms if hasattr(cut, "geoms") else [cut]):
        ys += [c[1] for c in getattr(part, "coords", [])]
    if len(ys) < 2:
        return None
    return (min(ys), max(ys))


def _ground_elevation(ground, x):
    """The ground surface elevation at ``x``, or ``None`` outside the section."""
    if x < ground[0][0] - 1e-9 or x > ground[-1][0] + 1e-9:
        return None
    best = None
    for (x1, y1), (x2, y2) in zip(ground, ground[1:]):
        lo, hi = (x1, x2) if x1 <= x2 else (x2, x1)
        if not (lo - 1e-9 <= x <= hi + 1e-9):
            continue
        y = y1 if x2 == x1 else y1 + (y2 - y1) * (x - x1) / (x2 - x1)
        best = y if best is None else max(best, y)
    return best


def _ground_snap(slope_data, ground, x):
    """The ground elevation at ``x``, computed the way the *slicer* computes it.

    An end point has to sit on the ground surface exactly enough that
    ``failure_surface.intersection(ground_surface)`` returns a point, and "exactly
    enough" is not a matter of decimal places -- it is a matter of which
    arithmetic. Interpolating the segment by hand and asking shapely whether the
    result lies on it are two different calculations, and on a long steep ground
    segment they disagree in the last bits: the slicer then reports "expected at
    least 2 intersection points, but got 1" and refuses a surface whose geometry is
    visibly correct. Cutting the ground with a vertical line -- letting shapely
    do the arithmetic -- gives an answer shapely agrees with.
    """
    surface = (slope_data or {}).get("ground_surface")
    if surface is not None and not getattr(surface, "is_empty", True):
        from shapely.geometry import LineString
        try:
            from .slice import get_y_from_intersection
            lo = min(p[1] for p in ground) - 1.0
            hi = max(p[1] for p in ground) + 1.0
            y = get_y_from_intersection(
                surface.intersection(LineString([(x, lo), (x, hi)])))
            if y is not None:
                return float(y)
        except Exception:                                # pragma: no cover
            pass
    return _ground_elevation(ground, x)


def _column_stress(slope_data, ground, x, y):
    """``sum(gamma_i * h_i)`` for the soil column between ``y`` and the ground.

    The same quantity ``generate_slices`` calls ``sum_gam_h`` and divides into the
    slice weight -- one ``gamma`` per material, the water table never entering the
    weight, exactly as the solver does it. Agreeing with the solver matters more
    here than any refinement the solver does not itself make.
    """
    from shapely.geometry import LineString
    gy = _ground_elevation(ground, x)
    if gy is None or gy <= y:
        return 0.0
    column = LineString([(x, y), (x, gy)])
    materials = (slope_data or {}).get("materials") or []
    total = 0.0
    for poly in (slope_data or {}).get("polygons") or []:
        shape = poly.get("polygon") if isinstance(poly, dict) else poly
        mat_id = poly.get("mat_id") if isinstance(poly, dict) else None
        if shape is None or mat_id is None:
            continue
        try:
            piece = column.intersection(shape)
        except Exception:                               # pragma: no cover
            continue
        if piece.is_empty:
            continue
        try:
            total += piece.length * float(materials[int(mat_id)]["gamma"])
        except (IndexError, KeyError, TypeError, ValueError):
            continue
    return total


def _pore_pressure(slope_data, material, x, y, sigma_v):
    """Pore pressure at a point, under the material's own ``u`` option.

    Mirrors ``generate_slices``: a piezometric head, a pore-pressure ratio on the
    soil column, or an interpolated seepage field. A model with no water answers
    zero, which is the same answer the solver gives it.
    """
    option = str((material or {}).get("u") or "none").strip().lower()
    gamma_w = float((slope_data or {}).get("gamma_water") or 0.0)
    if option == "piezo":
        line = (slope_data or {}).get("piezo_line") or []
        if len(line) < 2:
            return 0.0
        xs = [p[0] for p in line]
        if x < xs[0] or x > xs[-1]:
            return 0.0
        for k in range(len(xs) - 1):
            if xs[k] <= x <= xs[k + 1] and xs[k + 1] > xs[k]:
                t = (x - xs[k]) / (xs[k + 1] - xs[k])
                py = line[k][1] + t * (line[k + 1][1] - line[k][1])
                return max(0.0, (py - y) * gamma_w)
        return max(0.0, (line[-1][1] - y) * gamma_w)
    if option == "ru":
        try:
            return float(material.get("ru") or 0.0) * sigma_v
        except (TypeError, ValueError):
            return 0.0
    if option == "seep":
        mesh = (slope_data or {}).get("mesh")
        field = (slope_data or {}).get("seep_u")
        if mesh is None or field is None:
            return 0.0
        try:
            from .mesh import interpolate_at_point
            value, found = interpolate_at_point(
                mesh["nodes"], mesh["elements"], mesh["element_types"], field,
                (x, y), return_found=True, signed=True)
        except Exception:                               # pragma: no cover
            return 0.0
        return max(0.0, float(value)) if found else 0.0
    return 0.0


def _comparable_strength(slope_data, ground, material, x, y):
    """Mobilisable shear strength of a material at a point, across every model.

    This is the one quantity that makes zones rankable. Every strength option in
    the ``mat`` sheet reduces to a shear strength at the normal stress the zone
    carries, and the reductions are the solver's own: the undrained ``cp`` form is
    already a strength and needs no stress, a Hoek-Brown rock mass is linearised
    by :func:`xslope.hoekbrown.hb_tangent` at that stress, and a power envelope is
    evaluated on it. An ``elastic`` material cannot fail and is excluded upstream.
    """
    option = str((material or {}).get("option") or "mc").strip().lower()
    if option == "elastic":
        return None
    if option == "cp":
        try:
            return (float(material["c"])
                    + max(0.0, float(material["r_elev"]) - y) * float(material["cp"]))
        except (KeyError, TypeError, ValueError):
            return None
    sigma_v = _column_stress(slope_data, ground, x, y)
    sigma_n = max(0.0, sigma_v - _pore_pressure(slope_data, material, x, y, sigma_v))
    try:
        if option == "hb":
            from .hoekbrown import hb_tangent
            c_t, phi_t = hb_tangent(sigma_n, material["hb_sci"], material["hb_gsi"],
                                    material["hb_mi"], material["hb_d"])
            return float(c_t) + sigma_n * math.tan(math.radians(float(phi_t)))
        if option == "pow":
            a, b = float(material["pow_a"]), float(material["pow_b"])
            c_p, d_p = float(material["pow_c"]), float(material["pow_d"])
            s = max(sigma_n + d_p, 1e-4 * max(1.0, sigma_n))
            return a * s ** b + c_p
        return (float(material["c"])
                + sigma_n * math.tan(math.radians(float(material["phi"]))))
    except (KeyError, TypeError, ValueError):
        return None


def _friction_angle(slope_data, ground, material, x, y):
    """The friction angle a wedge ramp through this material should use."""
    option = str((material or {}).get("option") or "mc").strip().lower()
    if option in ("cp", "elastic", ""):
        return 0.0
    sigma_v = _column_stress(slope_data, ground, x, y)
    sigma_n = max(0.0, sigma_v - _pore_pressure(slope_data, material, x, y, sigma_v))
    try:
        if option == "hb":
            from .hoekbrown import hb_tangent
            _c, phi_t = hb_tangent(sigma_n, material["hb_sci"], material["hb_gsi"],
                                   material["hb_mi"], material["hb_d"])
            return float(phi_t)
        if option == "pow":
            a, b = float(material["pow_a"]), float(material["pow_b"])
            d_p = float(material["pow_d"])
            s = max(sigma_n + d_p, 1e-4 * max(1.0, sigma_n))
            return math.degrees(math.atan(a * b * s ** (b - 1.0)))
        return float(material["phi"])
    except (KeyError, TypeError, ValueError):
        return 0.0


def _zone_stations(polygon, ground, n, tol):
    """``(x, base, top)`` at ``n`` stations where the zone has real thickness.

    Stations where the polygon pinches to nothing are dropped rather than averaged
    in. A zone can be a degenerate sliver over part of its own bounding box -- the
    middle band of ``vp063`` is one, zero-thickness for the whole run left of the
    toe -- and sampling that sliver would rank the zone on ground it does not
    occupy and then try to build a track along it.
    """
    x0, _y0, x1, _y1 = polygon.bounds
    if x1 - x0 <= tol:
        return []
    out = []
    for i in range(n):
        x = x0 + (x1 - x0) * (i + 0.5) / n
        span = _vertical_span(polygon, x)
        if span is None:
            continue
        base, top = span
        if top - base <= tol:
            continue
        if _ground_elevation(ground, x) is None:
            continue
        out.append((x, base, top))
    return out


def rank_weak_zones(slope_data, samples=_ZONE_SAMPLES):
    """Every material zone, ranked weakest first by comparable shear strength.

    The ranking is over :class:`WeakZone` records, so a caller that has to *ask*
    which zone to seed on -- Studio's picker, or a script deciding for itself --
    shows the same numbers the automatic pick would have used.

    Strength is sampled just above the zone's own base, because that is where a
    surface tracking the zone would run. ``elastic`` materials are excluded: they
    cannot fail, so they cannot be the weakest anything.

    Returns
    -------
    tuple of WeakZone
        Weakest first. Ties break on polygon order, so the result is stable.
    """
    try:
        geom = slope_geometry(slope_data)
    except ValueError:
        return ()
    ground = list(geom.ground)
    tol = 1e-6 * max(geom.height, 1.0)
    materials = (slope_data or {}).get("materials") or []
    out = []
    for index, poly in enumerate((slope_data or {}).get("polygons") or []):
        shape = poly.get("polygon") if isinstance(poly, dict) else poly
        mat_id = poly.get("mat_id") if isinstance(poly, dict) else None
        if shape is None or mat_id is None or not hasattr(shape, "bounds"):
            continue
        try:
            material = materials[int(mat_id)]
        except (IndexError, TypeError, ValueError):
            continue
        if str(material.get("option") or "").strip().lower() == "elastic":
            continue
        stations = _zone_stations(shape, ground, samples, tol)
        taus = []
        for (x, base, top) in stations:
            tau = _comparable_strength(slope_data, ground, material,
                                       x, base + _OFFSET_FRAC * (top - base))
            if tau is not None:
                taus.append(tau)
        if not taus:
            continue
        xs = [s[0] for s in stations]
        out.append(WeakZone(index=index, mat_id=int(mat_id),
                            name=str(material.get("name") or f"zone {index + 1}"),
                            option=str(material.get("option") or "mc"),
                            tau=sum(taus) / len(taus), polygon=shape,
                            x_range=(min(xs), max(xs)),
                            y_range=(float(shape.bounds[1]), float(shape.bounds[3]))))
    out.sort(key=lambda z: (z.tau, z.index))
    return tuple(out)


def _ray_to_ground(ground, x0, y0, sx, tan_theta):
    """Where a ray rising at ``tan_theta`` in the ``sx`` direction meets the ground.

    Solved against the ground polyline segment by segment and the nearest crossing
    kept, which is what makes a ramp land on the first ground it reaches rather
    than on a bench further along.
    """
    best = None
    for (x1, y1), (x2, y2) in zip(ground, ground[1:]):
        dx, dy = x2 - x1, y2 - y1
        det = -sx * dy + dx * tan_theta
        if abs(det) < 1e-15:
            continue
        t = (-(x1 - x0) * dy + dx * (y1 - y0)) / det
        u = (sx * (y1 - y0) - tan_theta * (x1 - x0)) / det
        if t <= 1e-9 or not (-1e-9 <= u <= 1 + 1e-9):
            continue
        if best is None or t < best[0]:
            best = (t, (x0 + sx * t, y0 + tan_theta * t))
    return best[1] if best else None


#: How close the base of a zone has to come to the ground before the zone counts
#: as outcropping, as a fraction of the slope height. Below this the remaining
#: cover is thinner than any slice the surface would be cut into, so insisting on
#: a ramp there would build one inside a sliver.
_OUTCROP_COVER = 0.01

#: Stations used to hunt an outcrop between the track's end and the edge of the
#: zone. Fine, because the answer is an entry point on the ground surface and a
#: coarse scan would put it tens of feet from where the zone actually reaches day.
_OUTCROP_SCAN = 400


def _outcrop(zone, ground, x_from, x_to, cover):
    """Where the zone's base first comes within ``cover`` of the ground surface.

    A weak zone that daylights needs no ramp at that end -- the zone *is* the ramp.
    Griffiths & Lane's dipping band outcrops at both the crest and the toe, and its
    published critical surface is nothing but the band; a generator that insisted
    on building wedge ramps there would replace the mechanism with an approximation
    of it. Scanned outward from the track's end, so the answer is the *first*
    ground the zone reaches rather than the far side of the section.
    """
    if abs(x_to - x_from) <= 0:
        return None
    for i in range(1, _OUTCROP_SCAN + 1):
        x = x_from + (x_to - x_from) * i / _OUTCROP_SCAN
        span = _vertical_span(zone.polygon, x)
        gy = _ground_elevation(ground, x)
        if span is None or gy is None:
            continue
        if gy - span[0] <= cover:
            return x
    return None


def generate_noncircular_surface(slope_data, zone=None, separation=_SEPARATION,
                                 offset=_OFFSET_FRAC, points=_TRACK_POINTS,
                                 report=False):
    """A starting non-circular surface that tracks the model's weak zone.

    For a slope whose mechanism is controlled by a weak layer rather than by its
    own geometry, the useful starting surface is a polyline that runs just above
    the base of that layer along its continuous extent and reaches the ground
    surface at each end. A circular search cannot find this shape at all, and
    until now the only way to get one was to read it off a drawing by hand.

    How the zone is chosen
    ----------------------
    Every zone is ranked by :func:`rank_weak_zones` on the shear strength it can
    mobilise at the stress it actually carries, which is the only quantity
    comparable across an undrained clay, a frictional sand and a Hoek-Brown rock
    mass. Then:

    * **Clearly weakest** -- the weakest zone is at or below ``separation`` times
      the next weakest -- and the generator picks it and says so.
    * **No clear winner**, and it returns the ranked candidates instead of a guess.
      Two comparable candidate seams is a real situation, not a failure, and the
      answer to it is a question rather than a coin toss.
    * ``zone`` given explicitly, and that zone is used. Same call, same result
      shape, so a script overrides exactly the way the dialog does.

    The shape it builds
    -------------------
    * **The track** runs at ``offset`` of the zone's local thickness above its
      base, following the base through its own dip and kinks.
    * **Its ends** sit under the toe and under the crest, pulled in to where the
      zone actually exists -- so a seam that pinches out mid-slope ends the track
      at the pinch-out rather than at a place it does not reach.
    * **Where the zone daylights**, the track simply runs out to the ground and
      that is the entry or exit point. Where it does not, a straight ramp at the
      wedge angle (``45 - phi/2`` leaving at the toe, ``45 + phi/2`` entering at
      the crest) carries the surface up to the ground.
    * **No point is spare.** Where the track runs level the stations between its
      ends are dropped: an interior point is a ``Horiz`` one, its only move is a
      slide along x, and on a level run that slide leaves the surface exactly
      where it was. A dipping run keeps its stations, because there the same
      slide bends the surface and the point is a move the search can make.
    * **Every point carries an explicit Y and an explicit Movement.** A blank Y
      reaches the slicer as a ``TypeError`` and a blank Movement silently means
      ``Fixed``, which would freeze the search on the surface it was handed.

    Parameters
    ----------
    slope_data : dict
        A loaded model. Nothing is mutated.
    zone : int or str, optional
        The zone to seed on -- a polygon index or a material name. Omitted means
        pick automatically under the rule above.
    separation : float
        How much weaker the weakest zone must be than the next before it is picked
        without asking. Tuned against the corpus; see :data:`_SEPARATION`.
    offset : float
        Where the track runs inside the zone, as a fraction of local thickness
        above the base.
    points : int
        Interior track points.
    report : bool
        Return a dict carrying the surface, the chosen zone, the ranked
        candidates, a one-line summary and -- when nothing could be built -- the
        reason, instead of the bare surface. This is what a picker dialog needs,
        since it has to say *which* zone and *why* before anything is committed.

    Returns
    -------
    list of dict
        Points in the sheet's own form -- ``X``, ``Y``, ``Movement`` -- ordered
        left to right. Empty when no zone was clearly weakest and none was named,
        or when no surface could be built. With ``report=True``, a dict with
        ``surface``, ``zone``, ``candidates``, ``confident``, ``summary`` and
        ``reason``.

    Examples
    --------
    >>> out = generate_noncircular_surface(slope_data, report=True)   # doctest: +SKIP
    >>> if out["surface"]:                                            # doctest: +SKIP
    ...     slope_data["non_circ"] = out["surface"]                   # doctest: +SKIP
    ... else:                                                         # doctest: +SKIP
    ...     ask_the_user(out["candidates"])                           # doctest: +SKIP
    """
    def done(surface=(), picked=None, candidates=(), confident=False,
             summary="", reason=""):
        if not report:
            return list(surface)
        return {"surface": list(surface), "zone": picked,
                "candidates": list(candidates), "confident": bool(confident),
                "summary": summary, "reason": reason}

    try:
        geom = slope_geometry(slope_data)
    except ValueError as exc:
        return done(reason=str(exc))

    zones = rank_weak_zones(slope_data)
    if not zones:
        return done(reason="this model defines no material zone a failure surface "
                           "could be tracked along")

    picked, confident, why = None, False, ""
    if zone is not None:
        picked = _named_zone(zones, zone)
        if picked is None:
            return done(candidates=zones,
                        reason=f"no material zone matches {zone!r}")
        confident = True
        why = f"seeding on '{picked.name}' as asked"
    elif len(zones) == 1:
        return done(candidates=zones,
                    reason=("this model has one material zone, so it has no weak "
                            "layer for a surface to track. A non-circular surface "
                            "is for a slope whose mechanism follows a weak zone; "
                            "with one material the mechanism follows the slope's "
                            "own geometry, which is what a circular search finds"))
    else:
        weakest, runner_up = zones[0], zones[1]
        ratio = (weakest.tau / runner_up.tau) if runner_up.tau > 0 else 1.0
        if ratio <= separation:
            picked, confident = weakest, True
            why = (f"seeding on '{weakest.name}' -- mobilisable strength "
                   f"{weakest.tau:.3g} against {runner_up.tau:.3g} for the next "
                   f"weakest ('{runner_up.name}')")
        else:
            return done(candidates=zones, confident=False,
                        reason=(f"no zone is clearly the weakest: '{weakest.name}' "
                                f"at {weakest.tau:.3g} and '{runner_up.name}' at "
                                f"{runner_up.tau:.3g} are comparable. Choose which "
                                f"one the surface should track"))

    surface, note = _track_surface(slope_data, geom, picked, offset, points)
    if not surface:
        return done(picked=picked, candidates=zones, confident=confident,
                    reason=note)
    return done(surface=surface, picked=picked, candidates=zones,
                confident=confident, summary=f"{why}; {note}")


def _named_zone(zones, zone):
    """Resolve an explicit ``zone`` -- a polygon index or a material name."""
    if isinstance(zone, WeakZone):
        return zone
    if isinstance(zone, bool):
        return None
    if isinstance(zone, int):
        for z in zones:
            if z.index == zone:
                return z
        return None
    wanted = str(zone).strip().lower()
    matches = [z for z in zones if z.name.strip().lower() == wanted]
    return matches[0] if matches else None       # zones are already weakest-first


def _track_surface(slope_data, geom, zone, offset_frac, points):
    """Build the polyline for a chosen zone. Returns ``(surface, note)``."""
    ground = list(geom.ground)
    scale = max(geom.height, 1.0)
    tol = 1e-6 * scale
    face = _face_over(geom, zone)
    downhill = 1 if face.right_facing else -1

    cap = _OFFSET_CAP * scale

    def track_y(x):
        span = _vertical_span(zone.polygon, x)
        if span is None:
            return None
        base, top = span
        if top - base <= tol:
            return None
        return base + min(offset_frac * (top - base), cap)

    if abs(face.crest[0] - face.toe[0]) <= 0.1 * scale:
        return [], ("this section's slope face is vertical, so it has no horizontal "
                    "run for a track to follow beneath. A weak-zone surface is for a "
                    "slope, not a wall")

    window = _usable_window(zone, face, ground, track_y, tol)
    if window is None:
        return [], ("this zone never sits below the ground surface with enough cover "
                    "for a failure surface to run inside it")
    x_lo, x_hi = window

    # The track's own ends: under the toe and under the crest, pulled in to where
    # the zone actually exists. A zone that pinches out mid-slope ends the track at
    # the pinch-out; one that runs off the section is stopped at the slope instead,
    # because a surface reaching a vertical model edge is clipped, not critical.
    x_down = min(max(face.toe[0], x_lo), x_hi)
    x_up = min(max(face.crest[0], x_lo), x_hi)
    if abs(x_up - x_down) <= tol:
        return [], ("this zone is too narrow beneath the slope to lay a track along "
                    "-- it exists over less than the width of one slice there")

    ends, ramps = {}, []
    zx0, zx1 = float(zone.polygon.bounds[0]), float(zone.polygon.bounds[2])
    for name, anchor, sign in (("exit", x_down, downhill), ("entry", x_up, -downhill)):
        edge = min(zx1, ground[-1][0]) if sign > 0 else max(zx0, ground[0][0])
        y_anchor = track_y(anchor)
        if y_anchor is None:
            return [], "the zone has no thickness beneath the slope"
        cut = _outcrop(zone, ground, anchor, edge, _OUTCROP_COVER * scale)
        if cut is not None and abs(cut - anchor) > tol:
            gy = _ground_snap(slope_data, ground, cut)
            ends[name] = (cut, gy, cut)          # daylight point, track end == it
            continue
        theta = _wedge_angle(slope_data, ground, zone, anchor, y_anchor, name)
        if name == "exit":
            # The passive wedge emerges AT the toe. Cutting back from the toe into
            # the slope, rather than ramping outward from under it, is worth the
            # extra construction: a surface daylighting out on the flat has to shear
            # a block of ground that resists it and drives nothing, and on vp009 --
            # where that flat carries a 20 kPa surcharge -- the two constructions
            # are the difference between finding the mechanism and finding a local
            # minimum a third above it.
            back = _cut_back(face.toe, track_y, -sign, theta, x_lo, x_hi, tol)
            if back is not None:
                ends[name] = (face.toe[0], face.toe[1], back)
                ramps.append((name, theta))
                continue
        landing = _ramp_landing(slope_data, ground, anchor, y_anchor, sign, theta,
                                tol * scale)
        if landing is None:
            return [], (f"no ramp from this zone reaches the ground surface inside "
                        f"the model at the {'toe' if name == 'exit' else 'crest'} "
                        f"end -- the section has no room beyond the slope for the "
                        f"surface to daylight on the ground")
        ends[name] = (landing[0], landing[1], anchor)
        ramps.append((name, landing[2]))

    x_a, x_b = sorted((ends["exit"][2], ends["entry"][2]))
    n = max(2, points)
    stations = [x_a + (x_b - x_a) * i / (n - 1) for i in range(n)]
    # The zone's own vertices join the even stations. A seam that dips, flattens and
    # dips again has its bends AT those vertices, and a track sampled only on a
    # uniform grid rounds them off -- on Griffiths & Lane's three-segment band that
    # rounding is worth 3% of the factor of safety, because the search has to spend
    # its budget putting back a kink the geometry already knew about.
    for vx in _zone_vertices(zone.polygon):
        if x_a + tol < vx < x_b - tol:
            stations.append(vx)
    interior = []
    for x in sorted(stations):
        y = track_y(x)
        if y is not None and (not interior or x > interior[-1][0] + tol):
            interior.append((x, y))
    if len(interior) < 2:
        return [], "the zone has no continuous run beneath the slope to track along"

    left, right = sorted((ends["exit"], ends["entry"]), key=lambda e: e[0])
    pts = [(left[0], left[1])]
    for (x, y) in interior:
        if x > pts[-1][0] + tol and x < right[0] - tol:
            pts.append((x, y))
    pts.append((right[0], right[1]))
    pts = _drop_inert_points(pts, scale)
    if len(pts) < 3:
        return [], "the generated surface collapsed to a single segment"

    pts, nudge = _ensure_crossings(slope_data, pts, scale)
    if pts is None:
        return [], ("the generated surface does not cut the ground surface twice, "
                    "so the slicer would refuse it")

    problem = _surface_problem(pts, ground, geom, tol, nudge)
    if problem:
        return [], problem

    surface = [{"X": float(x), "Y": float(y),
                "Movement": "Free" if i in (0, len(pts) - 1) else "Horiz"}
               for i, (x, y) in enumerate(pts)]
    where = "; ".join(f"a {ang:.0f} degree ramp to the ground at the "
                      f"{'toe' if n == 'exit' else 'crest'}" for n, ang in ramps)
    note = (f"a {len(surface)}-point surface tracking {offset_frac:.0%} of the "
            f"zone's thickness above its base from x = {pts[0][0]:g} to "
            f"{pts[-1][0]:g}")
    note += f", with {where}" if ramps else ", daylighting where the zone itself does"
    return surface, note


def _drop_inert_points(pts, scale):
    """Drop interior points that neither shape the surface nor give the search
    anything to move.

    The even stations are laid across the whole track, so wherever the zone's base
    runs straight they land in the middle of a straight segment. A point there is
    worth keeping only if something can still use it, and there are only two
    things that could: the polyline, and the search that moves it. Removing a
    point that lies on the line through its neighbours leaves the polyline
    unchanged, so the whole case for it rests on the search -- which slides an
    interior point horizontally at a fixed elevation, because that is what a
    ``Horiz`` point is.

    On a **dipping** straight run that slide bends the surface, and the point
    earns its place: it is a degree of freedom the search can spend, which is why
    a dipping seam keeps its subdivision. On a **level** run the slide moves the
    point along the segment it already sits on, and the surface after the move is
    the surface before it. Such a point is inert -- it costs a row in the table
    and a slice boundary and buys nothing at all -- so the level run becomes what
    it always was, its two ends.

    Applied repeatedly, so a level run carrying three stations collapses to those
    two ends rather than to a shorter subdivision of itself.
    """
    tol = _COLLINEAR_TOL * scale
    out = list(pts)
    i = 1
    while i <= len(out) - 2:
        (xa, ya), (xb, yb), (xc, yc) = out[i - 1], out[i], out[i + 1]
        span = math.hypot(xc - xa, yc - ya)
        off = (abs((xc - xa) * (yb - ya) - (yc - ya) * (xb - xa)) / span
               if span > 0 else 0.0)
        if off <= tol and abs(yc - ya) <= tol:
            del out[i]                     # removal changes nothing; nor would a move
            continue
        i += 1
    return out


def _zone_vertices(polygon):
    """The x values of a zone's own boundary vertices, single- or multi-part.

    A zone read off the polygon sheet is one ring, but a zone that pinches out and
    reappears is two, and asking a MultiPolygon for its exterior raises rather than
    answering."""
    parts = getattr(polygon, "geoms", None)
    rings = list(parts) if parts is not None else [polygon]
    out = []
    for ring in rings:
        exterior = getattr(ring, "exterior", None)
        if exterior is None:
            continue
        out += [float(x) for x, _y in exterior.coords]
    return out


def _face_over(geom, zone):
    """The slope face this zone runs under -- the one a track should span."""
    best, best_overlap = geom.primary, -1.0
    for face in geom.faces:
        lo, hi = sorted((face.toe[0], face.crest[0]))
        overlap = min(hi, zone.x_range[1]) - max(lo, zone.x_range[0])
        if overlap > best_overlap:
            best, best_overlap = face, overlap
    return best


def _usable_window(zone, face, ground, track_y, tol):
    """The x range over which a failure surface could run inside this zone.

    Lateral continuity, resolved: a zone is used over its longest **continuous**
    run, not over its bounding box. A seam that pinches to nothing mid-section, or
    that surfaces and re-enters, is two zones as far as a sliding mass is
    concerned, and a track laid across the gap would leave the material it is
    supposed to be tracking. Where more than one run survives, the one that covers
    the most of the slope's own horizontal extent wins -- a fragment out on the
    flat ground beyond the toe carries no mechanism.
    """
    x0, x1 = zone.x_range
    runs, run = [], []
    for i in range(_ZONE_SAMPLES + 1):
        x = x0 + (x1 - x0) * i / _ZONE_SAMPLES
        y = track_y(x)
        gy = _ground_elevation(ground, x)
        if y is None or gy is None or gy < y - tol:
            if len(run) >= 2:
                runs.append((run[0], run[-1]))
            run = []
            continue
        run.append(x)
    if len(run) >= 2:
        runs.append((run[0], run[-1]))
    if not runs:
        return None
    lo, hi = sorted((face.toe[0], face.crest[0]))
    return max(runs, key=lambda r: (min(hi, r[1]) - max(lo, r[0]), r[1] - r[0]))


def _wedge_angle(slope_data, ground, zone, x, y, kind):
    """The Rankine wedge inclination for a ramp leaving the track at ``(x, y)``.

    ``45 + phi/2`` for the active wedge that drops in behind the crest and
    ``45 - phi/2`` for the passive one that pushes out at the toe -- the textbook
    choice, and one that costs nothing because ``phi`` has already been read. It is
    taken from the material the ramp actually cuts, which is the ground *above* the
    weak zone, not the zone itself: the ramp is a scarp through the overburden.

    Capped below the search's own ``max_base_angle``, because a starting surface
    the search rejects on sight is worse than none.
    """
    span = _vertical_span(zone.polygon, x)
    gy = _ground_elevation(ground, x)
    top = span[1] if span else y
    above = 0.5 * (top + gy) if (gy is not None and gy > top) else top
    material = (_material_at(slope_data, (x, above))
                or _material_at(slope_data, (x, y)))
    phi = _friction_angle(slope_data, ground, material, x, y)
    base = (45.0 + 0.5 * phi) if kind == "entry" else (45.0 - 0.5 * phi)
    return min(_RAMP_MAX_DEG, max(_RAMP_MIN_DEG, base))


def _cut_back(toe, track_y, sign, theta, x_lo, x_hi, tol):
    """Where a ramp daylighting at the toe meets the track, as an x value.

    Run from the ground down rather than from the track up, so the exit point is
    the toe itself. Solved by scanning the gap between the ramp and the track for a
    sign change and bisecting it -- the track is not a straight line, so there is
    no closed form.
    """
    tan_t = math.tan(math.radians(theta))

    def gap(t):
        y = track_y(toe[0] + sign * t)
        return None if y is None else (toe[1] - tan_t * t) - y

    reach = abs((x_hi - x_lo))
    if reach <= tol:
        return None
    lo_t = None
    for i in range(1, 201):
        t = reach * i / 200.0
        x = toe[0] + sign * t
        if not (x_lo - tol <= x <= x_hi + tol):
            break
        g = gap(t)
        if g is None:
            continue
        if g <= 0:
            if lo_t is None:
                return None                 # the track is already above the toe
            a, b = lo_t, t
            for _ in range(60):
                mid = 0.5 * (a + b)
                gm = gap(mid)
                if gm is None or gm > 0:
                    a = mid
                else:
                    b = mid
            return toe[0] + sign * 0.5 * (a + b)
        lo_t = t
    return None


def _ramp_landing(slope_data, ground, x_end, y_end, sign, theta, tol):
    """Where a wedge ramp from the track meets the ground. ``(x, y, degrees)``."""
    gx0, gx1 = ground[0][0], ground[-1][0]
    inset = 1e-3 * (gx1 - gx0)
    for factor in _RAMP_LADDER:
        angle = min(_RAMP_MAX_DEG, max(_RAMP_MIN_DEG, theta * factor))
        hit = _ray_to_ground(ground, x_end, y_end, sign,
                             math.tan(math.radians(angle)))
        if hit is None:
            continue
        if not (gx0 + inset < hit[0] < gx1 - inset):
            continue
        if abs(hit[0] - x_end) <= tol:
            continue
        # Snap the landing onto the ground polyline rather than keeping the ray's
        # own arithmetic. The two agree to within rounding, and rounding is exactly
        # what matters: an end point a hair BELOW the ground never crosses it, and
        # the slicer refuses the surface with "expected at least 2 intersection
        # points, but got 1" -- a failure whose cause is invisible in the geometry.
        landed = _ground_snap(slope_data, ground, hit[0])
        return (hit[0], hit[1] if landed is None else landed, angle)
    return None


#: Nudges tried, as fractions of the section's width, when the slicer will not see
#: an end point that is sitting on the ground surface.
_NUDGE_LADDER = (0.0, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5)


def _crossings_seen(surface, pts):
    """How many ground crossings the *slicer* will find on this polyline."""
    from shapely.geometry import LineString
    try:
        hit = surface.intersection(LineString(pts))
    except Exception:                                    # pragma: no cover
        return 0
    geoms = getattr(hit, "geoms", None)
    if geoms is not None:
        return sum(1 for g in geoms if g.geom_type == "Point")
    return 1 if hit.geom_type == "Point" else 0


def _extend(inner, end, distance):
    """``end``, moved away from ``inner`` by ``distance`` along their own line."""
    dx, dy = end[0] - inner[0], end[1] - inner[1]
    length = math.hypot(dx, dy)
    if length <= 0:
        return end
    return (end[0] + dx / length * distance, end[1] + dy / length * distance)


def _ensure_crossings(slope_data, pts, scale):
    """Return the polyline in a form ``generate_slices`` will accept, and by how
    much its ends had to be pushed to get there.

    An end point that lies exactly *on* the ground surface is a tangency, and
    whether GEOS reports a tangency as an intersection depends on the last bits of
    two independently-computed coordinates. When it does not, the slicer refuses
    the surface with "expected at least 2 intersection points, but got 1" -- and
    the geometry that provoked it is visibly, measurably correct, which makes the
    failure very hard to read. So the generator asks the question itself, and where
    the answer is no, extends the end segments outward by the smallest nudge that
    turns the tangency into an ordinary crossing. The largest nudge on offer is a
    hundred-thousandth of the section, orders below anything a plot or a slice
    width can resolve, and the search's first accepted move of that end point
    snaps it back onto the ground exactly.
    """
    surface = (slope_data or {}).get("ground_surface")
    if surface is None or getattr(surface, "is_empty", True):
        return list(pts), 0.0
    for frac in _NUDGE_LADDER:
        nudge = frac * scale
        out = list(pts)
        if nudge > 0:
            out[0] = _extend(pts[1], pts[0], nudge)
            out[-1] = _extend(pts[-2], pts[-1], nudge)
        if _crossings_seen(surface, out) >= 2:
            return out, nudge
    return None, 0.0


def _surface_problem(pts, ground, geom, tol, nudge=0.0):
    """Why a candidate polyline is not usable, or ``""`` when it is."""
    if any(b[0] <= a[0] + tol for a, b in zip(pts, pts[1:])):
        return "the generated points are not strictly ordered left to right"
    gx0, gx1 = ground[0][0], ground[-1][0]
    if pts[0][0] <= gx0 + tol or pts[-1][0] >= gx1 - tol:
        return ("the surface reaches a vertical edge of the section rather than "
                "daylighting on the ground surface. The section needs flat ground "
                "beyond the toe and the crest for a failure surface to exit on it")
    for i, (x, y) in enumerate(pts):
        gy = _ground_elevation(ground, x)
        if gy is None:
            return "a generated point falls outside the section"
        if i in (0, len(pts) - 1):
            if abs(gy - y) > max(1e-6 * max(1.0, abs(gy)), 2.0 * nudge):
                return "an end point does not sit on the ground surface"
        elif y >= gy - tol:
            return "a track point is not below the ground surface"
        if y < geom.floor - tol:
            return "the surface runs below the bottom of the model"
    return ""


def _stress_scale(slope_data, geom):
    """``gamma * H`` for the model -- the stress a cohesion is compared against."""
    gammas = []
    for m in (slope_data or {}).get("materials") or []:
        try:
            g = float(m.get("gamma"))
        except (TypeError, ValueError):
            continue
        if g > 0:
            gammas.append(g)
    gamma = max(gammas) if gammas else 1.0
    return max(1e-9, gamma * max(geom.height, 1e-9))
