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

* **Centre.** ``Xo`` halfway between toe and crest, ``Yo`` at the toe elevation
  plus twice the slope height.
* **A toe circle** -- one passing *through* the toe, ``R = dist(centre, toe)``.
  This is deliberately not a circle whose bottom sits at the toe *elevation*; the
  two are different surfaces and the distinction is easy to lose.
* **A circle tangent to the base of each distinct material layer**, with ``Depth``
  set to that layer's base elevation.
* **A skimming circle on a cohesionless face.** Where the material exposed on the
  face has ``c = 0``, the critical surface is an arbitrarily shallow face-parallel
  slide, and a toe-and-base seeded search converges to a deep local minimum with a
  non-conservatively high answer. A large-radius circle approximates the plane;
  its centre lands far outside the model, which is expected and correct.
* **Daylighting.** Every generated circle must cut the ground surface twice,
  inside the model rather than at a vertical edge, and must not dip below the
  bottom of the domain -- so the set can be handed straight to a search.

``Depth`` is an **elevation**, not a depth below ground: it is the elevation of
the circle's lowest point, and the radius is ``Yo - Depth``. That is the circles
sheet's own convention and a persistent source of confusion.

The generated circles are a *starting set*, not an answer. They exist to give the
search a seed in every family that could win.
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

#: Centre elevations to try, as multiples of the slope height above the toe. The
#: first is the rule; the rest are a fallback for a section transcribed too narrow
#: for it, where a standard circle would daylight past the model's own edge. A
#: lowered centre is reported, because a section that needs one is cropped -- the
#: real repair is to widen the geometry.
_CENTRE_LADDER = (2.0, 1.5, 1.0, 0.75)

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

    Only crossings **below the circle's centre** count, which is the slicer's own
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

    The centre sits on the outward side of the face, so the arc sags into the
    slope. ``k = R / L``; the centre lands far outside the model, which is what
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
        for mult in _CENTRE_LADDER:
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
            if mult != _CENTRE_LADDER[0]:
                note += (f", with the centre at {mult:g} x the slope height rather than "
                         f"the usual {_CENTRE_LADDER[0]:g} -- the section is not wide "
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
