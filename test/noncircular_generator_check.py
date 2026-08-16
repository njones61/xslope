"""Standing checks on the non-circular starting-surface generator -- the weak-zone
ranking, the threshold that decides whether it picks or asks, the geometry of the
surface it builds, and the Studio button and dialog that carry it.

A generated starting surface is not checked by looking at it. It is checked by the
properties that make it usable, and each of these is a property something could
break in silence:

  A. THE METRIC -- zones are ranked on the shear strength each can mobilize at the
     stress it actually carries, which is the only quantity comparable across a
     frictional sand, an undrained clay and a Hoek-Brown rock mass. Ranking on
     cohesion alone or friction alone picks a different zone, and the check asserts
     it picks a different one: a metric that agreed with `c` would not be doing any
     work. An `elastic` zone cannot fail and is excluded outright.
  B. THE THRESHOLD, FROM BOTH SIDES -- a separation just wide enough auto-picks and
     says which zone and why; a separation one step too narrow returns the ranked
     candidates instead. Asserted against a measured ratio rather than a literal,
     so the check tests the comparison rather than restating it. The SHIPPED value
     is pinned separately, inside the bracket the corpus measured it into.
  C. THE GEOMETRY -- the track lies inside the weak zone and above its base; the
     end points sit on the ground surface carrying an explicit float Y (never None
     and never NaN -- a blank Y reaches the slicer as a TypeError); every point
     carries an explicit Movement (a blank one silently means Fixed, which would
     freeze the search on the surface it was handed); and x increases strictly.
     Checked on a left-facing slope and on its mirror image, because a generator
     that only works downhill in one direction is a generator that works on half
     the corpus. No point on it is inert: an interior point that lies on the line
     through its neighbours on a LEVEL run neither shapes the surface nor gives
     the search anything to move, since the only move it has is a slide along the
     segment it already sits on. Checked from both sides -- a level run collapses
     to its two ends, a dipping run keeps its stations, because there the same
     slide bends the surface and the point is a degree of freedom worth having.
  D. THE SLICER -- the surface a corpus weak-seam model generates is one
     `generate_slices` accepts, cut into slices whose base runs inside the seam.
     This is the one check that would have caught the tangency failure, where a
     visibly correct polyline was refused because an end point sat exactly ON the
     ground surface.
  E. STUDIO -- the Generate button reaches the generator and its rows land in the
     table; it is DIMMED with the reason in its tooltip on a model that has no weak
     layer; and an ambiguous model raises the zone picker, whose choice the
     generator honours.

Skips cleanly (exit 0) when PySide6 is not installed; checks A-D run either way.
"""
import math
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from shapely.geometry import LineString, Polygon

from xslope import generators as G

#: A corpus model whose mechanism is controlled by a weak seam: ACADS problem 3(b),
#: a 0.5-unit weak layer under a 12-unit slope.
ACADS = os.path.join(_REPO, "docs/lem/files/xslope_acads_weak_layer.xlsx")
#: A corpus model with two comparable candidate zones -- Loukidis' three dipping
#: bands, where the ranking prefers the upper one and the mechanism rides the middle.
AMBIGUOUS = os.path.join(_REPO, "docs/verification/files/rocscience/vp063.xlsx")
#: A corpus model with a single material zone: nothing to track, button dimmed.
ONE_ZONE = os.path.join(_REPO, "docs/verification/files/rocscience/vp047.xlsx")
#: The tutorial model LEM-5 teaches: a flat soft-clay seam under a sand fill, whose
#: generated surface the page quotes point for point.
SAMPLE7 = os.path.join(_REPO, "docs/lem/files/xslope_noncircular.xlsx")

#: The bracket the corpus measured the separation threshold into: every ranking at
#: or under 0.56 picks a zone whose seeded search reaches or beats the row's own
#: answer, and the first ranking that picks the WRONG zone is 0.625 (Guo &
#: Griffiths' embankment over foundation). See generators._SEPARATION.
_BRACKET = (0.56, 0.625)


# ---------------------------------------------------------------------------
# A synthetic three-layer slope, and its mirror image
# ---------------------------------------------------------------------------

def _model(right_facing=True, seam=None, extra=None):
    """A slope over a thin weak seam over a stiff base.

    Ground runs from elevation 10 down to 0 across a 20-unit face. The seam is
    0.5 thick at y = -2.0 to -1.5, so a track inside it is unambiguous: any y
    strictly between those two is in the seam and nothing else is.
    """
    ground = [(0.0, 10.0), (20.0, 10.0), (40.0, 0.0), (60.0, 0.0)]
    bands = [
        # (polygon vertices, material)
        ([(0, -1.5), (0, 10), (20, 10), (40, 0), (60, 0), (60, -1.5)],
         dict(name="Stiff fill", option="mc", gamma=20.0, c=50.0, phi=30.0,
              cp=0.0, r_elev=0.0, u="none", ru=0.0)),
        ([(0, -2.0), (0, -1.5), (60, -1.5), (60, -2.0)],
         seam or dict(name="Weak seam", option="mc", gamma=20.0, c=5.0, phi=10.0,
                      cp=0.0, r_elev=0.0, u="none", ru=0.0)),
        # c = 0 and phi = 45: the LOWEST cohesion in the model and by a wide margin
        # the strongest zone at the stress it carries. A ranking on cohesion picks
        # this; the ranking that matters picks the seam.
        ([(0, -10.0), (0, -2.0), (60, -2.0), (60, -10.0)],
         dict(name="Dense sand", option="mc", gamma=20.0, c=0.0, phi=45.0,
              cp=0.0, r_elev=0.0, u="none", ru=0.0)),
    ]
    if extra:
        bands.append(extra)
    flip = (lambda x: x) if right_facing else (lambda x: 60.0 - x)
    def fix(pts):
        out = [(flip(x), y) for x, y in pts]
        return out if right_facing else out
    polygons, materials = [], []
    for i, (pts, mat) in enumerate(bands):
        polygons.append({"polygon": Polygon(fix(pts)), "mat_id": i})
        materials.append(mat)
    gcoords = sorted([(flip(x), y) for x, y in ground])
    domain = Polygon(fix([(0, -10), (0, 10), (20, 10), (40, 0), (60, 0), (60, -10)]))
    return {"ground_surface": LineString(gcoords), "polygons": polygons,
            "materials": materials, "domain_polygon": domain,
            "gamma_water": 9.81, "piezo_line": [], "non_circ": [], "circles": []}


def _outcrop_model():
    """The same slope, but with a seam that reaches the ground at both ends.

    A lens pinching to nothing where it meets the crest flat and again at the toe.
    Nothing here needs a ramp: the zone IS the ramp, which is the shape Griffiths &
    Lane's dipping band has and the reason their published critical surface is
    nothing but the band.
    """
    sd = _model()
    lens = Polygon([(10.0, 10.0), (25.0, 2.0), (40.0, 0.0), (25.0, 1.4)])
    sd["polygons"][1] = {"polygon": lens, "mat_id": 1}
    return sd


def _dipping_model():
    """The same slope, but the seam dips across the section instead of lying flat.

    Its base is still one straight line, so a track along it is still collinear --
    but not level, and that is the whole difference: sliding a ``Horiz`` point
    along x at fixed elevation bends this surface, so every station on it is a
    degree of freedom the search can spend.
    """
    sd = _model()
    fill = Polygon([(0, -1.5), (0, 10), (20, 10), (40, 0), (60, 0), (60, -7.5)])
    seam = Polygon([(0, -2.0), (0, -1.5), (60, -7.5), (60, -8.0)])
    sand = Polygon([(0, -10.0), (0, -2.0), (60, -8.0), (60, -10.0)])
    for i, poly in enumerate((fill, seam, sand)):
        sd["polygons"][i] = {"polygon": poly, "mat_id": i}
    return sd


def _pinchout_model():
    """The same slope, but the seam pinches out before the toe and reappears on the
    flat beyond it.

    Two lobes with a gap between them. As far as a sliding mass is concerned that
    is two zones, and a track laid straight across the gap would leave the material
    it is supposed to be tracking. The lobe out on the flat carries no mechanism at
    all, so the one under the slope is the one to use.
    """
    from shapely.geometry import MultiPolygon
    sd = _model()
    under_slope = Polygon([(5.0, -2.0), (5.0, -1.5), (38.0, -1.5), (38.0, -2.0)])
    beyond_toe = Polygon([(45.0, -2.0), (45.0, -1.5), (58.0, -1.5), (58.0, -2.0)])
    sd["polygons"][1] = {"polygon": MultiPolygon([under_slope, beyond_toe]),
                         "mat_id": 1}
    return sd


def _fail(cond, message):
    return [] if cond else [message]


#: How far off the chord through its neighbours a point may sit and still count as
#: on it, as a fraction of the slope height. Loose against the generator's own
#: `_COLLINEAR_TOL`, deliberately: the check has to catch a point the generator
#: left behind, not agree with it about the last bit.
_ON_THE_LINE = 1e-4


def _inert_points(surface):
    """The indices of interior points that do nothing.

    A point is inert when removing it would leave the polyline unchanged AND the
    move it is allowed cannot change the polyline either. Interior points are
    ``Horiz``, so their move is a slide along x at fixed elevation: on a level run
    that lands them back on the segment they started on. Both conditions are
    measured here rather than assumed, so a generator that starts emitting a
    different Movement is not silently waved through.
    """
    pts = [(p["X"], p["Y"]) for p in surface]
    out = []
    for i in range(1, len(pts) - 1):
        if surface[i]["Movement"] != "Horiz":
            continue
        (xa, ya), (xb, yb), (xc, yc) = pts[i - 1], pts[i], pts[i + 1]
        scale = max(abs(xc - xa), abs(yc - ya), 1e-12)
        span = math.hypot(xc - xa, yc - ya)
        off = (abs((xc - xa) * (yb - ya) - (yc - ya) * (xb - xa)) / span
               if span > 0 else 0.0)
        if off <= _ON_THE_LINE * scale and abs(yc - ya) <= _ON_THE_LINE * scale:
            out.append(i)
    return out


# ---------------------------------------------------------------------------
# A. The metric
# ---------------------------------------------------------------------------

def test_metric_ranks_on_mobilizable_strength():
    """The weakest zone by mobilizable strength, not by c and not by phi.

    The synthetic model is built so the three ways of ranking disagree: the seam
    has the LOWEST cohesion and the LOWEST friction angle, but the stiff base has
    the HIGHEST of both while carrying the most stress. The check asserts the
    ranking picks the seam AND that a cohesion-only ranking would have picked a
    different zone -- otherwise the metric could be a no-op.
    """
    out = []
    sd = _model()
    zones = G.rank_weak_zones(sd)
    out += _fail(len(zones) == 3, f"expected 3 rankable zones, got {len(zones)}")
    if not zones:
        return out
    out += _fail(zones[0].name == "Weak seam",
                 f"weakest zone is {zones[0].name!r}, expected 'Weak seam'")
    out += _fail(all(zones[i].tau <= zones[i + 1].tau for i in range(len(zones) - 1)),
                 "rank_weak_zones did not return the zones weakest first")
    # The model is built so a cohesion-only ranking disagrees: 'Dense sand' has the
    # lowest c in the model and is the strongest zone in it. If the ranking ever
    # agrees with cohesion here, it has stopped reading the stress.
    by_c = min(sd["materials"], key=lambda m: float(m["c"]))
    out += _fail(by_c["name"] == "Dense sand" and zones[0].name != by_c["name"],
                 "the ranking agrees with a ranking on cohesion alone, so it is no "
                 "longer evaluating strength at the stress each zone carries")

    # A frictional zone with no cohesion at low stress must still rank BELOW an
    # undrained zone with real cohesion -- which is the whole point of evaluating
    # at the stress the zone carries. Swap the seam for a cp material strong
    # enough to beat the fill and check the ranking moves.
    strong_seam = dict(name="Weak seam", option="cp", gamma=20.0, c=400.0,
                       cp=0.0, r_elev=0.0, phi=0.0, u="none", ru=0.0)
    zones2 = G.rank_weak_zones(_model(seam=strong_seam))
    out += _fail(zones2 and zones2[0].name != "Weak seam",
                 "a seam given 400 kPa of undrained strength still ranked weakest, "
                 "so the metric is not reading the strength model")

    # An elastic zone cannot fail, so it is not a candidate at any strength.
    elastic = ([(0, -12.0), (0, -10.0), (60, -10.0), (60, -12.0)],
               dict(name="Bedrock", option="elastic", gamma=20.0, c=0.0, phi=0.0,
                    cp=0.0, r_elev=0.0, u="none", ru=0.0))
    zones3 = G.rank_weak_zones(_model(extra=elastic))
    out += _fail(all(z.name != "Bedrock" for z in zones3),
                 "an option='elastic' zone was offered as a candidate weak zone")
    return out


def test_metric_spans_every_strength_model():
    """Each strength option reduces to a finite, comparable shear strength."""
    out = []
    sd = _model()
    ground = list(sd["ground_surface"].coords)
    cases = {
        "mc": dict(option="mc", c=10.0, phi=25.0),
        "cp": dict(option="cp", c=30.0, cp=1.5, r_elev=0.0),
        "pow": dict(option="pow", pow_a=2.0, pow_b=0.8, pow_c=1.0, pow_d=5.0),
        "hb": dict(option="hb", hb_sci=50000.0, hb_gsi=45.0, hb_mi=10.0, hb_d=0.0),
    }
    for name, mat in cases.items():
        tau = G._comparable_strength(sd, ground, mat, 30.0, -1.75)
        out += _fail(tau is not None and math.isfinite(tau) and tau > 0,
                     f"option={name!r} produced no comparable strength (got {tau!r})")
    out += _fail(G._comparable_strength(sd, ground, {"option": "elastic"},
                                        30.0, -1.75) is None,
                 "an elastic material returned a finite strength")
    return out


# ---------------------------------------------------------------------------
# B. The threshold, from both sides
# ---------------------------------------------------------------------------

def test_threshold_both_sides():
    """One step wider than the ranking and it picks; one step narrower and it asks."""
    out = []
    sd = _model()
    zones = G.rank_weak_zones(sd)
    if len(zones) < 2:
        return ["fewer than two zones to separate"]
    ratio = zones[0].tau / zones[1].tau

    wide = G.generate_noncircular_surface(sd, separation=ratio + 1e-6, report=True)
    out += _fail(bool(wide["surface"]),
                 f"a separation just above the ranking ({ratio:.4f}) did not pick")
    out += _fail(wide["confident"] is True,
                 "the automatic pick did not report itself as confident")
    out += _fail(wide["zone"] is not None and wide["zone"].name == "Weak seam",
                 "the automatic pick chose a zone other than the weak seam")
    out += _fail("Weak seam" in wide["summary"] and "Stiff" in wide["summary"],
                 f"the pick did not state the zone and its runner-up: "
                 f"{wide['summary'][:120]!r}")

    narrow = G.generate_noncircular_surface(sd, separation=ratio - 1e-6, report=True)
    out += _fail(not narrow["surface"],
                 f"a separation just below the ranking ({ratio:.4f}) still picked")
    out += _fail(narrow["confident"] is False,
                 "an ambiguous outcome reported itself as confident")
    out += _fail(len(narrow["candidates"]) == len(zones),
                 f"the ambiguous outcome offered {len(narrow['candidates'])} "
                 f"candidates, expected all {len(zones)}")
    out += _fail(bool(narrow["reason"]),
                 "the ambiguous outcome carried no reason to show the user")

    # An explicitly named zone bypasses the threshold entirely, by index and by name.
    for handle in (zones[1].index, zones[1].name):
        forced = G.generate_noncircular_surface(sd, zone=handle,
                                                separation=ratio - 1e-6, report=True)
        out += _fail(bool(forced["surface"]) and forced["zone"].index == zones[1].index,
                     f"an explicit zone={handle!r} was not honoured under an "
                     f"ambiguous separation")
    return out


def test_shipped_threshold_is_the_tuned_one():
    """The shipped separation sits inside the bracket the corpus measured."""
    lo, hi = _BRACKET
    return _fail(lo < G._SEPARATION < hi,
                 f"the shipped separation {G._SEPARATION} is outside the corpus "
                 f"bracket ({lo}, {hi}); re-tune it against the weak-seam rows "
                 f"rather than moving it by hand")


# ---------------------------------------------------------------------------
# C. Geometry invariants, both facings
# ---------------------------------------------------------------------------

def _geometry_failures(sd, label):
    """Every invariant, plus the toe-side anchoring that applies to a buried seam."""
    result = G.generate_noncircular_surface(sd, report=True)
    if not result["surface"]:
        return [f"{label}: no surface was generated ({result['reason'][:90]})"]
    out = _geometry_failures_common(sd, result, label)
    geom = G.slope_geometry(sd)
    # This seam does not reach the ground, so the passive end must daylight AT the
    # toe: a surface emerging out on the flat beyond it has to shear a block that
    # resists and drives nothing.
    downhill_end = result["surface"][-1] if geom.right_facing else result["surface"][0]
    out += _fail(abs(downhill_end["X"] - geom.toe[0]) < 1e-6 * max(1.0, geom.height),
                 f"{label}: the toe-side end is at x={downhill_end['X']:.6g}, not at "
                 f"the toe (x={geom.toe[0]:g})")
    return out


def _geometry_failures_common(sd, result, label):
    """The invariants that hold for every generated surface."""
    out = []
    surface = result["surface"]
    zone = result["zone"]

    xs = [p["X"] for p in surface]
    out += _fail(all(b > a for a, b in zip(xs, xs[1:])),
                 f"{label}: x is not strictly increasing: {xs}")

    for i, p in enumerate(surface):
        y = p["Y"]
        out += _fail(isinstance(y, float) and not math.isnan(y),
                     f"{label}: point {i} has Y={y!r} -- every point needs an "
                     f"explicit, finite Y")
        out += _fail(p["Movement"] in ("Free", "Horiz", "Fixed"),
                     f"{label}: point {i} has Movement={p['Movement']!r}; a blank "
                     f"or mistyped Movement silently means Fixed")
    out += _fail(surface[0]["Movement"] == "Free" and surface[-1]["Movement"] == "Free",
                 f"{label}: the end points are not Free, so the search cannot slide "
                 f"them along the ground surface")
    out += _fail(all(p["Movement"] == "Horiz" for p in surface[1:-1]),
                 f"{label}: an interior track point is not Horiz")

    for i in _inert_points(surface):
        out.append(f"{label}: point {i} at x = {surface[i]['X']:.6g} is inert -- it "
                   f"sits on the line through its neighbours on a level run, so "
                   f"neither removing it nor the only move it has (a horizontal "
                   f"slide) changes the surface at all")

    ground = list(sd["ground_surface"].coords)
    scale = G.slope_geometry(sd).height
    for i in (0, len(surface) - 1):
        gy = G._ground_elevation(ground, surface[i]["X"])
        out += _fail(gy is not None and abs(gy - surface[i]["Y"]) <= 1e-4 * scale,
                     f"{label}: end point {i} at y={surface[i]['Y']:.6g} is not on "
                     f"the ground surface (y={gy})")

    for p in surface[1:-1]:
        span = G._vertical_span(zone.polygon, p["X"])
        out += _fail(span is not None,
                     f"{label}: track point at x={p['X']:.6g} is outside the zone")
        if span is None:
            continue
        base, top = span
        out += _fail(p["Y"] > base + 1e-12,
                     f"{label}: track point at x={p['X']:.6g} sits at or below the "
                     f"zone base ({p['Y']:.6g} vs {base:.6g}) -- a surface on a "
                     f"material boundary is a slicing hazard")
        out += _fail(p["Y"] < top,
                     f"{label}: track point at x={p['X']:.6g} is above the zone top "
                     f"({p['Y']:.6g} vs {top:.6g}), so it is not in the zone at all")
        gy = G._ground_elevation(ground, p["X"])
        out += _fail(gy is not None and p["Y"] < gy,
                     f"{label}: track point at x={p['X']:.6g} is not below the ground")
    return out


def test_geometry_both_facings():
    """The invariants hold on a right-facing slope and on its mirror image."""
    out = []
    for right in (True, False):
        label = "right-facing" if right else "left-facing"
        sd = _model(right_facing=right)
        out += _geometry_failures(sd, label)
    # The two are mirror images, so they must produce mirror-image surfaces.
    a = G.generate_noncircular_surface(_model(True))
    b = G.generate_noncircular_surface(_model(False))
    if a and b:
        mirrored = sorted(60.0 - p["X"] for p in b)
        out += _fail(len(a) == len(b) and all(
            abs(p["X"] - m) < 1e-6 for p, m in zip(a, mirrored)),
            "the left-facing surface is not the mirror of the right-facing one")
    return out


def test_outcropping_zone_needs_no_ramp():
    """Where the zone reaches the ground, the track runs out to it -- no ramp.

    The ends must land where the seam meets the ground, not where a wedge ramp
    from under the slope would have put them: a ramp built through the overburden
    at an outcrop replaces the mechanism with an approximation of it.
    """
    sd = _outcrop_model()
    result = G.generate_noncircular_surface(sd, report=True)
    if not result["surface"]:
        return [f"outcropping seam: nothing generated ({result['reason'][:100]})"]
    out = _geometry_failures_common(sd, result, "outcropping seam")
    out += _fail("daylighting where the zone itself does" in result["summary"],
                 f"the summary claims a ramp on a zone that outcrops: "
                 f"{result['summary'][-90:]!r}")
    out += _fail("ramp" not in result["summary"],
                 "a ramp was built at an end where the zone reaches the ground")
    surface = result["surface"]
    zx0, zx1 = result["zone"].polygon.bounds[0], result["zone"].polygon.bounds[2]
    out += _fail(abs(surface[0]["X"] - zx0) < 0.5 and abs(surface[-1]["X"] - zx1) < 0.5,
                 f"the ends are at x={surface[0]['X']:.4g} and "
                 f"{surface[-1]['X']:.4g}, not at the seam's own outcrops "
                 f"({zx0:g}, {zx1:g})")
    # The lens bends at x = 25. The track has to bend there too: an evenly sampled
    # polyline rounds the corner off, and the search then has to spend its budget
    # putting back a kink the geometry already knew about.
    out += _fail(any(abs(p["X"] - 25.0) < 1e-6 for p in surface),
                 f"no track point sits on the seam's own kink at x = 25: "
                 f"{[round(p['X'], 3) for p in surface]}")
    return out


def test_pinched_out_zone_uses_its_continuous_run():
    """A zone in two pieces is tracked over the piece that lies under the slope."""
    sd = _pinchout_model()
    result = G.generate_noncircular_surface(sd, report=True)
    if not result["surface"]:
        return [f"pinched-out seam: nothing generated ({result['reason'][:100]})"]
    out = _geometry_failures_common(sd, result, "pinched-out seam")
    xs = [p["X"] for p in result["surface"]]
    out += _fail(max(xs) < 42.0,
                 f"the surface reaches x={max(xs):.4g}, crossing the gap into the "
                 f"lobe beyond the toe that carries no mechanism")
    out += _fail(min(xs) >= 5.0 - 1e-9,
                 f"the surface starts at x={min(xs):.4g}, outside the zone")
    return out


def test_wedge_ramps_are_active_and_passive():
    """The back scarp is steeper than the toe ramp, which is what makes them a
    Rankine pair -- 45 + phi/2 dropping in behind the crest, 45 - phi/2 pushing out
    at the toe. Swapped, the surface has a shallow scarp and a steep toe, which is
    neither construction."""
    out = []
    for right in (True, False):
        sd = _model(right_facing=right)
        surface = G.generate_noncircular_surface(sd)
        if len(surface) < 3:
            out.append(f"{'right' if right else 'left'}-facing: too few points")
            continue
        geom = G.slope_geometry(sd)
        pts = [(p["X"], p["Y"]) for p in surface]
        first, last = (pts[0], pts[1]), (pts[-2], pts[-1])
        uphill, downhill = (first, last) if geom.right_facing else (last, first)
        def incline(seg):
            (xa, ya), (xb, yb) = seg
            return abs(yb - ya) / max(abs(xb - xa), 1e-12)
        # phi is the OVERBURDEN's, not the seam's: the ramp is a scarp cut through
        # the ground above the weak layer, not through the layer. Reading it off the
        # seam here would give 50 and 40 degrees instead of 60 and 30.
        phi_over = 30.0                        # 'Stiff fill' in the synthetic model
        for label, seg, want in (("crest-side", uphill, 45.0 + 0.5 * phi_over),
                                 ("toe-side", downhill, 45.0 - 0.5 * phi_over)):
            got = math.degrees(math.atan(incline(seg)))
            out += _fail(abs(got - want) < 3.0,
                         f"{'right' if right else 'left'}-facing: the {label} ramp is "
                         f"{got:.1f} deg, not the {want:.0f} deg Rankine wedge angle "
                         f"for the phi = {phi_over:g} overburden it cuts")
    return out


def test_level_run_collapses_but_dipping_run_does_not():
    """A straight run keeps exactly the points something can use.

    The two halves are the same geometry differing only in dip, which is what
    makes this a test of the rule rather than of one model:

    * **Level.** The flat seam's track is one horizontal segment, and the stations
      the generator lays across it are points on the middle of a straight line
      whose only move slides them along that same line. It must emit the run's two
      ends and nothing between them.
    * **Dipping.** The same seam tilted. Its track is still one straight segment
      -- asserted here, so the check cannot pass by the run having quietly
      acquired a bend -- but now a horizontal slide lifts a station off the line,
      so the stations are degrees of freedom and must survive.
    """
    out = []
    flat = G.generate_noncircular_surface(_model())
    if len(flat) < 3:
        return ["the flat seam generated no track at all"]
    interior = flat[1:-1]
    out += _fail(len(interior) == 2,
                 f"the flat seam's track carries {len(interior)} interior points; "
                 f"a level run has two ends and nothing worth putting between them: "
                 f"{[round(p['X'], 3) for p in flat]}")
    out += _fail(len({round(p["Y"], 9) for p in interior}) == 1,
                 "the flat seam's interior points are not at one elevation, so this "
                 "is no longer the level run the check is about")

    dipping = G.generate_noncircular_surface(_dipping_model())
    if len(dipping) < 3:
        return out + ["the dipping seam generated no track at all"]
    track = [(p["X"], p["Y"]) for p in dipping[1:-1]]
    out += _fail(len(track) >= 3,
                 f"the dipping seam's track was cut to {len(track)} interior points; "
                 f"its stations are moves the search can make and must be kept: "
                 f"{[round(x, 3) for x, _y in track]}")
    if len(track) >= 3:
        (xa, ya), (xc, yc) = track[0], track[-1]
        span = math.hypot(xc - xa, yc - ya)
        worst = max(abs((xc - xa) * (yb - ya) - (yc - ya) * (xb - xa)) / span
                    for xb, yb in track[1:-1])
        out += _fail(worst <= 1e-6 * span,
                     f"the dipping seam's track bends by {worst:.3g} over its own "
                     f"length, so it is not the straight run this half is about")
        rise = abs(yc - ya)
        level_tol = _ON_THE_LINE * max(abs(xc - xa), rise)
        out += _fail(rise > 100 * level_tol,
                     f"the dipping seam's track rises {rise:.4g} over its length, "
                     f"which the level test ({level_tol:.3g}) can barely tell from "
                     f"flat -- both halves of the check are testing the same thing")
    out += _fail(not _inert_points(dipping),
                 f"the dipping seam's track carries an inert point at "
                 f"{_inert_points(dipping)}")
    return out


def test_corpus_surface_carries_no_spare_points():
    """The tutorial model's generated surface is the four points it needs.

    Sample 7's seam is flat and its track is one horizontal segment between two
    end ramps, so the whole surface is the two ground contacts and the two ramp
    junctions. Pinned on the shipped file because this is the surface the LEM-5
    page quotes and a reader compares against.
    """
    if not os.path.exists(SAMPLE7):
        return [f"missing {SAMPLE7}"]
    import contextlib
    import io

    from xslope.fileio import load_slope_data

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        sd = load_slope_data(SAMPLE7)
    result = G.generate_noncircular_surface(sd, report=True)
    if not result["surface"]:
        return [f"sample 7: nothing generated ({result['reason'][:90]})"]
    surface = result["surface"]
    out = _geometry_failures_common(sd, result, "sample 7")
    out += _fail(len(surface) == 4,
                 f"sample 7 generated {len(surface)} points; its flat seam admits "
                 f"four -- two ground contacts and two ramp junctions: "
                 f"{[(round(p['X'], 3), round(p['Y'], 3)) for p in surface]}")
    out += _fail(f"a {len(surface)}-point surface" in result["summary"],
                 f"the summary does not state the point count it built: "
                 f"{result['summary'][-140:]!r}")
    return out


def test_generator_is_deterministic():
    """The same model twice gives the same surface, coordinate for coordinate."""
    first = G.generate_noncircular_surface(_model())
    second = G.generate_noncircular_surface(_model())
    return _fail(first == second,
                 "two runs on the same model produced different surfaces")


def test_declines_state_a_reason():
    """Every refusal carries a sentence, and never a half-built surface."""
    out = []
    flat = _model()
    flat["ground_surface"] = LineString([(0.0, 0.0), (60.0, 0.0)])
    for label, sd in (("flat ground", flat),
                      ("no geometry", {"ground_surface": LineString([]),
                                       "polygons": [], "materials": []})):
        res = G.generate_noncircular_surface(sd, report=True)
        out += _fail(not res["surface"], f"{label}: a surface was built anyway")
        out += _fail(bool(res["reason"]), f"{label}: refused with no reason given")
        out += _fail(G.generate_noncircular_surface(sd) == [],
                     f"{label}: the bare form did not return an empty surface")
    return out


# ---------------------------------------------------------------------------
# D. The slicer accepts what the generator builds
# ---------------------------------------------------------------------------

def test_corpus_surface_slices():
    """A corpus weak-seam model's generated surface is one the slicer accepts.

    The tangency check: an end point sitting exactly ON the ground surface is a
    coordinate agreement between two independent calculations, and when they
    disagree in the last bits the slicer refuses a surface whose geometry is
    correct. Only running it finds that.
    """
    if not os.path.exists(ACADS):
        return [f"missing {ACADS}"]
    import contextlib
    import io

    from xslope.fileio import load_slope_data
    from xslope.slice import generate_slices

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        sd = load_slope_data(ACADS)
    out = []
    result = G.generate_noncircular_surface(sd, report=True)
    if not result["surface"]:
        return [f"ACADS weak layer: nothing generated ({result['reason'][:90]})"]
    out += _fail(result["zone"].name == "Weak Layer",
                 f"ACADS weak layer: seeded on {result['zone'].name!r}")
    with contextlib.redirect_stdout(buf):
        ok, res = generate_slices(sd, non_circ=result["surface"], num_slices=50,
                                  check_inputs=False)
    if not ok:
        return out + [f"ACADS weak layer: generate_slices refused the generated "
                      f"surface -- {res}"]
    slices, _surface = res
    base = result["zone"].polygon.bounds
    inside = sum(1 for y in slices["y_cb"] if base[1] - 1e-9 <= y <= base[3] + 1e-9)
    out += _fail(inside >= 0.5 * len(slices),
                 f"ACADS weak layer: only {inside} of {len(slices)} slice bases fall "
                 f"inside the weak zone's elevation band -- the track is not "
                 f"tracking the seam")
    return out


# ---------------------------------------------------------------------------
# E. Studio: the button and the picker
# ---------------------------------------------------------------------------

def test_studio_button():
    """The Generate button reaches the generator, and dims with a reason when it
    cannot run."""
    import contextlib
    import io

    from xslope.fileio import load_slope_data
    from studio.editors import CATEGORY_EDITORS

    editor = CATEGORY_EDITORS["non_circ"]
    buf = io.StringIO()
    out = []

    with contextlib.redirect_stdout(buf):
        sd = load_slope_data(ACADS)
    dlg = editor.build(sd, None)
    button = getattr(dlg, "generate_button", None)
    if button is None:
        return ["the non-circular editor has no Generate button"]
    out += _fail(button.isEnabled(),
                 "the Generate button is dimmed on a model with a clear weak zone")
    # A model with a clear weak zone must be built WITHOUT asking. Any dialog
    # raised here would block a headless run forever, which is also exactly what it
    # would do to a user who did not need to be asked, so refuse rather than wait.
    from PySide6.QtWidgets import QDialog
    real_exec, raised = QDialog.exec, []

    def no_dialog(self):
        raised.append(type(self).__name__)
        return 0
    QDialog.exec = no_dialog
    try:
        rows, message, reason = dlg._generate["propose"](None)
    finally:
        QDialog.exec = real_exec
    out += _fail(not raised,
                 f"a model with a clear weak zone raised {raised} instead of just "
                 f"generating")
    out += _fail(bool(rows), f"the button's generator produced nothing ({reason})")
    out += _fail("Weak Layer" in message,
                 f"the confirmation does not name the zone: {message[:100]!r}")
    if rows:
        # Press the button for real -- propose, confirm, replace -- rather than
        # calling the halves separately. The wiring between them is the thing this
        # check exists for.
        from PySide6.QtWidgets import QMessageBox
        real_box_exec, confirmed = QMessageBox.exec, []

        def apply_it(self):
            confirmed.append(self.informativeText())
            return QMessageBox.Apply
        QMessageBox.exec = apply_it
        QDialog.exec = no_dialog
        try:
            dlg._run_generate()
        finally:
            QMessageBox.exec = real_box_exec
            QDialog.exec = real_exec
        out += _fail(len(confirmed) == 1,
                     f"pressing Generate raised {len(confirmed)} confirmations; it "
                     f"must say what it will do before replacing the surface")
        out += _fail(confirmed and "Weak Layer" in confirmed[0],
                     "the confirmation does not name the zone it will seed on")
        after = dlg.result_rows()
        out += _fail(len(after) == len(rows),
                     f"pressing Generate left {len(after)} rows in the table, not the "
                     f"{len(rows)} the generator proposed")
        out += _fail(all(abs(a["X"] - b["X"]) < 1e-9 and abs(a["Y"] - b["Y"]) < 1e-9
                         and a["Movement"] == b["Movement"]
                         for a, b in zip(after, rows)),
                     "the rows in the table are not the ones the generator built")
        out += _fail(all(isinstance(r.get("Y"), float) and not math.isnan(r["Y"])
                         for r in after),
                     "a row reached the table without an explicit numeric Y")
        out += _fail(editor.apply.__name__ == "apply", "editor lost its apply hook")
        editor.apply(sd, dlg)
        out += _fail(len(sd["non_circ"]) == len(rows),
                     "apply() did not write the generated surface into the model")
    dlg.deleteLater()

    with contextlib.redirect_stdout(buf):
        sd1 = load_slope_data(ONE_ZONE)
    dlg1 = editor.build(sd1, None)
    out += _fail(not dlg1.generate_button.isEnabled(),
                 "the Generate button is live on a single-material model")
    tip = dlg1.generate_button.toolTip()
    out += _fail("one material zone" in tip,
                 f"the dimmed button's tooltip does not say why: {tip[:100]!r}")
    dlg1.deleteLater()
    return out


def test_studio_zone_picker():
    """An ambiguous model raises the picker, and its choice is what gets built."""
    import contextlib
    import io

    from PySide6.QtWidgets import QDialog

    from xslope.fileio import load_slope_data
    from studio.editors import CATEGORY_EDITORS, WeakZoneDialog

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        sd = load_slope_data(AMBIGUOUS)
    out = []
    zones = G.rank_weak_zones(sd)
    if len(zones) < 2:
        return [f"{os.path.basename(AMBIGUOUS)} is not an ambiguous model"]
    auto = G.generate_noncircular_surface(sd, report=True)
    out += _fail(not auto["surface"],
                 f"{os.path.basename(AMBIGUOUS)} auto-picked; it is meant to be "
                 f"ambiguous at the shipped separation")

    dlg = WeakZoneDialog(zones, None, headline="pick one", preselect=1)
    out += _fail(dlg.table.rowCount() == len(zones),
                 f"the picker listed {dlg.table.rowCount()} zones, expected "
                 f"{len(zones)}")
    out += _fail(dlg.chosen_zone() == zones[1].index,
                 "the picker's pre-selection is not the zone it was given")
    for row, zone in enumerate(zones):
        out += _fail(dlg.table.item(row, 0).text() == zone.name,
                     f"picker row {row} does not name its zone")
        out += _fail(not dlg.table.item(row, 0).icon().isNull(),
                     f"picker row {row} carries no material colour swatch")
        out += _fail(f"{zone.tau:.4g}" == dlg.table.item(row, 2).text(),
                     f"picker row {row} does not show the computed strength")
    dlg.deleteLater()

    # The button's propose() must RAISE the picker on this model and build what it
    # returns. Accept the dialog with the second zone selected, without a user.
    editor = CATEGORY_EDITORS["non_circ"]
    editor_dlg = editor.build(sd, None)
    seen = {}
    real_exec = QDialog.exec

    def fake_exec(self):
        if isinstance(self, WeakZoneDialog):
            seen["shown"] = True
            self.table.selectRow(1)
            return 1
        return real_exec(self)                       # pragma: no cover

    QDialog.exec = fake_exec
    try:
        rows, message, reason = editor_dlg._generate["propose"](None)
    finally:
        QDialog.exec = real_exec
    out += _fail(seen.get("shown"),
                 "the ambiguous model did not raise the zone picker")
    out += _fail(bool(rows),
                 f"nothing was built after the picker chose a zone ({reason})")
    expected = G.generate_noncircular_surface(sd, zone=zones[1].index)
    out += _fail(rows == expected,
                 "the surface built after the picker is not the one for the zone "
                 "the picker chose")
    editor_dlg.deleteLater()
    return out


CHECKS = [
    ("metric ranks on mobilizable strength", test_metric_ranks_on_mobilizable_strength),
    ("metric spans every strength model", test_metric_spans_every_strength_model),
    ("threshold picks and asks on both sides", test_threshold_both_sides),
    ("shipped threshold is the tuned one", test_shipped_threshold_is_the_tuned_one),
    ("geometry invariants, both facings", test_geometry_both_facings),
    ("outcropping zone needs no ramp", test_outcropping_zone_needs_no_ramp),
    ("pinched-out zone uses its continuous run", test_pinched_out_zone_uses_its_continuous_run),
    ("ramps are an active/passive pair", test_wedge_ramps_are_active_and_passive),
    ("level run collapses, dipping run does not",
     test_level_run_collapses_but_dipping_run_does_not),
    ("corpus surface carries no spare points",
     test_corpus_surface_carries_no_spare_points),
    ("generator is deterministic", test_generator_is_deterministic),
    ("refusals state a reason", test_declines_state_a_reason),
    ("corpus surface slices in the seam", test_corpus_surface_slices),
    ("studio button reaches the generator", test_studio_button),
    ("studio zone picker fires and is honoured", test_studio_zone_picker),
]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    try:
        import PySide6                                    # noqa: F401
        checks = CHECKS
    except Exception:
        print("non-circular generator: PySide6 not installed — Studio checks skipped.")
        checks = [c for c in CHECKS
                  if c[1] not in (test_studio_button, test_studio_zone_picker)]
    failures = []
    for name, fn in checks:
        try:
            fs = fn()
        except Exception as exc:
            import traceback
            traceback.print_exc()
            fs = [f"{name} raised: {exc!r}"]
        print(f"  {name:44s} {'ok' if not fs else f'FAIL ({len(fs)})'}")
        failures += fs
    return failures


def main():
    print("non-circular starting surface:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll non-circular generator checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication
        _app = QApplication.instance() or QApplication([])
    except Exception:
        pass
    main()
