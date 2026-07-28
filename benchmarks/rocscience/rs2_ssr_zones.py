"""Extract RS2 SSR constraint polygons ("SSR Search Area" / "SSR Exclusion Area")
from a vendor ``.fez`` / ``.fea``.

RS2's slope-stability verification models carry their per-case strength-reduction
constraint as an exact vertex polygon in the ``.fea`` text, in a
``SSR_polygonal_zones:`` block::

    SSR_polygonal_zones:
    1                 <- number of polygons
    10 0 1            <- vertex count N, then two flags (see below)
    8.516 12.255      <- vertex 1
    ...               <- ... N-1 more "x y" lines ...
    8.516 12.255      <- vertex N == vertex 1 (the ring is closed)

**The THIRD header field is the zone kind**, and it is the one that matters:

* ``0`` — the ring is a **Search Area**: strength reduction is applied only to
  material INSIDE it, and everything outside keeps full strength.
* ``1`` — the ring is an **Exclusion Area**: material inside it is NOT reduced,
  and everything outside is.

That decode is validated across all 77 polygons in the RS2 verification archive,
and corroborated by the vendor's own published figures: RS2 Part 4 Fig. 5.3 labels
the kind-1 ring of ``#005.fez`` "SSR Exclusion Area", RS2 Part 4 Fig. 68.3 labels
the kind-0 ring of ``#068.fez`` "SSR Search Area", and ``#067.fez``'s kind-1 ring
is the documented "SSR Exclusion Area below El. 81".

The SECOND field's meaning is **not established**, and it must not be read as a
"whole-model default / inert ring" marker: all six polygons in the archive that
carry ``1`` there are live constraints that genuinely clip their models (the
``#067_04`` / ``#067_06`` rectangles cut the dam at x = 102.3). It is preserved
verbatim on :attr:`SSRZone.header_flag` and carried no further.

The polygons are RS2's own constraints, read verbatim from the vendor model
files — they are NOT authored here.

**The two kinds are not interchangeable, and this module refuses to mix them.**
``solve_ssrm``'s ``ssr_zone`` argument (and the ``ssr_zone=`` test tag that feeds
it) means *search area*: elements whose centroid lies outside the polygon are held
at full strength. An exclusion area is the opposite mask, so it is served to that
argument only through :func:`exclusion_to_search_area`, which returns the
polygon's **complement** inside the model domain. Handing an exclusion ring
straight to a search-area consumer would silently invert the constraint, so
:func:`search_areas`, :func:`exclusion_areas` and :func:`zone_to_tag` raise
:class:`SSRZoneKindError` rather than serve the wrong kind.

This helper is used by the RS2-4 / RS2-61 / RS2-64 / Part IV VP68 corpus builds and
to record each case's constraint polygon in docs/verification/rs2.md. It also
re-exports the material strengths (via xslope.rs2.read_fez) so a corpus .xlsx case
can be matched to its vendor .fez BY CONTENT rather than by filename ordering.

Run from the repo root, e.g.::

    python benchmarks/rocscience/rs2_ssr_zones.py path/to/'slope stability #061_02.fez'
"""

import os
import re
import sys
import zipfile
from collections import namedtuple

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from xslope.rs2 import read_fez  # noqa: E402


def _fea_text(path):
    """Return the .fea model text from a ``.fez`` zip or a raw ``.fea`` file."""
    if zipfile.is_zipfile(path):
        zf = zipfile.ZipFile(path)
        fea = next((n for n in zf.namelist() if n.lower().endswith(".fea")), None)
        if fea is None:
            raise ValueError(f"{os.path.basename(path)} has no .fea model inside.")
        return zf.read(fea).decode("latin-1")
    with open(path, "rb") as fh:
        return fh.read().decode("latin-1")


class SSRZoneKindError(ValueError):
    """Raised when an SSR *exclusion* polygon is asked to stand in for a *search*
    area, or the reverse. The two rings mask opposite halves of the model, so
    substituting one for the other inverts the constraint silently — this module
    refuses instead."""


class SSRZone(namedtuple("SSRZone", "verts header_flag is_exclusion")):
    """One ``SSR_polygonal_zones`` ring, with its zone kind decoded.

    ``verts`` is the ring exactly as written in the vendor file (the closing
    vertex, a repeat of the first, is kept so the polygon is recorded verbatim).
    ``is_exclusion`` selects the mask sense, from the header's third field.
    ``header_flag`` is the header's second field, preserved verbatim because its
    meaning is not established — in particular it does NOT mark an inert or
    whole-model-default ring, since every polygon in the archive carrying it is a
    live constraint that clips its model.
    """

    __slots__ = ()

    @property
    def kind(self):
        """``'exclusion'`` or ``'search'`` — the vendor's own two zone kinds."""
        return "exclusion" if self.is_exclusion else "search"

    def __len__(self):
        return len(self.verts)

    def __repr__(self):
        return f"<SSRZone {self.kind} {len(self.verts)} verts flag={self.header_flag:d}>"


def read_ssr_zones(path):
    """Parse the ``SSR_polygonal_zones:`` block of a ``.fez`` / ``.fea``.

    Returns a list of :class:`SSRZone` — one per ring, carrying the ring's
    vertices AND its decoded zone kind, so a caller can tell a Search Area from an
    Exclusion Area. Returns ``[]`` when the model declares no SSR zone at all
    (unconstrained SSR over the whole domain).

    The typed return is deliberate: an earlier version of this reader discarded
    the header fields and returned bare vertex lists, which would have served an
    exclusion ring to a search-area consumer as though the two meant the same
    thing. Use :func:`search_areas` / :func:`exclusion_areas` to consume one kind.
    """
    lines = _fea_text(path).splitlines()
    i = next((k for k, ln in enumerate(lines)
              if ln.strip() == "SSR_polygonal_zones:"), -1)
    if i < 0:
        return []
    j = i + 1
    try:
        n_poly = int(lines[j].strip())
    except (ValueError, IndexError):
        return []
    j += 1
    zones = []
    for _ in range(n_poly):
        # Header line: "<n_vertices> <undecoded flag> <is_exclusion>"
        header = lines[j].split()
        n_vertices = int(header[0])
        flag = int(header[1]) if len(header) > 1 else 0
        is_excl = bool(int(header[2])) if len(header) > 2 else False
        j += 1
        verts = []
        for _ in range(n_vertices):
            xs, ys = lines[j].split()[:2]
            verts.append((float(xs), float(ys)))
            j += 1
        zones.append(SSRZone(verts, flag, is_excl))
    return zones


def _zones_of_kind(path, want_exclusion, label):
    zones = read_ssr_zones(path)
    if not zones:
        return []
    wanted = [z for z in zones if bool(z.is_exclusion) == want_exclusion]
    if not wanted:
        other = ", ".join(sorted({z.kind for z in zones}))
        raise SSRZoneKindError(
            f"{os.path.basename(str(path))} declares no SSR {label} — its "
            f"{len(zones)} zone(s) are {other} area(s). An exclusion area masks the "
            f"complement of a search area; convert it with exclusion_to_search_area() "
            f"instead of consuming it as one.")
    return wanted


def search_areas(path):
    """The model's SSR **Search Area** rings as :class:`SSRZone` objects.

    Raises :class:`SSRZoneKindError` when the model declares SSR zones but none of
    them is a search area — the caller is about to use an exclusion ring as though
    reduction ran inside it. Returns ``[]`` for an unconstrained model."""
    return _zones_of_kind(path, False, "Search Area")


def exclusion_areas(path):
    """The model's SSR **Exclusion Area** rings as :class:`SSRZone` objects.

    Raises :class:`SSRZoneKindError` when the model declares SSR zones but none of
    them is an exclusion area. Returns ``[]`` for an unconstrained model."""
    return _zones_of_kind(path, True, "Exclusion Area")


def exclusion_to_search_area(zone, domain):
    """Convert an SSR **Exclusion Area** into the equivalent **search area** —
    its complement inside the model domain.

    ``solve_ssrm``'s ``ssr_zone`` argument has one sense only: strength reduction
    is applied inside the polygon and everything outside is held at full strength.
    RS2's exclusion areas are the opposite mask, so the faithful translation of
    "do not reduce inside this ring" is "reduce inside everything else" — i.e.
    ``domain - ring``. This is the conversion the RS2-4 companion lock uses to
    reproduce RS2's own Part 4 run of the Talbingo dam, whose ``#005.fez`` holds
    the downstream benched shell at full strength.

    ``zone`` is an :class:`SSRZone` (an exclusion one — a search ring raises) or a
    bare vertex list; ``domain`` is the model outline, a shapely polygon or a
    vertex list. Returns the complement's exterior ring as ``[(x, y), ...]``,
    ready for :func:`zone_to_tag`. Raises when the complement is empty or breaks
    into disjoint pieces, since ``ssr_zone`` takes a single ring.
    """
    from shapely.geometry import Polygon as _Polygon
    if isinstance(zone, SSRZone):
        if not zone.is_exclusion:
            raise SSRZoneKindError(
                "exclusion_to_search_area() was handed a Search Area; a search "
                "area is already ssr_zone's own sense and needs no complement.")
        verts = zone.verts
    else:
        verts = list(zone)
    dom = domain if hasattr(domain, "exterior") else _Polygon(domain)
    rest = dom.difference(_Polygon(verts).buffer(0)).buffer(0)
    if rest.is_empty:
        raise ValueError("exclusion area covers the whole domain: nothing left to reduce")
    if rest.geom_type != "Polygon":
        raise ValueError(
            f"the complement is a {rest.geom_type} of {len(rest.geoms)} pieces; "
            "ssr_zone takes a single ring, so this case needs a hand-authored zone")
    return [(float(x), float(y)) for x, y in rest.exterior.coords]


def read_mc_footprint(path, mc_names=None, simplify_tol=0.05):
    """Union polygon of the vendor mesh elements whose material CAN yield (Mohr-
    Coulomb), i.e. the "MC corridor" of the material partition.

    RS2's #64 verification models constrain the SSR two ways: the SSR_polygonal_zones
    polygon AND a material partition — the Mohr-Coulomb material is placed only in a
    corridor/wedge hugging the proposed slip surface, while the rest of the domain is
    a second material with "Plasticity Specifications: Non" (linear-elastic, cannot
    yield). This reads the per-element ``material:`` tags in the ``elements:`` block,
    unions the MC elements' polygons, and returns that footprint — the authoritative
    "where a mechanism may form" region.

    Returns (shapely (Multi)Polygon or None, counts dict). ``mc_names`` defaults to
    the materials whose ``model`` is 'MohrCoulomb' (from read_fez). None is returned
    when every element is MC (no partition placed — the whole domain can yield).
    """
    from shapely.geometry import Polygon as _Polygon
    from shapely.ops import unary_union
    if mc_names is None:
        mats = read_fez(path)["materials"]
        mc_names = [m["name"] for m in mats if m["model"] == "MohrCoulomb"]
    mc_set = set(mc_names)
    lines = _fea_text(path).splitlines()

    # nodes: "<id> x: <x> y: <y>"
    coord = {}
    i = next((k for k, ln in enumerate(lines) if ln.startswith("nodes:")), -1)
    j = i + 1
    while j < len(lines):
        m = re.match(r"\s*(\d+)\s+x:\s*(\S+)\s+y:\s*(\S+)", lines[j])
        if m:
            coord[int(m.group(1))] = (float(m.group(2)), float(m.group(3)))
        elif lines[j].startswith("elements:"):
            break
        j += 1

    # elements: "<id> type: <t> nodes: [n1,n2,...] material: <name> materialID: <id>"
    e = next((k for k, ln in enumerate(lines) if ln.startswith("elements:")), -1)
    polys = []
    counts = {}
    for k in range(e + 1, len(lines)):
        ln = lines[k]
        mm = re.match(r"\s*\d+\s+type:\s*\S+\s+nodes:\s*\[([\d,\s]+)\]\s+material:\s*(\S+)", ln)
        if not mm:
            if ln.strip() and not ln[0].isdigit() and ":" in ln and not ln.startswith(" "):
                break  # left the elements block
            continue
        ids = [int(x) for x in mm.group(1).split(",")]
        name = mm.group(2)
        counts[name] = counts.get(name, 0) + 1
        if name in mc_set:
            corners = ids[0::2]  # even indices = element corner nodes (tri6/quad8/9)
            pts = [coord[c] for c in corners if c in coord]
            if len(pts) >= 3:
                p = _Polygon(pts)
                if p.is_valid and p.area > 0:
                    polys.append(p)
    total = sum(counts.values())
    n_mc = sum(counts.get(n, 0) for n in mc_set)
    if not polys or n_mc == total:
        return None, counts          # no partition placed: whole domain can yield
    foot = unary_union(polys)
    if simplify_tol:
        foot = foot.simplify(simplify_tol)
    return foot, counts


def zone_to_tag(zone):
    """Flatten a polygon ``[(x, y), ...]`` to the ``x1;y1;x2;y2;...`` tag string
    consumed by run_tests.py's ``ssr_zone=`` field.

    Accepts an :class:`SSRZone` or a bare vertex list. An exclusion :class:`SSRZone`
    raises :class:`SSRZoneKindError`: the ``ssr_zone=`` field means *search area*,
    so writing an exclusion ring into it would invert the constraint. Convert it
    with :func:`exclusion_to_search_area` first."""
    if isinstance(zone, SSRZone):
        if zone.is_exclusion:
            raise SSRZoneKindError(
                "ssr_zone= means SEARCH area (reduce inside); this ring is an SSR "
                "Exclusion Area (do NOT reduce inside). Pass it through "
                "exclusion_to_search_area(zone, domain) first.")
        verts = zone.verts
    else:
        verts = zone
    return ";".join(f"{v:g}" for xy in verts for v in xy)


def corridor_centerline(verts, n=40):
    """Medial line between the two long edges of a thin corridor polygon.

    RS2's SSR Search-Area / material corridor is a narrow band hugging Teoman's
    (unpublished) digitized slip surface: a closed ring that runs down one long
    edge and back along the other. This returns the band's centerline — a simple,
    reproducible medial line usable as a non-circular LEM surface.

    Method (deliberately simple, no medial-axis library):

    1. Drop any closing duplicate vertex; take the ring vertices.
    2. Principal axis (PCA/SVD) of the vertex cloud = the band's long direction.
       Project every vertex onto it ("station" along the band).
    3. The two band tips are the vertices at min / max station.
    4. Split the ring into its two chains between the tips (the two long edges),
       orient both tip_lo -> tip_hi.
    5. At ``n`` equal stations, interpolate each edge's (x, y) vs station and
       average the two -> the medial point at that station.

    Returns ``[(x, y), ...]`` ordered left-to-right (ascending x), the slip
    surface approximation. Valid for a monotone (crest-to-toe) band, which is
    exactly the RS2 #64 case; a strongly re-entrant band would need a true
    medial axis and is out of scope here.
    """
    import numpy as np
    P = np.asarray([(float(x), float(y)) for x, y in verts], float)
    if len(P) > 1 and np.allclose(P[0], P[-1]):
        P = P[:-1]
    m = P.mean(axis=0)
    _, _, vt = np.linalg.svd(P - m, full_matrices=False)
    axis = vt[0]
    proj = (P - m) @ axis
    i_lo, i_hi = int(proj.argmin()), int(proj.argmax())
    k = len(P)

    def _chain(a, b):
        out, i = [a], a
        while i != b:
            i = (i + 1) % k
            out.append(i)
        return out

    A = P[_chain(i_lo, i_hi)]              # one long edge, tip_lo -> tip_hi
    B = P[_chain(i_hi, i_lo)][::-1]        # other long edge, tip_lo -> tip_hi

    def _interp(edge, stations):
        s = (edge - m) @ axis
        order = np.argsort(s)
        return np.column_stack([np.interp(stations, s[order], edge[order, 0]),
                                np.interp(stations, s[order], edge[order, 1])])

    stations = np.linspace(proj[i_lo], proj[i_hi], n)
    center = 0.5 * (_interp(A, stations) + _interp(B, stations))
    if center[0, 0] > center[-1, 0]:      # order left-to-right
        center = center[::-1]
    return [(float(x), float(y)) for x, y in center]


def summarize(path):
    """One-line-per-entity summary of a vendor file: materials (name/c/phi/gamma)
    and each SSR polygon's vertex count and bounding box. For case matching."""
    d = read_fez(path)
    zones = read_ssr_zones(path)
    xs = [v[0] for z in zones for v in z.verts]
    ys = [v[1] for z in zones for v in z.verts]
    out = {
        "file": os.path.basename(path),
        "materials": [(m["name"], m["c"], m["phi"], m.get("gamma", 0.0),
                       m["model"]) for m in d["materials"]],
        "n_zones": len(zones),
        "zones": zones,
        "zone_bbox": (min(xs), min(ys), max(xs), max(ys)) if xs else None,
        "gw_type": d.get("gw_type", ""),
        "srf": d.get("srf", {}),
    }
    return out


if __name__ == "__main__":
    for p in sys.argv[1:]:
        s = summarize(p)
        print(f"=== {s['file']} ===")
        print(f"  materials: {s['materials']}")
        print(f"  gw_type: {s['gw_type']!r}   srf keys: {list(s['srf'])}")
        print(f"  n_zones: {s['n_zones']}   bbox: {s['zone_bbox']}")
        for k, z in enumerate(s["zones"]):
            tag = ("(exclusion — complement it before use as ssr_zone)"
                   if z.is_exclusion else zone_to_tag(z))
            print(f"  zone {k}: {z.kind}, flag {z.header_flag}, "
                  f"{len(z)} verts -> {tag}")
