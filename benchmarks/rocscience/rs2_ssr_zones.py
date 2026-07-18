"""Extract RS2 "SSR Search Area" polygons from a vendor ``.fez`` / ``.fea``.

RS2's native slope-stability verification models carry their per-case
strength-reduction search constraint as an exact vertex polygon in the ``.fea``
text, in a ``SSR_polygonal_zones:`` block::

    SSR_polygonal_zones:
    1                 <- number of polygons
    10 0 0            <- vertex count N (of the first polygon), then "0 0"
    8.516 12.255      <- vertex 1
    ...               <- ... N-1 more "x y" lines ...
    8.516 12.255      <- vertex N == vertex 1 (the ring is closed)

The polygon is RS2's own constraint, read verbatim from the vendor model file —
it is NOT authored here. Its semantics (confirmed against the geometry and the
manual Part 3 figures): strength reduction is applied only to material INSIDE the
polygon; everything outside keeps full strength. That is the mask ``solve_ssrm``'s
``ssr_zone`` argument builds — elements whose centroid lies outside the polygon
are excluded from reduction.

This helper is used by the RS2-61 / RS2-64 corpus builds and to record each case's
constraint polygon in docs/verification/rs2.md. It also re-exports the material
strengths (via xslope.rs2.read_fez) so a corpus .xlsx case can be matched to its
vendor .fez BY CONTENT rather than by filename ordering.

Run from the repo root, e.g.::

    python benchmarks/rocscience/rs2_ssr_zones.py path/to/'slope stability #061_02.fez'
"""

import os
import re
import sys
import zipfile

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


def read_ssr_zones(path):
    """Parse the ``SSR_polygonal_zones:`` block of a ``.fez`` / ``.fea``.

    Returns a list of polygons, one per SSR search area, each a list of ``(x, y)``
    vertices exactly as written in the file (the closing vertex, a repeat of the
    first, is kept so the polygon is recorded verbatim). Returns ``[]`` when the
    model declares no SSR search area (unconstrained SSR over the whole domain).
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
    polys = []
    for _ in range(n_poly):
        # Header line: "<n_vertices> 0 0"
        header = lines[j].split()
        n_vertices = int(header[0])
        j += 1
        verts = []
        for _ in range(n_vertices):
            xs, ys = lines[j].split()[:2]
            verts.append((float(xs), float(ys)))
            j += 1
        polys.append(verts)
    return polys


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


def zone_to_tag(verts):
    """Flatten a polygon ``[(x, y), ...]`` to the ``x1;y1;x2;y2;...`` tag string
    consumed by run_tests.py's ``ssr_zone=`` field."""
    return ";".join(f"{v:g}" for xy in verts for v in xy)


def summarize(path):
    """One-line-per-entity summary of a vendor file: materials (name/c/phi/gamma)
    and each SSR polygon's vertex count and bounding box. For case matching."""
    d = read_fez(path)
    zones = read_ssr_zones(path)
    xs = [v[0] for z in zones for v in z]
    ys = [v[1] for z in zones for v in z]
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
            print(f"  zone {k}: {len(z)} verts -> {zone_to_tag(z)}")
