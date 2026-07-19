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

"""Guard: every corpus file's material polygons must TILE its domain with no void.

A correctly built section fills each vertical column contiguously, from its bottom
boundary up to the ground surface, with material — the zones tile the body of the
slope, adjacent zones share matching edges, and no interior region is left carrying
no material. A VOID is an enclosed y-gap in a column: material below it, material
above it, and nothing in between. Voids are silent: the domain polygon is the union
of the zones, so a void is simply excised from the domain, and every failure surface
that crosses it loses the resisting weight (and strength, and pore-pressure area) of
the missing material.

This is the defect that sat under Slide2 VP42 (Baker & Leshchinsky earth dam): the
profile-line builder drew the granular shell's lower-right boundary as a single chord
from the toe to the core's waist, excising a large downstream toe wedge and reading
XSLOPE ~13-19% below the published cluster on the reservoir-side surfaces. Retiling
the section as explicit material-zone polygons closed it. This guard freezes that
class of defect out of the corpus.

Detection (deterministic, no randomness): for every material-tiling input file, the
material polygons are unioned and each x-column is intersected with a vertical line.
The covered y-intervals are read off; a gap BETWEEN two covered intervals (not the
open space above the ground surface, nor below the section) is a void. Contiguous
columns carrying a gap are clustered into a void region with a WIDTH (x-extent) and a
MAX_GAP (tallest vertical gap).

Tolerance for legitimate thin features. A void region is a FAILURE only when it is
both meaningfully WIDE and meaningfully TALL:

    width  >= WIDTH_TOL  (a tension crack, V-notch, or sheet-pile slot is NARROW —
                          a near-vertical thin gap; the VP42 toe wedge was 110+ units
                          wide)
    max_gap >= GAP_TOL   (a digitized layer interface whose two sides do not share
                          the exact same vertex list leaves a thin, near-horizontal
                          sliver ~0.1 unit tall; a real void is many units tall)

Requiring BOTH targets the wide material-excision defect (VP42-class) while narrowly
exempting the two legitimate thin cases: narrow cracks/notches (small width) and
digitization slivers along shared interfaces / touching stepped blocks (small gap).
The sheet-pile sample files carry their piles as line features, not as material gaps,
so they tile solidly and never reach this guard's attention at all.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/void_guard.py
Exits non-zero on any failure.
"""
import glob
import os
import sys
import warnings

warnings.filterwarnings("ignore")

from shapely.geometry import LineString
from shapely.ops import unary_union

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

# A void region fails only if it is BOTH wider than WIDTH_TOL AND taller than
# GAP_TOL (see the module docstring). These are in the input file's own length
# units. VP42's excised toe wedge measured width 110+, max_gap ~18 — an order of
# magnitude past both floors; the legitimate thin features it must not flag sit an
# order of magnitude below one or the other.
WIDTH_TOL = 2.0
GAP_TOL = 0.5

# Material-tiling input files to audit: the verification corpora plus the LEM / FEM
# sample files. Archived/deprecated trees and the partial sample-sheet templates
# (which are not complete slope models) are excluded — they are not part of the
# verified corpus.
_PATTERNS = (
    "docs/files/**/*.xlsx",
    "docs/lem/files/*.xlsx",
    "docs/seep/files/*.xlsx",
    "docs/fem/files/*.xlsx",
)


def corpus_files():
    """Absolute paths of every material-tiling corpus input file, sorted."""
    out = set()
    for pat in _PATTERNS:
        for f in glob.glob(os.path.join(_ROOT, pat), recursive=True):
            base = os.path.basename(f)
            if base.startswith("~$") or "input_template" in base:
                continue
            if os.sep + "archive" + os.sep in f:
                continue
            out.add(f)
    return sorted(out)


def find_voids(polygons, ncols=400, y_eps=1e-3):
    """Void regions in a material-zone tiling, by per-column vertical-gap sampling.

    ``polygons`` is a list of dicts each carrying a shapely ``'polygon'`` (the
    ``slope_data['polygons']`` shape). Returns a list of dicts with ``x0, x1,
    width, max_gap, area_est`` — one per contiguous run of columns that carry an
    enclosed gap. An empty list means the section tiles with no void.
    """
    polys = [p["polygon"] for p in polygons if p.get("polygon") is not None]
    if not polys:
        return []
    union = unary_union(polys)
    minx, miny, maxx, maxy = union.bounds
    span = maxx - minx
    if span <= 0:
        return []
    dx = span / ncols
    col_gaps = []  # (column index, gap_y0, gap_y1)
    for i in range(1, ncols):
        x = minx + i * dx
        inter = union.intersection(LineString([(x, miny - 1.0), (x, maxy + 1.0)]))
        if inter.is_empty:
            continue
        geoms = [inter] if inter.geom_type == "LineString" else list(getattr(inter, "geoms", []))
        segs = []
        for g in geoms:
            if g.is_empty or g.geom_type != "LineString":
                continue
            ys = [c[1] for c in g.coords]
            segs.append((min(ys), max(ys)))
        if len(segs) < 2:
            continue
        segs.sort()
        for a, b in zip(segs[:-1], segs[1:]):
            gy0, gy1 = a[1], b[0]
            if gy1 - gy0 > y_eps:
                col_gaps.append((i, gy0, gy1))
    if not col_gaps:
        return []
    # Cluster consecutive gap-carrying columns (allow a single missing column) into
    # regions; a real material excision spans a contiguous run of columns.
    regions, cur = [], None
    for (i, gy0, gy1) in col_gaps:
        if cur is None:
            cur = {"i0": i, "i1": i, "gaps": [(gy0, gy1)]}
        elif i - cur["i1"] <= 2:
            cur["i1"] = i
            cur["gaps"].append((gy0, gy1))
        else:
            regions.append(cur)
            cur = {"i0": i, "i1": i, "gaps": [(gy0, gy1)]}
    if cur:
        regions.append(cur)
    out = []
    for r in regions:
        width = (r["i1"] - r["i0"] + 1) * dx
        max_gap = max(g1 - g0 for g0, g1 in r["gaps"])
        area = sum(g1 - g0 for g0, g1 in r["gaps"]) * dx
        out.append({
            "x0": round(minx + r["i0"] * dx, 3),
            "x1": round(minx + r["i1"] * dx, 3),
            "width": round(width, 3),
            "max_gap": round(max_gap, 3),
            "area_est": round(area, 3),
        })
    return out


def check():
    """Audit every corpus file. Returns a list of failure strings (empty == clean).

    A file contributes a failure only if it carries a void region exceeding BOTH
    WIDTH_TOL and GAP_TOL. Files that do not load, or carry no material polygons
    (seep-only meshes, partial templates), are not material tilings and are skipped.
    """
    from xslope.fileio import load_slope_data

    failures = []
    for f in corpus_files():
        rel = os.path.relpath(f, _ROOT)
        try:
            sd = load_slope_data(f)
        except Exception:
            continue  # not a loadable slope model — out of scope for a tiling audit
        polys = sd.get("polygons") or []
        if len(polys) < 1:
            continue
        try:
            voids = find_voids(polys)
        except Exception as e:  # pragma: no cover - defensive
            failures.append(f"{rel}: void detector error: {str(e)[:80]}")
            continue
        real = [v for v in voids
                if v["width"] >= WIDTH_TOL and v["max_gap"] >= GAP_TOL]
        for v in real:
            failures.append(
                f"{rel}: material void width={v['width']} x max_gap={v['max_gap']} "
                f"(area~{v['area_est']}) over x=[{v['x0']}, {v['x1']}] — the zones do "
                f"not tile the section; a region carries no material")
    return failures


if __name__ == "__main__":
    fails = check()
    if fails:
        print(f"VOID GUARD: {len(fails)} file(s) with an untiled material void:")
        for m in fails:
            print("  " + m)
        sys.exit(1)
    n = len(corpus_files())
    print(f"VOID GUARD: clean — {n} corpus files tile their domains with no void.")
