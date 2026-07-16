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

"""Rocscience RS2 (.fez) import.

RS2 is a 2D finite-element program; its slope-stability answer comes from a
Shear-Strength-Reduction (SSR) analysis, not a limit-equilibrium search. A
``.fez`` is a ZIP holding one sectioned ``key: value`` ASCII model file
(``.fea``) plus binary sidecars (``.p2m`` / ``.p2mrv`` / ``.config``). The
``.fea`` grammar is the same family as Slide2's ``.sli`` (see slide2.py), with
two dialect differences: it uses nested ``start:``/``end:`` blocks for geometry,
and a few keys carry their value on the FOLLOWING line (``version:\\n11.028``).

This module reads that text and builds an xslope ``slope_data``:

  - read_fez(path): parse the model into its raw entities.
  - fez_to_slope_data(d): build slope_data for the model, plus a list of caveats
    for anything that did not map cleanly.
  - import_fez(path, template, out_path): the script route — parse and write an
    xslope .xlsx input file.

Only what RS2 and xslope have in common is imported. Everything else is reported
as a caveat rather than silently dropped: the SSR/strength-reduction settings
(xslope's analog is the FEM ``solve_ssrm`` path — the LEM import records them as
metadata, it does not fabricate an equivalent), strength models other than
Mohr-Coulomb, finite-element pore-pressure fields and transient seepage, joints,
liners, bolts/piles, field stress and distributed/line loads.

RS2 defines NO limit-equilibrium slip surface (its SSR result is a field), so an
imported model arrives without a failure surface and must be given circles or a
non-circular surface before it will solve — the same as a search-only Slide2
model.

Format learned from the public RS2 verification model files (Rocscience help
site); version 11.028 seen.
"""

import os
import re
import zipfile

from shapely.geometry import LineString, Polygon
from shapely.ops import polygonize, unary_union

from .fileio import build_ground_surface_from_polygons


# Unit weight of water, used to price pore pressure and to name the unit system.
_GAMMA_W_METRIC = 9.81       # kN/m3   (RFCunits "Metric kPa")
_GAMMA_W_IMPERIAL = 62.4     # lb/ft3  (RFCunits "Imperial psf")

# RS2 plasticity model that maps 1:1 onto xslope's Mohr-Coulomb. Everything else
# is imported with whatever c/phi it carries and flagged, the same discipline
# geostudio.py and slide2.py use. "Non" is RS2's word for an elastic material with
# no failure criterion — it has no strength at all.
_MC_MODEL = "MohrCoulomb"
_MODEL_NAMES = {
    "Non": "elastic (no failure criterion)",
    "GeneralizedHoekBrown": "generalized Hoek-Brown",
    "DruckerPrager": "Drucker-Prager",
    "PowerCurve": "power-curve (nonlinear) strength",
}

# The SSR (Shear-Strength-Reduction) settings we surface as metadata. RS2 reduces
# c and tan(phi) by a factor and searches for the factor at which the FE model
# stops converging — its "critical SRF" is the factor of safety. xslope's analog
# is the FEM solve_ssrm path; the LEM importer records these numbers so the user
# can reproduce the run, and never invents an LEM search to stand in for them.
_SRF_KEYS = ("strength_reduction_analysis", "auto_SRF", "initial_SRF", "final_SRF",
             "change_in_SRF", "delta_FS", "maxiter_SRF", "tolerance_SRF",
             "tensilestrength_SRF")


# --------------------------------------------------------------------------------------
# Text grammar
# --------------------------------------------------------------------------------------

def _find(lines, header, start=0):
    """Index of the first line whose stripped text equals ``header``, or -1."""
    for i in range(start, len(lines)):
        if lines[i].strip() == header:
            return i
    return -1


def _block(lines, header):
    """The [start, end] line indices of a ``<name> start:`` / ``<name> end:`` block.

    ``header`` is the start line (e.g. ``'v6 geometry start:'``); the matching end
    swaps 'start' for 'end'. Returns (-1, -1) if the block is absent.
    """
    s = _find(lines, header)
    if s < 0:
        return -1, -1
    end_header = header.replace("start:", "end:")
    e = _find(lines, end_header, s + 1)
    return s, e


def _model_description(lines):
    """Parse the ``model description:`` section into a flat ``{key: value}`` dict.

    Handles the .fea dialect where a key's value can sit on the following line
    (``version:`` then ``11.028``): a ``key:`` with an empty value takes the next
    body line as its value when that line is itself not a ``key: value`` pair.
    """
    s = _find(lines, "model description:")
    if s < 0:
        return {}
    # The section runs to the first blank line (every entry is a key line, and the
    # dialect where a value sits on its own line — ``version:`` then ``11.028`` —
    # rules out breaking on "a line that ends in ':'", since ``version:`` does).
    body = []
    for i in range(s + 1, len(lines)):
        if lines[i].strip() == "":
            break
        body.append(lines[i])
    md = {}
    for i, ln in enumerate(body):
        m = re.match(r"\s*([\w ]+?):\s*(.*)$", ln)
        if not m:
            continue
        key, val = m.group(1).strip(), m.group(2).strip()
        if val == "" and i + 1 < len(body):
            nxt = body[i + 1]
            if nxt.strip() and ":" not in nxt:
                val = nxt.strip()
        md[key] = val
    return md


def _fnum(v, default=0.0):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def _read_dpoint_array(lines, i):
    """Read a ``dpoint array start:`` block beginning at/after line ``i``.

    Returns (list of (x, y), index just past the array end). The array is
    ``num points: N`` followed by ``k: x, y`` lines.
    """
    while i < len(lines) and "dpoint array start:" not in lines[i]:
        i += 1
    n = 0
    j = i + 1
    while j < len(lines):
        s = lines[j].strip()
        if s.startswith("num points:"):
            n = int(s.split(":", 1)[1])
            j += 1
            break
        j += 1
    pts = []
    for k in range(j, j + n):
        after = lines[k].split(":", 1)[1]
        x, y = after.split(",")
        pts.append((float(x), float(y)))
    return pts, j + n


def _parse_boundaries(lines):
    """The v6-geometry boundaries: ``[{'type': str, 'pts': [(x, y), ...]}, ...]``.

    Each boundary carries only its defining vertices (its discretization — the FE
    meshing subdivision — is skipped). ``type`` is 'external', 'material' or
    'pile_geotextile'.
    """
    s, e = _block(lines, "v6 geometry start:")
    if s < 0:
        return []
    boundaries = []
    i = s
    while i < e:
        m = re.match(r'\s*type:\s*"([^"]+)"', lines[i])
        if m:
            cur = {"type": m.group(1), "pts": []}
            boundaries.append(cur)
            i += 1
            continue
        if lines[i].strip() == "vertices start:":
            pts, i = _read_dpoint_array(lines, i)
            if boundaries:
                boundaries[-1]["pts"] = pts
            continue
        i += 1
    return boundaries


def _parse_material_mesh_seeds(lines):
    """The coarse material-mesh triangles: ``[(Polygon, material_number), ...]``.

    RS2 tags material zones with a sparse seed triangulation — not a full tiling,
    so it cannot build the zones by itself (unlike Slide2's cells), but each
    triangle lies inside exactly one zone and names its material. That is enough to
    label the polygons the boundaries polygonize into.
    """
    s, e = _block(lines, "materials mesh start:")
    if s < 0:
        return []
    seeds = []
    i = s
    while i < e:
        if lines[i].strip().startswith("P1:"):
            def pt(k):
                x, y = lines[k].split(":", 1)[1].split(",")
                return (float(x), float(y))
            tri = Polygon([pt(i), pt(i + 1), pt(i + 2)])
            mm = next((k for k in range(i, min(i + 8, e))
                       if lines[k].strip().startswith("stage 1:")), None)
            if mm is not None and tri.is_valid and tri.area > 0:
                seeds.append((tri, int(lines[mm].split(":", 1)[1])))
            i += 3
            continue
        i += 1
    return seeds


def _parse_material_types(lines):
    """The defined materials: ``[{'name','model','c','phi'}, ...]``, 1-based order.

    Read from the labelled ``material types:`` section (``material N: name`` /
    ``Plasticity Specifications: MODEL`` / ``C: c phi: phi``). Unit weight is NOT
    here — it is filled in from ``material properties:`` by the caller.
    """
    s = _find(lines, "material types:")
    if s < 0:
        return []
    # Section ends at the next column-0 header.
    e = s + 1
    while e < len(lines) and not (lines[e] and not lines[e][0].isspace()
                                  and lines[e].rstrip().endswith(":")):
        e += 1
    mats = []
    cur = None
    for i in range(s + 1, e):
        ln = lines[i]
        m = re.match(r"material\s+(\d+):\s*(.*)$", ln.strip())
        if m:
            cur = {"num": int(m.group(1)), "name": m.group(2).strip() or f"material {m.group(1)}",
                   "model": "", "c": 0.0, "phi": 0.0}
            mats.append(cur)
            continue
        if cur is None:
            continue
        m = re.match(r"Plasticity Specifications:\s*(\S+)", ln.strip())
        if m:
            cur["model"] = m.group(1)
            continue
        m = re.match(r"C:\s*(\S+)\s+phi:\s*(\S+)", ln.strip())
        if m:
            cur["c"] = _fnum(m.group(1))
            cur["phi"] = _fnum(m.group(2))
    return mats


def _parse_material_gammas(lines):
    """Unit weight per material number, from the ``material properties:`` table.

    The section is name/numeric-row pairs in material-number order; the unit weight
    is the SECOND field of the numeric row (verified against known ACADS values:
    field 1 = 20, 19.5). Returns ``{material_number: gamma}`` (1-based).
    """
    s = _find(lines, "material properties:")
    if s < 0:
        return {}
    e = s + 1
    while e < len(lines) and not (lines[e] and not lines[e][0].isspace()
                                  and lines[e].rstrip().endswith(":")):
        e += 1
    gammas = {}
    num = 0
    i = s + 1
    while i < e:
        row = lines[i].strip()
        # A numeric row starts with a number; a name row does not.
        if row and re.match(r"^-?\d", row):
            fields = row.split()
            if len(fields) >= 2:
                num += 1
                gammas[num] = _fnum(fields[1])
        i += 1
    return gammas


def _parse_piezos(lines):
    """The piezometric lines: ``{line_id: [(x, y), ...]}``.

    The ``piezos:`` section is: line count, then per line an id, a point count, and
    that many ``x y`` rows. Absent section (a dry model) -> empty dict.
    """
    s = _find(lines, "piezos:")
    if s < 0:
        return {}
    i = s + 1
    if i >= len(lines):
        return {}
    try:
        nlines = int(lines[i].strip())
    except (ValueError, IndexError):
        return {}
    i += 1
    out = {}
    for _ in range(nlines):
        if i + 1 >= len(lines):
            break
        pid = int(lines[i].strip())
        npts = int(lines[i + 1].strip())
        i += 2
        pts = []
        for k in range(i, i + npts):
            parts = lines[k].strip().split()
            if len(parts) >= 2:
                pts.append((float(parts[0]), float(parts[1])))
        out[pid] = pts
        i += npts
    return out


def _parse_material_piezos(lines):
    """Per-material piezo-line id (0 = none): ``{material_number: piezo_id}``."""
    s = _find(lines, "material piezos:")
    if s < 0:
        return {}
    out = {}
    num = 0
    for i in range(s + 1, len(lines)):
        row = lines[i].strip()
        if not row or not re.match(r"^-?\d", row):
            break
        num += 1
        out[num] = int(_fnum(row))
    return out


def _count_block_children(lines, header, child_prefix):
    """Count children (lines beginning with ``child_prefix``) inside a start/end block."""
    s, e = _block(lines, header)
    if s < 0:
        return 0
    return sum(1 for i in range(s + 1, e) if lines[i].strip().startswith(child_prefix))


# --------------------------------------------------------------------------------------
# Reader
# --------------------------------------------------------------------------------------

def read_fez(path):
    """Parse an RS2 ``.fez`` (or a raw ``.fea``) into its raw entities.

    Returns a dict with the model's version, its ``model description`` settings and
    the SSR settings pulled from them, the defined materials (name/model/c/phi/gamma),
    the geometry boundaries, the material-mesh seed triangles, the piezometric lines
    and their per-material assignment, and counts of the features that will be
    reported but not imported (joints, liners, bolts/piles, loads).

    Coordinates are as-authored; no unit conversion is applied here.

    Raises ValueError if the file is not a readable RS2 model.
    """
    if zipfile.is_zipfile(path):
        zf = zipfile.ZipFile(path)
        fea = next((n for n in zf.namelist() if n.lower().endswith(".fea")), None)
        if fea is None:
            raise ValueError(f"{os.path.basename(path)} is not an RS2 .fez "
                             f"(no .fea model inside).")
        text = zf.read(fea).decode("latin-1")
    else:
        with open(path, "rb") as fh:
            text = fh.read().decode("latin-1")
    if "model description:" not in text:
        raise ValueError(f"{os.path.basename(path)} is not an RS2 model "
                         f"(no 'model description:' section).")
    lines = text.splitlines()

    md = _model_description(lines)
    srf = {k: md[k] for k in _SRF_KEYS if k in md}

    # gw_type lives in the (indented) 'groundwater setup:' section, not the model
    # description — e.g. 'Static Analysis' | 'Finite Element Analysis' | 'Transient'.
    gw_type = ""
    for ln in lines:
        m = re.match(r'\s*gw_type:\s*"([^"]+)"', ln)
        if m:
            gw_type = m.group(1)
            break

    mats = _parse_material_types(lines)
    gammas = _parse_material_gammas(lines)
    for m in mats:
        m["gamma"] = gammas.get(m["num"], 0.0)

    return {
        "path": str(path),
        "version": md.get("version", ""),
        "title": md.get("title", ""),
        "model_description": md,
        "srf": srf,
        "units": md.get("RFCunits", ""),
        "materials": mats,                                   # 1-based order
        "boundaries": _parse_boundaries(lines),
        "seeds": _parse_material_mesh_seeds(lines),
        "piezos": _parse_piezos(lines),
        "material_piezos": _parse_material_piezos(lines),
        "gw_type": gw_type,
        "counts": {
            "distributed_loads": _count_block_children(
                lines, "new distributed loads start:", "distributed load"),
            "line_loads": _count_block_children(
                lines, "new line loads start:", "line load"),
            "bolts": _count_block_children(lines, "bolts start:", "bolt"),
            "joint_networks": _count_block_children(
                lines, "joint networks start:", "joint"),
        },
    }


# --------------------------------------------------------------------------------------
# slope_data construction
# --------------------------------------------------------------------------------------

def _blank_material(name):
    """A material with xslope's full key set, everything zeroed but the name —
    the same shape ``load_slope_data`` produces (mirrors slide2.py/geostudio.py)."""
    return {
        "name": str(name), "gamma": 0.0, "gamma_sat": None, "option": "", "c": 0.0,
        "phi": 0.0, "cp": 0.0, "r_elev": 0.0, "d": 0, "psi": 0, "u": "none", "ru": 0.0,
        "sigma_gamma": 0.0, "sigma_c": 0.0, "sigma_phi": 0.0, "sigma_cp": 0.0,
        "sigma_d": 0.0, "sigma_psi": 0.0, "k1": 0.0, "k2": 0.0, "alpha": 0.0,
        "unsat": "lf", "kr0": 0.0, "h0": 0.0, "vg_a": 0.0, "vg_n": 0.0,
        "E": 0.0, "nu": 0.0,
    }


def _detect_units(rfcunits):
    """(unit_system, gamma_water) from the RFCunits string. Fail loudly on the rest."""
    u = (rfcunits or "").lower()
    if "metric" in u:
        return "metric", _GAMMA_W_METRIC
    if "imperial" in u:
        return "imperial", _GAMMA_W_IMPERIAL
    # No explicit unit tag: default to metric but say so.
    return "unknown", _GAMMA_W_METRIC


def _polygonize_zones(boundaries, seeds, caveats):
    """Build ``polygons`` (each ``{'polygon', 'mat_id'}``) from the RS2 geometry.

    RS2 stores an external boundary polygon plus internal *material* boundaries —
    open polylines that cut the domain into zones. Node them together and polygonize
    to recover the zone faces, then label each face with the material of whichever
    seed triangle falls inside it. Returns (polygons, used_material_numbers).
    """
    ext = next((b for b in boundaries if b["type"] == "external"), None)
    if ext is None or len(ext["pts"]) < 3:
        raise ValueError("This model has no external boundary to import.")
    mat_bounds = [b for b in boundaries if b["type"] == "material"]

    edges = [LineString(ext["pts"] + [ext["pts"][0]])]
    for b in mat_bounds:
        if len(b["pts"]) >= 2:
            edges.append(LineString(b["pts"]))
    faces = list(polygonize(unary_union(edges)))
    if not faces:
        raise ValueError("This model's boundaries did not enclose any zone.")

    polygons = []
    used = []
    unlabelled = 0
    for f in faces:
        mat_num = None
        for tri, num in seeds:
            rp = tri.representative_point()
            if f.contains(rp) or f.touches(rp):
                mat_num = num
                break
        if mat_num is None:
            unlabelled += 1
            continue
        polygons.append({"polygon": Polygon(f.exterior.coords), "_mat_num": mat_num})
        if mat_num not in used:
            used.append(mat_num)
    if unlabelled:
        caveats.append(
            f"{unlabelled} geometry zone(s) had no RS2 material seed inside them and "
            f"were dropped — check the imported geometry is complete")
    if not polygons:
        raise ValueError("No geometry zone could be matched to a material.")

    used.sort()
    remap = {num: i for i, num in enumerate(used)}
    for p in polygons:
        p["mat_id"] = remap[p.pop("_mat_num")]
    return polygons, used


def fez_to_slope_data(d):
    """Build an xslope ``slope_data`` from a parsed RS2 model.

    Parameters
    ----------
    d : dict
        Output of :func:`read_fez`.

    Returns
    -------
    (slope_data, caveats) : (dict, list of str)
        ``caveats`` names everything RS2 defined that xslope could not take
        faithfully. It is never silently empty when something was dropped.

    Notes
    -----
    RS2's slope-stability result is a finite-element SSR field, not a limit-
    equilibrium surface, so NO failure surface is imported: the model arrives with
    geometry, materials and water, and you must define circles or a non-circular
    surface before it will solve (or reload — an input file with no surface does not
    validate). The SSR settings are recorded as a caveat, because xslope's analog is
    the separate FEM ``solve_ssrm`` path and inventing an LEM search would answer a
    different question than the file asks.
    """
    caveats = []
    unit_system, gamma_water = _detect_units(d.get("units"))
    if unit_system == "unknown":
        caveats.append(
            f"this RS2 model declares no unit system (RFCunits={d.get('units')!r}); "
            f"imported assuming metric with unit weight of water {gamma_water:g} — "
            f"check it before solving")

    polygons, used = _polygonize_zones(d["boundaries"], d["seeds"], caveats)

    # Materials, in the order the zones use them, renumbered 0-based to match mat_id.
    by_num = {m["num"]: m for m in d["materials"]}
    materials = []
    for num in used:
        src = by_num.get(num)
        name = src["name"] if src else f"material {num}"
        mat = _blank_material(name)
        if src is None:
            caveats.append(f"a zone uses RS2 material {num}, which the file does not "
                           f"define — imported with no strength; set it before solving")
            mat["option"] = "mc"
            materials.append(mat)
            continue
        model = src.get("model") or ""
        c, phi, gamma = src["c"], src["phi"], src["gamma"]
        if model and model != _MC_MODEL:
            label = _MODEL_NAMES.get(model, f"the {model} model")
            caveats.append(
                f"material '{name}' uses RS2's {label}, which xslope does not have — "
                f"imported as Mohr-Coulomb with c = {c:g}, phi = {phi:g}, which is NOT "
                f"the same strength; check it before solving")
        mat.update(option="mc", gamma=gamma, c=c, phi=phi)
        materials.append(mat)
        if not gamma:
            caveats.append(f"material '{name}' has no unit weight in the file — set "
                           f"gamma before solving")
        if not c and not phi:
            caveats.append(
                f"material '{name}' came across with NO STRENGTH (c = 0, phi = 0) — set "
                f"its strength before solving or the factor of safety is meaningless")

    ground_surface, domain_polygon = build_ground_surface_from_polygons(polygons)

    # --- water -------------------------------------------------------------------
    piezo_line = []
    mat_piezos = d.get("material_piezos", {})
    piezos = d.get("piezos", {})
    gw_type = d.get("gw_type", "")
    used_piezo_ids = sorted({mat_piezos.get(num, 0) for num in used} - {0})
    if used_piezo_ids and piezos:
        pid = used_piezo_ids[0]
        piezo_line = list(piezos.get(pid, []))
        for num, mat in zip(used, materials):
            if mat_piezos.get(num, 0):
                mat["u"] = "piezo"
        dry = [m["name"] for num, m in zip(used, materials) if not mat_piezos.get(num, 0)]
        if piezo_line and dry:
            caveats.append(
                f"material(s) {', '.join(repr(n) for n in dry)} are not connected to a "
                f"water table in RS2 — imported with no pore pressure")
        if len(used_piezo_ids) > 1:
            caveats.append(
                f"the model uses {len(used_piezo_ids)} piezometric lines; xslope takes "
                f"one, so the first was imported")
    elif gw_type and gw_type != "Static Analysis":
        caveats.append(
            f"THE FACTOR OF SAFETY WILL BE WRONG: pore pressure comes from RS2's "
            f"'{gw_type}' finite-element groundwater analysis, which xslope cannot read "
            f"from this file — imported as zero. Run the seepage analysis in xslope, or "
            f"set a piezo line, before solving")
    elif any(mat_piezos.get(num, 0) for num in used) and not piezos:
        caveats.append(
            "materials reference a water table but the file defines no piezometric line "
            "— pore pressure imported as zero")

    # --- SSR settings (metadata, not fabricated) ---------------------------------
    srf = d.get("srf", {})
    if (srf.get("strength_reduction_analysis") or "").upper() == "ON":
        parts = []
        for k in ("initial_SRF", "final_SRF", "change_in_SRF", "tolerance_SRF"):
            if k in srf:
                parts.append(f"{k}={srf[k]}")
        caveats.append(
            "this is an RS2 Shear-Strength-Reduction (SSR) finite-element analysis"
            + (" (" + ", ".join(parts) + ")" if parts else "")
            + ". Its critical SRF is the factor of safety, found by an FE run — not a "
            "limit-equilibrium search. xslope's analog is the FEM solve_ssrm path; "
            "these settings are recorded but no LEM search was created from them")

    # --- things reported but not imported ----------------------------------------
    counts = d.get("counts", {})
    if counts.get("joint_networks"):
        caveats.append(f"{counts['joint_networks']} joint network(s) are defined in RS2 "
                       f"but were NOT imported — xslope has no jointed-rock model")
    n_pile = sum(1 for b in d["boundaries"] if b["type"] == "pile_geotextile")
    if n_pile or counts.get("bolts"):
        caveats.append(
            f"{n_pile + counts.get('bolts', 0)} reinforcement element(s) (bolts / "
            f"piles / geotextiles) are defined in RS2 but were NOT imported — RS2's "
            f"structural support does not map onto xslope's reinforcement, and guessing "
            f"the force would be worse than omitting it; the factor of safety will be "
            f"lower than RS2's. Add reinforcement by hand if needed")
    if counts.get("distributed_loads") or counts.get("line_loads"):
        caveats.append(
            f"{counts.get('distributed_loads', 0)} distributed and "
            f"{counts.get('line_loads', 0)} line load(s) are defined in RS2 but were NOT "
            f"imported — add them by hand; the factor of safety will differ")

    # Pseudo-static seismic: ``seismic: kh kv f1 f2 f3`` — horizontal and vertical
    # coefficients first, then load factors (surveyed across the 298 public
    # verification models; every one is static, kh = kv = ±0). The horizontal
    # coefficient maps onto xslope's k; the vertical has no xslope analog.
    k_seismic = 0.0
    _seis = str(d.get("model_description", {}).get("seismic", "")).split()
    if _seis:
        kh = _fnum(_seis[0])
        kv = _fnum(_seis[1]) if len(_seis) > 1 else 0.0
        if kh:
            k_seismic = kh
            caveats.append(
                f"a horizontal seismic coefficient ({kh:g}) was imported as k — RS2's "
                f"sign/direction convention has not been verified against a solved "
                f"seismic RS2 model (none exist in the public verification set), so "
                f"check the pseudo-static force direction before trusting the result")
        if kv:
            caveats.append(
                f"a vertical seismic coefficient ({kv:g}) is set in RS2 — xslope "
                f"models only the horizontal one, so it was dropped")

    caveats.append(
        "no failure surface was imported — an RS2 slope-stability model defines an SSR "
        "finite-element analysis, which has no xslope slip surface to take. Define "
        "circles or a non-circular surface before solving (an input file with no "
        "surface will not re-load)")

    caveats.append(
        f"this model is in {unit_system} units — xslope is unit-agnostic and imports "
        f"the numbers as they are, so results come back in the same system")

    slope_data = {
        "template_version": None,
        "gamma_water": gamma_water,
        "tcrack_depth": 0.0,
        "tcrack_water": 0.0,
        "k_seismic": k_seismic,
        "max_depth": None,
        "profile_lines": [],
        "polygons": polygons,
        "domain_polygon": domain_polygon,
        "ground_surface": ground_surface,
        "tcrack_surface": None,
        "materials": materials,
        "piezo_line": piezo_line,
        "piezo_line2": [],
        "piezo_phreatic": True,
        "piezo_phreatic2": True,
        "circular": False,
        "circles": [],
        "non_circ": [],
        "dloads": [],
        "dloads2": [],
        "reinforce_lines": [],
        "reinforcement_lines": [],
        "pile_lines": [],
        "line_loads": [],
        "seepage_bc": {"specified_heads": [], "specified_fluxes": [], "exit_face": []},
        "seepage_bc2": {"specified_heads": [], "specified_fluxes": [], "exit_face": []},
        "has_seepage_bc2": False,
        "mesh": None,
    }
    return slope_data, caveats


def import_fez(path, template, out_path):
    """Read an RS2 ``.fez`` and write it to an xslope .xlsx input file.

    The script-level route, mirroring :func:`xslope.slide2.import_slmd`: parse the
    RS2 model, convert it to ``slope_data``, and write it through a copy of an input
    template.

    Returns the caveat list from :func:`fez_to_slope_data`. Read it — a clean return
    does not mean a complete model: an RS2 import never carries a failure surface,
    because an SSR analysis has none.
    """
    from .fileio import save_slope_data_to_xlsx

    d = read_fez(path)
    slope_data, caveats = fez_to_slope_data(d)
    save_slope_data_to_xlsx(slope_data, out_path, template=template)
    return caveats
