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

"""Rocscience Slide2 (.sli / .slim / .slmd) import.

Slide2 stores a model as sectioned ``key: value`` ASCII text. Three containers
wrap that same grammar:

  - ``.sli``  — a single, uncompressed model (all sections inline).
  - ``.slim`` — a ZIP of one ``.sli`` plus view/theme sidecars.
  - ``.slmd`` — a ZIP of a shared model (``.slsp``, carrying the material
    strengths) and one scenario file (``.slsc``) per scenario, the way a
    GeoStudio ``.gsz`` holds several analyses. A ``.slmm`` index names them.

This module reads that text and builds an xslope ``slope_data``:

  - read_sli(path) / read_slim(path) / read_slmd(path): parse a file into its
    raw sections (and, for a ``.slmd``, its list of scenarios).
  - list_scenarios(d): the scenarios a file contains, for the caller to pick.
  - slide2_to_slope_data(d, scenario): build slope_data for ONE scenario, plus
    a list of caveats for anything that did not map cleanly.
  - import_slmd(path, template, out_path, scenario): the script route — parse
    and write an xslope .xlsx input file.

Only what Slide2 and xslope have in common is imported. Everything else is
reported as a caveat rather than silently dropped: strength models other than
Mohr-Coulomb, water-pressure grids and finite-element seepage, search
definitions (Slide2 defines a search, not a surface), support/anchors,
distributed and seismic loads, and design-standard partial factors.

Format learned from the Slide2 tutorial files (public downloads, Rocscience
help site); two format generations are handled — v7.010 (.slim) and v9.027
(.slmd).
"""

import os
import re
import zipfile

from shapely.geometry import Polygon
from shapely.ops import unary_union

from .fileio import build_ground_surface_from_polygons


# Slide2 material "type" codes and whether xslope can represent them. Only the
# Mohr-Coulomb type (0) maps 1:1; every other type is imported with whatever c/phi
# it carries and flagged, the same discipline geostudio.py uses. Correlated with the
# tutorials that use each: 8 = generalized Hoek-Brown (Damage Regions), 15/16 =
# anisotropic (Generalized Anisotropic, Sarma), 18 = SHANSEP (Staged Embankment),
# 1 = Mohr-Coulomb with a tensile-strength cutoff (Tensile Strength tutorial).
_MC_TYPE = 0
_TYPE_NAMES = {
    1: "Mohr-Coulomb with a tensile-strength cutoff",
    2: "undrained (phi = 0)",
    5: "infinite strength / impenetrable",
    8: "generalized Hoek-Brown",
    15: "anisotropic (c/phi with a second orientation)",
    16: "generalized anisotropic (function-based)",
    18: "SHANSEP",
}

# Slide2 water-pressure method (model description "water:"). Only 'hu' — a
# hydrostatic pressure below a water table — maps onto xslope's piezometric line.
_WATER_HU = "hu"

# Slide2 LEM method bitmask (model description "methods:"). Bits we recognise; the
# value is informational only — xslope solves with its own solvers regardless.
_METHOD_BITS = [
    (1, "Ordinary/Fellenius"), (2, "Bishop simplified"), (4, "Janbu simplified"),
    (8, "Spencer"), (128, "GLE/Morgenstern-Price"),
]

_GAMMA_W_METRIC = 9.81


# --------------------------------------------------------------------------------------
# Text grammar
# --------------------------------------------------------------------------------------

def _parse_sections(text):
    """Split a Slide2 payload into ``{section header: [body lines]}``.

    A section header is a line at column 0 whose (right-stripped) text ends in ':'.
    Everything indented under it, up to the next header, is its body. This is the
    whole grammar — every entity below is pulled out of one section's body.
    """
    out, cur = {}, None
    for line in text.splitlines():
        if line and not line[0].isspace() and line.rstrip().endswith(":"):
            cur = line.rstrip()[:-1]
            out.setdefault(cur, [])
        elif cur is not None:
            out[cur].append(line.rstrip())
    return out


def _kv(s):
    """Parse a run of ``key: value`` tokens into a dict.

    Values may be a bracketed vertex list ``[1,2,3]``, a GUID ``{...}``, or a bare
    token; a trailing comma (Slide2 writes ``wtable: 0,``) is stripped.
    """
    d = {}
    for m in re.finditer(r"(\w+):\s*(\[[^\]]*\]|\{[^}]*\}|[^\s]+)", s):
        val = m.group(2).rstrip(",")
        d[m.group(1)] = val
    return d


def _fnum(v, default=0.0):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def _vlist(token):
    """The integer vertex ids in a ``[a,b,c]`` token."""
    if not token:
        return []
    inner = token.strip().lstrip("[").rstrip("]")
    return [int(x) for x in inner.split(",") if x.strip().lstrip("-").isdigit()]


def _parse_materials(body):
    """``material types`` body -> ``{soil name: {prop: value}}``, order preserved.

    Works for both generations: a ``.sli`` / ``.slsp`` line carries the strengths
    inline (``c``, ``phi``, ``uw``, ``type``); a ``.slsc`` line carries only the
    per-scenario flags (``water``, ``wtable``). The caller merges the two.
    """
    mats = {}
    for ln in body:
        m = re.match(r"\s*(soil\d+)\s*=\s*(.*)", ln)
        if m:
            mats[m.group(1)] = _kv(m.group(2))
    return mats


def _parse_anchor_types(body):
    types = {}
    for ln in body:
        m = re.match(r"\s*(anchor\d+)\s*=\s*(.*)", ln)
        if m:
            types[m.group(1)] = _kv(m.group(2))
    return types


# --------------------------------------------------------------------------------------
# Readers
# --------------------------------------------------------------------------------------

def _model_from_text(text, name):
    """One parsed scenario: its sections plus the entities pulled from them."""
    sections = _parse_sections(text)

    md = {}
    for ln in sections.get("model description", []):
        m = re.match(r"\s*([\w ]+):\s*(.*)", ln)
        if m:
            md[m.group(1).strip()] = m.group(2).strip()

    verts = {}
    for ln in sections.get("vertices", []):
        m = re.match(r"\s*(\d+)\s+x:\s*(\S+)\s+y:\s*(\S+)", ln)
        if m:
            verts[int(m.group(1))] = (float(m.group(2)), float(m.group(3)))

    cells = []
    for ln in sections.get("cells", []):
        m = re.match(r"\s*(\d+)\s+vertices:\s*\[([^\]]*)\]\s*material:\s*(\S+)", ln)
        if m:
            ids = [int(x) for x in m.group(2).split(",") if x.strip()]
            cells.append((ids, m.group(3)))

    return {
        "name": name,
        "sections": sections,
        "model_description": md,
        "vertices": verts,
        "cells": cells,
        "materials": _parse_materials(sections.get("material types", [])),
    }


def read_sli(path):
    """Parse a Slide2 ``.sli`` (or the ``.sli`` text) into its raw sections.

    Returns a dict with a single scenario. Material strengths are inline in this
    generation, so the scenario's own materials are also the shared ones.

    Raises ValueError if the file is not readable Slide2 text.
    """
    with open(path, "rb") as fh:
        text = fh.read().decode("latin-1")
    return _from_sli_text(text, os.path.splitext(os.path.basename(path))[0], str(path))


def _from_sli_text(text, name, path):
    if "material types:" not in text:
        raise ValueError(f"{os.path.basename(path)} is not a Slide2 model "
                         f"(no 'material types:' section).")
    scen = _model_from_text(text, name)
    return {
        "path": path, "kind": "sli",
        "version": scen["model_description"].get("version", ""),
        "shared_materials": scen["materials"],
        "anchor_types": _parse_anchor_types(
            _parse_sections(text).get("anchor types", [])),
        "scenarios": [scen],
    }


def read_slim(path):
    """Parse a Slide2 ``.slim`` — a ZIP wrapping one ``.sli`` — into its sections.

    Same shape as :func:`read_sli`; the ZIP's view/theme sidecars are ignored.
    """
    if not zipfile.is_zipfile(path):
        raise ValueError(f"{os.path.basename(path)} is not a Slide2 .slim "
                         f"(expected a ZIP archive).")
    zf = zipfile.ZipFile(path)
    sli = next((n for n in zf.namelist() if n.lower().endswith(".sli")), None)
    if sli is None:
        raise ValueError(f"{os.path.basename(path)} contains no .sli model.")
    text = zf.read(sli).decode("latin-1")
    d = _from_sli_text(text, os.path.splitext(os.path.basename(path))[0], str(path))
    d["kind"] = "slim"
    return d


def read_slmd(path):
    """Parse a Slide2 ``.slmd`` — a multi-scenario ZIP — into its scenarios.

    The shared model (``.slsp``) carries the material strengths; each scenario
    (``.slsc``) carries its own geometry and settings and references those
    materials by name. The ``.slmm`` index names the scenarios and gives their
    order; the reader merges each scenario's material flags over the shared
    strengths so downstream code sees one complete material set per scenario.

    Returns a dict with a ``scenarios`` list. Raises ValueError if the archive is
    not a readable Slide2 multi-scenario file.
    """
    if not zipfile.is_zipfile(path):
        raise ValueError(f"{os.path.basename(path)} is not a Slide2 .slmd "
                         f"(expected a ZIP archive).")
    zf = zipfile.ZipFile(path)
    names = zf.namelist()

    slsp = next((n for n in names if n.lower().endswith(".slsp")), None)
    shared_materials, anchor_types, version = {}, {}, ""
    if slsp is not None:
        stext = zf.read(slsp).decode("latin-1")
        ssec = _parse_sections(stext)
        shared_materials = _parse_materials(ssec.get("material types", []))
        anchor_types = _parse_anchor_types(ssec.get("anchor types", []))
        for ln in ssec.get("shared model description", []):
            m = re.match(r"\s*version:\s*(\S+)", ln)
            if m:
                version = m.group(1)

    # The .slmm index gives each scenario a human name and an order. Map the
    # scenario file (a GUID .slsc) to its name; fall back to the file stem.
    slmm = next((n for n in names if n.lower().endswith(".slmm")), None)
    order, sc_names = [], {}
    if slmm is not None:
        mtext = zf.read(slmm).decode("latin-1")
        cur_file = cur_name = None
        for ln in mtext.splitlines():
            s = ln.strip()
            m = re.match(r'sc_filename:\s*"([^"]+)"', s)
            if m:
                cur_file = os.path.basename(m.group(1))
            m = re.match(r'sc_name:\s*"([^"]*)"', s)
            if m:
                cur_name = m.group(1)
            if cur_file and cur_name is not None:
                sc_names[cur_file] = cur_name
                order.append(cur_file)
                cur_file = cur_name = None

    slsc_files = [n for n in names if n.lower().endswith(".slsc")]
    # Honour the .slmm order, then append any scenario the index did not list.
    ordered = [f for f in order if any(n.endswith(f) for n in slsc_files)]
    ordered += [os.path.basename(n) for n in slsc_files
                if os.path.basename(n) not in ordered]

    scenarios = []
    for base in ordered:
        full = next((n for n in slsc_files if os.path.basename(n) == base), None)
        if full is None:
            continue
        text = zf.read(full).decode("latin-1")
        scen = _model_from_text(text, sc_names.get(base, os.path.splitext(base)[0]))
        scen["file"] = base
        scenarios.append(scen)

    if not scenarios:
        raise ValueError(f"{os.path.basename(path)} defines no scenarios (.slsc).")

    return {
        "path": str(path), "kind": "slmd", "version": version,
        "shared_materials": shared_materials, "anchor_types": anchor_types,
        "scenarios": scenarios,
    }


def read_slide2(path):
    """Read any Slide2 model file, dispatching on its extension."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".slmd":
        return read_slmd(path)
    if ext == ".slim":
        return read_slim(path)
    if ext == ".sli":
        return read_sli(path)
    raise ValueError(f"{os.path.basename(path)} is not a Slide2 model "
                     f"(.sli/.slim/.slmd).")


def list_scenarios(d):
    """The scenarios in a parsed Slide2 file: ``[{'index','name'}, ...]``.

    A ``.slmd`` bundles several (a base case plus variants); a ``.sli``/``.slim``
    has exactly one. The caller passes the chosen index to
    :func:`slide2_to_slope_data`.
    """
    return [{"index": i, "name": s["name"]} for i, s in enumerate(d["scenarios"])]


# --------------------------------------------------------------------------------------
# slope_data construction
# --------------------------------------------------------------------------------------

def _blank_material(name):
    """A material with xslope's full key set, everything zeroed but the name —
    the same shape ``load_slope_data`` produces (mirrors geostudio.py)."""
    return {
        "name": str(name), "gamma": 0.0, "gamma_sat": None, "option": "", "c": 0.0,
        "phi": 0.0, "cp": 0.0, "r_elev": 0.0, "d": 0, "psi": 0, "u": "none", "ru": 0.0,
        "sigma_gamma": 0.0, "sigma_c": 0.0, "sigma_phi": 0.0, "sigma_cp": 0.0,
        "sigma_d": 0.0, "sigma_psi": 0.0, "k1": 0.0, "k2": 0.0, "alpha": 0.0,
        "unsat": "lf", "kr0": 0.0, "h0": 0.0, "vg_a": 0.0, "vg_n": 0.0,
        "E": 0.0, "nu": 0.0,
    }


def _merge_material(shared, scen):
    """Shared strengths (from the .slsp) overlaid with per-scenario flags (.slsc).

    The scenario line carries only ``water``/``wtable``/``withru`` and never the
    strengths, so updating shared with the scenario dict fills in the flags without
    clobbering c/phi/uw. For a .sli both are the same inline line, so the merge is
    a no-op.
    """
    out = dict(shared or {})
    out.update(scen or {})
    return out


def _strength(props):
    """(c, phi, gamma, type_code, caveat) for one Slide2 material.

    Only ``type: 0`` (Mohr-Coulomb) maps 1:1. Any other type is imported with
    whatever c/phi it lists and flagged — the same behaviour geostudio.py gives an
    unknown SLOPE/W strength model: never silently a zero-strength soil, never
    silently the wrong strength.
    """
    tcode = int(_fnum(props.get("type", 0)))
    c = _fnum(props.get("c", 0.0))
    phi = _fnum(props.get("phi", 0.0))
    gamma = _fnum(props.get("uw", 0.0))
    if tcode == _MC_TYPE:
        return c, phi, gamma, tcode, None
    label = _TYPE_NAMES.get(tcode, f"strength type {tcode}")
    return c, phi, gamma, tcode, (
        f"material uses Slide2's {label} model, which xslope does not have — "
        f"imported as Mohr-Coulomb with c = {c:g}, phi = {phi:g}, which is NOT the "
        f"same strength; check it before solving")


def slide2_to_slope_data(d, scenario=None):
    """Build an xslope ``slope_data`` for one scenario of a parsed Slide2 file.

    Parameters
    ----------
    d : dict
        Output of :func:`read_slmd`, :func:`read_slim` or :func:`read_sli`.
    scenario : int or str, optional
        Which scenario to import — an index into :func:`list_scenarios`, or a
        scenario name. Defaults to the first. A ``.slmd``'s scenarios can differ
        in geometry, water and search, so this choice is not cosmetic.

    Returns
    -------
    (slope_data, caveats) : (dict, list of str)
        ``caveats`` names everything Slide2 defined that xslope could not take
        faithfully. It is never silently empty when something was dropped.

    Notes
    -----
    Slide2's *search* definition (a grid, block, path or metaheuristic search) is
    never imported: none is an xslope surface, and a wrong search is worse than
    none. A *specified* circle (``centers``) or a *specified* non-circular surface
    (``surfaces``) is imported, so a model built that way arrives complete. A
    search-only model arrives without a surface and must be given one before it
    will solve.
    """
    caveats = []
    scen = _pick_scenario(d, scenario)
    sections = scen["sections"]
    md = scen["model_description"]
    verts = scen["vertices"]
    cells = scen["cells"]

    if not cells:
        raise ValueError("This scenario has no material cells — there is no geometry "
                         "to import. (It may be a template/definition-only scenario.)")

    units = md.get("units", "metric")
    gamma_water = _fnum(md.get("gammaw", _GAMMA_W_METRIC), _GAMMA_W_METRIC)

    # --- materials & zones -------------------------------------------------------
    # Union the triangular cells of each material into its zone(s). Slide2 tiles the
    # domain with triangles tagged by material; xslope wants one polygon per zone.
    # Order the materials by their Slide2 index (soil1, soil2, ...) rather than by the
    # order cells happen to appear, so mat_id lines up with Slide2's own numbering.
    def _soil_num(name):
        m = re.search(r"(\d+)", name)
        return int(m.group(1)) if m else 0
    used_names = sorted({mname for _, mname in cells}, key=_soil_num)
    mat_index = {name: i for i, name in enumerate(used_names)}

    materials = []
    for name in used_names:
        props = _merge_material(d.get("shared_materials", {}).get(name), scen["materials"].get(name))
        c, phi, gamma, tcode, note = _strength(props)
        if note:
            caveats.append(f"material '{name}': " + note)
        mat = _blank_material(name)
        mat.update(option="mc", gamma=gamma, c=c, phi=phi)
        mat["_props"] = props        # kept for the water pass, dropped before return
        materials.append(mat)
        if not gamma:
            caveats.append(f"material '{name}' has no unit weight in the file — set "
                           f"gamma before solving")
        if not c and not phi:
            caveats.append(
                f"material '{name}' came across with NO STRENGTH (c = 0, phi = 0) — set "
                f"its strength before solving or the factor of safety is meaningless")

    polygons = []
    for name in used_names:
        tris = [Polygon([verts[i] for i in ids if i in verts])
                for ids, mname in cells if mname == name]
        tris = [t for t in tris if t.is_valid and t.area > 0]
        if not tris:
            continue
        merged = unary_union(tris)
        geoms = list(merged.geoms) if merged.geom_type == "MultiPolygon" else [merged]
        for g in geoms:
            polygons.append({"polygon": Polygon(g.exterior.coords),
                             "mat_id": mat_index[name]})
    if not polygons:
        raise ValueError("This scenario has no material zones to import.")

    ground_surface, domain_polygon = build_ground_surface_from_polygons(polygons)

    # --- water -------------------------------------------------------------------
    piezo_line = []
    water_mode = md.get("water", "none")
    wt_ids = []
    for ln in sections.get("water table", []):
        m = re.search(r"vertices:\s*\[([^\]]*)\]", ln)
        if m:
            wt_ids = [int(x) for x in m.group(1).split(",") if x.strip()]
    if water_mode == _WATER_HU and wt_ids:
        piezo_line = [verts[i] for i in wt_ids if i in verts]
        for mat in materials:
            props = mat["_props"]
            # A material draws pore pressure from the water table only if wtable = 1
            # (Slide2's per-material switch); one left off it is dry even with a table.
            if int(_fnum(props.get("wtable", 0))):
                mat["u"] = "piezo"
        dry = [m["name"] for m in materials if m["u"] != "piezo"]
        if piezo_line and dry:
            caveats.append(
                f"material(s) {', '.join(repr(n) for n in dry)} are not connected to the "
                f"water table in Slide2 — imported with no pore pressure")
    elif water_mode == "groundwater":
        caveats.append(
            "THE FACTOR OF SAFETY WILL BE WRONG: pore pressure comes from a Slide2 "
            "finite-element groundwater analysis, which xslope cannot read from this "
            "file — imported as zero. Run the seepage analysis in xslope, or set a "
            "piezo line, before solving")
    elif water_mode not in ("none", "", _WATER_HU):
        caveats.append(
            f"THE FACTOR OF SAFETY MAY BE WRONG: Slide2's water method '{water_mode}' has "
            f"no xslope equivalent — pore pressure imported as zero")

    # ru: a material with withru = 1 carries a pore-pressure ratio instead.
    for mat in materials:
        props = mat["_props"]
        if int(_fnum(props.get("withru", 0))) and mat["u"] == "none":
            mat["u"] = "ru"
            mat["ru"] = _fnum(props.get("ru", 0.0))
            caveats.append(f"material '{mat['name']}' uses an ru pore-pressure ratio "
                           f"({mat['ru']:g})")

    # A water-pressure grid (Water Pressure Grid / FE seepage tutorials) is a table of
    # (x, y, pressure) points, not a mesh xslope reads — flag it if present.
    if any(sections.get("pore pressures", [])):
        pts = [ln for ln in sections["pore pressures"] if re.search(r"\bv:", ln)]
        if pts:
            caveats.append(
                f"the model carries a Slide2 water-pressure grid ({len(pts)} points), "
                f"which xslope does not import — pore pressure from it is NOT applied; "
                f"set a piezo line or a seepage solution instead")

    hu_bad = [m["name"] for m in materials
              if m["u"] == "piezo" and int(_fnum(m["_props"].get("hutype", 0))) != 0]
    if hu_bad:
        caveats.append(
            f"material(s) {', '.join(repr(n) for n in hu_bad)} use a non-default Hu "
            f"pressure type in Slide2; xslope applies full hydrostatic pressure below "
            f"the water table — the pore pressure may differ")

    for mat in materials:
        mat.pop("_props", None)

    # --- failure surface ---------------------------------------------------------
    circles, non_circ = [], []
    for ln in sections.get("centers", []):
        m = re.match(r"\s*\d+\s+x:\s*(\S+)\s+y:\s*(\S+)\s+r:\s*(\S+)", ln)
        if m:
            xo, yo, r = float(m.group(1)), float(m.group(2)), float(m.group(3))
            circles.append({"Xo": xo, "Yo": yo, "R": r, "Depth": yo - r})
    for ln in sections.get("surfaces", []):
        m = re.search(r"vertices:\s*\[([^\]]*)\]", ln)
        if m:
            ids = [int(x) for x in m.group(1).split(",") if x.strip()]
            pts = [verts[i] for i in ids if i in verts]
            if len(pts) >= 2:
                non_circ = [{"X": x, "Y": y, "Movement": None} for x, y in pts]

    if circles:
        if len(circles) > 1:
            caveats.append(f"the file specifies {len(circles)} circles; all were "
                           f"imported as trial circles")
        caveats.append("a specified failure circle was imported; run a search to find "
                       "xslope's own critical surface")
    elif non_circ:
        caveats.append(
            "a specified non-circular surface was imported; Slide2 records no slip "
            "direction on it, so xslope infers movement from the geometry — check it")
    else:
        caveats.append(
            "no failure surface was imported — this scenario defines a SEARCH (grid, "
            "block, path or metaheuristic), which has no xslope equivalent, and no "
            "specified circle or surface to take one from. Define circles or a "
            "non-circular surface before solving (an input file with no surface will "
            "not re-load)")

    # --- tension crack -----------------------------------------------------------
    tcrack_depth = tcrack_water = 0.0
    tc_ids = []
    for ln in sections.get("tension crack", []):
        m = re.search(r"vertices:\s*\[([^\]]*)\]", ln)
        if m:
            tc_ids = [int(x) for x in m.group(1).split(",") if x.strip()]
    if tc_ids:
        tc_props = _kv(" ".join(sections.get("tension crack properties", [])))
        pct = {"filled": 1.0, "wet": 1.0, "dry": 0.0}.get(tc_props.get("type", "dry"), 0.0)
        depths = []
        for i in tc_ids:
            if i in verts:
                x, y = verts[i]
                yg = _y_on(ground_surface, x)
                if yg is not None and yg - y > 0:
                    depths.append(yg - y)
        if depths:
            tcrack_depth = max(depths)
            tcrack_water = pct * tcrack_depth
            caveats.append(
                f"tension crack imported at depth {tcrack_depth:.2f}"
                + (f", {100 * pct:.0f}% full of water" if pct else "")
                + " (xslope models a crack of constant depth; the deepest point was taken)")

    # --- things reported but not imported ----------------------------------------
    k_seismic = _fnum(md.get("seismic", 0.0))
    if _fnum(md.get("seismicv", 0.0)):
        caveats.append(
            f"a vertical seismic coefficient ({_fnum(md.get('seismicv')):g}) is set in "
            f"Slide2 — xslope models only the horizontal one, so it was dropped")

    n_anchors = sum(1 for ln in sections.get("anchors", []) if re.match(r"\s*\d+\s+x1:", ln))
    if n_anchors:
        caveats.append(
            f"{n_anchors} support/anchor(s) are defined in Slide2 but were NOT imported "
            f"— Slide2's support model (end-anchored/grouted, with its own capacity, "
            f"bond and orientation rules) does not map onto xslope's reinforcement, and "
            f"guessing the force would be worse than omitting it. THE FACTOR OF SAFETY "
            f"WILL BE LOWER than Slide2's; add reinforcement by hand if needed")

    n_forces = sum(1 for ln in sections.get("forces", []) if re.search(r"\bload:", ln))
    if n_forces:
        caveats.append(
            f"{n_forces} external load(s) (distributed/line loads) are defined in Slide2 "
            f"but were NOT imported — add them by hand; the factor of safety will differ")

    if md.get("design_selection", "0") not in ("0", ""):
        caveats.append(
            "a design standard (e.g. Eurocode 7) with partial factors is selected in "
            "Slide2 — xslope applies no partial factors, so the factor of safety is an "
            "unfactored value")

    if md.get("transient", "no") not in ("no", "", "0"):
        caveats.append(
            "this is a transient (time-dependent) analysis in Slide2; xslope imports a "
            "single state — check the water condition")

    # Methods Slide2 was set to run (informational; xslope solves with its own).
    methods_val = int(_fnum(md.get("methods", 0)))
    named = [name for bit, name in _METHOD_BITS if methods_val & bit]
    if named:
        caveats.append("Slide2 was set to run: " + ", ".join(named)
                       + " — solve with the matching xslope method for comparison")

    caveats.append(
        f"this model is in {units} units — xslope is unit-agnostic and imports the "
        f"numbers as they are, so results come back in the same system")

    slope_data = {
        "template_version": None,
        "gamma_water": gamma_water,
        "tcrack_depth": tcrack_depth,
        "tcrack_water": tcrack_water,
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
        "circular": bool(circles),
        "circles": circles,
        "non_circ": non_circ,
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


def _pick_scenario(d, scenario):
    scens = d["scenarios"]
    if scenario is None:
        return scens[0]
    if isinstance(scenario, int):
        if not 0 <= scenario < len(scens):
            raise ValueError(f"This file has no scenario index {scenario} "
                             f"(it has {len(scens)}).")
        return scens[scenario]
    for s in scens:
        if s["name"] == scenario:
            return s
    raise ValueError(f"This file has no scenario named {scenario!r}.")


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


def import_slmd(path, template, out_path, scenario=None):
    """Read a Slide2 file and write it to an xslope .xlsx input file.

    The script-level route, mirroring :func:`xslope.geostudio.import_gsz`: parse
    the Slide2 model, convert one scenario to ``slope_data``, and write it through
    a copy of an input template. Works for ``.slmd``, ``.slim`` and ``.sli``.

    Returns the caveat list from :func:`slide2_to_slope_data`. Read it — a clean
    return does not mean a complete model (a search-only scenario imports no
    failure surface).
    """
    from .fileio import save_slope_data_to_xlsx

    d = read_slide2(path)
    slope_data, caveats = slide2_to_slope_data(d, scenario)
    save_slope_data_to_xlsx(slope_data, out_path, template=template)
    return caveats
