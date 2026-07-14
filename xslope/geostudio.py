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

"""GeoStudio SLOPE/W (.gsz) import.

A .gsz is a ZIP holding the model as a single XML document (``GSIData`` root),
plus a mesh and — if the file was saved after solving — SLOPE/W's own results.
This module reads that XML and builds an xslope ``slope_data``:

  - read_gsz(path): parse the file into its raw entities (all analyses).
  - list_analyses(gsz): the analyses it contains, for the caller to choose from.
  - gsz_to_slope_data(gsz, analysis_id): build slope_data for ONE analysis,
    plus a list of caveats for anything that did not map cleanly.
  - import_gsz(path, template, out_path): the script route — parse and write an
    xslope .xlsx input file.
  - read_gsz_results(gsz, analysis_name): SLOPE/W's solved trial surfaces, when
    the file carries them. Not part of the model; used to compare answers.

Unlike DXF, a .gsz is *semantically complete* — regions, materials and water
conditions already know what they are, so no layer-mapping wizard is needed. The
one thing the caller must choose is which analysis to import, because a single
.gsz commonly holds several (and they can differ in materials as well as in
slip-surface definition — see ``gsz_to_slope_data``).

Only what SLOPE/W and xslope have in common is imported. Everything else is
reported as a caveat rather than silently dropped: unsupported strength models,
unsupported pore-pressure options, reinforcement, and distributed loads.
"""

import csv
import io
import os
import re
import zipfile
import xml.etree.ElementTree as ET

from shapely.geometry import Polygon

from .fileio import build_ground_surface_from_polygons


# Strength models SLOPE/W offers that xslope can represent. Anything else is
# imported as a placeholder and flagged.
_SUPPORTED_MODELS = {"MohrCoulomb"}

# Unit weight of water, used to infer the file's unit system (see _detect_units).
_GAMMA_W_METRIC = 9.807      # kN/m3
_GAMMA_W_IMPERIAL = 62.4     # lb/ft3


def _detect_units(gamma_water):
    """Infer the unit system from the unit weight of water.

    A .gsz carries no explicit unit-system field — the numbers are just numbers,
    and the same file format holds both metric and imperial models. The unit
    weight of water is the one quantity whose value is fixed by physics, so it
    identifies the system. Refuse to guess on anything else: silently assuming
    the wrong system would rescale an entire model.
    """
    if abs(gamma_water - _GAMMA_W_METRIC) < 0.15:
        return "metric"                       # kN/m3, m, kPa
    if abs(gamma_water - _GAMMA_W_IMPERIAL) < 1.0:
        return "imperial"                     # lb/ft3, ft, psf
    raise ValueError(
        f"Cannot determine the unit system of this file: it gives the unit weight "
        f"of water as {gamma_water:g}, which is neither {_GAMMA_W_METRIC:g} (kN/m3) "
        f"nor {_GAMMA_W_IMPERIAL:g} (lb/ft3)."
    )


def _text(node, tag, default=None):
    """Child element's text, treating GeoStudio's ``Missing="true"`` as absent."""
    if node is None:
        return default
    el = node.find(tag)
    if el is None or el.get("Missing") == "true" or el.text is None:
        return default
    return el.text.strip()


def _num(node, tag, default=None):
    v = _text(node, tag)
    if v is None or v == "":
        return default
    try:
        return float(v)
    except ValueError:
        return default


def read_gsz(path):
    """Parse a .gsz into its raw entities, without interpreting them.

    Returns a dict with the file's points, regions, materials, per-analysis
    material assignments, water conditions, piezometric surfaces and analyses.
    Coordinates are as-authored; no unit conversion is applied here.

    Raises ValueError if the file is not a readable GeoStudio project.
    """
    if not zipfile.is_zipfile(path):
        raise ValueError(f"{os.path.basename(path)} is not a GeoStudio .gsz "
                         f"(expected a ZIP archive).")
    zf = zipfile.ZipFile(path)
    names = zf.namelist()
    # The model lives in the one XML at the archive root; the per-analysis copies
    # in subfolders are snapshots written alongside results.
    root_xml = next((n for n in names if "/" not in n and n.lower().endswith(".xml")), None)
    if root_xml is None:
        raise ValueError(f"{os.path.basename(path)} contains no GeoStudio model XML.")
    try:
        root = ET.fromstring(zf.read(root_xml))
    except ET.ParseError as exc:
        raise ValueError(f"Could not parse the model XML in "
                         f"{os.path.basename(path)}: {exc}") from exc
    if root.tag != "GSIData":
        raise ValueError(f"{os.path.basename(path)} is not a GeoStudio project "
                         f"(XML root is <{root.tag}>, expected <GSIData>).")

    analyses = []
    for a in root.findall("./Analyses/Analysis"):
        aid = _text(a, "ID")
        if aid is None:
            continue
        analyses.append({
            "id": int(aid),
            "name": _text(a, "Name", f"Analysis {aid}"),
            "kind": _text(a, "Kind", ""),
            "method": _text(a, "Method", ""),
        })

    points, regions = {}, {}
    for g in root.findall("./Geometries/Geometry"):
        for p in g.findall("./Points/Point"):
            points[int(p.get("ID"))] = (float(p.get("X")), float(p.get("Y")))
        for r in g.findall("./Regions/Region"):
            ids = _text(r, "PointIDs")
            if ids:
                regions[int(_text(r, "ID"))] = [int(i) for i in ids.split(",")]

    materials = {}
    for m in root.findall("./Materials/Material"):
        ss = m.find("StressStrain")
        materials[int(_text(m, "ID"))] = {
            "name": _text(m, "Name", "material"),
            "model": _text(m, "SlopeModel", ""),
            "gamma": _num(ss, "UnitWeight", 0.0),
            "c": _num(ss, "CohesionPrime", 0.0),
            "phi": _num(ss, "PhiPrime", 0.0),
        }

    # Region -> material is defined PER ANALYSIS, in <Contexts>, not on the region.
    # The same geometry can therefore carry different soils in different analyses.
    contexts = {}
    for c in root.findall("./Contexts/Context"):
        aid = int(_text(c, "AnalysisID"))
        assign = {}
        for gum in c.findall("./GeometryUsesMaterials/GeometryUsesMaterial"):
            m = re.match(r"Regions-(\d+)$", gum.get("ID", ""))
            if m and gum.get("Entry"):
                assign[int(m.group(1))] = int(gum.get("Entry"))
        contexts[aid] = assign

    water = {}
    for w in root.findall("./WaterItems/WaterItem"):
        entry = w.find("Entry")
        water[int(_text(w, "AnalysisID"))] = {
            "gamma_water": _num(entry, "UnitWaterWeight", _GAMMA_W_METRIC),
            "option": _text(entry.find("ResultInputInfo"), "Option", "None"),
        }

    # Piezometric surfaces are polylines through the shared Points table — they
    # store point IDs, not coordinates.
    piezo = {}
    for ps in root.iter("PiezometricSurface"):
        pid = _text(ps, "ID")
        if pid is None:
            continue
        piezo[int(pid)] = [int(dp.text) for dp in ps.findall("./DataPoints/DataPoint")
                           if dp.text]

    stability = {}
    for s in root.findall("./StabilityItems/StabilityItem"):
        entry = s.find("Entry")
        seis = entry.find("Seismic") if entry is not None else None
        kh = seis.get("Horizontal") if seis is not None else ""
        stability[int(_text(s, "AnalysisID"))] = {
            "k_seismic": float(kh) if kh else 0.0,
        }

    return {
        "path": str(path), "zip": zf, "analyses": analyses, "points": points,
        "regions": regions, "materials": materials, "contexts": contexts,
        "water": water, "piezo": piezo, "stability": stability,
    }


def list_analyses(gsz):
    """The analyses in a parsed .gsz: ``[{'id','name','kind','method'}, ...]``.

    A .gsz routinely bundles several (a base case plus variants). They can differ
    in materials, not just in slip surface, so the caller must pick one.
    """
    return list(gsz["analyses"])


def _blank_material(name):
    """A material with xslope's full key set, everything zeroed but the name —
    the same shape ``load_slope_data`` produces."""
    return {
        "name": str(name), "gamma": 0.0, "gamma_sat": None, "option": "", "c": 0.0,
        "phi": 0.0, "cp": 0.0, "r_elev": 0.0, "d": 0, "psi": 0, "u": "none", "ru": 0.0,
        "sigma_gamma": 0.0, "sigma_c": 0.0, "sigma_phi": 0.0, "sigma_cp": 0.0,
        "sigma_d": 0.0, "sigma_psi": 0.0, "k1": 0.0, "k2": 0.0, "alpha": 0.0,
        "unsat": "lf", "kr0": 0.0, "h0": 0.0, "vg_a": 0.0, "vg_n": 0.0,
        "E": 0.0, "nu": 0.0,
    }


def gsz_to_slope_data(gsz, analysis_id=None, critical_surface=True):
    """Build an xslope ``slope_data`` for one analysis of a parsed .gsz.

    Parameters
    ----------
    gsz : dict
        Output of :func:`read_gsz`.
    analysis_id : int, optional
        Which analysis to import. Defaults to the first one. A .gsz's analyses
        can differ in *materials* as well as in slip surface, so this choice is
        not cosmetic — see :func:`list_analyses`.
    critical_surface : bool
        If the file was saved *solved*, import SLOPE/W's critical circle as the
        single trial circle (default). This is SLOPE/W's answer, not a search
        definition — it is imported because it makes the model complete and lets
        the two programs be compared on identical geometry. Set False to import
        no surface and define your own search.

    Returns
    -------
    (slope_data, caveats) : (dict, list of str)
        ``caveats`` names everything SLOPE/W defined that xslope could not take
        faithfully. It is never silently empty when something was dropped.

    Notes
    -----
    SLOPE/W's *search* definition (an entry/exit range or a grid-and-radius) is
    never imported: neither is an xslope search, and a wrong search is worse than
    none. If the file carries no results and ``critical_surface`` finds nothing,
    the model arrives without a failure surface and you must define one before it
    will solve — or reload, since an input file with no surface does not validate.
    """
    caveats = []
    pts = gsz["points"]

    if not gsz["analyses"]:
        raise ValueError("This .gsz defines no analyses.")
    if analysis_id is None:
        analysis_id = gsz["analyses"][0]["id"]
    if not any(a["id"] == analysis_id for a in gsz["analyses"]):
        raise ValueError(f"This .gsz has no analysis with ID {analysis_id}.")

    assign = gsz["contexts"].get(analysis_id, {})
    if not assign:
        raise ValueError(
            "This analysis assigns no materials to its regions — there is no "
            "model to import. (It may be a parent/definition-only analysis; try "
            "one of the others in the file.)")

    water = gsz["water"].get(analysis_id, {})
    gamma_water = water.get("gamma_water", _GAMMA_W_METRIC)
    units = _detect_units(gamma_water)

    # Import only the materials this analysis actually uses, renumbered 0-based.
    used = sorted(set(assign.values()))
    mat_index = {mid: i for i, mid in enumerate(used)}
    materials = []
    for mid in used:
        src = gsz["materials"].get(mid)
        if src is None:
            raise ValueError(f"This analysis references material {mid}, which the "
                             f"file does not define.")
        mat = _blank_material(src["name"])
        if src["model"] not in _SUPPORTED_MODELS:
            caveats.append(
                f"material '{src['name']}' uses SLOPE/W's {src['model'] or 'unnamed'} "
                f"strength model, which xslope does not have — imported as "
                f"Mohr-Coulomb with its c and phi; check it before solving")
        mat.update(option="mc", gamma=src["gamma"], c=src["c"], phi=src["phi"])
        if not mat["gamma"]:
            caveats.append(f"material '{src['name']}' has no unit weight in the "
                           f"file — set gamma before solving")
        materials.append(mat)

    polygons = []
    for region_id, mid in sorted(assign.items()):
        ring = gsz["regions"].get(region_id)
        if not ring:
            continue
        polygons.append({"polygon": Polygon([pts[i] for i in ring]),
                         "mat_id": mat_index[mid]})
    if not polygons:
        raise ValueError("This analysis has no material regions to import.")

    # Pore pressure. SLOPE/W names the option on the analysis; only a piezometric
    # surface maps onto xslope's piezo line.
    piezo_line = []
    option = water.get("option") or "None"
    if option == "PiezoSurface":
        if gsz["piezo"]:
            ids = next(iter(gsz["piezo"].values()))
            piezo_line = [pts[i] for i in ids if i in pts]
            for mat in materials:
                mat["u"] = "piezo"
            if len(gsz["piezo"]) > 1:
                caveats.append(
                    f"the file defines {len(gsz['piezo'])} piezometric surfaces; "
                    f"xslope takes one, so the first was imported")
        else:
            caveats.append("the analysis asks for a piezometric surface but the file "
                           "defines none — pore pressure imported as zero")
    elif option not in ("None", "", None):
        caveats.append(
            f"pore pressure comes from SLOPE/W's '{option}' option, which xslope "
            f"cannot read — imported as zero; set a piezo line, ru, or a seepage "
            f"solution before solving")

    ground_surface, domain_polygon = build_ground_surface_from_polygons(polygons)

    slope_data = {
        "template_version": None,
        "gamma_water": gamma_water,
        "tcrack_depth": 0.0,
        "tcrack_water": 0.0,
        "k_seismic": gsz["stability"].get(analysis_id, {}).get("k_seismic", 0.0),
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

    # Failure surface. SLOPE/W's search definition has no xslope equivalent, but a
    # solved file records every trial surface it evaluated — take the critical one,
    # so the model arrives complete and directly comparable.
    analysis = next(a for a in gsz["analyses"] if a["id"] == analysis_id)
    trials = read_gsz_results(gsz, analysis["name"]) if critical_surface else []
    if trials:
        best = min(trials, key=lambda t: t["fs"])
        slope_data["circles"] = [{"Xo": best["xo"], "Yo": best["yo"], "R": best["r"],
                                  "Depth": best["yo"] - best["r"]}]
        slope_data["circular"] = True
        caveats.append(
            f"the failure surface is SLOPE/W's own critical circle (its FS = "
            f"{best['fs']:.3f} by {analysis['method'] or 'its method'}), not a search "
            f"— it is imported so the model is complete and the two programs can be "
            f"compared on the same circle; run a search to find xslope's own critical "
            f"surface")
    else:
        caveats.append(
            "no failure surface was imported — SLOPE/W's search definition has no "
            "xslope equivalent, and this file carries no solved results to take a "
            "critical circle from; define circles or a non-circular surface before "
            "solving (an input file with no surface will not re-load)")

    caveats.append(f"model read as {units} units "
                   f"({'kN/m3, m, kPa' if units == 'metric' else 'lb/ft3, ft, psf'}) "
                   f"— confirm the units match your template")
    return slope_data, caveats


def read_gsz_results(gsz, analysis_name):
    """SLOPE/W's own solved trial slip surfaces for an analysis, if the file has them.

    Returns ``[{'fs','xo','yo','r','weight'}, ...]`` — one entry per trial surface
    SLOPE/W evaluated, with the factor of safety it computed. Empty if the file was
    saved unsolved.

    This is not part of the model. It exists so the same circles can be re-solved
    in xslope and the two programs' answers compared directly, with no difference
    in search to explain away.
    """
    zf = gsz["zip"]
    folder = analysis_name.lower() + "/"
    name = next((n for n in zf.namelist()
                 if n.lower().startswith(folder) and n.endswith("slip_surface.csv")), None)
    if name is None:
        return []

    out = []
    text = zf.read(name).decode("utf-8-sig")
    for row in csv.DictReader(io.StringIO(text)):
        try:
            fs = float(row["SlipFOS"])
            xo, yo = float(row["SlipCenterX"]), float(row["SlipCenterY"])
            radius = float(row["SlipRadiusX"])
        except (KeyError, TypeError, ValueError):
            continue                     # a trial SLOPE/W could not evaluate
        if not 0.0 < fs < 100.0:
            continue                     # SLOPE/W writes a sentinel FOS for failed trials
        try:
            weight = float(row.get("SlipWeight") or 0.0)
        except ValueError:
            weight = 0.0
        out.append({"fs": fs, "xo": xo, "yo": yo, "r": radius, "weight": weight})
    return out


def _xml_escape(text):
    return (str(text).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


def export_gsz(slope_data, gsz_path, analysis_name="xslope", method="Morgenstern-Price"):
    """Write a model to a GeoStudio SLOPE/W project (.gsz).

    The inverse of :func:`gsz_to_slope_data`, and the counterpart of
    :func:`xslope.cad.export_dxf`: material zones become regions, materials become
    Mohr-Coulomb materials, and a piezo line becomes a piezometric surface. The
    result is a single-analysis GeoStudio project.

    Parameters
    ----------
    slope_data : dict
        The model to write. Must be polygon-based (material zones); profile-line
        models have no region structure to map onto and are rejected.
    gsz_path : str
        Output .gsz path.
    analysis_name : str
        Name for the single analysis written into the file.
    method : str
        SLOPE/W method name recorded on the analysis.

    Returns
    -------
    list of str
        Caveats — what did not survive the trip. Read them: a .gsz cannot carry
        everything xslope models.

    Notes
    -----
    What is NOT written: failure surfaces (SLOPE/W defines a search, not a circle
    list), reinforcement, piles, distributed and line loads, seepage meshes, and
    non-Mohr-Coulomb strengths (c/phi models, power curves, Hoek-Brown). Each is
    reported rather than dropped in silence.

    Verified to round-trip through this module's own reader. Whether GeoStudio
    itself opens the file has not been verified against the application — check
    before relying on it.
    """
    caveats = []
    polygons = slope_data.get("polygons") or []
    materials = slope_data.get("materials") or []
    if not polygons:
        raise ValueError(
            "Only polygon-based models can be exported to .gsz — this model has no "
            "material zones. (A profile-line model has no regions to map onto; build "
            "the zones first.)")
    if not materials:
        raise ValueError("This model has no materials to export.")

    gamma_water = float(slope_data.get("gamma_water") or _GAMMA_W_METRIC)
    try:
        _detect_units(gamma_water)
    except ValueError:
        caveats.append(
            f"the model's unit weight of water ({gamma_water:g}) matches neither "
            f"metric nor imperial — GeoStudio may read the units wrongly")

    # One shared point table, as GeoStudio expects: zones index into it, and so does
    # the piezometric surface.
    points, point_id = {}, {}
    def _point(xy):
        key = (round(xy[0], 9), round(xy[1], 9))
        if key not in point_id:
            point_id[key] = len(point_id) + 1
            points[point_id[key]] = key
        return point_id[key]

    regions = []
    for poly in polygons:
        ring = list(poly["polygon"].exterior.coords)[:-1]      # drop the repeat
        regions.append({"ids": [_point(xy) for xy in ring],
                        "mat": int(poly["mat_id"]) + 1})       # GeoStudio is 1-based

    piezo_ids = []
    piezo_line = slope_data.get("piezo_line") or []
    uses_piezo = any((m.get("u") or "none") == "piezo" for m in materials)
    if piezo_line and uses_piezo:
        piezo_ids = [_point(xy) for xy in piezo_line]
    elif piezo_line:
        caveats.append("a piezo line is defined but no material uses it — written as "
                       "a piezometric surface anyway")
        piezo_ids = [_point(xy) for xy in piezo_line]

    other_u = sorted({(m.get("u") or "none") for m in materials} - {"none", "piezo"})
    if other_u:
        caveats.append(
            f"pore pressure from {', '.join(repr(u) for u in other_u)} is not written "
            f"— GeoStudio will open the model with no pore pressure from that source")

    non_mc = sorted({(m.get("option") or "") for m in materials} - {"mc", ""})
    if non_mc:
        caveats.append(
            f"strength model(s) {', '.join(repr(o) for o in non_mc)} have no SLOPE/W "
            f"equivalent — those materials are written with their c and phi, which for "
            f"a non-Mohr-Coulomb model is not the same strength")

    for key, label in [("circles", "failure circles"), ("non_circ", "non-circular surface"),
                       ("reinforce_lines", "reinforcement"), ("pile_lines", "piles"),
                       ("dloads", "distributed loads"), ("line_loads", "line loads")]:
        if slope_data.get(key):
            caveats.append(f"{label} not written — no .gsz mapping; re-create in GeoStudio")

    lines = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
             '<GSIData Version="11.11" AppVersion="25.2.1.4">',
             f'  <FileInfo Title="{_xml_escape(analysis_name)}" '
             f'Comments="Exported by xslope." />',
             '  <Analyses Len="1">',
             '    <Analysis><ID>1</ID>'
             f'<Name>{_xml_escape(analysis_name)}</Name>'
             '<Kind>SLOPE/W</Kind>'
             f'<Method>{_xml_escape(method)}</Method>'
             '<GeometryId>1</GeometryId></Analysis>',
             '  </Analyses>',
             '  <Geometries Len="1">',
             '    <Geometry><Name>2D Geometry</Name>',
             f'      <Points Len="{len(points)}">']
    for pid in sorted(points):
        x, y = points[pid]
        lines.append(f'        <Point ID="{pid}" X="{x:.10g}" Y="{y:.10g}" Pinned="true" />')
    lines.append('      </Points>')
    lines.append(f'      <Regions Len="{len(regions)}">')
    for i, r in enumerate(regions, start=1):
        lines.append(f'        <Region><ID>{i}</ID>'
                     f'<PointIDs>{",".join(str(j) for j in r["ids"])}</PointIDs></Region>')
    lines.append('      </Regions>')
    if piezo_ids:
        lines.append('      <PiezometricSurfaces Len="1">')
        lines.append('        <PiezometricSurface><ID>1</ID>'
                     f'<DataPoints Len="{len(piezo_ids)}">'
                     + "".join(f'<DataPoint>{i}</DataPoint>' for i in piezo_ids)
                     + '</DataPoints></PiezometricSurface>')
        lines.append('      </PiezometricSurfaces>')
    lines.append('    </Geometry>')
    lines.append('  </Geometries>')

    lines.append(f'  <Materials Len="{len(materials)}">')
    for i, m in enumerate(materials, start=1):
        lines.append(
            f'    <Material><ID>{i}</ID><Name>{_xml_escape(m.get("name") or f"material {i}")}</Name>'
            f'<SlopeModel>MohrCoulomb</SlopeModel><StressStrain>'
            f'<UnitWeight>{float(m.get("gamma") or 0.0):.10g}</UnitWeight>'
            f'<CohesionPrime>{float(m.get("c") or 0.0):.10g}</CohesionPrime>'
            f'<PhiPrime>{float(m.get("phi") or 0.0):.10g}</PhiPrime>'
            f'</StressStrain></Material>')
    lines.append('  </Materials>')

    # Region -> material, per analysis (GeoStudio keeps it here, not on the region).
    lines.append('  <Contexts Len="1">')
    lines.append('    <Context><AnalysisID>1</AnalysisID>')
    lines.append(f'      <GeometryUsesMaterials Len="{len(regions)}">')
    for i, r in enumerate(regions, start=1):
        lines.append(f'        <GeometryUsesMaterial ID="Regions-{i}" Entry="{r["mat"]}" />')
    lines.append('      </GeometryUsesMaterials>')
    lines.append('      <IsDefined>true</IsDefined></Context>')
    lines.append('  </Contexts>')

    lines.append('  <WaterItems Len="1">')
    lines.append('    <WaterItem><AnalysisID>1</AnalysisID><Entry>'
                 + ('<ResultInputInfo><Option>PiezoSurface</Option></ResultInputInfo>'
                    if piezo_ids else '')
                 + f'<UnitWaterWeight>{gamma_water:.10g}</UnitWaterWeight>'
                 '</Entry></WaterItem>')
    lines.append('  </WaterItems>')

    k = float(slope_data.get("k_seismic") or 0.0)
    lines.append('  <StabilityItems Len="1">')
    lines.append('    <StabilityItem><AnalysisID>1</AnalysisID><Entry>'
                 f'<Seismic Horizontal="{k:.10g}" Vertical="" />' if k else
                 '    <StabilityItem><AnalysisID>1</AnalysisID><Entry>'
                 '<Seismic Horizontal="" Vertical="" />')
    lines.append('    </Entry></StabilityItem>')
    lines.append('  </StabilityItems>')
    lines.append('</GSIData>')

    stem = os.path.splitext(os.path.basename(gsz_path))[0]
    with zipfile.ZipFile(gsz_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(f"{stem}.xml", "\n".join(lines) + "\n")
    return caveats


def import_gsz(gsz_path, template, out_path, analysis_id=None):
    """Read a .gsz and write it to an xslope .xlsx input file.

    The script-level route, mirroring :func:`xslope.cad.import_dxf`: parse the
    GeoStudio model, convert one analysis to ``slope_data``, and write it through
    a copy of an input template.

    Returns the caveat list from :func:`gsz_to_slope_data`. Read it — a clean
    return does not mean a complete model (no failure surface is ever imported).
    """
    from .fileio import save_slope_data_to_xlsx

    gsz = read_gsz(gsz_path)
    slope_data, caveats = gsz_to_slope_data(gsz, analysis_id)
    save_slope_data_to_xlsx(slope_data, out_path, template=template)
    return caveats
