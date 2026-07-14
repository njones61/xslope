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


def _detect_units(gamma_water, declared=None):
    """The file's unit system.

    GeoStudio states it outright, in ``<Coordinates><EngCoords UnitSystem="...">``,
    so use that when it is there. Fall back to inferring it from the unit weight of
    water — the one quantity whose value is fixed by physics — for files that omit
    the attribute, and refuse to guess on anything else: silently assuming the wrong
    system would rescale an entire model.
    """
    if declared:
        d = declared.strip().lower()
        if d == "metric":
            return "metric"                   # kN/m3, m, kPa
        if d == "imperial":
            return "imperial"                 # lb/ft3, ft, psf
        raise ValueError(f"Unrecognized GeoStudio unit system {declared!r}.")
    if abs(gamma_water - _GAMMA_W_METRIC) < 0.15:
        return "metric"
    if abs(gamma_water - _GAMMA_W_IMPERIAL) < 1.0:
        return "imperial"
    raise ValueError(
        f"Cannot determine the unit system of this file: it declares none, and gives "
        f"the unit weight of water as {gamma_water:g}, which is neither "
        f"{_GAMMA_W_METRIC:g} (kN/m3) nor {_GAMMA_W_IMPERIAL:g} (lb/ft3)."
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

    points, regions, glines = {}, {}, {}
    for g in root.findall("./Geometries/Geometry"):
        for p in g.findall("./Points/Point"):
            points[int(p.get("ID"))] = (float(p.get("X")), float(p.get("Y")))
        for r in g.findall("./Regions/Region"):
            ids = _text(r, "PointIDs")
            if ids:
                regions[int(_text(r, "ID"))] = [int(i) for i in ids.split(",")]
        for ln in g.findall("./Lines/Line"):
            p1, p2 = _text(ln, "PointID1"), _text(ln, "PointID2")
            if p1 and p2:
                glines[int(_text(ln, "ID"))] = (int(p1), int(p2))

    materials = {}
    for m in root.findall("./Materials/Material"):
        ss = m.find("StressStrain")
        hyd = m.find("Hydraulic")
        # Keep the whole StressStrain block. WHICH field holds the strength depends on
        # the model: a MohrCoulomb material uses CohesionPrime/PhiPrime, but an
        # UndrainedPhiZero one puts its undrained shear strength in <Cohesion> and has
        # no PhiPrime at all. Reading only the drained fields silently yields a
        # zero-strength soil. The model-specific mapping happens in gsz_to_slope_data.
        strength = {c.tag: c.text for c in (ss if ss is not None else [])
                    if c.get("Missing") != "true" and c.text}
        materials[int(_text(m, "ID"))] = {
            "name": _text(m, "Name", "material"),
            "model": _text(m, "SlopeModel", ""),
            "color": _text(m, "Color", ""),
            "gamma": _num(ss, "UnitWeight", 0.0),
            "strength": strength,
            # SEEP/W side. The hydraulic model is a reference into <Functions>, not
            # inline: <Hydraulic KFnNum="1" VolWCFnNum="1"/>.
            "seep_model": _text(m, "SeepModel", ""),
            "kfn": int(hyd.get("KFnNum")) if (hyd is not None and hyd.get("KFnNum"))
                   else None,
        }

    # Conductivity functions: a shared library keyed by ID, each a table of
    # (matric suction, k) points. See _fit_vg for how they become xslope's model.
    kfns = {}
    for kfn in root.findall("./Functions/Material/Hydraulic/KFns/KFn"):
        kid = _text(kfn, "ID")
        if kid is None:
            continue
        pts = [(float(p.get("X")), float(p.get("Y")))
               for p in kfn.findall("./Points/Point")]
        if pts:
            kfns[int(kid)] = {"name": _text(kfn, "Name", ""), "points": pts}

    # Boundary-condition library. The hydraulic ones encode their type in an
    # expression string, e.g. HydHeadvsTime(Variability=Constant,Value=8) or
    # HydTotalFlux(Variability=Constant,Value=0,Review=true).
    bcs = {}
    for bc in root.findall("./BCs/BC"):
        expr = bc.get("Hydraulic")
        if not expr:
            continue                       # displacement/air/contaminant BC — not ours
        bcs[int(bc.get("ID"))] = {"name": bc.get("Name", ""), "expr": expr}

    # Region -> material is defined PER ANALYSIS, in <Contexts>, not on the region.
    # The same geometry can therefore carry different soils in different analyses.
    # Hydraulic BCs are attached the same way, keyed by geometry item.
    contexts, hyd_bcs = {}, {}
    for c in root.findall("./Contexts/Context"):
        aid = int(_text(c, "AnalysisID"))
        assign = {}
        for gum in c.findall("./GeometryUsesMaterials/GeometryUsesMaterial"):
            m = re.match(r"Regions-(\d+)$", gum.get("ID", ""))
            if m and gum.get("Entry"):
                assign[int(m.group(1))] = int(gum.get("Entry"))
        contexts[aid] = assign
        applied = []
        for g in c.findall("./GeometryUsesHydraulicBCs/GeometryUsesHydraulicBC"):
            gid, entry = g.get("ID", ""), g.get("Entry")
            if entry:
                applied.append((gid, int(entry)))
        hyd_bcs[aid] = applied

    water = {}
    for w in root.findall("./WaterItems/WaterItem"):
        entry = w.find("Entry")
        water[int(_text(w, "AnalysisID"))] = {
            "gamma_water": _num(entry, "UnitWaterWeight", _GAMMA_W_METRIC),
            "option": _text(entry.find("ResultInputInfo"), "Option", "None"),
        }

    # Piezometric surfaces live on the StabilityItem, not the geometry, and are
    # polylines through the shared Points table — they store point IDs, not
    # coordinates. <MaterialUsesPiezs> says which materials actually draw pore
    # pressure from which surface.
    piezo = {}
    for ps in root.iter("PiezometricSurface"):
        pid = _text(ps, "ID")
        if pid is None:
            continue
        piezo[int(pid)] = [int(dp.text) for dp in ps.findall("./DataPoints/DataPoint")
                           if dp.text]

    # Reinforcement product library (top level, shared across analyses).
    reinforcements = {}
    for rf in root.findall("./Reinforcements/Reinforcement"):
        rid = _text(rf, "ID")
        if rid is None:
            continue
        reinforcements[int(rid)] = {
            "name": _text(rf, "Name", ""),
            "type": _text(rf, "Type", ""),                 # Nail | Geosynthetic | Anchor | ...
            "tensile": _num(rf, "Tensile", 0.0),
            "plate": _num(rf, "PlateCapacity", 0.0),
            "pullout": _num(rf, "PulloutResistance", 0.0),
            "spacing": _num(rf, "Spacing", None),          # out-of-plane, None = per metre
            "distribution": _text(rf, "ForceDistribution", ""),
        }

    stability = {}
    for s in root.findall("./StabilityItems/StabilityItem"):
        entry = s.find("Entry")
        if entry is None:
            continue
        seis = entry.find("Seismic")
        kh = seis.get("Horizontal") if seis is not None else ""
        kv = seis.get("Vertical") if seis is not None else ""
        uses = {}
        for mu in entry.findall("./MaterialUsesPiezs/MaterialUsesPiez"):
            if mu.get("ID") and mu.get("UsesID"):
                uses[int(mu.get("ID"))] = int(mu.get("UsesID"))

        # Entry/DataPoints is a LOCAL point list for this analysis, numbered from 1.
        # The tension crack, surcharges and reinforcement lines all index into it —
        # NOT into the geometry's point table.
        local = {}
        for d in entry.findall("./DataPoints/DataPoint"):
            if d.get("Number") and d.get("X") is not None:
                local[int(d.get("Number"))] = (float(d.get("X")), float(d.get("Y")))

        tc = entry.find("TensionCrack")
        tcrack = None
        if tc is not None:
            ids = [int(d.text) for d in tc.findall("./DataPoints/DataPoint") if d.text]
            if ids:                       # a bare <TensionCrack><UnitWaterWeight/> is not a crack
                tcrack = {
                    "option": _text(tc, "TensionOption", ""),
                    "pct_water": _num(tc, "PctFilledWithWater", 0.0),
                    "angle": _num(tc, "Angle", None),
                    "coords": [local[i] for i in ids if i in local],
                }

        surcharges = []
        for sc in entry.findall("./Surcharges/Surcharge"):
            ids = [int(d.text) for d in sc.findall("./DataPoints/DataPoint") if d.text]
            coords = [local[i] for i in ids if i in local]
            if len(coords) >= 2:
                surcharges.append({"pressure": _num(sc, "Pressure", 0.0),
                                   "coords": coords})

        reinf_lines = []
        for rl in entry.findall("./ReinforcementLines/ReinforcementLine"):
            p1, p2 = rl.get("Point1Id"), rl.get("Point2Id")
            if not (p1 and p2):
                continue
            a, b = local.get(int(p1)), local.get(int(p2))
            if a and b:
                reinf_lines.append({"reinforcement": int(rl.get("Reinforcement") or 0),
                                    "p1": a, "p2": b,
                                    "on_ground": rl.get("OnGroundSurface") == "true"})

        stability[int(_text(s, "AnalysisID"))] = {
            "k_seismic": float(kh) if kh else 0.0,
            "k_seismic_v": float(kv) if kv else 0.0,
            "material_piezo": uses,          # material ID -> piezometric surface ID
            "points": local,
            "tcrack": tcrack,
            "surcharges": surcharges,
            "reinf_lines": reinf_lines,
            # Everything GeoStudio put on this analysis, so gsz_to_slope_data can
            # report anything it does not import instead of dropping it in silence.
            "elements": [c.tag for c in entry],
        }

    # GeoStudio declares its unit system rather than leaving it to be inferred.
    eng = root.find("./Coordinates/EngCoords")
    unit_system = eng.get("UnitSystem") if eng is not None else None

    return {
        "path": str(path), "zip": zf, "analyses": analyses, "points": points,
        "regions": regions, "lines": glines, "materials": materials,
        "contexts": contexts, "water": water, "piezo": piezo,
        "stability": stability, "kfns": kfns, "bcs": bcs, "hyd_bcs": hyd_bcs,
        "reinforcements": reinforcements, "unit_system": unit_system,
    }


# Every element GeoStudio can hang on a stability analysis, and what we do with it.
# The value is None when the element carries nothing xslope needs (a solver setting, or
# a list the other entries index into); otherwise it is the sentence the user is told
# when we cannot take it. THE POINT of the table is that it is exhaustive: an element
# that is not in it at all is reported as unrecognized rather than ignored, so a
# GeoStudio feature we have never seen can never silently change someone's answer.
_ENTRY_ELEMENTS = {
    # --- imported ---
    "Seismic": None,
    "PiezometricSurfaces": None,
    "MaterialUsesPiezs": None,
    "TensionCrack": None,
    "Surcharges": None,
    # --- carries nothing xslope needs ---
    "DataPoints": None,            # the local point list the others index into
    "LambdaSettings": None,        # solver setting
    "IterationControls": None,
    "UnderRelaxationCriteria": None,
    "Convergence": None,
    "InitialInputInfo": None,
    "ResultInputInfo": None,
    "SlipSurface": None,           # the search — reported separately, deliberately skipped
    # --- present, load-bearing, and NOT imported ---
    # Each of these changes the factor of safety. Reinforcement is deliberately not
    # guessed at: GeoStudio's nails/geosynthetics carry a pullout law, an out-of-plane
    # spacing and a force distribution that do not line up with xslope's model 1:1, and
    # a wrong reinforcement force is worse than an obviously absent one.
    "ReinforcementLines": "reinforcement (nails / geosynthetics / anchors)",
    "ReinforcementSets": "reinforcement sets",
    "LineLoadPoints": "line loads",
}


def list_analyses(gsz):
    """The analyses in a parsed .gsz: ``[{'id','name','kind','method'}, ...]``.

    A .gsz routinely bundles several (a base case plus variants). They can differ
    in materials, not just in slip surface, so the caller must pick one.
    """
    return list(gsz["analyses"])


def _strength(src):
    """(c, phi, caveat) for one GeoStudio material, according to its SlopeModel.

    Which field holds the strength depends on the model, and getting this wrong is
    silent: a MohrCoulomb material carries CohesionPrime/PhiPrime, but an
    UndrainedPhiZero one carries its undrained shear strength in <Cohesion> and has no
    PhiPrime at all. Reading the drained fields regardless hands back c = 0, phi = 0 —
    a soil with no strength, and a factor of safety near zero, with nothing to say so.
    """
    s = src["strength"]
    model = src["model"] or ""

    def f(key, default=0.0):
        try:
            return float(s[key])
        except (KeyError, TypeError, ValueError):
            return default

    if model == "MohrCoulomb":
        return f("CohesionPrime"), f("PhiPrime"), None

    if model == "UndrainedPhiZero":
        # phi = 0 by definition; Su lives in Cohesion (CohesionPrime is sometimes
        # written alongside it, but Cohesion is the authoritative one).
        su = f("Cohesion", f("CohesionPrime"))
        return su, 0.0, (
            f"material '{src['name']}' is undrained (phi = 0) in GeoStudio — imported "
            f"as Mohr-Coulomb with c = Su = {su:g} and phi = 0")

    if model == "Bedrock":
        # Handled by excluding the region from the domain (see gsz_to_slope_data) —
        # xslope's impenetrable boundary is geometric, not a material property.
        return f("CohesionPrime"), f("PhiPrime"), None

    # An unknown model. Take whatever strength fields exist, and say so plainly.
    c = f("CohesionPrime", f("Cohesion"))
    phi = f("PhiPrime")
    return c, phi, (
        f"material '{src['name']}' uses GeoStudio's {model or 'unnamed'} strength model, "
        f"which xslope does not have — imported as Mohr-Coulomb with c = {c:g}, "
        f"phi = {phi:g}, which is NOT the same strength; check it before solving")


def _circle_usable(slope_data, circle):
    """Can xslope actually build a slip surface from this circle on this model?

    SLOPE/W reports the circle it searched, but the surface it *scored* may have been
    truncated against an impenetrable boundary (a composite surface). That circle does
    not fit a domain whose floor is the bedrock, and importing it would hand the user a
    model that will not solve. Ask xslope's own slice generator rather than re-deriving
    the geometry here — it is the authority on what it will accept, and hand-rolling the
    check is how you end up rejecting circles that were fine.
    """
    from .slice import generate_slices          # local: slice imports fileio, we import both

    probe = dict(slope_data)
    probe["circles"] = [circle]
    probe["circular"] = True
    try:
        ok, _ = generate_slices(probe, circle=circle, num_slices=12, debug=False)
        return bool(ok)
    except Exception:
        return False


def _y_on(line, x):
    """Elevation of a left-to-right polyline at x, or None if x is off its ends."""
    coords = list(line.coords)
    if not coords or x < coords[0][0] or x > coords[-1][0]:
        return None
    for (x1, y1), (x2, y2) in zip(coords, coords[1:]):
        if x1 <= x <= x2:
            if x2 == x1:
                return max(y1, y2)
            return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
    return None


def ponded_water_dload(ground_surface, piezo_line, gamma_water):
    """Water standing above the ground surface, as a distributed load.

    GeoStudio has no ponded-water object: where the piezometric line rises above the
    ground, SLOPE/W simply takes the water's weight as a surcharge. xslope needs that
    surcharge to exist explicitly, so importing the piezo line alone would quietly lose
    the weight of the reservoir.

    Returns a list of dload blocks — one per submerged stretch — each a list of
    ``{'X','Y','Normal'}`` along the ground surface, with ``Normal`` the water pressure
    gamma_w * depth (zero at the waterline, so the load tapers out correctly).
    """
    from shapely.geometry import LineString

    if not piezo_line or ground_surface is None or ground_surface.is_empty:
        return []
    piezo = LineString(piezo_line)

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
    units = _detect_units(gamma_water, gsz.get("unit_system"))

    # GeoStudio's Bedrock is IMPENETRABLE: slip surfaces cannot enter it. xslope has no
    # such material — its impenetrable boundary is the domain itself. So drop bedrock
    # regions from the model and let the domain floor land on their top surface, which
    # reproduces SLOPE/W's behaviour exactly. Importing them as ordinary soil would let
    # trial surfaces cut straight through the bedrock.
    bedrock = {mid for mid, m in gsz["materials"].items() if m["model"] == "Bedrock"}
    dropped_bedrock = sorted({gsz["materials"][m]["name"]
                              for m in assign.values() if m in bedrock})
    if dropped_bedrock:
        keep = {r: m for r, m in assign.items() if m not in bedrock}
        if keep:
            assign = keep
            caveats.append(
                f"region(s) of {', '.join(repr(b) for b in dropped_bedrock)} are "
                f"GeoStudio's impenetrable Bedrock — excluded from the model, so the "
                f"domain now ends at the top of the bedrock. That is how xslope makes a "
                f"boundary impenetrable, and it matches SLOPE/W: slip surfaces cannot "
                f"enter it")
        else:
            caveats.append(
                f"every region in this analysis is Bedrock — imported as ordinary soil, "
                f"because excluding it would leave no model at all")

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
        c, phi, note = _strength(src)
        if note:
            caveats.append(note)
        mat.update(option="mc", gamma=src["gamma"], c=c, phi=phi)
        if not mat["gamma"]:
            caveats.append(f"material '{src['name']}' has no unit weight in the "
                           f"file — set gamma before solving")
        if not c and not phi:
            caveats.append(
                f"material '{src['name']}' came across with NO STRENGTH (c = 0, "
                f"phi = 0) — it will behave as a liquid; set its strength before "
                f"solving or the factor of safety will be meaningless")
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
    # surface maps onto xslope's piezo line. <MaterialUsesPiezs> says WHICH materials
    # draw from it — a material not listed there is dry even when a surface exists.
    piezo_line = []
    option = water.get("option") or "None"
    uses = gsz["stability"].get(analysis_id, {}).get("material_piezo", {})
    if option == "PiezoSurface":
        if gsz["piezo"]:
            ids = next(iter(gsz["piezo"].values()))
            piezo_line = [pts[i] for i in ids if i in pts]
            for mid, mat in zip(used, materials):
                # No MaterialUsesPiezs at all -> the surface applies to everything.
                if not uses or mid in uses:
                    mat["u"] = "piezo"
            dry = [gsz["materials"][m]["name"] for m in used if uses and m not in uses]
            if dry:
                caveats.append(
                    f"material(s) {', '.join(repr(d) for d in dry)} are not connected "
                    f"to the piezometric surface in GeoStudio — imported with no pore "
                    f"pressure")
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
    best = min(trials, key=lambda t: t["fs"]) if trials else None

    # Does SLOPE/W's critical circle actually fit this model? It may not: with an
    # impenetrable bedrock, SLOPE/W truncates the circle ALONG the bedrock and reports
    # the FS of that composite surface. The circle itself then cuts below our domain
    # floor, and importing it would hand the user a model that cannot be solved.
    if best is not None:
        trial = {"Xo": best["xo"], "Yo": best["yo"], "R": best["r"],
                 "Depth": best["yo"] - best["r"]}
        if not _circle_usable(slope_data, trial):
            caveats.append(
                f"xslope could not build a slip surface from SLOPE/W's critical circle "
                f"(SLOPE/W's FS = {best['fs']:.3f}, centre ({best['xo']:.1f}, "
                f"{best['yo']:.1f}), R = {best['r']:.1f}), so NO failure surface was "
                f"imported — define one, or run a search. Common causes: SLOPE/W scored a "
                f"COMPOSITE surface (the circle truncated along an impenetrable "
                f"boundary), or the mechanism is one xslope's slice generator rejects")
            best = None

    if best is not None:
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

    stab = gsz["stability"].get(analysis_id, {})

    # --- distributed loads -------------------------------------------------------
    # GeoStudio surcharges map straight onto xslope dloads: a pressure over a run of
    # ground, normal to the surface.
    dloads = []
    for sc in stab.get("surcharges", []):
        dloads.append([{"X": x, "Y": y, "Normal": sc["pressure"]}
                       for x, y in sc["coords"]])
    if dloads:
        caveats.append(f"{len(dloads)} surcharge load(s) imported as distributed loads")

    # Ponded water. GeoStudio stores no such object — where the piezometric line rises
    # above the ground, SLOPE/W just takes the water's weight as a surcharge. xslope
    # needs it explicitly, so synthesize it, or the weight of the reservoir is lost.
    if piezo_line:
        ponded = ponded_water_dload(ground_surface, piezo_line, gamma_water)
        if ponded:
            deepest = max(p["Normal"] for b in ponded for p in b) / (gamma_water or 1.0)
            dloads.extend(ponded)
            caveats.append(
                f"the piezometric line runs above the ground surface — that is ponded "
                f"water, and GeoStudio carries its weight implicitly. It has been added "
                f"as {len(ponded)} distributed load(s), up to {deepest:.2f} deep "
                f"({gamma_water * deepest:.1f} pressure at the deepest point)")
    slope_data["dloads"] = dloads

    # --- tension crack -----------------------------------------------------------
    # GeoStudio's crack is an arbitrary polyline; xslope's is a uniform depth below the
    # ground surface. Convert by depth, and say so when the two cannot agree.
    tc = stab.get("tcrack")
    if tc and tc.get("coords"):
        depths = []
        for x, y in tc["coords"]:
            yg = _y_on(ground_surface, x)
            if yg is not None and yg - y > 0:
                depths.append(yg - y)
        if depths:
            depth = max(depths)
            slope_data["tcrack_depth"] = depth
            slope_data["tcrack_water"] = (tc.get("pct_water") or 0.0) * depth
            if max(depths) - min(depths) > 0.05 * max(depths):
                caveats.append(
                    f"GeoStudio's tension crack is a line at varying depth "
                    f"({min(depths):.2f} to {max(depths):.2f} below ground); xslope "
                    f"models a crack of constant depth, so the deepest was taken "
                    f"({depth:.2f}) — conservative, but not the same crack")
            else:
                caveats.append(f"tension crack imported at depth {depth:.2f}"
                               + (f", {100*(tc.get('pct_water') or 0):.0f}% full of water"
                                  if tc.get("pct_water") else ""))
            if tc.get("angle") is not None:
                caveats.append(
                    f"the tension crack is inclined ({tc['angle']:g} degrees) in "
                    f"GeoStudio — xslope's crack is vertical, so the inclination was "
                    f"dropped")
        else:
            caveats.append("GeoStudio defines a tension crack, but it does not sit "
                           "below the ground surface — it was NOT imported")

    # Anything GeoStudio attached to this analysis that we did NOT take. Driven off the
    # exhaustive _ENTRY_ELEMENTS table, so an element we have never seen is reported as
    # unrecognized rather than passed over — silence here means a wrong factor of safety
    # with nothing to warn the user.
    for tag in stab.get("elements", []):
        if tag not in _ENTRY_ELEMENTS:
            caveats.append(
                f"this analysis has a GeoStudio '{tag}' that xslope does not recognise "
                f"and did NOT import — if it affects the model, the imported factor of "
                f"safety will be wrong")
        elif _ENTRY_ELEMENTS[tag]:
            caveats.append(
                f"{_ENTRY_ELEMENTS[tag]} defined in GeoStudio were NOT imported — "
                f"xslope has no mapping for them; re-enter them or the factor of safety "
                f"will be wrong")

    if stab.get("k_seismic_v"):
        caveats.append(
            f"a vertical seismic coefficient ({stab['k_seismic_v']:g}) is set in "
            f"GeoStudio — xslope models only the horizontal one, so it was dropped")

    caveats.append(f"model read as {units} units "
                   f"({'kN/m3, m, kPa' if units == 'metric' else 'lb/ft3, ft, psf'}) "
                   f"— confirm the units match your template; xslope does not convert "
                   f"between unit systems")
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

    Checked two ways: the file round-trips through this module's reader, and every
    tag path it writes is one that GeoStudio's own files use. The round-trip alone
    proves nothing — a reader and writer sharing the same wrong idea of the schema
    agree perfectly with each other — so the conformance check is the load-bearing
    one. GeoStudio remains the only authority on its own format.
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

    # The piezometric surface gets its OWN points, allocated in the order the line
    # runs, rather than reusing coincident region vertices. GeoStudio expects a
    # surface's point IDs to ascend along the polyline; reusing a region corner hands
    # it a low ID out of sequence, and it reports the surface as corrupt.
    piezo_ids = []
    piezo_line = slope_data.get("piezo_line") or []
    uses_piezo = any((m.get("u") or "none") == "piezo" for m in materials)
    if piezo_line and not uses_piezo:
        caveats.append("a piezo line is defined but no material uses it — written as "
                       "a piezometric surface anyway")
    for xy in piezo_line:
        pid = len(point_id) + 1
        point_id[("piezo", len(piezo_ids))] = pid     # keyed so it never dedups
        points[pid] = (round(xy[0], 9), round(xy[1], 9))
        piezo_ids.append(pid)

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

    # Edges. GeoStudio draws a region from its <Lines>, not from its point list, so a
    # file without them opens with nothing visible. Each region ring contributes its
    # consecutive pairs; an edge shared by two regions is written once.
    edges, edge_id = {}, {}
    for r in regions:
        ring = r["ids"]
        for a, b in zip(ring, ring[1:] + ring[:1]):
            key = (min(a, b), max(a, b))
            if key not in edge_id:
                edge_id[key] = len(edge_id) + 1
                edges[edge_id[key]] = (a, b)

    xs = [p[0] for p in points.values()]
    ys = [p[1] for p in points.values()]
    mx = max((max(xs) - min(xs)) * 0.1, 1.0)
    my = max((max(ys) - min(ys)) * 0.1, 1.0)
    x0, x1 = min(xs) - mx, max(xs) + mx
    y0, y1 = min(ys) - my, max(ys) + my
    unit_system = "Metric" if abs(gamma_water - _GAMMA_W_METRIC) < 1.0 else "Imperial"
    # The drawing scale relates model units to the printed page. GeoStudio's page is in
    # inches, so the model unit -> inch factor differs by system (m: 39.37, ft: 12).
    per_inch = 39.3701 if unit_system == "Metric" else 12.0
    # Pick a round scale that lands the model on a page about 10 inches wide.
    scale = max(1.0, round((x1 - x0) * per_inch / 10.0, -1)) if (x1 - x0) else 200.0

    L = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
         '<GSIData Version="11.11" AppVersion="25.2.1.4">',
         f'  <FileInfo Title="{_xml_escape(analysis_name)}" '
         f'Comments="Exported by xslope." />',
         '  <Analyses Len="1">',
         f'    <Analysis><ID>1</ID><Name>{_xml_escape(analysis_name)}</Name>'
         f'<Kind>SLOPE/W</Kind><Method>{_xml_escape(method)}</Method>'
         '<GeometryId>1</GeometryId></Analysis>',
         '  </Analyses>',
         # The view window, in engineering coordinates, and the declared unit system.
         # Without these GeoStudio has no canvas extent to draw the model into.
         '  <Coordinates>',
         f'    <EngCoords HorzScale="{scale:.6g}" VertScale="{scale:.6g}" '
         f'XPageLeft="{x0:.6g}" XPageRight="{x1:.6g}" '
         f'YPageBottom="{y0:.6g}" YPageTop="{y1:.6g}" '
         f'XPageOrg="{-x0:.6g}" YPageOrg="{-y0:.6g}" '
         f'MaxSnapDist="20" UnitSystem="{unit_system}" LockScales="false" />',
         f'    <PageCoords Units="in" PageWidth="{(x1-x0)/scale*per_inch:.6g}" '
         f'PageHeight="{(y1-y0)/scale*per_inch:.6g}" PageXOrg="0" PageYOrg="0" />',
         '  </Coordinates>',
         '  <Geometries Len="1">',
         '    <Geometry><Name>2D Geometry</Name>',
         f'      <Points Len="{len(points)}">']
    for pid in sorted(points):
        x, y = points[pid]
        L.append(f'        <Point ID="{pid}" X="{x:.10g}" Y="{y:.10g}" Pinned="true" />')
    L.append('      </Points>')
    L.append(f'      <Lines Len="{len(edges)}">')
    for lid in sorted(edges):
        a, b = edges[lid]
        L.append(f'        <Line><ID>{lid}</ID><PointID1>{a}</PointID1>'
                 f'<PointID2>{b}</PointID2></Line>')
    L.append('      </Lines>')
    L.append(f'      <Regions Len="{len(regions)}">')
    for i, r in enumerate(regions, start=1):
        L.append(f'        <Region><ID>{i}</ID>'
                 f'<PointIDs>{",".join(str(j) for j in r["ids"])}</PointIDs>'
                 f'<Mesh Pattern="Unstructured Mixed Elements" /></Region>')
    L.append('      </Regions>')
    # No <Window>: Base/Zoom are GeoStudio's saved scroll position in its own screen
    # units, which we cannot compute from a model. Omitting it lets GeoStudio pick its
    # own view; writing a guess left the model off-centre with the axes in shot.
    L.append('    </Geometry>')
    L.append('  </Geometries>')

    # Carry xslope's own material colours across, so the zones are visually
    # distinguishable in GeoStudio instead of all landing on its default yellow.
    from .plot import get_material_color

    L.append(f'  <Materials Len="{len(materials)}">')
    for i, m in enumerate(materials, start=1):
        try:
            rgb = get_material_color(i - 1)                 # tab10, 0-based mat_id
            r8, g8, b8 = (int(round(255 * c)) for c in rgb[:3])
        except Exception:
            r8 = g8 = b8 = 200
        L.append(
            f'    <Material><ID>{i}</ID><Color>RGB=({r8},{g8},{b8})</Color>'
            f'<Name>{_xml_escape(m.get("name") or f"material {i}")}</Name>'
            f'<SlopeModel>MohrCoulomb</SlopeModel><StressStrain>'
            f'<UnitWeight>{float(m.get("gamma") or 0.0):.10g}</UnitWeight>'
            f'<CohesionPrime>{float(m.get("c") or 0.0):.10g}</CohesionPrime>'
            f'<PhiPrime>{float(m.get("phi") or 0.0):.10g}</PhiPrime>'
            f'</StressStrain></Material>')
    L.append('  </Materials>')

    # Region -> material, per analysis (GeoStudio keeps it here, not on the region).
    L.append('  <Contexts Len="1">')
    L.append('    <Context><AnalysisID>1</AnalysisID>')
    L.append(f'      <GeometryUsesMaterials Len="{len(regions)}">')
    for i, r in enumerate(regions, start=1):
        L.append(f'        <GeometryUsesMaterial ID="Regions-{i}" Entry="{r["mat"]}" />')
    L.append('      </GeometryUsesMaterials>')
    L.append('      <IsDefined>true</IsDefined></Context>')
    L.append('  </Contexts>')

    L.append('  <WaterItems Len="1">')
    L.append('    <WaterItem><AnalysisID>1</AnalysisID><Entry>'
             + ('<ResultInputInfo><Option>PiezoSurface</Option></ResultInputInfo>'
                if piezo_ids else '')
             + f'<UnitWaterWeight>{gamma_water:.10g}</UnitWaterWeight>'
             '</Entry></WaterItem>')
    L.append('  </WaterItems>')

    # The piezometric surface belongs to the StabilityItem, NOT the geometry, and the
    # materials that draw pore pressure from it are named in <MaterialUsesPiezs>.
    k = float(slope_data.get("k_seismic") or 0.0)
    L.append('  <StabilityItems Len="1">')
    L.append('    <StabilityItem><AnalysisID>1</AnalysisID><Entry>')
    L.append(f'      <Seismic Horizontal="{k:.10g}" Vertical="" />' if k else
             '      <Seismic Horizontal="" Vertical="" />')
    if piezo_ids:
        L.append('      <PiezometricSurfaces Len="1">')
        L.append('        <PiezometricSurface><ID>1</ID>'
                 f'<DataPoints Len="{len(piezo_ids)}">'
                 + "".join(f'<DataPoint>{i}</DataPoint>' for i in piezo_ids)
                 + '</DataPoints><CapSuction>false</CapSuction>'
                 '<MaxSuction>0</MaxSuction></PiezometricSurface>')
        L.append('      </PiezometricSurfaces>')
        wet = [i for i, m in enumerate(materials, start=1)
               if (m.get("u") or "none") == "piezo"] or list(range(1, len(materials) + 1))
        L.append(f'      <MaterialUsesPiezs Len="{len(wet)}">')
        for i in wet:
            L.append(f'        <MaterialUsesPiez ID="{i}" UsesID="1" />')
        L.append('      </MaterialUsesPiezs>')
    L.append('    </Entry></StabilityItem>')
    L.append('  </StabilityItems>')
    L.append('</GSIData>')

    stem = os.path.splitext(os.path.basename(gsz_path))[0]
    with zipfile.ZipFile(gsz_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(f"{stem}.xml", "\n".join(L) + "\n")
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
