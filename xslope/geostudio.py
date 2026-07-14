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
import datetime
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
            # An ALTERNATIVE pullout law: bond from soil adhesion + friction on the
            # interface, which grows with the normal stress on the bar. xslope's bond is a
            # constant rate, so this one cannot be reproduced -- it is reported instead.
            "adhesion": _num(rf, "InterfaceAdhesion", None),
            "shear_angle": _num(rf, "InterfaceShearAngle", None),
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

        # Piezometric surfaces. They index this analysis's LOCAL DataPoints -- the same
        # list the tension crack, surcharges and reinforcement index -- NOT the shared
        # geometry <Points> table. Reading them against the geometry is silent nonsense:
        # in one verification model it produced a water table that doubled back on
        # itself, and pore pressures low enough to lift the factor of safety by 10%.
        piezo = {}
        for ps in entry.findall("./PiezometricSurfaces/PiezometricSurface"):
            pid = _text(ps, "ID")
            if pid is None:
                continue
            ids = [int(dp.text) for dp in ps.findall("./DataPoints/DataPoint") if dp.text]
            piezo[int(pid)] = [local[i] for i in ids if i in local]

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

        # Line loads. The point is the local DataPoint with the same ID; the direction is
        # a trend/plunge pair — plunge is the angle BELOW horizontal, trend picks which
        # way along x. xslope wants an angle measured counter-clockwise from +x, so a
        # plunge of 90 (straight down) becomes -90, which is xslope's own default.
        line_loads = []
        for lp_ in entry.findall("./LineLoadPoints/LineLoadPoint"):
            pid = _text(lp_, "ID")
            ll = lp_.find("LineLoad")
            if pid is None or ll is None or int(pid) not in local:
                continue
            value = _num(ll, "Value", 0.0)
            d = ll.find("Direction")
            plunge = float(d.get("Plunge")) if (d is not None and d.get("Plunge")) else 90.0
            trend = float(d.get("Trend")) if (d is not None and d.get("Trend")) else 0.0
            angle = -plunge if abs(trend) < 90 else -(180.0 - plunge)
            x, y = local[int(pid)]
            if value > 0:
                line_loads.append({"x": x, "y": y, "P": value, "angle": angle,
                                   "label": f"Line load {pid}"})

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
            "piezo": piezo,                  # piezometric surface ID -> [(x, y), ...]
            "points": local,
            "tcrack": tcrack,
            "surcharges": surcharges,
            "reinf_lines": reinf_lines,
            "line_loads": line_loads,
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
        "contexts": contexts, "water": water,
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
    "ReinforcementLines": None,
    "LineLoadPoints": None,
    # --- present, load-bearing, and NOT imported ---
    "ReinforcementSets": "reinforcement sets (out-of-plane staging groups)",
}

# How GeoStudio says the tension crack is defined. 'Surface' means "by the line whose
# points follow", which is the only form xslope can convert. A crack that is switched
# OFF has NO <TensionOption> at all -- but GeoStudio keeps its geometry, so the presence
# of a crack line means nothing on its own. Anything else here is reported, not guessed.
_TCRACK_OPTIONS = {"Surface"}

# SLOPE/W's SlipOrigSearchMethod, in the solved results: which searches produce a real
# circle. 5 is entry-and-exit; 6 inherits the surface from a parent analysis, which is a
# circle only if the parent's was. 2 (fully specified) and 4 (block) are not circles at
# all -- though SLOPE/W still writes a centre and radius for them, which is the trap.
# Callers must confirm a 6 really is a circle; read_gsz_results cannot know.
_CIRCULAR_SEARCHES = {"5", "6"}


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


# GeoStudio reinforcement Type -> xslope (type, dir, appl). These mirror fileio's
# _TYPE_PRESETS, which is local to load_slope_data and so cannot be imported: 'axial'
# acts along the bar, 'tangent' reorients to the slip surface; 'passive' means the force
# is divided by the factor of safety (an ultimate capacity), 'active' means it is not.
#
# The file does NOT record a direction or an F-of-S dependence: GeoStudio implies both
# from the Type. So does this table, and each row is measured against SLOPE/W's own
# factor of safety rather than reasoned about -- see run_gsz_corpus().
_REINF_TYPE = {
    "Nail":         ("nail",         "axial",   "active"),
    "Geosynthetic": ("geosynthetic", "axial",   "active"),
    "Anchor":       ("anchor",       "axial",   "active"),
    "Tieback":      ("tieback",      "axial",   "active"),
}


def _reinforcement(product, line, ground_surface):
    """One GeoStudio reinforcement line -> one xslope ``reinforcement_lines`` entry.

    The two models line up more closely than they first appear. xslope's capacity at a
    point is ``cap_i = min(t_max, tend_i + t_max * d_i / lp_i)`` from each end i, taking
    the smaller — i.e. the force available at the end itself, plus bond resistance earned
    over the length d_i back to that end. ``t_max / lp_i`` IS a pullout rate, so:

        lp = t_max / (pullout resistance per unit length)

    reproduces GeoStudio's "pullout resistance x embedded length" exactly, and the plate
    (facing) capacity is the end capacity ``tend`` at whichever end is at the face.

    Everything is stored PER UNIT WIDTH: GeoStudio quotes per-element values with an
    out-of-plane Spacing, and xslope's ``reinforcement_lines`` are already divided by it.
    """
    spacing = product.get("spacing") or 1.0
    t_max = (product.get("tensile") or 0.0) / spacing
    plate = (product.get("plate") or 0.0) / spacing
    rate = (product.get("pullout") or 0.0) / spacing      # force per unit length
    lp = (t_max / rate) if rate > 0 else 0.0              # 0 => that end is fully anchored

    p1, p2 = line["p1"], line["p2"]
    # The plate/facing sits at the end that is AT the face — the endpoint on (or nearest)
    # the ground surface. The other end is buried and develops force only by bond.
    if ground_surface is not None and not ground_surface.is_empty:
        from shapely.geometry import Point
        d1 = ground_surface.distance(Point(*p1))
        d2 = ground_surface.distance(Point(*p2))
        face_is_1 = d1 <= d2
    else:
        face_is_1 = True

    xtype, xdir, xappl = _REINF_TYPE.get(product.get("type") or "",
                                         ("", "tangent", "active"))
    return {
        "x1": p1[0], "y1": p1[1], "x2": p2[0], "y2": p2[1],
        "t_max": t_max,
        "t_res": float("nan"),          # GeoStudio has no post-peak drop: hold capacity
        "lp1": lp, "lp2": lp,
        "tend1": plate if face_is_1 else 0.0,
        "tend2": 0.0 if face_is_1 else plate,
        "E": float("nan"), "area": float("nan"),
        "spacing": 1.0,                 # already per unit width — never divide twice
        "label": product.get("name") or "reinforcement",
        "type": xtype, "dir": xdir, "appl": xappl,
    }


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


def material_above_ground_dload(ground_surface, upper_line, unit_weight):
    """The weight of whatever fills the gap between a line and the ground, as dloads.

    Two GeoStudio things are this same object, and both would otherwise be lost:

    * **Ponded water.** GeoStudio has no ponded-water object -- where the piezometric
      line rises above the ground, SLOPE/W simply carries the water's weight. Pass the
      piezo line and gamma_w.
    * **A surcharge.** GeoStudio's ``<Surcharge>`` is a body of fill between a drawn
      line and the ground, and its ``<Pressure>`` is the fill's UNIT WEIGHT, not a
      pressure -- so the load varies with how deep the fill is. (Verified against
      SLOPE/W's own per-slice surcharge forces: a 20 kN/m3 fill under a line running
      1 m to 2 m above flat ground gives 300 kN/m, not 20 x 10 = 200.)

    Returns a list of dload blocks -- one per stretch where the line is above ground --
    each a list of ``{'X','Y','Normal'}``, with ``Normal`` the pressure
    ``unit_weight * depth``, zero at the crossing so the load tapers out correctly.
    """
    from shapely.geometry import LineString

    if not upper_line or ground_surface is None or ground_surface.is_empty:
        return []
    piezo = LineString(upper_line)
    gamma_water = unit_weight

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


def ponded_water_dload(ground_surface, piezo_line, gamma_water):
    """Water standing above the ground surface, as a distributed load.

    A thin name for :func:`material_above_ground_dload` -- ponded water is just the
    material above the ground being water.
    """
    return material_above_ground_dload(ground_surface, piezo_line, gamma_water)


def _dload_pressure_at(blocks, x):
    """The pressure a set of dload blocks applies at x (0 outside them)."""
    for b in blocks:
        xs = [p["X"] for p in b]
        if not xs or x < min(xs) - 1e-9 or x > max(xs) + 1e-9:
            continue
        for p, q in zip(b, b[1:]):
            lo, hi = sorted((p["X"], q["X"]))
            if lo - 1e-9 <= x <= hi + 1e-9:
                if abs(q["X"] - p["X"]) < 1e-12:
                    return max(p["Normal"], q["Normal"])
                t = (x - p["X"]) / (q["X"] - p["X"])
                return p["Normal"] + t * (q["Normal"] - p["Normal"])
    return 0.0


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
    stab = gsz["stability"].get(analysis_id, {})
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
    uses = stab.get("material_piezo", {})
    surfaces = stab.get("piezo", {})
    if option == "PiezoSurface":
        if surfaces:
            piezo_line = list(next(iter(surfaces.values())))
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
            if len(surfaces) > 1:
                caveats.append(
                    f"the file defines {len(surfaces)} piezometric surfaces; "
                    f"xslope takes one, so the first was imported")
        else:
            caveats.append("the analysis asks for a piezometric surface but the file "
                           "defines none — pore pressure imported as zero")
    elif option == "Parent":
        # The water condition is a whole head field computed by a parent SEEP/W run.
        # It is not just pore pressure: SLOPE/W also raises the water standing against
        # the slope from that same field, and that load is often the bigger half. In
        # GeoStudio's own rapid-drawdown example, dropping both costs 13% of the FS.
        caveats.append(
            "THE FACTOR OF SAFETY WILL BE WRONG: this analysis takes its water from a "
            "parent SEEP/W analysis, and xslope imports neither the pore pressure nor "
            "the weight of the water standing against the slope. The soil and geometry "
            "are right; the water is simply absent. Set a piezo line or a seepage "
            "solution before you trust any number this model produces")
    elif option not in ("None", "", None):
        caveats.append(
            f"THE FACTOR OF SAFETY WILL BE WRONG: pore pressure comes from SLOPE/W's "
            f"'{option}' option, which xslope cannot read — imported as zero. Set a "
            f"piezo line, ru, or a seepage solution before solving")

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
    # Only circular trials. SLOPE/W writes a centre and radius even for block and
    # fully-specified surfaces, but they describe a circle FITTED to the surface, not
    # the surface it solved -- importing one would silently swap the mechanism.
    circular = [t for t in trials if t["circular"]]
    best = min(circular, key=lambda t: t["fs"]) if circular else None
    if trials and not circular:
        caveats.append(
            "SLOPE/W's critical surface here is non-circular (a block or fully "
            "specified surface), which is not a trial circle xslope can take — no "
            "failure surface was imported; define one, or run a search")

    # Does SLOPE/W's critical circle actually fit this model? It may not: with an
    # impenetrable bedrock, SLOPE/W truncates the circle ALONG the bedrock and reports
    # the FS of that composite surface. The circle itself then cuts below our domain
    # floor, and importing it would hand the user a model that cannot be solved.
    # SLOPE/W's own critical circle may not be one xslope can build -- with an
    # impenetrable bedrock SLOPE/W truncates the circle ALONG the bedrock and scores
    # that COMPOSITE surface, and the circle itself then cuts below our domain floor. So
    # walk its trials from the most critical up and take the first that does build: the
    # model still arrives complete and comparable, just on a slightly less severe circle.
    chosen = rejected = None
    for t in sorted(circular, key=lambda t: t["fs"]):
        trial = {"Xo": t["xo"], "Yo": t["yo"], "R": t["r"], "Depth": t["yo"] - t["r"]}
        if _circle_usable(slope_data, trial):
            chosen = t
            break
        if rejected is None:
            rejected = t

    if chosen is not None:
        slope_data["circles"] = [{"Xo": chosen["xo"], "Yo": chosen["yo"], "R": chosen["r"],
                                  "Depth": chosen["yo"] - chosen["r"]}]
        slope_data["circular"] = True
        caveats.append(
            f"the failure surface is a trial circle of SLOPE/W's own (its FS = "
            f"{chosen['fs']:.3f} by {analysis['method'] or 'its method'}), not a search "
            f"— it is imported so the model is complete and the two programs can be "
            f"compared on the same circle; run a search to find xslope's own critical "
            f"surface")
        if rejected is not None:
            caveats.append(
                f"note that SLOPE/W's MOST critical surface (its FS = "
                f"{rejected['fs']:.3f}) is not a circle xslope can build — most likely "
                f"a COMPOSITE surface, truncated along the impenetrable bedrock. The "
                f"circle imported above is the most critical one that does build, so "
                f"its FS is higher than SLOPE/W's")
    elif circular:
        caveats.append(
            f"none of SLOPE/W's {len(circular)} trial circles can be built on this "
            f"model (its critical FS = {rejected['fs']:.3f}), so NO failure surface was "
            f"imported — define one, or run a search. Usually this means SLOPE/W scored "
            f"COMPOSITE surfaces, truncated along an impenetrable boundary")
    else:
        caveats.append(
            "no failure surface was imported — SLOPE/W's search definition has no "
            "xslope equivalent, and this file carries no solved circles to take one "
            "from; define circles or a non-circular surface before solving (an input "
            "file with no surface will not re-load)")

    # --- distributed loads -------------------------------------------------------
    # A GeoStudio surcharge is NOT a uniform pressure, despite its <Pressure> tag: it is
    # a body of fill between the line you draw and the ground, and <Pressure> is that
    # fill's UNIT WEIGHT. So the load grows with the depth of the fill, and importing it
    # as a constant pressure is wrong everywhere the wedge is not exactly 1 deep.
    dloads = []
    n_sc = 0
    for sc in stab.get("surcharges", []):
        blocks = material_above_ground_dload(ground_surface, sc["coords"], sc["pressure"])
        dloads.extend(blocks)
        n_sc += bool(blocks)
        if not blocks:
            caveats.append(
                f"a surcharge of unit weight {sc['pressure']:g} does not sit above the "
                f"ground surface anywhere, so it carries no load — NOT imported")
    if n_sc:
        caveats.append(f"{n_sc} surcharge load(s) imported as distributed loads — "
                       f"GeoStudio's surcharge is a wedge of fill above the ground, so "
                       f"the load varies with its depth")

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
    opt = (tc.get("option") or "") if tc else ""
    if tc and not opt:
        # A crack that is switched OFF still keeps its geometry: GeoStudio drops the
        # <TensionOption> element and leaves the DataPoints behind. Applying that
        # leftover line as a real crack is silently wrong -- in a c'=0 soil a
        # water-filled crack is pure driving force.
        tc = None
    elif tc and opt not in _TCRACK_OPTIONS:
        caveats.append(f"GeoStudio defines the tension crack by '{opt}', which xslope "
                       f"does not recognise — the crack was NOT imported, so the factor "
                       f"of safety will be too high")
        tc = None

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
        else:
            caveats.append("GeoStudio defines a tension crack, but it does not sit "
                           "below the ground surface — it was NOT imported")

    # --- reinforcement -----------------------------------------------------------
    from .fileio import build_reinforce_lines

    rlines = []
    for rl in stab.get("reinf_lines", []):
        product = gsz["reinforcements"].get(rl["reinforcement"])
        if product is None:
            caveats.append("a reinforcement line references a product the file does not "
                           "define — skipped")
            continue
        rlines.append(_reinforcement(product, rl, ground_surface))
    if rlines:
        slope_data["reinforcement_lines"] = rlines
        slope_data["reinforce_lines"] = build_reinforce_lines(rlines)
        kinds = sorted({r["type"] or "generic" for r in rlines})
        caveats.append(
            f"{len(rlines)} reinforcement line(s) imported ({', '.join(kinds)}), acting "
            f"along the bar as a known load — the convention that reproduces SLOPE/W's "
            f"own factor of safety. GeoStudio's pullout resistance became xslope's bond "
            f"length (Lp = Tmax / pullout rate), and its plate capacity the end capacity "
            f"at the face")

        # A pullout law xslope cannot reproduce: bond that grows with the normal stress on
        # the bar. Silence here would hand back a plausible force computed the wrong way.
        used = [gsz["reinforcements"][rl["reinforcement"]] for rl in stab["reinf_lines"]
                if rl["reinforcement"] in gsz["reinforcements"]]
        interface = sorted({p["name"] for p in used
                            if p.get("adhesion") or p.get("shear_angle")})
        if interface:
            caveats.append(
                f"reinforcement {', '.join(repr(n) for n in interface)} takes its pullout "
                f"from the soil INTERFACE (adhesion + friction), which grows with the "
                f"normal stress on the bar. xslope's bond is a constant rate, so the "
                f"imported bond length is an approximation — check the developed forces "
                f"against GeoStudio if the bars are not fully anchored")
        spread = sorted({p["name"] for p in used
                         if (p.get("distribution") or "Concentrated") != "Concentrated"})
        if spread:
            caveats.append(
                f"reinforcement {', '.join(repr(n) for n in spread)} spreads its force "
                f"over the bar in GeoStudio rather than applying it at the slip surface; "
                f"xslope applies it where the bar crosses the surface")

    # --- line loads --------------------------------------------------------------
    lloads = []
    for ll in stab.get("line_loads", []):
        lloads.append({"x": ll["x"], "y": ll["y"], "P": ll["P"],
                       "angle": ll["angle"], "label": ll["label"]})
    if lloads:
        slope_data["line_loads"] = lloads
        caveats.append(f"{len(lloads)} line load(s) imported")

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

    # A statement of fact, not a warning. There is nothing to convert and nothing to
    # confirm: xslope has no unit setting anywhere, and an import builds a whole new
    # model whose geometry, strengths and unit weight of water all come out of this one
    # file. They are consistent with each other by construction. All the user needs to
    # know is which system the numbers they are about to look at are in.
    caveats.append(
        f"this model is in {units} units "
        f"({'kN/m3, m, kPa' if units == 'metric' else 'lb/ft3, ft, psf'}) — xslope is "
        f"unit-agnostic and imports them as they are, so results come back in the same "
        f"system")
    return slope_data, caveats


def gsz_style(gsz, analysis_id=None):
    """The material colours from a .gsz, as an xslope style sidecar.

    GeoStudio gives every material a colour, and it is how people recognise their own
    model at a glance -- so a model that comes back in xslope's default palette looks
    like someone else's. Returns ``{'materials': {'<mat_id>': {'color': '#rrggbb'}}}``,
    keyed the way :mod:`xslope.style` expects, or ``{}`` if the file names no colours.

    Kept out of ``gsz_to_slope_data`` on purpose: colour is presentation, and
    ``slope_data`` is the model that gets written to the input file.
    """
    if analysis_id is None:
        analysis_id = gsz["analyses"][0]["id"] if gsz["analyses"] else None
    assign = gsz["contexts"].get(analysis_id, {})
    bedrock = {mid for mid, m in gsz["materials"].items() if m["model"] == "Bedrock"}

    # Same ordering rule gsz_to_slope_data uses, so mat_id lines up with the polygons.
    used = sorted({m for m in assign.values() if m not in bedrock})
    colors = {}
    for i, mid in enumerate(used):
        rgb = _parse_rgb(gsz["materials"][mid].get("color"))
        if rgb:
            colors[str(i)] = {"color": "#%02x%02x%02x" % rgb}
    return {"materials": colors} if colors else {}


def _parse_rgb(text):
    """GeoStudio writes colour as the string ``RGB=(211,201,137)``."""
    m = re.search(r"RGB=\((\d+),\s*(\d+),\s*(\d+)\)", text or "")
    if not m:
        return None
    r, g, b = (min(255, max(0, int(v))) for v in m.groups())
    return r, g, b


def read_gsz_results(gsz, analysis_name):
    """SLOPE/W's own solved trial slip surfaces for an analysis, if the file has them.

    Returns ``[{'fs','xo','yo','r','weight','circular'}, ...]`` — one entry per trial
    surface SLOPE/W evaluated, with the factor of safety it computed. Empty if the file
    was saved unsolved.

    ``weight`` is SLOPE/W's own weight of the sliding mass, which is worth as much as the
    factor of safety: it checks the imported geometry and unit weights on their own,
    without the solver in the way.

    ``circular`` is False for a surface SLOPE/W did not generate as a circle (a block or
    fully-specified one). Those rows still carry a centre and a radius, but they describe
    a *fitted* circle, not the surface that was solved -- rebuilding them as circles
    silently scores the wrong geometry.

    This is not part of the model. It exists so the same surfaces can be re-solved
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
        out.append({"fs": fs, "xo": xo, "yo": yo, "r": radius, "weight": weight,
                    "circular": row.get("SlipOrigSearchMethod") in _CIRCULAR_SEARCHES})
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

    # The piezometric surface does NOT live in the geometry. Its points go in the
    # analysis's own <DataPoints> list -- the local, self-contained one that carries X
    # and Y and is numbered from 1 -- and the surface indexes THAT. Writing geometry
    # point IDs here instead makes GeoStudio resolve them against an empty local list,
    # and the water table silently disappears.
    piezo_line = slope_data.get("piezo_line") or []
    piezo_ids = list(range(1, len(piezo_line) + 1))
    uses_piezo = any((m.get("u") or "none") == "piezo" for m in materials)
    if piezo_line and not uses_piezo:
        caveats.append("a piezo line is defined but no material uses it — written as "
                       "a piezometric surface anyway")

    # Needed to tell ponded water (which the piezo line already carries) from a real
    # surcharge, and to sit a surcharge's fill wedge on top of the ground.
    ground_surface = slope_data.get("ground_surface")
    if ground_surface is None or getattr(ground_surface, "is_empty", True):
        ground_surface = build_ground_surface_from_polygons(polygons)

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
                       ("line_loads", "line loads")]:
        if slope_data.get(key):
            caveats.append(f"{label} not written — no .gsz mapping; re-create in GeoStudio")

    # Distributed loads. Two very different cases hide behind one xslope key, and telling
    # the user to "re-create them in GeoStudio" was wrong for the common one:
    #
    #  * PONDED WATER. GeoStudio has no ponded-water object -- it DERIVES the reservoir,
    #    and the load it puts on the slope, from the piezometric surface itself. So the
    #    load is already carried by the piezo line we wrote. Re-creating it by hand would
    #    DOUBLE-COUNT the water.
    #  * A REAL SURCHARGE. That we can now write: GeoStudio's <Surcharge> is a wedge of
    #    fill between a drawn line and the ground, and its <Pressure> is the fill's unit
    #    weight. Any pressure profile p(x) is therefore exactly a wedge of unit weight
    #    gamma whose upper line sits at y_ground(x) + p(x)/gamma.
    ponded = material_above_ground_dload(ground_surface, piezo_line, gamma_water)
    surcharges = []
    for block in (slope_data.get("dloads") or []):
        if all(abs(p["Normal"] - _dload_pressure_at(ponded, p["X"])) <= 1e-6 *
               max(1.0, abs(p["Normal"])) for p in block):
            continue                      # the piezo surface already carries this
        surcharges.append(block)
    if ponded and len(surcharges) < len(slope_data.get("dloads") or []):
        caveats.append(
            "the ponded water was NOT written as a load, and must not be re-created: "
            "GeoStudio derives the reservoir and its pressure on the slope from the "
            "piezometric surface, which was written. Adding it by hand would count the "
            "water twice")

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

    now = datetime.datetime.now()

    L = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
         '<GSIData Version="11.11" AppVersion="25.2.1.4">',
         f'  <FileInfo FileVersion="11.11" AppVersion="25.2.1.4" '
         f'Title="{_xml_escape(analysis_name)}" Comments="Exported by xslope." '
         f'Author="xslope" LastAuthor="xslope" RevNumber="1" '
         f'Date="{now:%d/%m/%Y}" Time="{now:%I:%M:%S %p}" />',
         '  <Analyses Len="1">',
         '    <Analysis>',
         '      <ID>1</ID>',
         f'      <Name>{_xml_escape(analysis_name)}</Name>',
         '      <Kind>SLOPE/W</Kind>',
         f'      <Method>{_xml_escape(method)}</Method>',
         '      <GeometryId>1</GeometryId>',
         # GeoStudio writes all three of these on every analysis. They are not decoration:
         # <ComputedPhysics> is how it knows this analysis SOLVES slope stability. Without
         # it the file opens, the geometry draws and the materials are named and coloured
         # -- but there is no physics attached to them, so the material Properties button
         # is greyed out and the strengths cannot be seen or edited. It looks like the
         # strengths were lost, and nothing says otherwise.
         '      <TimeIncrements>',
         '        <Duration Missing="true" />',
         '        <TimeSteps Len="1"><TimeStep Save="true" /></TimeSteps>',
         '      </TimeIncrements>',
         '      <IterationControls Len="1">',
         '        <IterationControl>',
         '          <Key>SlopeStability</Key>',
         '          <Entry MaxIterations="100" MaxReviewIterations="10" />',
         '        </IterationControl>',
         '      </IterationControls>',
         '      <ComputedPhysics SlopeStability="true" />',
         '    </Analysis>',
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
        # Shaped exactly like one GeoStudio writes. The Missing="true" placeholders are
        # how it records "this parameter is not set" -- they are part of a well-formed
        # StressStrain block, not noise, and it writes them on every material.
        L.append(
            f'    <Material><ID>{i}</ID><Color>RGB=({r8},{g8},{b8})</Color>'
            f'<Name>{_xml_escape(m.get("name") or f"material {i}")}</Name>'
            f'<SlopeModel>MohrCoulomb</SlopeModel>'
            f'<StressModel>MohrCoulomb</StressModel>'
            f'<StressStrain>'
            f'<JointEffectiveCohesion Missing="true" />'
            f'<IntactRockParam Missing="true" />'
            f'<GeologicalStrengthIndex Missing="true" />'
            f'<DisturbanceFactor Missing="true" />'
            f'<MaxConfiningStress Missing="true" />'
            f'<UseTensileStrength>false</UseTensileStrength>'
            f'<UnitWeight>{float(m.get("gamma") or 0.0):.10g}</UnitWeight>'
            f'<CohesionPrime>{float(m.get("c") or 0.0):.10g}</CohesionPrime>'
            f'<PhiPrime>{float(m.get("phi") or 0.0):.10g}</PhiPrime>'
            f'<TensileStrength Missing="true" />'
            f'<OCRatio Missing="true" />'
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

    # GeoStudio writes the bulk modulus of water on every WaterItem. It is a constant
    # (SLOPE/W never uses it -- it is there for the coupled products), but it is quoted in
    # the file's own pressure unit, so it has to follow the unit system.
    bulk = 2.08333333e6 if unit_system == "Metric" else 4.351132120308411e7
    L.append('  <WaterItems Len="1">')
    L.append('    <WaterItem><AnalysisID>1</AnalysisID><Entry>'
             + ('<ResultInputInfo><Option>PiezoSurface</Option></ResultInputInfo>'
                if piezo_ids else '')
             + f'<UnitWaterWeight>{gamma_water:.10g}</UnitWaterWeight>'
             + f'<WaterBulkMod>{bulk:.9g}</WaterBulkMod>'
             '</Entry></WaterItem>')
    L.append('  </WaterItems>')

    # The piezometric surface belongs to the StabilityItem, NOT the geometry, and the
    # materials that draw pore pressure from it are named in <MaterialUsesPiezs>.
    k = float(slope_data.get("k_seismic") or 0.0)

    # The analysis-local point list. The piezometric surface AND the surcharges both index
    # into it, so it has to be built before either is written. A surcharge is a wedge of
    # fill of unit weight gamma between its line and the ground, so a pressure profile
    # p(x) becomes the line y_ground(x) + p(x)/gamma. Any gamma reproduces p exactly; the
    # unit weight of water is used because it is already in the file and unit-consistent.
    local_pts = list(piezo_line)                        # numbered 1..n
    sc_ids = []
    gamma_s = gamma_water or 1.0
    for block in surcharges:
        ids = []
        for p in block:
            yg = _y_on(ground_surface, p["X"])
            if yg is None:
                continue
            local_pts.append((p["X"], yg + p["Normal"] / gamma_s))
            ids.append(len(local_pts))
        if len(ids) >= 2:
            sc_ids.append(ids)

    L.append('  <StabilityItems Len="1">')
    L.append('    <StabilityItem><AnalysisID>1</AnalysisID><Entry>')
    if local_pts:
        L.append(f'      <DataPoints Len="{len(local_pts)}">')
        for i, (x, y) in enumerate(local_pts, start=1):
            L.append(f'        <DataPoint Number="{i}" X="{x:.10g}" Y="{y:.10g}" />')
        L.append('      </DataPoints>')
    L.append(f'      <Seismic Horizontal="{k:.10g}" Vertical="" />' if k else
             '      <Seismic Horizontal="" Vertical="" />')
    if sc_ids:
        L.append(f'      <Surcharges Len="{len(sc_ids)}">')
        for i, ids in enumerate(sc_ids, start=1):
            L.append(f'        <Surcharge><ID>{i}</ID>'
                     f'<DataPoints Len="{len(ids)}">'
                     + "".join(f'<DataPoint>{j}</DataPoint>' for j in ids)
                     + f'</DataPoints><Pressure>{gamma_s:.10g}</Pressure></Surcharge>')
        L.append('      </Surcharges>')
        caveats.append(
            f"{len(sc_ids)} distributed load(s) written as GeoStudio surcharges — a "
            f"surcharge there is a wedge of fill, so they appear as fill of unit weight "
            f"{gamma_s:g} whose depth reproduces the pressure exactly")
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
