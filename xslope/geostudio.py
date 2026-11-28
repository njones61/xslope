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
  - read_gsz_seep(gsz, analysis_name): the pore-pressure field a parent SEEP/W
    analysis computed, as an xslope mesh + nodal u.

Unlike DXF, a .gsz is *semantically complete* — regions, materials and water
conditions already know what they are, so no layer-mapping wizard is needed. The
one thing the caller must choose is which analysis to import, because a single
.gsz commonly holds several (and they can differ in materials as well as in
slip-surface definition — see ``gsz_to_slope_data``).

Only what SLOPE/W and xslope have in common is imported -- regions, materials
(Mohr-Coulomb and undrained; SHANSEP as ``ratio``; S=f(depth) as ``cp``; bedrock as
``elastic``), the piezometric surface
(with SLOPE/W's phreatic correction where the analysis applies it), Ru, water
above ground, a solved SEEP/W field, a staged rapid drawdown (both piezometric
surfaces and the drawdown strengths), a separate unit weight above the water table,
surcharges and line loads, reinforcement, the tension crack, the seismic coefficient
and the solved critical circle. Everything else is reported as a caveat rather than
silently dropped: unsupported strength models, unsupported pore-pressure options
(B-bar and spatial functions), partial factors, SLOPE/W's search definition and an
optimized or non-circular critical surface, reinforcement sets, inclined tension
cracks and the vertical seismic coefficient.
"""

import csv
import datetime
import io
import math
import os
import re
import struct
import zipfile
import xml.etree.ElementTree as ET

import numpy as np
from shapely.geometry import Polygon

from .fileio import (build_ground_surface_from_polygons, cp_from_layer_top,
                     main_sheet_cell)
from .mesh import ensure_ccw_elements, interpolate_at_point
from .units import GAMMA_W, require_gamma_water, units_check, normalize_unit_system
# Water geometry lives in a neutral module. material_above_ground_dload is used here
# for GeoStudio's SURCHARGE (a wedge of fill, priced the same way); derive_water_loads
# is what the import reports before the engine applies it. ponded_water_dload is no
# longer called from this module -- the reservoir is derived rather than synthesised --
# and is re-exported so `from xslope.geostudio import ponded_water_dload` keeps working.
from .water import (_y_on, derive_water_loads, material_above_ground_dload,
                    ponded_water_dload)


# Strength models SLOPE/W offers that xslope can represent (see _strength). Anything
# else is imported as a placeholder and flagged.
_SUPPORTED_MODELS = {"MohrCoulomb", "UndrainedPhiZero", "Bedrock", "SHANSEP", "SFnDepth"}

# Unit weight of water, used to infer the file's unit system (see _detect_units) and
# as the value assigned when a Metric GeoStudio file omits <UnitWaterWeight>.
# NOTE: GeoStudio's canonical metric gamma_w is 9.807, NOT xslope's SI 9.81
# (GAMMA_W["si"]). This vendor-specific value is kept deliberately so importing a
# GeoStudio file reproduces its own number bit-for-bit; only the Imperial value
# (62.4) coincides with xslope's canonical and is sourced from the single-source-of-
# truth. Detection tolerates the 0.003 gap either way (bands are +/- 0.15).
_GAMMA_W_METRIC = 9.807                   # kN/m3 -- GeoStudio vendor canonical
_GAMMA_W_IMPERIAL = GAMMA_W["imperial"]   # lb/ft3


def _water_loads_cell():
    """How a message names the Water loads cell: by its label, and by its address on
    the current template, which is the one an import is written into."""
    return f"Water loads cell ({main_sheet_cell('water_loads')})"


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


#: A probabilistic modifier's key, e.g. ``Materials[3].Applied.EffectiveCohesion``.
#: The bracketed number is the material's <ID>, not its position in the file.
_PROB_KEY_RE = re.compile(r"^\s*Materials\[(\d+)\]\.(.+?)\s*$")

#: A distribution expression, e.g. ``Normal(Mean=0,SD=1150,Min=-2000,Max=2000)``.
_PROB_FN_RE = re.compile(r"^\s*([A-Za-z]\w*)\s*\((.*)\)\s*$")


def _parse_prob_function(expr):
    """``("Normal", {"Mean": 0.0, "SD": 1150.0, ...})`` from a distribution expression.

    Returns ``(None, {})`` for anything that is not a ``Name(k=v,...)`` expression.
    Values that are not numbers are kept as their raw text, so an unrecognised
    distribution can still be NAMED in a caveat rather than silently dropped.
    """
    m = _PROB_FN_RE.match(expr or "")
    if not m:
        return None, {}
    params = {}
    for part in m.group(2).split(","):
        if "=" not in part:
            continue
        k, _, v = part.partition("=")
        v = v.strip()
        try:
            params[k.strip()] = float(v)
        except ValueError:
            params[k.strip()] = v
    return m.group(1), params


def _prob_modifiers(ss_node):
    """The <ProbabilisticModifiers> of one <SampleSettings>, parsed but uninterpreted.

    Each entry is ``{"key", "mat_id", "prop", "dist", "params", "corr_tag", "corr"}``.
    ``mat_id``/``prop`` are None for a key that does not name a material parameter (the
    schema also modifies loads and reinforcement), so the caller can report those rather
    than mistake them for a material's sigma.
    """
    out = []
    if ss_node is None:
        return out
    for mod in ss_node.findall("./ProbabilisticModifiers/Modifier"):
        key = _text(mod, "Key")
        entry = mod.find("Entry")
        if key is None or entry is None:
            continue
        dist, params = _parse_prob_function(entry.get("Function"))
        km = _PROB_KEY_RE.match(key)
        try:
            corr = float(entry.get("CorrelationCoefficient"))
        except (TypeError, ValueError):
            corr = None
        out.append({
            "key": key,
            "mat_id": int(km.group(1)) if km else None,
            "prop": km.group(2) if km else None,
            "dist": dist,
            "params": params,
            "corr_tag": entry.get("CorrelationTag"),
            "corr": corr,
        })
    return out


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
        parent = _text(a, "ParentID")
        # A transient analysis saves a result set per time step. Its SLOPE/W child
        # solves the slope at every one of them, so an imported model has to say WHICH.
        # Step 0 is the initial condition (inherited, t = 0); the listed steps follow it.
        times = [0.0]
        for ts in a.findall("./TimeIncrements/TimeSteps/TimeStep"):
            if ts.get("Save") == "true" and ts.get("ElapsedTime"):
                times.append(float(ts.get("ElapsedTime")))
        # An analysis names its geometry by <GeometryId>. GeoStudio puts no <ID> on the
        # <Geometry> element itself — the reference is POSITIONAL, 1-based, in document
        # order. A single-geometry file (the vast majority) always says GeometryId 1.
        gid = _text(a, "GeometryId")
        # Probabilistic sampling. SLOPE/W keeps the uncertainty on the ANALYSIS, not on
        # the material: <SampleSettings> names the sampling type and holds one <Modifier>
        # per uncertain parameter. The same material can therefore be deterministic in one
        # analysis of a file and probabilistic in another, which is why this is parsed here
        # and not alongside <Materials>.
        ss_node = a.find("SampleSettings")
        analyses.append({
            "id": int(aid),
            "name": _text(a, "Name", f"Analysis {aid}"),
            "kind": _text(a, "Kind", ""),
            "method": _text(a, "Method", ""),
            "parent": int(parent) if parent else None,
            "geometry_id": int(gid) if gid else None,
            "times": times,
            "sampling": _text(ss_node, "Type", "") or "",
            "trials": _num(ss_node, "NumMonteCarloTrials", None),
            "modifiers": _prob_modifiers(ss_node),
        })

    # <Geometries> holds one <Geometry> PER geometry, in document order, keyed 1-based to
    # match the analyses' <GeometryId>. They must be kept SEPARATE, not merged: a file with
    # two geometries (e.g. the same embankment at two heights) reuses point and region IDs
    # across them — both number from 1 — so a single merged table lets the last geometry
    # parsed silently overwrite the first, and every analysis then resolves to whichever
    # geometry came last. gsz_to_slope_data picks the analysis's own geometry.
    geometries = {}
    for gidx, g in enumerate(root.findall("./Geometries/Geometry"), start=1):
        gpoints, gregions, gglines = {}, {}, {}
        for p in g.findall("./Points/Point"):
            # A deleted point keeps its ID slot but loses its coordinates
            # (<Point ID="4" /> — seen in the staged rapid-drawdown examples).
            # Skip it: nothing may reference it, and if something does, the
            # plain dict lookup downstream raises KeyError on the real culprit
            # instead of this parse dying on every orphaned slot.
            if p.get("X") is None or p.get("Y") is None:
                continue
            gpoints[int(p.get("ID"))] = (float(p.get("X")), float(p.get("Y")))
        for r in g.findall("./Regions/Region"):
            ids = _text(r, "PointIDs")
            if ids:
                gregions[int(_text(r, "ID"))] = [int(i) for i in ids.split(",")]
        for ln in g.findall("./Lines/Line"):
            p1, p2 = _text(ln, "PointID1"), _text(ln, "PointID2")
            if p1 and p2:
                gglines[int(_text(ln, "ID"))] = (int(p1), int(p2))
        geometries[gidx] = {"points": gpoints, "regions": gregions, "lines": gglines}

    # Flat, merged views for any caller that reads gsz["points"]/["regions"]/["lines"]
    # directly. For a single-geometry file these are the whole model; for a multi-geometry
    # one they are last-geometry-wins, which is exactly why gsz_to_slope_data does NOT use
    # them — it goes through gsz["geometries"] keyed by the analysis's GeometryId instead.
    points, regions, glines = {}, {}, {}
    for geom in geometries.values():
        points.update(geom["points"])
        regions.update(geom["regions"])
        glines.update(geom["lines"])

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

    # Material functions: a shared library of (x, y) tables keyed by ID, which a
    # material names by number -- SHANSEP's OCR against elevation (OCRFnNum) and its
    # undrained strength against vertical effective stress (ShearVerticalFnNum).
    # <Function> holds the fit and its input/output, e.g. Spline(OutputParam=
    # OverconsolidationRatio,InputParam=Y,...).
    mfns = {}
    for fn in root.findall("./Functions/Material/Fn2Ds/Fn2D"):
        fid = _text(fn, "ID")
        if fid is None:
            continue
        pts = [(float(p.get("X")), float(p.get("Y")))
               for p in fn.findall("./Points/Point")]
        mfns[int(fid)] = {"name": _text(fn, "Name", ""),
                          "function": _text(fn, "Function", ""), "points": pts}

    # Spatial Mohr-Coulomb's functions. A material whose <SlopeCohesionModel> (or Phi,
    # or Weight) is "Function" names one of these curves by <CohesionFnNum>: the value
    # against one coordinate, InputParam=X or Y. "Spatial Function" names a Fn3D by
    # <CohesionSpatialFnNum>: the value over (x, y), its points stored as X, Y, Z.
    sfns = {}
    for kind in ("CohesionFn", "PhiFn", "WeightFn"):
        for fn in root.findall(f"./Functions/Material/StressStrain/{kind}s/{kind}"):
            fid = _text(fn, "ID")
            if fid is None:
                continue
            expr = _text(fn, "Function", "") or ""
            inp = re.search(r"InputParam=(\w+)", expr)
            sfns[(kind, int(fid))] = {
                "name": _text(fn, "Name", ""), "function": expr,
                "input": inp.group(1) if inp else None,
                "points": [(float(p.get("X")), float(p.get("Y")))
                           for p in fn.findall("./Points/Point")]}
    fn3ds = {}
    for fn in root.findall("./Functions/Func3Ds/Fn3D"):
        fid = _text(fn, "ID")
        if fid is None:
            continue
        fn3ds[int(fid)] = {
            "name": _text(fn, "Name", ""), "output": _text(fn, "OutputType", ""),
            "points": [(float(p.get("X")), float(p.get("Y")), float(p.get("Z")))
                       for p in fn.findall("./Points/Points_")]}

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
        # The tensile capacity is <Tensile> in GeoStudio 2025.2 files, which write the
        # plate capacity beside it as <PlateCapacity>. An earlier 2025 build wrote one
        # <TensileCapacity> and no plate capacity: SLOPE/W's own Tie-Back Wall example
        # comes both ways, and its 2025.2 copy carries <Tensile> and <PlateCapacity>
        # both equal to the older file's <TensileCapacity>, with the same solved factor
        # of safety. So the older tag is the same capacity, and stands for the plate
        # capacity too when the file gives none.
        tensile = _num(rf, "Tensile", None)
        old_tag = tensile is None and _num(rf, "TensileCapacity", None) is not None
        if old_tag:
            tensile = _num(rf, "TensileCapacity", 0.0)
        plate = _num(rf, "PlateCapacity", None)
        plate_from_tensile = old_tag and plate is None
        reinforcements[int(rid)] = {
            "name": _text(rf, "Name", ""),
            "type": _text(rf, "Type", ""),                 # Nail | Geosynthetic | Anchor | ...
            "tensile": tensile or 0.0,
            "plate": (tensile if plate_from_tensile else plate) or 0.0,
            "capacity_tag": "TensileCapacity" if old_tag else "Tensile",
            "plate_from_tensile": plate_from_tensile,
            "pullout": _num(rf, "PulloutResistance", 0.0),
            "spacing": _num(rf, "Spacing", None),          # out-of-plane, None = per metre
            # GeoStudio writes a field only when it differs from the default, and
            # these four defaults are measured from SLOPE/W's own files: no <Type> is
            # an Anchor ("Reinforcement with Anchors" has no <Type> on any of its four
            # products, which carry the anchor-only BondLength and BondDiameter); no
            # <ForceDistribution> is Distributed (that file's "Distributed load"
            # analysis uses the product with no tag, its "Concentrated load" analysis
            # the one tagged Concentrated); no <FofSDependent> is No.
            "distribution": _text(rf, "ForceDistribution", "") or "Distributed",
            "fofs_dependent": (_text(rf, "FofSDependent", "") or "").strip().lower()
                              == "true",
            "shear": _num(rf, "Shear", None),              # anchor/nail/pile shear force
            "shear_option": _text(rf, "ShearOption", ""),
            # The ALTERNATIVE pullout law, and xslope's Adhesion/Delta columns:
            # bond from soil adhesion + friction on the interface, growing with the
            # effective overburden. GeoStudio states the interface strength per unit
            # AREA, tau = a + sigma'v*tan(delta) (Slope Stability Modeling, Eq. 67),
            # and multiplies it by a length measure xi to get force per unit length
            # (Eq. 65): for a geosynthetic xi = SurfaceAreaFactor (default 2, both
            # faces), for an anchor or nail xi = pi*BondDiameter (Eq. 69). Which law
            # is in force is CalculatedPulloutResistance -- a product may carry both
            # a constant PulloutResistance and the interface pair.
            "adhesion": _num(rf, "InterfaceAdhesion", None),
            "shear_angle": _num(rf, "InterfaceShearAngle", None),
            "calculated_pullout": (_text(rf, "CalculatedPulloutResistance", "")
                                   or "").strip().lower() == "true",
            "surface_area_factor": _num(rf, "SurfaceAreaFactor", None),
            "bond_diameter": _num(rf, "BondDiameter", None),
            "bond_length": _num(rf, "BondLength", None),
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
        # SLOPE/W parameterises unsaturated (negative pore-pressure) shear strength on
        # the piezometric surface, not the material: CapSuction (whether to limit the
        # suction magnitude) and MaxSuction (the ceiling). This is a threshold model on
        # the pore pressure, NOT xslope's phi_b/s_cap apparent-cohesion pair, so the
        # importer reports it rather than converting it — the values are kept here so it
        # can say what SLOPE/W was set to.
        piezo = {}
        suction = None
        for ps in entry.findall("./PiezometricSurfaces/PiezometricSurface"):
            pid = _text(ps, "ID")
            if pid is None:
                continue
            ids = [int(dp.text) for dp in ps.findall("./DataPoints/DataPoint") if dp.text]
            piezo[int(pid)] = [local[i] for i in ids if i in local]
            cap = (_text(ps, "CapSuction", "") or "").lower() == "true"
            mx = _num(ps, "MaxSuction", None)
            if cap or mx:
                suction = {"cap": cap, "max": mx}

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

        # Analysis switches, e.g. <IsPiezPhreatic> (SLOPE/W's phreatic correction on
        # the piezometric surfaces) and <IsRapidDrawdown> (its staged rapid drawdown).
        # GeoStudio writes a switch only when it is on.
        control = {c.tag: (c.text or "").strip() for c in entry.findall("./Control/*")}

        # Ru and B-bar, per material ID. Each is read by SLOPE/W only under its own
        # water option (Ru; PiezoSurfaceBbar or Bbar), and GeoStudio keeps the table
        # on analyses that use neither -- often with empty values. B-bar's 'IsIncluded'
        # marks a material whose weight drives the excess pressure ("Add weight").
        ru = {}
        for p in entry.findall("./MaterialUsesRuParams/MaterialUsesRuParam"):
            pid, e = _text(p, "ID"), p.find("Entry")
            if pid and e is not None:
                ru[int(pid)] = (e.get("Function") or "").strip()
        bbar = {}
        for p in entry.findall("./MaterialUsesBBarParams/MaterialUsesBBarParam"):
            pid, e = _text(p, "ID"), p.find("Entry")
            if pid and e is not None:
                bbar[int(pid)] = {"bbar": (e.get("Function") or "").strip(),
                                  "add_weight": e.get("IsIncluded") == "true"}

        stability[int(_text(s, "AnalysisID"))] = {
            "k_seismic": float(kh) if kh else 0.0,
            "k_seismic_v": float(kv) if kv else 0.0,
            "material_piezo": uses,          # material ID -> piezometric surface ID
            "piezo": piezo,                  # piezometric surface ID -> [(x, y), ...]
            "suction": suction,              # {'cap','max'} unsaturated-strength setting
            "points": local,
            "tcrack": tcrack,
            "surcharges": surcharges,
            "reinf_lines": reinf_lines,
            "line_loads": line_loads,
            "control": control,              # analysis switches: tag -> text
            "phreatic": control.get("IsPiezPhreatic", "").lower() == "true",
            "rapid_drawdown": control.get("IsRapidDrawdown", "").lower() == "true",
            "ru": ru,                        # material ID -> Ru (as written)
            "bbar": bbar,                    # material ID -> {'bbar', 'add_weight'}
            # Everything GeoStudio put on this analysis, so gsz_to_slope_data can
            # report anything it does not import instead of dropping it in silence.
            "elements": [c.tag for c in entry],
        }

    # GeoStudio declares its unit system rather than leaving it to be inferred.
    eng = root.find("./Coordinates/EngCoords")
    unit_system = eng.get("UnitSystem") if eng is not None else None

    # Partial factor sets (Eurocode 7 design approaches), shared by the analyses, which
    # name the one they apply by <Control><DesignFactorID>. Each factor is a key and a
    # favorable and an unfavorable value, either of which may be blank.
    design_factors = {}
    for ds in root.findall("./DesignFactors/DesignFactors_"):
        did = _text(ds, "ID")
        if did is None:
            continue
        factors = []
        for f in ds.findall("./Factors/Factors_"):
            e = f.find("Entry")
            factors.append((_text(f, "Key", ""),
                            (e.get("Favorable") or "").strip() if e is not None else "",
                            (e.get("Unfavorable") or "").strip() if e is not None else ""))
        design_factors[int(did)] = {"name": _text(ds, "Name", f"design factors {did}"),
                                    "factors": factors}

    return {
        "path": str(path), "zip": zf, "analyses": analyses, "points": points,
        "regions": regions, "lines": glines, "geometries": geometries,
        "materials": materials,
        "contexts": contexts, "water": water,
        "stability": stability, "kfns": kfns, "mfns": mfns, "sfns": sfns,
        "fn3ds": fn3ds, "bcs": bcs,
        "hyd_bcs": hyd_bcs, "reinforcements": reinforcements,
        "unit_system": unit_system, "design_factors": design_factors,
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
    "Control": None,               # switches; each is checked against _CONTROL_HANDLED
    "MaterialUsesRuParams": None,  # read with the water option that uses them
    "MaterialUsesBBarParams": None,
    # --- present, load-bearing, and NOT imported ---
    "ReinforcementSets": "reinforcement sets (out-of-plane staging groups)",
}

# The <Control> switches gsz_to_slope_data acts on: the phreatic correction and a
# staged rapid drawdown are both imported, and a set of partial factors is named in a
# caveat. Any other switch is reported as unrecognized, the same as an unknown element.
_CONTROL_HANDLED = {"IsPiezPhreatic", "IsRapidDrawdown", "DesignFactorID"}

# How GeoStudio says the tension crack is defined. 'Surface' means "by the line whose
# points follow", which is the only form xslope can convert. A crack that is switched
# OFF has NO <TensionOption> at all -- but GeoStudio keeps its geometry, so the presence
# of a crack line means nothing on its own. Anything else here is reported, not guessed.
_TCRACK_OPTIONS = {"Surface"}

# SLOPE/W's SlipOrigSearchMethod, in the solved results: the search that produced each
# trial surface. Read off SLOPE/W's own verification files by pairing each analysis's
# <SlipSurface><Option> with the codes in its slip_surface.csv. Every analysis carries
# a single code on all its rows, and each option always the same code:
#
#   1  GridAndRadius     4  Block                      7  Automatic (auto-locate)
#   2  FullySpecified    5  entry and exit, the        6  FromOther: the parent
#                           default (no <Option>)         analysis's surface
#
# 1 and 5 search circles. 2 and 4 are not circles at all -- though SLOPE/W still writes
# a center and radius for them, a circle FITTED to the surface, which is the trap. 7 is
# not taken either. 6 inherits the parent's surface, so it is whatever the parent's
# search was: read_gsz_results looks the parent up (_search_option) rather than count
# a 6 as a circle -- Tandjiria's reinforced block analyses are 6s with a Block parent.
_SEARCH_KIND = {"1": "circle", "5": "circle", "2": "specified", "4": "block",
                "7": "automatic"}

# The same, by the analysis's <SlipSurface><Option>, for a 6 (by its parent's) and for
# older files (AppVersion 11.4 in the vendor's set), which write no SlipOrigSearchMethod
# column at all and name the radius column SlipRadius rather than SlipRadiusX. None is
# the default entry-and-exit search, which GeoStudio writes as no <Option>.
_OPTION_KIND = {"GridAndRadius": "circle", None: "circle", "FullySpecified": "specified",
                "Block": "block", "Automatic": "automatic"}

# How a caveat names a surface that is not a circle, by its kind. 'optimized' is the
# surface SLOPE/W's optimization makes out of the critical one (see read_gsz_results).
_SURFACE_WORDS = {
    "block": "a block surface",
    "specified": "a fully specified surface",
    "automatic": "a surface from SLOPE/W's automatic search",
    "optimized": "an optimized surface",
    "unknown": "a surface from a search XSLOPE does not recognize",
}


def list_analyses(gsz):
    """The analyses in a parsed .gsz: ``[{'id','name','kind','method'}, ...]``.

    A .gsz routinely bundles several (a base case plus variants). They can differ
    in materials, not just in slip surface, so the caller must pick one.
    """
    return list(gsz["analyses"])


# StressStrain fields that belong to SIGMA/W's stress-strain models (stiffness, Cam
# clay, Hoek-Brown inputs GeoStudio writes as placeholders) rather than to a SLOPE/W
# strength. Every material carries some of them whatever its SlopeModel.
_STRESS_STRAIN_FIELDS = {
    "UnitWeight", "ResponseType", "Rf", "OCRatio", "ConsolLambda", "OCKappa",
    "InitVoidRatio", "YieldSurfaceShape", "LimitOfAnisotropy",
    "RateOfAnisotropyEvolution", "YoungsPrimeModulus", "EffectivePoissonsRatio",
    "EPrimeModulusXPrime", "EffectivePoissonsRatioXPrime", "EPrimeModulusYPrime",
    "EffectivePoissonsRatioYPrime", "ShearModulusXYPrime", "JointEffectiveCohesion",
    "IntactRockParam", "GeologicalStrengthIndex", "DisturbanceFactor",
    "MaxConfiningStress", "UseTensileStrength", "SlopeSuctionModel",
    "ShearAnisoFnNum",          # reported on its own, for every model
    "DryWeight", "SlopeUseUnitWeightAboveWT",   # the unit weight above the water table
}

# The fields each mapped strength model reads, plus the ones GeoStudio leaves behind
# from a material's earlier model (a SHANSEP material keeps its old CohesionPrime and
# PhiPrime; its manual and its results show they play no part). Anything else on one
# of these materials is reported by name rather than ignored -- SLOPE/W's S=f(depth)
# has an optional maximum strength, and no file we have seen writes it.
_MODEL_FIELDS = {
    "SHANSEP": {"TauSigmaRatio", "SlopeUndrainedStrengthModel",
                "NCUndrainedStrengthRatio", "SHANSEPExponent", "OCRFnNum",
                "ShearVerticalFnNum", "CohesionPrime", "PhiPrime", "Cohesion"},
    "SFnDepth": {"CTopOfLayer", "CRateOfIncrease", "CohesionPrime", "PhiPrime"},
    # TauSigmaRatio is SHANSEP's, left behind on a material switched to Spatial
    # Mohr-Coulomb: SLOPE/W's own example carries it and its results ignore it.
    "MohrCoulombSpatial": {"TauSigmaRatio", "CohesionPrime", "PhiPrime",
                           "SlopeCohesionModel",
                           "SlopePhiModel", "SlopeWeightModel", "CohesionFnNum",
                           "PhiFnNum", "WeightFnNum", "CohesionSpatialFnNum",
                           "PhiSpatialFnNum", "WeightSpatialFnNum"},
}


def _strength(src, mfns=None, sfns=None, fn3ds=None):
    """``(fields, notes, from_top)`` for one GeoStudio material, by its SlopeModel.

    ``fields`` updates the xslope material (option, c, phi, cp); ``notes`` are its
    caveats; ``from_top`` is True when r_elev is the top of the material's own layer,
    which the caller sets once the zones are built.

    Which field holds the strength depends on the model, and getting this wrong is
    silent: a MohrCoulomb material carries CohesionPrime/PhiPrime, but an
    UndrainedPhiZero one carries its undrained shear strength in <Cohesion> and has no
    PhiPrime at all. Reading the drained fields regardless hands back c = 0, phi = 0 —
    a soil with no strength, and a factor of safety near zero, with nothing to say so.

    Two SLOPE/W models map onto XSLOPE's undrained options (SLOPE/W's Slope Stability
    Modeling book, sections 4.9 and 4.14):

    - SHANSEP with a constant ratio, s_u = (tau/sigma) x sigma'_v, is ``ratio`` with
      c/p = the ratio and no minimum (c = 0). Its older name is S=f(overburden).
      SLOPE/W's sigma'_v is the overburden, surcharges and ponded water included,
      point loads not, less the pore pressure: the same as XSLOPE's.
    - S=f(depth), s_u = CTopOfLayer + CRateOfIncrease x (depth below the top of the
      layer), is ``cp`` measured from that top.
    """
    s = src["strength"]
    model = src["model"] or ""
    name = src["name"]
    mfns = mfns or {}

    def f(key, default=0.0):
        try:
            return float(s[key])
        except (KeyError, TypeError, ValueError):
            return default

    def unread():
        extra = sorted(set(s) - _STRESS_STRAIN_FIELDS - _MODEL_FIELDS.get(model, set()))
        if not extra:
            return []
        return [f"material '{name}' also sets GeoStudio field(s) {', '.join(extra)}, "
                f"which XSLOPE does not read. Check whether they change its strength."]

    if model == "MohrCoulomb":
        return {"option": "mc", "c": f("CohesionPrime"), "phi": f("PhiPrime")}, [], False

    if model == "UndrainedPhiZero":
        # phi = 0 by definition; Su lives in Cohesion (CohesionPrime is sometimes
        # written alongside it, but Cohesion is the authoritative one).
        su = f("Cohesion", f("CohesionPrime"))
        return {"option": "mc", "c": su, "phi": 0.0}, [
            f"material '{name}' is undrained (phi = 0) in GeoStudio — imported "
            f"as Mohr-Coulomb with c = Su = {su:g} and phi = 0"], False

    if model == "Bedrock":
        # Handled by excluding the region from the domain (see gsz_to_slope_data) —
        # xslope's impenetrable boundary is geometric, not a material property.
        return {"option": "mc", "c": f("CohesionPrime"), "phi": f("PhiPrime")}, [], False

    if model == "SHANSEP":
        how = s.get("SlopeUndrainedStrengthModel") or ""
        if not how:
            # The constant ratio (GeoStudio writes no <SlopeUndrainedStrengthModel>).
            r = f("TauSigmaRatio")
            return {"option": "ratio", "c": 0.0, "phi": 0.0, "cp": r}, [
                f"material '{name}' uses SLOPE/W's SHANSEP model with a constant "
                f"ratio — imported as option ratio with c/p = {r:g} and no minimum "
                f"strength (c = 0)"] + unread(), False
        if how == "Custom":
            # s_u = sigma'_v x S x OCR^m, OCR a function of elevation.
            S, m = f("NCUndrainedStrengthRatio"), f("SHANSEPExponent")
            fn = mfns.get(int(f("OCRFnNum", -1)))
            ocrs = [y for _, y in fn["points"]] if fn and fn["points"] else []
            if not ocrs:
                return {"option": "ratio", "c": 0.0, "phi": 0.0, "cp": S}, [
                    f"material '{name}' uses SLOPE/W's SHANSEP model with the ratio "
                    f"S x OCR^m (S = {S:g}, m = {m:g}), but the file gives no OCR "
                    f"function XSLOPE can read. It was imported as option ratio with "
                    f"c/p = S = {S:g}, as if OCR = 1. Check it before solving."
                ] + unread(), False
            rs = sorted(S * max(o, 0.0) ** m for o in (min(ocrs), max(ocrs)))
            r_lo, r_hi = rs
            if r_hi - r_lo <= 1e-12 * max(1.0, r_hi):
                return {"option": "ratio", "c": 0.0, "phi": 0.0, "cp": r_lo}, [
                    f"material '{name}' uses SLOPE/W's SHANSEP model with the ratio "
                    f"S x OCR^m and a constant OCR of {ocrs[0]:g} — imported as option "
                    f"ratio with c/p = {r_lo:g} and no minimum strength (c = 0)"
                ] + unread(), False
            fname = fn["name"] or f"function {int(f('OCRFnNum'))}"
            return {"option": "ratio", "c": 0.0, "phi": 0.0, "cp": r_lo}, [
                f"material '{name}' uses SLOPE/W's SHANSEP model with the ratio "
                f"S x OCR^m (S = {S:g}, m = {m:g}), and OCR varies with elevation "
                f"('{fname}': {min(ocrs):g} to {max(ocrs):g} at its points). XSLOPE has "
                f"one ratio per material, so it was imported as option ratio with the "
                f"lowest, c/p = {r_lo:g}. Where OCR is higher, XSLOPE's strength is "
                f"lower than SLOPE/W's, by up to a factor of {r_hi / r_lo:.3g}."
                if r_lo > 0 else
                f"material '{name}' uses SLOPE/W's SHANSEP model with the ratio "
                f"S x OCR^m (S = {S:g}, m = {m:g}), and OCR varies with elevation "
                f"('{fname}': {min(ocrs):g} to {max(ocrs):g} at its points). XSLOPE has "
                f"one ratio per material, so it was imported as option ratio with the "
                f"lowest, c/p = {r_lo:g}, which leaves it no strength. Set its strength "
                f"before solving."] + unread(), False
        # 'Function': s_u read off a curve against sigma'_v. No single ratio is that
        # curve, so the material is refused rather than given a made-up one.
        if how == "Function":
            fn = mfns.get(int(f("ShearVerticalFnNum", -1)))
            what = f"'{fn['name']}'" if fn and fn.get("name") else "a function"
            why = (f"with its strength read from a curve of vertical effective stress "
                   f"({what}). XSLOPE cannot represent that curve")
        else:
            why = (f"with the option '{how}', which XSLOPE does not know. It cannot "
                   f"tell what strength SLOPE/W gives it")
        return {"option": "mc", "c": 0.0, "phi": 0.0}, [
            f"material '{name}' uses SLOPE/W's SHANSEP model {why}, so the material "
            f"was imported with no strength. Set its strength before solving."
        ] + unread(), False

    if model == "SFnDepth":
        c_top, rate = f("CTopOfLayer"), f("CRateOfIncrease")
        return {"option": "cp", "c": c_top, "phi": 0.0, "cp": rate}, unread(), True

    if model == "MohrCoulombSpatial":
        fields, notes = _spatial_mc(name, s, f, sfns or {}, fn3ds or {})
        return fields, notes + unread(), False

    # An unknown model. Take whatever strength fields exist, and say so plainly.
    c = f("CohesionPrime", f("Cohesion"))
    phi = f("PhiPrime")
    return {"option": "mc", "c": c, "phi": phi}, [
        f"material '{name}' uses GeoStudio's {model or 'unnamed'} strength model, "
        f"which xslope does not have — imported as Mohr-Coulomb with c = {c:g}, "
        f"phi = {phi:g}, which is NOT the same strength; check it before solving"], False


def _linear_in_y(points):
    """``(c_top, rate, y_top, y_low, c_low)`` when a curve of value against elevation
    is one straight line falling with elevation (so rising with depth), else None.

    SLOPE/W's own Spatial Mohr-Coulomb example (in its SHANSEP example file) is one:
    c' = 27.52 at y = 0 and 5.4 at y = 9, and its saved column strengths are that line
    exactly, held at the end values beyond the two points. Above its top that is
    XSLOPE's ``cp`` (c, c/p and r_elev), which holds c above r_elev too.
    """
    pts = sorted(points, key=lambda p: p[0])
    if len(pts) < 2 or pts[-1][0] - pts[0][0] <= 0:
        return None
    (y0, c0), (y1, c1) = pts[0], pts[-1]
    rate = (c0 - c1) / (y1 - y0)                 # strength gained per unit depth
    if rate <= 0:
        return None
    tol = 1e-6 * max(1.0, abs(c0), abs(c1))
    if any(abs(c - (c1 + rate * (y1 - y))) > tol for y, c in pts):
        return None
    return c1, rate, y1, y0, c0


def _spatial_mc(name, s, f, sfns, fn3ds):
    """``(fields, notes)`` for SLOPE/W's Spatial Mohr-Coulomb material.

    Its cohesion, friction angle and unit weight can each be a constant, a curve
    against x or y ("Function", <CohesionFnNum>), or a value over (x, y) ("Spatial
    Function", <CohesionSpatialFnNum>). A cohesion that is a straight line rising
    with depth, on a friction angle of zero, is XSLOPE's ``cp`` and imports as it.
    Every other function is replaced by one value -- the lowest of a strength
    function, the average of a unit weight function -- and the notes name the
    material, the function, its range and the value used.
    """
    fields, notes = {}, []

    def curve(prop, tag):
        """(label, values, fn) for the property's function, or None if constant."""
        how = s.get(f"Slope{prop}Model") or ""
        if how == "Function":
            fn = sfns.get((f"{prop}Fn", int(f(f"{prop}FnNum", -1))))
            if fn is None or not fn["points"]:
                return (f"a {tag} function the file does not define", [], None)
            label = (f"'{fn['name']}', {tag} against {fn['input'] or 'a coordinate'}"
                     if fn["name"] else f"a {tag} function")
            return label, [v for _, v in fn["points"]], fn
        if how == "Spatial Function":
            fn = fn3ds.get(int(f(f"{prop}SpatialFnNum", -1)))
            if fn is None or not fn["points"]:
                return (f"a spatial {tag} function the file does not define", [], None)
            label = (f"'{fn['name']}', {tag} over x and y" if fn["name"]
                     else f"a spatial {tag} function over x and y")
            return label, [z for _, _, z in fn["points"]], fn
        return None

    c_fn, phi_fn, w_fn = (curve("Cohesion", "cohesion"), curve("Phi", "friction angle"),
                          curve("Weight", "unit weight"))
    c, phi = f("CohesionPrime"), f("PhiPrime")

    if phi_fn is not None:
        label, vals, _ = phi_fn
        if vals:
            phi = min(vals)
            notes.append(
                f"material '{name}' takes its friction angle from {label}, from "
                f"{min(vals):g} to {max(vals):g} degrees at its points. XSLOPE has one "
                f"friction angle per material, so it was imported with the lowest, "
                f"phi = {phi:g}. Where SLOPE/W's is higher, XSLOPE's strength is lower.")
        else:
            notes.append(f"material '{name}' takes its friction angle from {label}; "
                         f"it was imported with phi = {phi:g}. Check it before solving.")

    if c_fn is not None:
        label, vals, fn = c_fn
        # XSLOPE's cp is an undrained strength: it takes no friction angle. So only a
        # cohesion line on phi = 0 is cp exactly; with friction it falls to the
        # substitution below.
        line = (_linear_in_y(fn["points"])
                if fn is not None and fn.get("input") == "Y" and phi == 0
                and phi_fn is None else None)
        if line is not None:
            c_top, rate, y_top, y_low, c_low = line
            fields.update(option="cp", c=c_top, cp=rate, r_elev=y_top,
                          _cp_floor=(y_low, c_low))
            notes.append(
                f"material '{name}' takes its cohesion from {label}, a straight line "
                f"rising from {c_top:g} at y = {y_top:g} to {c_low:g} at y = {y_low:g}. "
                f"That is XSLOPE's cp option: c = {c_top:g}, c/p = {rate:g} per unit of "
                f"depth below r_elev = {y_top:g}.")
        elif vals:
            c = min(vals)
            notes.append(
                f"material '{name}' takes its cohesion from {label}, from {min(vals):g} "
                f"to {max(vals):g} at its points. XSLOPE cannot represent that function, "
                f"so the material was imported as Mohr-Coulomb with the lowest, c = "
                f"{c:g}. Where SLOPE/W's is higher, XSLOPE's strength is lower"
                + (f", by up to a factor of {max(vals) / c:.3g}" if c > 0 else "")
                + ". Check it before solving.")
        else:
            notes.append(f"material '{name}' takes its cohesion from {label}; it was "
                         f"imported with c = {c:g}. Check it before solving.")
    if "option" not in fields:
        fields.update(option="mc", c=c)
    fields["phi"] = phi

    if w_fn is not None:
        label, vals, _ = w_fn
        if vals:
            fields["gamma"] = sum(vals) / len(vals)
            notes.append(
                f"material '{name}' takes its unit weight from {label}, from "
                f"{min(vals):g} to {max(vals):g} at its points. XSLOPE has one unit "
                f"weight per material, so it was imported with the average of those "
                f"points, {fields['gamma']:.4g}. Check it before solving.")
        else:
            notes.append(f"material '{name}' takes its unit weight from {label}. "
                         f"Set gamma before solving.")
    return fields, notes


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
        # check_inputs=False: this is an internal viability probe on a partial dict
        # built during import, not a user's run.
        ok, _ = generate_slices(probe, circle=circle, num_slices=12, debug=False,
                                check_inputs=False)
        return bool(ok)
    except Exception:
        return False


# GeoStudio reinforcement Type -> xslope (type, dir, appl). Deliberately NOT
# fileio.REINFORCE_TYPE_PRESETS: these are GeoStudio's conventions for the same four
# words, and they differ (a GeoStudio geosynthetic acts axially, an xslope one
# tangentially). 'axial' acts along the bar, 'tangent' reorients to the slip surface;
# 'passive' means the force is divided by the factor of safety (an ultimate capacity),
# 'active' means it is not.
#
# The file does not record a direction: GeoStudio implies it from the Type. So does
# this table, and each row is measured against SLOPE/W's own factor of safety rather
# than reasoned about -- see run_gsz_corpus(). The F-of-S dependence IS recorded, as
# <FofSDependent>true</FofSDependent> (absent = No), and a dependent product imports
# as 'passive' whatever its Type.
_REINF_TYPE = {
    "Nail":         ("nail",         "axial",   "active"),
    "Geosynthetic": ("geosynthetic", "axial",   "active"),
    "Anchor":       ("anchor",       "axial",   "active"),
    "Tieback":      ("tieback",      "axial",   "active"),
}
# A product with no <Type> is GeoStudio's default type, the Anchor: GeoStudio writes a
# field only when it differs from the default. Every untyped product in SLOPE/W's own
# files is an anchor (they carry the anchor-only BondLength and BondDiameter, and the
# example "Reinforcement with Anchors" types none of its four), and on the Tie-Back
# Wall example an untyped product read as an anchor matches SLOPE/W's factor of safety.
_REINF_DEFAULT_TYPE = "Anchor"


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

    # An anchor develops its pullout only over its bond length, the grouted stretch at
    # the buried end; the free length in front of it carries no bond. SLOPE/W's anchor
    # force where a slip surface crosses is the smallest of the tensile capacity, the
    # bond behind the crossing (rate x the bonded length behind it, at most rate x
    # BondLength) and, at the face, the plate capacity plus any bond in front. Capping
    # the member at rate x BondLength and ramping from the buried end over Lp =
    # capacity / rate reproduces the first two exactly; on SLOPE/W's Reinforcement with
    # Anchors example its column forces are 300 kN at both anchors, rate x BondLength.
    # The face end is exact while the plate capacity is at least that capacity (every
    # anchor in SLOPE/W's files); a smaller plate also earns bond over the free length
    # here. (Nails are bonded over their whole length: their BondLength plays no part.)
    gtype = product.get("type") or _REINF_DEFAULT_TYPE
    p1, p2 = line["p1"], line["p2"]
    if gtype == "Anchor" and _bond_limited(product, line):
        t_max = min(t_max, rate * product["bond_length"])
    lp = (t_max / rate) if rate > 0 else 0.0              # 0 => that end is fully anchored

    # The plate/facing sits at the end that is AT the face — the endpoint on (or nearest)
    # the ground surface. The other end is buried and develops force only by bond.
    if ground_surface is not None and not ground_surface.is_empty:
        from shapely.geometry import Point
        d1 = ground_surface.distance(Point(*p1))
        d2 = ground_surface.distance(Point(*p2))
        face_is_1 = d1 <= d2
    else:
        face_is_1 = True

    xtype, xdir, xappl = _REINF_TYPE.get(gtype, ("", "tangent", "active"))
    if product.get("fofs_dependent"):
        xappl = "passive"
    adhesion, delta = _interface_pullout(product, spacing)
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
        "adhesion": adhesion, "delta": delta,
    }


def _bond_limited(product, line):
    """Does this anchor's bond length stop short of the bar, so that its pullout
    develops over the bond length only? (The test ``_reinforcement`` applies.)"""
    bond = product.get("bond_length")
    p1, p2 = line["p1"], line["p2"]
    return bool(bond and (product.get("pullout") or 0.0) > 0
                and bond < math.hypot(p2[0] - p1[0], p2[1] - p1[1]))


def _interface_pullout(product, spacing):
    """GeoStudio's calculated pullout law -> xslope's (Adhesion, Delta), or (NaN, NaN).

    GeoStudio states the interface strength per unit AREA and multiplies it by a
    length measure xi to reach force per unit length of bar::

        r = xi * (a + sigma'v * tan(delta_gs))

    with xi = SurfaceAreaFactor for a geosynthetic (default 2.0 — soil bearing on
    both faces) and xi = pi * BondDiameter for an anchor or nail. xslope writes the
    two-faced planar case into the law itself::

        r = 2 * (Adhesion + sigma'v * tan(Delta))

    so the two agree exactly when the entered pair carries the factor xi/2::

        Adhesion  = a * xi / 2
        Delta     = atan( tan(delta_gs) * xi / 2 )

    For the SLOPE/W default (a geosynthetic at SAF = 2) that factor is 1 and the
    numbers transfer unchanged; a product that sets SAF = 1, or an anchor whose
    xi = pi*D is anything but 2, is scaled here rather than silently mismatched.
    The out-of-plane Spacing divides it for the same reason it divides every other
    capacity term: xslope stores reinforcement per unit width of slope.
    """
    if not product.get("calculated_pullout"):
        return float("nan"), float("nan")
    a = product.get("adhesion")
    d = product.get("shear_angle")
    if a is None and d is None:
        return float("nan"), float("nan")
    a = float(a or 0.0)
    d = float(d or 0.0)
    if not (0.0 < d < 90.0):
        # An interface with no friction angle is pure adhesion, which xslope's law
        # cannot state on its own: Delta has to be an angle strictly inside (0, 90).
        # Left unmapped, and reported.
        return float("nan"), float("nan")
    gtype = (product.get("type") or "").strip()
    if gtype == "Geosynthetic":
        xi = product.get("surface_area_factor")
        xi = 2.0 if xi is None else float(xi)
    else:                                    # anchor / nail: the bonded perimeter
        dia = product.get("bond_diameter")
        xi = math.pi * float(dia) if dia else 2.0
    factor = xi / (2.0 * (spacing or 1.0))
    delta = math.degrees(math.atan(math.tan(math.radians(d)) * factor))
    return a * factor, delta


_REINF_TYPE_OUT = {x: g for g, (x, _, _) in _REINF_TYPE.items()}   # xslope type -> GeoStudio


def _reinforcement_out(r, gamma_water):
    """One xslope reinforcement line -> a GeoStudio product, as the inverse of the import.

    xslope stores everything PER UNIT WIDTH already, so the product is written with
    ``Spacing`` 1 and the values unscaled. The bond length runs the other way:
    ``lp = t_max / rate`` on the way in, so ``rate = t_max / lp`` on the way out. An
    ``lp`` of zero means "fully anchored at once", which GeoStudio expresses as a pullout
    resistance so large it is never the limit.

    Returns (product_dict, [caveat, ...]) -- the caveats being everything the GeoStudio
    model cannot hold.
    """
    notes = []
    t_max = float(r.get("t_max") or 0.0)
    lp1, lp2 = float(r.get("lp1") or 0.0), float(r.get("lp2") or 0.0)
    tend1, tend2 = float(r.get("tend1") or 0.0), float(r.get("tend2") or 0.0)

    lp = max(lp1, lp2)
    if lp1 > 0 and lp2 > 0 and abs(lp1 - lp2) > 1e-9 * max(lp1, lp2):
        notes.append(
            f"reinforcement {r.get('label') or ''!r} has a different bond length at each "
            f"end ({lp1:g} and {lp2:g}); GeoStudio has ONE pullout resistance per product, "
            f"so the longer was written — its capacity develops more slowly than xslope's")
    rate = (t_max / lp) if lp > 0 else t_max * 1.0e6      # 0 => never the limit

    plate = max(tend1, tend2)
    if tend1 > 0 and tend2 > 0:
        notes.append(
            f"reinforcement {r.get('label') or ''!r} has an end capacity at BOTH ends "
            f"({tend1:g} and {tend2:g}); GeoStudio has one plate capacity, so the larger "
            f"was written")

    xtype = (r.get("type") or "").lower()
    gtype = _REINF_TYPE_OUT.get(xtype)
    if gtype is None:
        gtype = "Nail"
        notes.append(
            f"reinforcement {r.get('label') or ''!r} has type {r.get('type')!r}, which "
            f"GeoStudio does not have — written as a Nail, which is axial and unfactored; "
            f"check the force it develops")

    # GeoStudio does not STORE a direction or an F-of-S dependence: it implies both from
    # the type. So a line whose dir/appl disagree with the convention for its type cannot
    # be expressed, and would behave differently in GeoStudio without saying so.
    want_dir, want_appl = _REINF_TYPE.get(gtype, (None, "axial", "active"))[1:]
    have_dir = (r.get("dir") or want_dir).lower()
    have_appl = (r.get("appl") or want_appl).lower()
    if (have_dir, have_appl) != (want_dir, want_appl):
        notes.append(
            f"reinforcement {r.get('label') or ''!r} is {have_dir}/{have_appl} in xslope, "
            f"but GeoStudio has no field for either — it infers {want_dir}/{want_appl} "
            f"from the type '{gtype}'. THE FORCE WILL DIFFER; re-check it in GeoStudio")

    # The overburden pullout law. GeoStudio states the same law per unit AREA and
    # multiplies it by a length measure -- SurfaceAreaFactor for a geosynthetic,
    # pi*BondDiameter for an anchor or nail -- and its manual restricts the
    # calculated law to anchors and geosynthetics. Writing SAF = 2 makes GeoStudio's
    # rate 2*(a + sigma'v*tan delta), which is xslope's, so the pair transfers
    # unchanged. A line on the law whose type is not a geosynthetic keeps its
    # constant rate and says so.
    adhesion, delta = r.get("adhesion"), r.get("delta")
    interface = None
    if (adhesion is not None and delta is not None
            and float(adhesion) == float(adhesion) and float(delta) == float(delta)):
        if gtype == "Geosynthetic":
            interface = (float(adhesion), float(delta))
        else:
            notes.append(
                f"reinforcement {r.get('label') or ''!r} develops its pullout from the "
                f"effective overburden (Adhesion/Delta), which GeoStudio offers for "
                f"geosynthetics and anchors but not for a '{gtype}'; the constant "
                f"pullout rate was written instead and the developed force will differ")

    return {"type": gtype, "name": str(r.get("label") or gtype), "spacing": 1.0,
            "tensile": t_max, "plate": plate, "pullout": rate,
            "interface": interface}, notes


def _line_load_out(ll):
    """xslope line-load angle (CCW from +x) -> GeoStudio trend/plunge.

    Plunge is the angle BELOW horizontal; trend picks which way along x. The inverse of
    what the reader does, and anchored to a real file: the Amherst wall's facing load is
    Trend=0 / Plunge=90, i.e. straight down, which is xslope's -90.
    """
    ang = math.radians(float(ll.get("angle", -90.0)))
    dx, dy = math.cos(ang), math.sin(ang)
    plunge = math.degrees(math.atan2(-dy, abs(dx)))     # >0 is downward
    trend = 0.0 if dx >= 0 else 180.0
    return trend, plunge


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
        # v16 tensile cutoff and v17 matric-suction strength (None = unset, the pre-
        # v16/v17 default), mirroring load_slope_data's material shape.
        "t_cut": None, "phi_b": None, "s_cap": None,
        "sigma_gamma": 0.0, "sigma_c": 0.0, "sigma_phi": 0.0, "sigma_cp": 0.0,
        "sigma_d": 0.0, "sigma_psi": 0.0, "k1": 0.0, "k2": 0.0, "alpha": 0.0,
        "unsat": "lf", "kr0": 0.0, "h0": 0.0, "vg_a": 0.0, "vg_n": 0.0,
        "vg_l": 0.5,
        "E": 0.0, "nu": 0.0,
    }


#: GeoStudio's probabilistic parameter names -> the xslope material sigma they are.
#: The Total* names belong to a total-stress material (the undrained one keeps its
#: strength in <Cohesion>), which _strength imports into the SAME c/phi slots, so its
#: standard deviation belongs in the same sigma slot.
_PROB_SIGMA = {
    "Applied.UnitWeight": "sigma_gamma",
    "Applied.EffectiveCohesion": "sigma_c",
    "Applied.TotalCohesion": "sigma_c",
    "Applied.EffectivePhi": "sigma_phi",
    "Applied.TotalPhi": "sigma_phi",
}


def _apply_prob_modifiers(analysis, mat_index, materials, bedrock, caveats):
    """Import a probabilistic analysis's standard deviations onto the materials.

    SLOPE/W states uncertainty as a distribution of the DEVIATION from the applied
    value -- ``Normal(Mean=0,SD=1150,...)`` -- so its SD is exactly xslope's sigma for
    that parameter, in the file's own units, and needs no conversion. Dropping these
    was silent and total: a probabilistic model arrived with every sigma zero, so
    ``reliability()`` refused it for having no standard deviations at all.

    Everything SLOPE/W can vary that xslope's sigma set has no slot for (Ru, load and
    reinforcement modifiers), and every distribution that is not normal, is reported
    rather than approximated.
    """
    if (analysis.get("sampling") or "") != "Probabilistic":
        return
    mods = analysis.get("modifiers") or []
    imported, dropped, truncated, correlated = [], [], [], []

    for mod in mods:
        mat_id, prop = mod.get("mat_id"), mod.get("prop")
        if mat_id is None:
            dropped.append(f"{mod['key']} (not a material parameter)")
            continue
        if mat_id not in mat_index:
            continue                       # a material this analysis does not use
        mat = materials[mat_index[mat_id]]
        name = mat["name"]
        field = _PROB_SIGMA.get(prop)
        if field is None:
            dropped.append(f"'{name}' {prop} (no xslope equivalent)")
            continue
        if mat_id in bedrock and field != "sigma_gamma":
            # A Bedrock material's c/phi are not imported at all -- it cannot fail --
            # so a standard deviation on them has nothing to attach to.
            dropped.append(f"'{name}' {prop} (impenetrable Bedrock has no strength)")
            continue
        if mod.get("dist") != "Normal":
            dropped.append(f"'{name}' {prop} ({mod.get('dist') or 'unreadable'} "
                           f"distribution; xslope's reliability is normal/lognormal "
                           f"about the mean)")
            continue
        sd = mod["params"].get("SD")
        if not isinstance(sd, float) or not sd > 0:
            dropped.append(f"'{name}' {prop} (no usable standard deviation)")
            continue
        mean = mod["params"].get("Mean", 0.0)
        if isinstance(mean, float) and mean != 0.0:
            dropped.append(f"'{name}' {prop} mean shift of {mean:g} (xslope varies "
                           f"about the entered value)")
        mat[field] = sd
        imported.append(f"'{name}' {prop} -> {field} = {sd:g}")
        if any(k in mod["params"] for k in ("Min", "Max")):
            truncated.append(f"'{name}' {prop}")
        if mod.get("corr_tag") and mod.get("corr"):
            correlated.append(f"'{name}' {prop} with {mod['corr_tag']} "
                              f"(r = {mod['corr']:g})")

    trials = analysis.get("trials")
    if imported:
        caveats.append(
            f"this is a PROBABILISTIC SLOPE/W analysis"
            + (f" ({int(trials)} Monte Carlo trials)" if trials else "")
            + f" — {len(imported)} standard deviation(s) were imported into the mat "
              f"sheet: {', '.join(imported)}. Run reliability() (or reliability_mc) "
              f"to use them; a deterministic solve ignores them")
    elif mods:
        caveats.append(
            "this is a PROBABILISTIC SLOPE/W analysis, but none of its standard "
            "deviations could be imported (see below) — the model arrives "
            "deterministic")
    if truncated:
        caveats.append(
            f"SLOPE/W truncates the sampled range for {', '.join(truncated)} (a Min/Max "
            f"on the distribution) — xslope has no truncation, so its sampled tails are "
            f"wider than SLOPE/W's")
    if correlated:
        caveats.append(
            f"correlated parameters {', '.join(correlated)} — xslope's reliability "
            f"treats material parameters as independent, so the correlation is not "
            f"reproduced")
    if dropped:
        caveats.append(
            f"probabilistic parameter(s) not imported: {', '.join(dropped)}")


def _import_seep(gsz, analysis, step, caveats, polygons):
    """Pick a time step and read its SEEP/W field. Returns (mesh, u, step) or None."""
    name = analysis["name"]
    steps = gsz_seep_steps(gsz, name)
    if not steps:
        return None

    if step is None and len(steps) > 1:
        # An xslope model is ONE state in time; SLOPE/W solved every step, and the
        # factor of safety differs at each. Default to the step SLOPE/W itself found
        # most critical -- the condition a drawdown analysis exists to find. Picking
        # first or last instead would be arbitrary, and would quietly answer a
        # different question than the one the file was built to ask.
        scored = []
        for s in steps:
            trials = read_gsz_results(gsz, name, s["step"])
            if trials:
                scored.append((min(t["fs"] for t in trials), s))
        if scored:
            worst_fs, worst = min(scored, key=lambda p: p[0])
            step = worst["step"]
            table = ", ".join(f"{s['step']}={fs:.3f}" for fs, s in scored)
            caveats.append(
                f"this analysis is transient: SLOPE/W solved it at {len(steps)} saved "
                f"time steps, with factors of safety {table}. An xslope model is one "
                f"state in time, so step {step} was imported — the GOVERNING one, where "
                f"SLOPE/W's own factor of safety is lowest ({worst_fs:.3f}) at t = "
                f"{worst['time']:g}. Re-import with step='NNN' for any of the others")

    got = read_gsz_seep(gsz, name, step, polygons)
    if got is None:
        return None
    mesh, u, chosen = got
    caveats.append(
        f"pore pressure came from the parent SEEP/W analysis as a finite element field "
        f"({len(mesh['nodes'])} nodes, step {chosen}) — imported onto the mesh SEEP/W "
        f"solved on, with every material set to the 'seep' option. It is SLOPE/W's own "
        f"pore pressure. xslope did not run the seepage analysis")
    return mesh, u, chosen


def _piezo_short_ends(line, ground_surface, domain_polygon, which="piezometric surface",
                      solved=True):
    """The import note for a piezometric line that stops short of the model's edges,
    or None when it spans them. The line itself is imported as drawn.

    SLOPE/W applies no pore pressure outside a line's x-range (not documented; its
    saved slice results read exactly zero there). XSLOPE has no such rule: its slicer
    refuses a slice past a line's end, and the preflight warns (piezo.extent_short).
    The import does not invent the missing water. The note says what SLOPE/W did and,
    for each short end, what the end's position suggests the user do: an end on the
    ground surface is a water table coming out on the ground; one above the ground at
    the model's edge is a pool, which continues level; one below the ground goes on
    as the water does, which only the user knows.
    """
    if len(line) < 2 or domain_polygon is None or domain_polygon.is_empty:
        return None
    x0m, y0m, x1m, y1m = domain_polygon.bounds
    pts = sorted(line)
    tol_x = 1e-6 * max(1.0, x1m - x0m)
    # "On the ground": within 0.1% of the model's size. Walter Bouldin's line ends
    # 0.08 ft above the ground on a model 216 ft wide; the other short ends in
    # SLOPE/W's files are tens of feet off it.
    tol_y = 1e-3 * max(1.0, x1m - x0m, y1m - y0m)
    ends = []
    for side, (x, y), edge in (("left", pts[0], x0m), ("right", pts[-1], x1m)):
        if (side == "left" and x <= edge + tol_x) or (side == "right" and x >= edge - tol_x):
            continue
        past = f"{side} of x = {x:g}"
        yg = _y_on(ground_surface, x)
        if yg is not None and abs(y - yg) <= tol_y:
            advice = (f"At its {side} end the line is on the ground surface: the water "
                      f"table comes out on the ground there. If the water continues, "
                      f"extend the line along the ground surface to x = {edge:g}")
        elif yg is not None and y > yg:
            advice = (f"At its {side} end the line is {y - yg:.3g} above the ground: a "
                      f"pool. Extend the line level to x = {edge:g}")
        else:
            below = "" if yg is None else f"{yg - y:.3g} "
            advice = (f"At its {side} end the line is {below}below the ground. Extend "
                      f"it to x = {edge:g} as the water continues")
        ends.append((past, advice))
    if not ends:
        return None
    pasts = " and ".join(p for p, _ in ends)
    return (
        f"the {which} stops short of the model: it runs from x = {pts[0][0]:g} to "
        f"{pts[-1][0]:g}, and the model from x = {x0m:g} to {x1m:g}. SLOPE/W applies "
        f"no pore pressure outside a line's x-range; XSLOPE refuses a slice there "
        f"instead. "
        + (f"This file's SLOPE/W results were computed with no water {pasts}. "
           if solved else f"SLOPE/W would compute this file with no water {pasts}. ")
        + ". ".join(a for _, a in ends))


def _circle_past_line(trial, lines, domain_polygon):
    """Does the sliding mass of SLOPE/W's trial circle ``trial`` (xo, yo, r) reach
    past either end of any of ``lines``, where XSLOPE's slicer refuses a slice?"""
    if not lines or domain_polygon is None or domain_polygon.is_empty:
        return False
    from shapely.geometry import Point
    mass = Point(trial["xo"], trial["yo"]).buffer(trial["r"], 256).intersection(
        domain_polygon)
    if mass.is_empty:
        return False
    x0, _, x1, _ = mass.bounds
    tol = 1e-6 * max(1.0, x1 - x0)
    return any(x0 < min(p[0] for p in ln) - tol or x1 > max(p[0] for p in ln) + tol
               for ln in lines if len(ln) >= 2)


def _line_on_ground(line, ground_surface, n=20):
    """The x ranges, as ``[(x0, x1), ...]``, where a SLOPING stretch of a piezometric
    line runs along the ground surface -- where SLOPE/W's phreatic correction leaves
    the head alone and an XSLOPE 'phreatic' line does not (Greco's weak-layer slope,
    the one verification file with the correction on). A level stretch is skipped:
    cos^2 of a zero slope is 1, so the two agree there anyway. No file shows what
    SLOPE/W does where a sloping line is above the ground, so that is not claimed."""
    if ground_surface is None or ground_surface.is_empty or len(line) < 2:
        return []
    gx0, _, gx1, _ = ground_surface.bounds
    tol = 1e-6 * max(1.0, gx1 - gx0)
    spans = []
    for (xa, ya), (xb, yb) in zip(line, line[1:]):
        if abs(xb - xa) < 1e-12 or abs(yb - ya) < 1e-12:
            continue
        on = []
        for k in range(n + 1):
            x = xa + (xb - xa) * k / n
            yg = _y_on(ground_surface, x)
            if yg is not None and abs(ya + (yb - ya) * k / n - yg) <= tol:
                on.append(x)
        if on:
            lo, hi = min(on), max(on)
            if spans and lo <= spans[-1][1] + tol:
                spans[-1] = (spans[-1][0], max(spans[-1][1], hi))
            else:
                spans.append((lo, hi))
    return [(a, b) for a, b in spans if b - a > tol]


def _kc1_envelope(c_r, phi_r, phi):
    """XSLOPE's K_c = 1 envelope ``(d, psi)`` from SLOPE/W's drawdown strengths.

    A SLOPE/W material's <DrawdownCohesion> and <DrawdownPhi> are c_R and phi_R, the
    intercept and slope of the R (total stress) envelope drawn tangent to the Mohr
    circles. SLOPE/W turns them into the isotropic-consolidation envelope by its
    Eqs. 77-78 (Slope Stability Modeling, sec. 9.3), the same as Duncan, Wright and
    Brandon's (2014) Eqs. 9.6-9.7:

        d   = c_R cos(phi_R) cos(phi') / (1 - sin(phi_R))
        psi = atan( sin(phi_R) cos(phi') / (1 - sin(phi_R)) )

    That envelope is XSLOPE's d and psi. On SLOPE/W's four rapid drawdown examples,
    scored on SLOPE/W's own critical circles, the converted pair lands within 1.9% of
    SLOPE/W's factor of safety; c_R and phi_R taken as d and psi unchanged land 4-10%
    low. ``None`` when phi_R or phi' is 90 degrees or more, where the equations fail.
    """
    if not (phi_r < 90.0 and phi < 90.0):
        return None
    sr, cr = math.sin(math.radians(phi_r)), math.cos(math.radians(phi_r))
    cp = math.cos(math.radians(phi))
    return (c_r * cr * cp / (1.0 - sr),
            math.degrees(math.atan(sr * cp / (1.0 - sr))))


def _r_envelope(d, psi, phi):
    """The inverse of :func:`_kc1_envelope`: the ``(c_R, phi_R)`` that SLOPE/W converts
    back to exactly this ``(d, psi)``. With t = tan(psi) / cos(phi'), Eq. 78 gives
    sin(phi_R) = t / (1 + t), and Eq. 77 then gives c_R. ``None`` when phi' is 90
    degrees or more, or psi is negative."""
    if not (0.0 <= psi < 90.0 and phi < 90.0):
        return None
    cp = math.cos(math.radians(phi))
    t = math.tan(math.radians(psi)) / cp
    sr = t / (1.0 + t)
    phi_r = math.degrees(math.asin(sr))
    return d * (1.0 - sr) / (math.cos(math.radians(phi_r)) * cp), phi_r


def _y_along(line, x):
    """The elevation of a polyline at ``x``, or None outside its x-extent."""
    for (xa, ya), (xb, yb) in zip(line, line[1:]):
        lo, hi = min(xa, xb), max(xa, xb)
        if lo <= x <= hi:
            return ya if xb == xa else ya + (yb - ya) * (x - xa) / (xb - xa)
    return None


def _before_after(a, b):
    """Which of two piezometric surfaces is the water before drawdown.

    Returns ``(high, low, crossing)``: ``high`` stands at or above ``low`` wherever
    both are defined unless ``crossing`` is the x of a point where the order flips, in
    which case ``high`` is the one that is higher on average over their common
    extent."""
    xs = sorted({x for x, _ in a} | {x for x, _ in b})
    diffs = []
    for x in xs:
        ya, yb = _y_along(a, x), _y_along(b, x)
        if ya is not None and yb is not None:
            diffs.append((x, ya - yb))
    if not diffs:                       # no common extent: compare their averages
        ma = sum(y for _, y in a) / len(a)
        mb = sum(y for _, y in b) / len(b)
        return (a, b, None) if ma >= mb else (b, a, None)
    ys = [y for _, y in a + b]
    tol = 1e-9 * max(1.0, max(ys) - min(ys))
    signs = [(x, 1 if dy > tol else -1) for x, dy in diffs if abs(dy) > tol]
    if all(s > 0 for _, s in signs):
        return a, b, None
    if all(s < 0 for _, s in signs):
        return b, a, None
    flip = next(x for (x, s), (_, s0) in zip(signs[1:], signs) if s != s0)
    mean = sum(dy for _, dy in diffs) / len(diffs)
    return (a, b, flip) if mean >= 0 else (b, a, flip)


def _rapid_drawdown(gsz, used, materials, bedrock, stab, option, piezo_line):
    """Map a SLOPE/W staged rapid drawdown onto XSLOPE's own.

    SLOPE/W (after Duncan, Wright and Wong 1990) reads two piezometric surfaces, the
    higher as the water before drawdown and the lower as the water after it, and an R
    envelope (c_R, phi_R) on each soil that does not drain; a free-draining soil has
    neither. XSLOPE's rapid drawdown takes the same three things: piezo line 1, piezo
    line 2, and the K_c = 1 envelope (d, psi) on each Mohr-Coulomb material that does
    not drain, converted from c_R and phi_R by SLOPE/W's own equations
    (:func:`_kc1_envelope`).

    Sets ``d`` and ``psi`` on ``materials`` in place and returns ``(line1, line2,
    caveats)``. Anything that is not simply two surfaces, one above the other, with
    drawdown strengths on Mohr-Coulomb soils, is reported.
    """
    caveats, soils, free = [], [], []
    for mid, mat in zip(used, materials):
        src = gsz["materials"][mid]
        s = src["strength"]
        try:
            c_r = float(s.get("DrawdownCohesion") or 0.0)
            phi_r = float(s.get("DrawdownPhi") or 0.0)
        except ValueError:
            caveats.append(
                f"material '{src['name']}' has drawdown strengths XSLOPE cannot read "
                f"(c_R = {s.get('DrawdownCohesion')!r}, phi_R = "
                f"{s.get('DrawdownPhi')!r}), so it was imported as free-draining. Set "
                f"its d and psi by hand")
            continue
        if mid in bedrock:
            continue
        if not (c_r or phi_r):
            free.append(repr(src["name"]))
            continue
        if src["model"] != "MohrCoulomb":
            caveats.append(
                f"material '{src['name']}' carries drawdown strengths (c_R = {c_r:g}, "
                f"phi_R = {phi_r:g}), but SLOPE/W's staged rapid drawdown reads them only "
                f"on a Mohr-Coulomb soil, and this one is {src['model'] or 'unnamed'}. "
                f"They were not imported")
            continue
        env = _kc1_envelope(c_r, phi_r, mat["phi"])
        if env is None:
            caveats.append(
                f"material '{src['name']}' has drawdown strengths c_R = {c_r:g}, "
                f"phi_R = {phi_r:g} with phi' = {mat['phi']:g}, which SLOPE/W's "
                f"conversion cannot take, so it was imported as free-draining. Set its "
                f"d and psi by hand")
            continue
        mat["d"], mat["psi"] = env
        soils.append(f"'{src['name']}' c_R = {c_r:g}, phi_R = {phi_r:g} -> "
                     f"d = {env[0]:.4g}, psi = {env[1]:.3g}")
        if (env[0] == 0.0) != (env[1] == 0.0):
            caveats.append(
                f"material '{src['name']}' has c_R = {c_r:g} and phi_R = {phi_r:g}, so "
                f"one of its d and psi is zero. XSLOPE reads a zero there as a blank "
                f"and refuses a rapid drawdown on a material with only one of the "
                f"two, so this model will not run as one until that is resolved")

    # The water before and after drawdown.
    line1, line2, how = list(piezo_line), [], None
    surfaces = stab.get("piezo", {})
    lines = [(sid, list(pts)) for sid, pts in surfaces.items() if len(pts) >= 2]
    if option not in ("PiezoSurface", "PiezoSurfaceBbar"):
        caveats.append(
            f"SLOPE/W's staged rapid drawdown here takes its pore pressure from its "
            f"'{option}' option. XSLOPE's takes the water before and after drawdown "
            f"as piezo lines 1 and 2, and this file states neither, so set both before "
            f"running XSLOPE's rapid drawdown")
    elif len(lines) < 2:
        caveats.append(
            f"SLOPE/W's staged rapid drawdown reads two piezometric surfaces, the water "
            f"before and after drawdown, and this analysis has {len(lines)}. "
            + ("It was imported as piezo line 1. " if lines else "")
            + "Draw piezo line 2, the water after drawdown, before running XSLOPE's "
              "rapid drawdown")
    else:
        if len(lines) > 2:
            mean = {sid: sum(y for _, y in pts) / len(pts) for sid, pts in lines}
            hi_id = max(mean, key=mean.get)
            lo_id = min(mean, key=mean.get)
            caveats.append(
                f"SLOPE/W's staged rapid drawdown reads two piezometric surfaces, and "
                f"this analysis has {len(lines)}. XSLOPE took the highest (ID {hi_id}) "
                f"as the water before drawdown and the lowest (ID {lo_id}) as the water "
                f"after it; the others were not imported")
            lines = [(sid, pts) for sid, pts in lines if sid in (hi_id, lo_id)]
        (ida, a), (idb, b) = lines
        high, low, crossing = _before_after(a, b)
        hi_id, lo_id = (ida, idb) if high is a else (idb, ida)
        line1, line2 = high, low
        if crossing is not None:
            caveats.append(
                f"SLOPE/W's two piezometric surfaces cross near x = {crossing:g}, so "
                f"neither is the water before drawdown everywhere. XSLOPE took surface "
                f"{hi_id}, the higher on average, as piezo line 1 (before drawdown) and "
                f"surface {lo_id} as piezo line 2. Check them against the GeoStudio "
                f"model")
        uses = stab.get("material_piezo", {})
        lowered = [repr(gsz["materials"][mid]["name"]) for mid in used
                   if uses.get(mid) == lo_id and mid not in bedrock]
        if lowered:
            caveats.append(
                f"SLOPE/W assigns {', '.join(lowered)} to the lowered piezometric "
                f"surface. In XSLOPE's rapid drawdown every material that takes its "
                f"pore pressure from a line reads line 1 before drawdown and line 2 "
                f"after it, so check that material's water against the GeoStudio model")
        how = ("the full-pool piezometric surface as piezo line 1, the lowered one "
               + (f"({', '.join(f'({x:g}, {y:g})' for x, y in low)}) "
                  if len(low) <= 12 else "")
               + "as piezo line 2")

    caveats.insert(0, (
        "SLOPE/W solved this analysis as a staged rapid drawdown (Duncan, Wright and "
        "Wong 1990), and it was imported as XSLOPE's rapid drawdown"
        + (f": {how}" if how else "")
        + (f"{', and' if how else ':'} the drawdown strengths as d and psi, converted "
           f"from SLOPE/W's R envelope by its own equations ({'; '.join(soils)})"
           if soils else
           ". No soil carries drawdown strengths, so every soil drains and XSLOPE's "
           "rapid drawdown refuses the model until one is given d and psi")
        + (f". {', '.join(free)} carr{'ies' if len(free) == 1 else 'y'} none and "
           f"drain{'s' if len(free) == 1 else ''} freely" if free and soils else "")
        + ". An XSLOPE file does not record which analysis to run, so run it as a rapid "
          "drawdown (in Studio, Run LEM with Rapid drawdown ticked; in a script, "
          "rapid=True). Run as an ordinary analysis it is the slope before drawdown, "
          "and its factor of safety is not SLOPE/W's"))
    return line1, line2, caveats


# SLOPE/W's partial factor keys, in words.
_FACTOR_WORDS = {
    "PermanentLoads": "permanent loads", "VariableLoads": "variable loads",
    "SoilWeight": "soil weight", "EffectiveCohesion": "c'",
    "EffectiveTanPhi": "tan phi'", "UndrainedStrength": "undrained strength",
    "ShearStrengthOther": "other shear strengths",
    "PulloutResistance": "pullout resistance",
    "StructuralCapacities": "structural capacities",
    "CompressiveStrength": "compressive strength",
    "SeismicCoefficients": "seismic coefficients", "EarthResistance": "earth resistance",
}


def _design_factor_caveat(gsz, dfid):
    """The caveat for an analysis SLOPE/W solves with a set of partial factors.

    A set (<DesignFactors>, named by the analysis's <Control><DesignFactorID>) gives
    each factor a favorable and an unfavorable value: the strengths are divided by
    theirs, the loads multiplied. The caveat names every factor that is not 1."""
    try:
        ds = gsz.get("design_factors", {}).get(int(dfid))
    except (TypeError, ValueError):
        ds = None
    if ds is None:
        return (f"SLOPE/W applies partial factor set {dfid} to this analysis, which the "
                f"file does not define. XSLOPE applies no partial factors, so its factor "
                f"of safety is not SLOPE/W's factored one")
    parts = []
    for key, fav, unfav in ds["factors"]:
        vals = []
        for v, side in ((fav, "favorable"), (unfav, "unfavorable")):
            try:
                x = float(v)
            except ValueError:
                continue
            if abs(x - 1.0) > 1e-12:
                vals.append((x, side))
        if not vals:
            continue
        word = _FACTOR_WORDS.get(key, key)
        if len(vals) == 1 and vals[0][1] == "favorable":
            parts.append(f"{word} {vals[0][0]:g}")
        else:
            parts.append(f"{word} " + " and ".join(f"{x:g} ({s})" for x, s in vals))
    if not parts:
        return (f"SLOPE/W applies the partial factor set '{ds['name']}' to this "
                f"analysis, and every factor in it is 1, so the model is unchanged by it")
    return (f"SLOPE/W applies the partial factors of '{ds['name']}' to this analysis: "
            f"{', '.join(parts)}. XSLOPE applies no partial factors, so the strengths "
            f"and loads were imported as the file states them, and XSLOPE's factor of "
            f"safety is not SLOPE/W's factored one. To compare the two, divide each "
            f"strength by its factor and scale each load by its own before solving")


def _unit_weight_above_wt(gsz, used, materials):
    """Apply SLOPE/W's unit weight above the water table to ``materials``; the caveats.

    A SLOPE/W material weighs <UnitWeight> everywhere unless its
    <SlopeUseUnitWeightAboveWT> switch says otherwise; a <DryWeight> without the
    switch is a value SLOPE/W keeps but does not use (its own pumped-storage rapid
    drawdown example weighs its rockfill at <UnitWeight> 142 in every saved column,
    its <DryWeight> 128 notwithstanding). The switch has two settings:

    - ``ConstantValue``: <DryWeight> above the water table and <UnitWeight> below it,
      which is XSLOPE's gamma and gamma_sat. On SLOPE/W's own Unit Weight Definition
      example (Case 2, 16 above and 20 below), the split puts XSLOPE's weight within
      0.07% of SLOPE/W's over the critical slip surface, column by column within 0.2%
      wherever up to 46% of a column is below the water table.
    - ``CalculateFromVolFn``: the weight above the water table from the soil's
      volumetric water content. XSLOPE has no such weight, so the material keeps
      <UnitWeight> throughout and is reported.
    """
    caveats, split = [], []
    for mid, mat in zip(used, materials):
        s = gsz["materials"][mid]["strength"]
        how = (s.get("SlopeUseUnitWeightAboveWT") or "").strip()
        if not how:
            continue
        if how == "ConstantValue":
            try:
                dry = float(s.get("DryWeight"))
            except (TypeError, ValueError):
                caveats.append(
                    f"material '{mat['name']}' takes a separate unit weight above the "
                    f"water table, but the file gives none XSLOPE can read, so it "
                    f"weighs {mat['gamma']:g} throughout")
                continue
            if abs(dry - mat["gamma"]) > 1e-12:
                mat["gamma_sat"], mat["gamma"] = mat["gamma"], dry
                split.append(f"'{mat['name']}' {dry:g} above and {mat['gamma_sat']:g} below")
        elif how == "CalculateFromVolFn":
            caveats.append(
                f"material '{mat['name']}': SLOPE/W computes its unit weight above the "
                f"water table from its volumetric water content. XSLOPE has one unit "
                f"weight above the water table, so the material weighs {mat['gamma']:g}, "
                f"its saturated unit weight, throughout, and XSLOPE's slices above the "
                f"water table are heavier than SLOPE/W's")
        else:
            caveats.append(
                f"material '{mat['name']}' sets its unit weight above the water table "
                f"by '{how}', which XSLOPE does not know, so it weighs "
                f"{mat['gamma']:g} throughout")
    if split:
        caveats.insert(0, (
            f"SLOPE/W gives {'this material' if len(split) == 1 else 'these materials'} "
            f"a separate unit weight above the water table: {'; '.join(split)}. It was "
            f"imported as gamma above the water table and gamma_sat below it, split where "
            f"XSLOPE finds the water table"))
    return caveats


def gsz_to_slope_data(gsz, analysis_id=None, critical_surface=True, step=None):
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
    step : str, optional
        For an analysis whose water comes from a *transient* SEEP/W run: which
        saved time step to import ('000', '001', …). An xslope model is one state
        in time, but SLOPE/W solves the slope at every step and the factor of
        safety differs at each — in a drawdown problem, by a lot. The default is
        the **governing step**, the one where SLOPE/W's own factor of safety is
        lowest, which is the condition the analysis exists to find. The choice is
        always reported. See :func:`gsz_seep_steps` for what is on offer.

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

    if not gsz["analyses"]:
        raise ValueError("This .gsz defines no analyses.")
    if analysis_id is None:
        analysis_id = gsz["analyses"][0]["id"]
    analysis = next((a for a in gsz["analyses"] if a["id"] == analysis_id), None)
    if analysis is None:
        raise ValueError(f"This .gsz has no analysis with ID {analysis_id}.")

    # Resolve THIS analysis's geometry. A multi-geometry .gsz (e.g. the same embankment
    # imported at two heights) stores one <Geometry> per analysis and the analysis names
    # its own by GeometryId; because the geometries reuse point and region IDs, only the
    # analysis's own geometry carries the right coordinates. Fall back to the sole/last
    # geometry for an older file that names none.
    geom = gsz["geometries"].get(analysis.get("geometry_id"))
    if geom is None:
        # No GeometryId named (older single-geometry files): fall back to the flat
        # merged tables, which are exactly what the importer used before it became
        # geometry-aware — so a single-geometry import is byte-for-byte unchanged.
        geom = {"points": gsz["points"], "regions": gsz["regions"]}
    pts = geom["points"]

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

    # GeoStudio's Bedrock is IMPENETRABLE: slip surfaces cannot enter it. xslope's v16
    # 'elastic' option is the exact analog — an impenetrable zone the LEM cannot cut and
    # the FEM holds out of plasticity — so a Bedrock material is imported as elastic and
    # its region KEPT (earlier imports dropped the region and floored the domain at the
    # bedrock top; keeping it as elastic reproduces SLOPE/W just as faithfully and lets
    # the material round-trip). Importing it as ordinary soil would let trial surfaces
    # cut straight through the bedrock.
    bedrock = {mid for mid, m in gsz["materials"].items() if m["model"] == "Bedrock"}

    # Import only the materials this analysis actually uses, renumbered 0-based.
    used = sorted(set(assign.values()))
    mat_index = {mid: i for i, mid in enumerate(used)}
    materials = []
    depth_from_top = []     # S=f(depth) materials, whose r_elev is their layer's top
    spatial_cp = []         # (index, (y, c)): a cp cohesion SLOPE/W holds below y
    for mid in used:
        src = gsz["materials"].get(mid)
        if src is None:
            raise ValueError(f"This analysis references material {mid}, which the "
                             f"file does not define.")
        mat = _blank_material(src["name"])
        if mid in bedrock:
            # Impenetrable Bedrock -> elastic. Keep the unit weight (an elastic zone
            # still has weight); its c/phi are meaningless because it cannot fail.
            mat.update(option="elastic", gamma=src["gamma"])
            caveats.append(
                f"material '{src['name']}' is GeoStudio's impenetrable Bedrock — "
                f"imported as an elastic material: it cannot fail and a slip surface "
                f"cannot cut through it, which matches SLOPE/W")
            materials.append(mat)
            if not mat["gamma"]:
                caveats.append(f"material '{src['name']}' has no unit weight in the "
                               f"file — set gamma before solving")
            continue
        fields, notes, from_top = _strength(src, gsz.get("mfns"), gsz.get("sfns"),
                                           gsz.get("fn3ds"))
        caveats.extend(notes)
        cp_floor = fields.pop("_cp_floor", None)
        mat.update(gamma=src["gamma"])
        mat.update(fields)
        if cp_floor is not None:
            spatial_cp.append((mat_index[mid], cp_floor))
        if from_top:
            depth_from_top.append(mat_index[mid])
        # SLOPE/W can scale any strength by a function of the base inclination (its
        # anisotropic modifier). XSLOPE has no such scaling.
        aniso = src["strength"].get("ShearAnisoFnNum")
        if aniso is not None:
            caveats.append(
                f"material '{src['name']}' scales its strength with the slope of the "
                f"slice base (SLOPE/W anisotropy function {aniso}). XSLOPE does not, "
                f"so its strength is the unscaled one.")

        # v17 matric-suction strength. SLOPE/W's own unsaturated model is a threshold on
        # the piezometric surface (CapSuction/MaxSuction), reported below — not a per-
        # material phi_b. A few SLOPE/W builds DO write a material <PhiB>, so read it if
        # present (none appears in the reference files, so this path is format-doc-only).
        phib = src["strength"].get("PhiB")
        if phib is not None:
            try:
                mat["phi_b"] = float(phib)
                caveats.append(
                    f"material '{src['name']}' carries a GeoStudio unsaturated friction "
                    f"angle (phi_b = {float(phib):g}) — imported as phi_b; s_cap left "
                    f"unset")
            except (TypeError, ValueError):
                pass

        if not mat["gamma"]:
            caveats.append(f"material '{src['name']}' has no unit weight in the "
                           f"file — set gamma before solving")
        # A ratio or cp material with c = 0 still has strength through c/p.
        if not mat["c"] and not mat["phi"] and not (
                mat["option"] in ("cp", "ratio") and mat["cp"]):
            caveats.append(
                f"material '{src['name']}' came across with NO STRENGTH (c = 0, "
                f"phi = 0) — it will behave as a liquid; set its strength before "
                f"solving or the factor of safety will be meaningless")
        materials.append(mat)

    # Standard deviations, if this analysis is probabilistic. Applied after the loop
    # because a modifier names its material by ID, not by position in `used`.
    _apply_prob_modifiers(analysis, mat_index, materials, bedrock, caveats)

    if used and all(mid in bedrock for mid in used):
        caveats.append(
            "every material in this analysis is impenetrable Bedrock (imported as "
            "elastic) — none can fail, so give the model at least one ordinary soil "
            "before solving or no slip surface is possible")

    polygons = []
    for region_id, mid in sorted(assign.items()):
        ring = geom["regions"].get(region_id)
        if not ring:
            continue
        polygons.append({"polygon": Polygon([pts[i] for i in ring]),
                         "mat_id": mat_index[mid]})
    if not polygons:
        raise ValueError("This analysis has no material regions to import.")

    # S=f(depth) measures depth from the top of the material's own layer, column by
    # column; cp measures it from one elevation. They agree when that top is level.
    for i in depth_from_top:
        caveats.extend(cp_from_layer_top(materials[i], polygons, i,
                                         "SLOPE/W", "S=f(depth)"))
    # A Spatial Mohr-Coulomb cohesion line imported as cp: SLOPE/W holds the strength
    # at the line's lowest point below it, and cp keeps rising.
    for i, (y_low, c_low) in spatial_cp:
        bottom = min((p["polygon"].bounds[1] for p in polygons if p["mat_id"] == i),
                     default=None)
        if bottom is not None and bottom < y_low - 1e-9:
            m = materials[i]
            c_bot = m["c"] + m["cp"] * (m["r_elev"] - bottom)
            caveats.append(
                f"material '{m['name']}' reaches y = {bottom:g}, below its cohesion "
                f"function's lowest point (y = {y_low:g}). SLOPE/W holds the cohesion "
                f"at {c_low:g} below that point; XSLOPE's cp keeps rising, to "
                f"{c_bot:g} at y = {bottom:g}, so XSLOPE's strength is higher there")

    # Pore pressure. SLOPE/W names the option on the analysis; only a piezometric
    # surface maps onto xslope's piezo line. <MaterialUsesPiezs> says WHICH materials
    # draw from it — a material not listed there is dry even when a surface exists.
    piezo_line = []
    seep = None
    option = water.get("option") or "None"
    uses = stab.get("material_piezo", {})
    surfaces = stab.get("piezo", {})
    rapid = bool(stab.get("rapid_drawdown"))
    # PiezoSurfaceBbar is a piezometric surface plus B-bar's excess pressure. The
    # surface part is plain head, the same as PiezoSurface (SLOPE/W's own B-bar
    # example: its per-column pressure is the head plus B-bar x the added load), so
    # it is imported; the B-bar part is reported below.
    if option in ("PiezoSurface", "PiezoSurfaceBbar"):
        if surfaces:
            first_id = next(iter(surfaces))
            piezo_line = list(surfaces[first_id])
            # SLOPE/W can give each material its own piezometric surface; XSLOPE has
            # one line per stage, so every connected material reads the first one.
            # Say which, and which materials SLOPE/W gave another surface.
            other = [gsz["materials"][m]["name"] for m in used
                     if uses.get(m) not in (None, first_id) and uses.get(m) in surfaces]
            if other and not rapid:
                caveats.append(
                    f"GeoStudio gives material(s) {', '.join(repr(o) for o in other)} "
                    f"a different piezometric surface from the rest; XSLOPE has one "
                    f"piezometric line, so every material was imported reading "
                    f"piezometric surface {first_id}. Check the water in those materials")
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
            if len(surfaces) > 1 and not rapid and not other:
                # (A staged rapid drawdown has two by design, the water before and
                # after drawdown; _rapid_drawdown imports both. When materials use
                # different surfaces, the note above already says which was taken.)
                caveats.append(
                    f"the file defines {len(surfaces)} piezometric surfaces; "
                    f"xslope takes one, so the first was imported")
        else:
            caveats.append("the analysis asks for a piezometric surface but the file "
                           "defines none — pore pressure imported as zero")
    if option in ("PiezoSurfaceBbar", "Bbar"):
        bb = stab.get("bbar") or {}
        coef = [f"{gsz['materials'][m]['name']!r} {bb[m]['bbar']}" for m in used
                if bb.get(m, {}).get("bbar")]
        added = [repr(gsz["materials"][m]["name"]) for m in used
                 if bb.get(m, {}).get("add_weight")]
        caveats.append(
            "THE FACTOR OF SAFETY WILL BE TOO HIGH: SLOPE/W adds excess pore pressure "
            "from B-bar"
            + (f" ({', '.join(coef)})" if coef else "")
            + (f", driven by the weight of {', '.join(added)}" if added else "")
            + ". XSLOPE has no B-bar, so that pressure was not imported. "
            + ("Only the piezometric line was. " if option == "PiezoSurfaceBbar"
               else "Pore pressure was imported as zero. ")
            + "Add the excess pressure by hand before you trust the result")
        # A ratio material's strength follows sigma'_v, which the missing pressure
        # also lowers. Leaving it out raises the strength as well as cutting u.
        ratio_bb = [repr(m["name"]) for m, mid in zip(materials, used)
                    if m["option"] == "ratio" and bb.get(mid, {}).get("bbar")]
        if ratio_bb:
            caveats.append(
                f"material(s) {', '.join(ratio_bb)} take their strength from the "
                f"vertical effective stress. In SLOPE/W, B-bar's excess pore pressure "
                f"keeps that stress, and so the strength, as it was before the added "
                f"weight. XSLOPE counts the added weight in full, so its strength there "
                f"is higher than SLOPE/W's")
    elif option == "Ru":
        # SLOPE/W's Ru: u = Ru x the weight of the soil column over the base, per
        # material (<MaterialUsesRuParams>). On SLOPE/W's own Ru files every column's
        # pressure is exactly Ru x its weight / width. XSLOPE's u = 'ru' is the same
        # definition, so the value carries over unchanged.
        ru, set_ru, bad = stab.get("ru") or {}, [], []
        for mid, mat in zip(used, materials):
            text = ru.get(mid) or ""
            if mid in bedrock or not text:
                continue
            try:
                value = float(text)
            except ValueError:
                bad.append(f"{mat['name']!r} ({text})")
                continue
            if value:
                mat["u"], mat["ru"] = "ru", value
                set_ru.append(f"{mat['name']!r} ru = {value:g}")
        if set_ru:
            caveats.append(
                f"pore pressure is SLOPE/W's Ru, imported as XSLOPE's ru with the same "
                f"value for each material: {', '.join(set_ru)}. Both programs take "
                f"u = ru times the weight of the soil above the slice base")
        elif not bad:
            caveats.append("pore pressure is SLOPE/W's Ru, but no material has a Ru "
                           "value, so the model is dry")
        if bad:
            caveats.append(
                f"THE FACTOR OF SAFETY WILL BE WRONG: SLOPE/W's Ru for "
                f"{', '.join(bad)} is not a number XSLOPE can read, so that material "
                f"was imported dry. Set its ru by hand")
        if set_ru and stab.get("surcharges"):
            caveats.append(
                "SLOPE/W counts the surcharge in the weight that Ru multiplies. "
                "XSLOPE's ru uses the soil weight alone, so its pore pressure under "
                "the surcharge is lower than SLOPE/W's")
    elif option == "Parent":
        # The water condition is a whole head field computed by a parent SEEP/W run.
        # SLOPE/W writes the TRANSFERRED field into this analysis's own result folder,
        # so the parent never has to be traversed or re-solved — the pore pressure
        # xslope needs is already in the file, on the mesh SEEP/W solved it on.
        #
        # It is not just pore pressure, either: SLOPE/W also raises the water standing
        # against the slope from that same field, and in a drawdown problem that load is
        # the bigger half of the answer. seep_ponded_dload recovers it. Dropping both
        # cost 13% of the factor of safety on GeoStudio's own rapid-drawdown example.
        seep = _import_seep(gsz, analysis, step, caveats, polygons)
        if seep:
            for mat in materials:
                mat["u"] = "seep"
        else:
            caveats.append(
                "THE FACTOR OF SAFETY WILL BE WRONG: this analysis takes its water from "
                "a parent SEEP/W analysis, but the file was saved UNSOLVED, so it holds "
                "no pore-pressure field to import — and xslope cannot re-solve the "
                "seepage problem from a .gsz. The soil and geometry are right; the water "
                "is simply absent. Solve the file in GeoStudio and re-import, or set a "
                "piezo line, before you trust any number this model produces")
    elif option not in ("None", "", None, "PiezoSurface"):
        caveats.append(
            f"THE FACTOR OF SAFETY WILL BE WRONG: pore pressure comes from SLOPE/W's "
            f"'{option}' option, which xslope cannot read — imported as zero. Set a "
            f"piezo line, ru, or a seepage solution before solving")

    # A staged rapid drawdown (Duncan, Wright and Wong 1990). SLOPE/W reads two
    # piezometric surfaces -- the higher one before drawdown, the lower one after -- and
    # drawdown strengths (c_R, phi_R) for the soils that do not drain. XSLOPE's rapid
    # drawdown takes the same three things: piezo lines 1 and 2, and d and psi.
    piezo_line2 = []
    if rapid:
        piezo_line, piezo_line2, notes = _rapid_drawdown(
            gsz, used, materials, bedrock, stab, option, piezo_line)
        caveats.extend(notes)
    caveats.extend(_unit_weight_above_wt(gsz, used, materials))

    # Partial factors. SLOPE/W divides the strengths and multiplies the loads by the set
    # the analysis names before it solves; XSLOPE has no partial factors, so the set is
    # named rather than dropped in silence.
    dfid = (stab.get("control") or {}).get("DesignFactorID")
    if dfid:
        caveats.append(_design_factor_caveat(gsz, dfid))

    # SLOPE/W's unsaturated shear-strength setting rides on the piezometric surface as
    # a suction ceiling (CapSuction/MaxSuction), NOT as a per-material phi_b. It is a
    # threshold on negative pore pressure, semantically different from xslope's Fredlund
    # phi_b/s_cap apparent-cohesion pair, so it is reported rather than converted: mapping
    # a suction cap onto s_cap without a phi_b would credit no strength at all, and
    # guessing a phi_b would invent one the file never set.
    suction = stab.get("suction")
    if suction and (suction.get("cap") or suction.get("max")):
        mx = suction.get("max")
        caveats.append(
            "GeoStudio limits the negative pore-pressure (suction) "
            + (f"to a maximum of {mx:g}" if mx else "with CapSuction on")
            + " for shear strength. That is SLOPE/W's threshold model. It is not xslope's "
            "phi_b/s_cap apparent-cohesion pair and was NOT converted; set phi_b (and "
            "s_cap) by hand if you want unsaturated strength credited")

    ground_surface, domain_polygon = build_ground_surface_from_polygons(polygons)

    # A line that stops short of the model's edges is imported as drawn, with a note
    # (see _piezo_short_ends): SLOPE/W has no water past its ends, XSLOPE refuses a
    # slice there.
    solved = bool(read_gsz_results(gsz, analysis["name"], step))
    short_lines = []
    for line, which in (
            (piezo_line, "piezometric surface" if not piezo_line2 else
             "piezometric surface before drawdown (line 1)"),
            (piezo_line2, "piezometric surface after drawdown (line 2)")):
        note = _piezo_short_ends(line or [], ground_surface, domain_polygon, which,
                                 solved)
        if note:
            caveats.append(note)
            short_lines.append(line)

    phreatic = bool(stab.get("phreatic"))
    if piezo_line and phreatic:
        # SLOPE/W's correction is head x cos^2 of the line's slope, the same as an
        # XSLOPE 'phreatic' line, but SLOPE/W skips it where the line runs along the
        # ground: on Greco's weak-layer slope its pressure there is the plain head,
        # while every column under the buried part of the line is corrected.
        spans = _line_on_ground(piezo_line, ground_surface)
        caveats.append(
            "SLOPE/W applies its phreatic correction to this analysis's piezometric "
            "line, so the line was imported as a 'phreatic' line: the head under it "
            "is reduced by cos² of the line's slope"
            + ("" if not spans else
               ". SLOPE/W leaves the head uncorrected where the line runs along the "
               "ground surface ("
               + ", ".join(f"x = {a:g} to {b:g}" for a, b in spans)
               + "), and XSLOPE corrects it there too, so XSLOPE's pore pressure is "
                 "lower than SLOPE/W's there"))

    slope_data = {
        "template_version": None,
        # Persist the unit system SLOPE/W declared (or that _detect_units inferred from
        # the unit weight of water). _detect_units returns GeoStudio's own "metric"; the
        # normalizer folds that to "si". It never returns "unknown" (it raises on an
        # ambiguous file), so this is always "si"/"imperial" for a file that imported.
        # time_unit stays None: SLOPE/W's Units Time= is a display preference only, and
        # xslope imports SOLVED pore pressures, not raw conductivity, so no time-bearing
        # quantity we ingest carries a declared time base (units plan phase 2, sec 2d).
        "unit_system": normalize_unit_system(units),
        "time_unit": None,
        "gamma_water": gamma_water,
        "tcrack_depth": 0.0,
        "tcrack_water": 0.0,
        "k_seismic": gsz["stability"].get(analysis_id, {}).get("k_seismic", 0.0),
        # SLOPE/W's vertical coefficient is signed, positive acting downward (it adds
        # to each slice's weight): xslope's own convention.
        "k_vertical": gsz["stability"].get(analysis_id, {}).get("k_seismic_v", 0.0),
        "max_depth": None,
        "profile_lines": [],
        "polygons": polygons,
        "domain_polygon": domain_polygon,
        "ground_surface": ground_surface,
        "tcrack_surface": None,
        "materials": materials,
        "piezo_line": piezo_line,
        "piezo_line2": piezo_line2,
        # SLOPE/W's phreatic correction is a per-analysis switch, <Control>
        # <IsPiezPhreatic>, written only when it is on. Without it SLOPE/W's pore
        # pressure is the plain head under the line: on every such analysis in
        # SLOPE/W's verification files with per-column output, its pressure is
        # gamma_w x head exactly. Importing the line as 'phreatic' regardless cut
        # the pressure by cos^2 of the line's slope (Zhu & Lee wet: 41.7 for
        # SLOPE/W's 52.2) and raised the factor of safety 3-4%.
        "piezo_phreatic": phreatic,
        "piezo_phreatic2": phreatic,
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

    # The SEEP/W head field, and the reservoir it implies. Both need the ground surface,
    # so they land here rather than up in the pore-pressure section.
    if seep:
        mesh, seep_u, seep_step = seep
        slope_data["mesh"] = mesh
        slope_data["seep_u"] = seep_u
        step = seep_step

        mx0, my0, mx1, my1 = (mesh["nodes"][:, 0].min(), mesh["nodes"][:, 1].min(),
                              mesh["nodes"][:, 0].max(), mesh["nodes"][:, 1].max())
        dx0, dy0, dx1, dy1 = domain_polygon.bounds
        if mx0 > dx0 + 1e-6 or mx1 < dx1 - 1e-6 or my0 > dy0 + 1e-6 or my1 < dy1 - 1e-6:
            caveats.append(
                f"the SEEP/W mesh ({mx0:g}..{mx1:g}, {my0:g}..{my1:g}) does not cover "
                f"the whole slope model ({dx0:g}..{dx1:g}, {dy0:g}..{dy1:g}) — a slice "
                f"base outside the mesh gets NO pore pressure, which is not the same as "
                f"zero pore pressure being correct there. Check the uncovered ground")

    # Failure surface. SLOPE/W's search definition has no xslope equivalent, but a
    # solved file records every trial surface it evaluated — take the critical one,
    # so the model arrives complete and directly comparable.
    trials = read_gsz_results(gsz, analysis["name"], step) if critical_surface else []
    # Only circular trials. SLOPE/W writes a center and radius even for block and
    # fully-specified surfaces, but they describe a circle FITTED to the surface, not
    # the surface it solved -- importing one would silently swap the mechanism.
    circular = [t for t in trials if t["circular"]]
    if trials and not circular:
        worst = min(trials, key=lambda t: t["fs"])
        caveats.append(
            f"SLOPE/W's critical surface here is non-circular "
            f"({_SURFACE_WORDS.get(worst['kind'], _SURFACE_WORDS['unknown'])}"
            f"{', taken from its parent analysis' if worst.get('inherited') else ''}), "
            f"which is not a trial circle xslope can take — no failure surface was "
            f"imported; define one, or run a search")
    # An optimized search's final surface is not a circle (see read_gsz_results), but
    # the circles it searched first are, each with SLOPE/W's own FS. Rebuilding the
    # optimized row's fitted circle scored Greco's layered slope 24.6% high; the most
    # critical searched circle is a surface both programs can solve, so that one is
    # imported, and the caveat says SLOPE/W's final answer is lower.
    optimized = [t for t in trials if t["kind"] == "optimized"]

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
            f"the failure surface is a trial circle of SLOPE/W's own (its "
            f"{'rapid drawdown ' if rapid else ''}FS = "
            f"{chosen['fs']:.3f} by {analysis['method'] or 'its method'}), not a search "
            f"— it is imported so the model is complete and the two programs can be "
            f"compared on the same circle; run a search to find xslope's own critical "
            f"surface")
        if rejected is not None:
            # The bedrock explanation is offered only where there is bedrock: on the
            # verification files without any, the circle that will not build runs
            # out of the model, or past the end of the piezometric line, instead.
            bedrock_note = any(m["option"] == "elastic" for m in materials)
            caveats.append(
                f"note that the MOST critical of SLOPE/W's circles"
                f"{' before the optimization' if optimized else ''} (its FS = "
                f"{rejected['fs']:.3f}) is one xslope can't build on this model"
                + (" — most likely a COMPOSITE surface, truncated along the "
                   "impenetrable bedrock" if bedrock_note else
                   " — it reaches past the end of the piezometric line"
                   if _circle_past_line(rejected, short_lines, domain_polygon)
                   else "")
                + ". The circle imported above is the most critical one that does "
                  "build, so its FS is higher than SLOPE/W's")
    elif circular:
        caveats.append(
            f"none of SLOPE/W's {len(circular)} trial circles can be built on this "
            f"model (its critical FS = {rejected['fs']:.3f}), so NO failure surface was "
            f"imported — define one, or run a search. "
            + ("They reach past the end of the piezometric line"
               if all(_circle_past_line(t, short_lines, domain_polygon)
                      for t in circular) else
               "Usually this means SLOPE/W scored COMPOSITE surfaces, truncated along "
               "an impenetrable boundary"))
    else:
        caveats.append(
            "no failure surface was imported — SLOPE/W's search definition has no "
            "xslope equivalent, and this file carries no solved circles to take one "
            "from; define circles or a non-circular surface before solving (an input "
            "file with no surface will not re-load)")
    if optimized and circular:
        opt = min(optimized, key=lambda t: t["fs"])
        caveats.append(
            f"SLOPE/W optimized its critical circle into a non-circular surface (its "
            f"FS = {opt['fs']:.3f}), which xslope can't take as a circle. "
            + ("The circle imported is the most critical one SLOPE/W scored before "
               "the optimization, so its FS is higher than SLOPE/W's final answer"
               if chosen is not None else
               "None of the circles SLOPE/W scored before the optimization could be "
               "imported either"))

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

    # --- ponded water ------------------------------------------------------------
    # GeoStudio stores no ponded-water object at all: where the water rises above the
    # ground, SLOPE/W simply carries its weight. From template v22 xslope models the
    # reservoir the same way — Water loads = auto, and the engine derives the standing
    # water from the model's own water definition at solve time. So where the file
    # states a water SURFACE, the import carries that definition and writes NO water
    # distributed load: the two programs now share one implicit model, and the import
    # is simpler and more faithful at the same time. What used to be synthesised here
    # was only ever a reconstruction of what SLOPE/W's solver already implied.
    #
    # The exception is a SEEP/W-coupled analysis. It records no water surface anywhere
    # — only a solved head field — and the engine's derivation reads a piezometric line
    # or seepage head boundaries, neither of which an imported field is. That case
    # keeps the synthesis (seep_ponded_dload recovers the surface from the field) and
    # is imported as MANUAL, because a load nothing downstream can re-derive has to be
    # carried explicitly or it is lost.
    water_mode = "auto"
    derived = derive_water_loads(slope_data)
    if derived["blocks"]:
        deepest = derived["peak"] / (gamma_water or 1.0)
        caveats.append(
            f"water stands above the ground surface — that is ponded water, and "
            f"GeoStudio carries its weight implicitly. So does xslope: water loads "
            f"derive automatically from the imported water definition (the main "
            f"sheet's {_water_loads_cell()} is set to auto, and the engine measures "
            f"the reservoir off {derived['source']} at solve time), which comes to "
            f"{len(derived['blocks'])} load(s) up to {deepest:.2f} deep "
            f"({derived['peak']:.1f} pressure at the deepest point)"
            + (f". In a rapid drawdown's later stages it measures the lowered pool off "
               f"{derive_water_loads(slope_data, stage=2)['source']} the same way"
               if piezo_line2 else "")
            + ". Nothing was written to the dloads sheet, and nothing should be added "
              "there by hand")
        for note in derived["cut"]:
            caveats.append(
                f"the derived water load stops where a boundary ends: {note}. Check "
                f"it against the GeoStudio model, and if the water goes on past that "
                f"point, set the main sheet's {_water_loads_cell()} to manual with "
                f"the load typed in")
    elif seep:
        ponded = seep_ponded_dload(ground_surface, seep[0], seep[1], gamma_water)
        if ponded:
            deepest = max(p["Normal"] for b in ponded for p in b) / (gamma_water or 1.0)
            dloads.extend(ponded)
            water_mode = "manual"
            caveats.append(
                f"water stands against the slope in the SEEP/W field — SLOPE/W stores "
                f"no water surface here at all, and derives the reservoir from the "
                f"head field itself. xslope's automatic water loads read a piezometric "
                f"line or seepage head boundaries, and an imported head field is "
                f"neither, so this reservoir was recovered from the field and written "
                f"as {len(ponded)} distributed load(s), up to {deepest:.2f} deep "
                f"({gamma_water * deepest:.1f} pressure at the deepest point), with "
                f"the main sheet's {_water_loads_cell()} set to manual. Do not "
                f"re-create it by hand")
    slope_data["dloads"] = dloads
    slope_data["water_loads"] = water_mode

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
                       f"does not recognize — the crack was NOT imported, so the factor "
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
    bonded = set()          # anchors whose pullout develops over their bond length only
    for rl in stab.get("reinf_lines", []):
        product = gsz["reinforcements"].get(rl["reinforcement"])
        if product is None:
            caveats.append("a reinforcement line references a product the file does not "
                           "define — skipped")
            continue
        rlines.append(_reinforcement(product, rl, ground_surface))
        if rlines[-1]["type"] == "anchor" and _bond_limited(product, rl):
            bonded.add(rlines[-1]["label"])
    if rlines:
        from .fileio import attach_reinforce_pullout
        slope_data["reinforcement_lines"] = rlines
        slope_data["reinforce_lines"] = build_reinforce_lines(rlines)
        # A line imported onto the overburden law needs its profile before the
        # envelope means anything — the geometry, materials and water are all in
        # slope_data by now.
        attach_reinforce_pullout(slope_data)
        kinds = sorted({r["type"] or "generic" for r in rlines})
        caveats.append(
            f"{len(rlines)} reinforcement line(s) imported ({', '.join(kinds)}), acting "
            f"along the bar as a known load — the convention that reproduces SLOPE/W's "
            f"own factor of safety. GeoStudio's pullout resistance became xslope's bond "
            f"length (Lp = Tmax / pullout rate), and its plate capacity the end capacity "
            f"at the face")

        # The interface pullout law. xslope has the same law, so a product using it
        # arrives as Adhesion/Delta rather than as an approximation -- but the length
        # measure GeoStudio multiplies it by is type-dependent, and a product whose
        # interface pair could not be expressed still has to say so.
        used = [gsz["reinforcements"][rl["reinforcement"]] for rl in stab["reinf_lines"]
                if rl["reinforcement"] in gsz["reinforcements"]]
        older = sorted({p["name"] for p in used if p.get("plate_from_tensile")})
        if older:
            caveats.append(
                f"reinforcement {', '.join(repr(n) for n in older)} states its capacity "
                f"as <TensileCapacity>, the tag an earlier GeoStudio 2025 build wrote, "
                f"with no plate capacity. GeoStudio 2025.2 writes the same model with a "
                f"tensile and a plate capacity both equal to it, so the plate capacity "
                f"was taken as the tensile capacity")
        mapped = sorted({r["label"] for r in rlines
                         if r.get("delta") == r.get("delta")})
        if mapped:
            caveats.append(
                f"reinforcement {', '.join(repr(n) for n in mapped)} takes its pullout "
                f"from the soil interface (adhesion + friction on the effective "
                f"overburden); xslope reads the same law, so it was imported as the "
                f"Adhesion and Delta columns and the bond length is not used there")
        unmapped = sorted({p["name"] for p in used
                           if p.get("calculated_pullout")
                           and (p.get("adhesion") or p.get("shear_angle"))
                           and not any(r["label"] == p["name"]
                                       and r.get("delta") == r.get("delta")
                                       for r in rlines)})
        if unmapped:
            caveats.append(
                f"reinforcement {', '.join(repr(n) for n in unmapped)} states its "
                f"interface pullout with no friction angle (pure adhesion), which "
                f"xslope's law cannot express -- it needs an interface friction angle "
                f"strictly between 0 and 90. The constant bond rate was imported "
                f"instead; check the developed forces against GeoStudio")
        spread = sorted({p["name"] for p in used
                         if (p.get("distribution") or "Concentrated") != "Concentrated"})
        if spread:
            caveats.append(
                f"reinforcement {', '.join(repr(n) for n in spread)} spreads its force "
                f"over the bar in GeoStudio (its Distributed option, GeoStudio's default) "
                f"rather than applying it at the slip surface; xslope applies it where "
                f"the bar crosses the surface")

        untyped = sorted({p["name"] for p in used if not p.get("type")})
        if untyped:
            caveats.append(
                f"reinforcement {', '.join(repr(n) for n in untyped)} has no type in the "
                f"file, which is GeoStudio's default type, an anchor. It was imported as "
                f"an anchor: axial, along the bar")
        if bonded:
            caveats.append(
                f"anchor {', '.join(repr(n) for n in sorted(bonded))} develops its pullout only "
                f"over its bond length at the buried end, so its capacity was set to "
                f"the smaller of the tensile capacity and the pullout resistance times "
                f"the bond length, reached over that length from the buried end, as in "
                f"SLOPE/W")
        passive = sorted({p["name"] for p in used if p.get("fofs_dependent")})
        if passive:
            caveats.append(
                f"reinforcement {', '.join(repr(n) for n in passive)} is F of S "
                f"dependent in GeoStudio, so it was imported with Appl = passive: its "
                f"force is divided by the factor of safety, as in SLOPE/W")
        sheared = sorted({p["name"] for p in used
                          if p.get("shear") and (p.get("type") or "Anchor") != "Pile"})
        if sheared:
            caveats.append(
                f"reinforcement {', '.join(repr(n) for n in sheared)} also carries a "
                f"shear force in GeoStudio, which XSLOPE's reinforcement lines do not "
                f"have. It was not imported, so XSLOPE's factor of safety will be lower "
                f"than SLOPE/W's")
        other = sorted({f"{p['name']!r} ({p['type']})" for p in used
                        if p.get("type") and p["type"] not in _REINF_TYPE})
        if other:
            caveats.append(
                f"reinforcement {', '.join(other)} is a GeoStudio type XSLOPE does not "
                f"map. It was imported as a tensile line with the capacities the file "
                f"gives, which may be none; re-create it by hand before solving")

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
        if tag == "Control" and set(stab.get("control") or {}) - _CONTROL_HANDLED:
            # A switch this import neither applies nor reports on its own.
            caveats.append(
                f"this analysis has a GeoStudio '{tag}' that xslope does not recognize "
                f"and did NOT import — if it affects the model, the imported factor of "
                f"safety will be wrong")
        elif tag not in _ENTRY_ELEMENTS:
            caveats.append(
                f"this analysis has a GeoStudio '{tag}' that xslope does not recognize "
                f"and did NOT import — if it affects the model, the imported factor of "
                f"safety will be wrong")
        elif _ENTRY_ELEMENTS[tag]:
            caveats.append(
                f"{_ENTRY_ELEMENTS[tag]} defined in GeoStudio were NOT imported — "
                f"xslope has no mapping for them; re-enter them or the factor of safety "
                f"will be wrong")


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
    # Unit sanity cross-checks (gamma_water vs the declared system, material gammas,
    # E magnitude) as caveats -- warnings, never errors (units plan phase 2, sec 2e).
    caveats.extend(units_check(slope_data))
    return slope_data, caveats


def gsz_style(gsz, analysis_id=None):
    """The material colors from a .gsz, as an xslope style sidecar.

    GeoStudio gives every material a colour, and it is how people recognize their own
    model at a glance -- so a model that comes back in xslope's default palette looks
    like someone else's. Returns ``{'materials': {'<mat_id>': {'color': '#rrggbb'}}}``,
    keyed the way :mod:`xslope.style` expects, or ``{}`` if the file names no colours.

    Kept out of ``gsz_to_slope_data`` on purpose: colour is presentation, and
    ``slope_data`` is the model that gets written to the input file.
    """
    if analysis_id is None:
        analysis_id = gsz["analyses"][0]["id"] if gsz["analyses"] else None
    assign = gsz["contexts"].get(analysis_id, {})

    # Same ordering rule gsz_to_slope_data uses, so mat_id lines up with the polygons.
    # Bedrock is now KEPT (imported as an elastic material), so it is included here too;
    # excluding it would shift every later material's colour by one.
    used = sorted(set(assign.values()))
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


def read_gsz_results(gsz, analysis_name, step=None):
    """SLOPE/W's own solved trial slip surfaces for an analysis, if the file has them.

    Returns ``[{'fs','xo','yo','r','weight','circular','kind','inherited'}, ...]`` — one
    entry per trial surface SLOPE/W evaluated, with the factor of safety it computed.
    Empty if the file was saved unsolved.

    ``step`` names a saved time step ('000', '001', …) for a transient analysis, which
    has a whole result set — and a different factor of safety on every surface — at each
    one. The default is the first. See :func:`gsz_seep_steps`.

    ``weight`` is SLOPE/W's own weight of the sliding mass, which is worth as much as the
    factor of safety: it checks the imported geometry and unit weights on their own,
    without the solver in the way.

    ``kind`` says what the surface is: 'circle', 'block', 'specified', 'automatic',
    'optimized', 'axis' (a circle of a search with a moment axis, whose center the
    file does not hold) or 'unknown', and ``circular`` is ``kind == 'circle'``. Every row
    carries a center and a radius, but for anything but a circle they describe a
    *fitted* circle, not the surface that was solved -- rebuilding them as circles
    silently scores the wrong geometry. Each row's SlipOrigSearchMethod decides it
    (``_SEARCH_KIND``); a 6, and every row of an older file with no such column, is
    decided by the search option of the analysis that made the surface
    (``_OPTION_KIND``, through :func:`_search_option`). ``inherited`` marks a surface
    taken from a parent analysis.

    An optimized search (``<IsOptimizationOn>``) saves its search rows as usual, then
    ONE more row, the last, for the surface the optimization made out of the critical
    one: no longer a circle, though it is written with the search's own code and a
    fitted circle, and it has the lowest factor of safety. That row is 'optimized'.
    An analysis that inherits an optimized parent's surface solves both of the
    parent's critical surfaces, and the file does not say which row is which, so none
    of its rows is taken as a circle.

    This is not part of the model. It exists so the same surfaces can be re-solved
    in xslope and the two programs' answers compared directly, with no difference
    in search to explain away.
    """
    zf = gsz["zip"]
    folder = _result_folder(gsz, analysis_name)
    if folder is None:
        return []
    prefix = folder.lower() + "/"
    if step is not None:
        prefix += step.lower() + "/"
    name = next((n for n in sorted(zf.namelist())
                 if n.lower().startswith(prefix) and n.endswith("slip_surface.csv")), None)
    if name is None:
        return []

    out = []
    text = zf.read(name).decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(text))
    fields = reader.fieldnames or []
    rows = list(reader)
    # An older file names the radius SlipRadius and records no search code per row.
    radius_col = "SlipRadiusX" if "SlipRadiusX" in fields else "SlipRadius"
    coded = "SlipOrigSearchMethod" in fields
    search = _search_option(gsz, analysis_name)
    by_option = (_OPTION_KIND.get(search["option"], "unknown") if search["found"]
                 else "unknown")
    # A search with a moment axis (<Axis>) writes the AXIS as every row's SlipCenter,
    # not the circle's center: on SLOPE/W's "Reinforcement with Anchors" every row
    # gives (19, 18), and the critical circle's column bases lie on a circle centred
    # at (187.64, 119.34) with the row's radius. Rebuilt from the row, the circle is
    # a different surface (a sliding mass of 2.6 times SLOPE/W's weight). The one
    # surface whose columns are saved gets its center back from them; the others are
    # 'axis' rows, not rebuildable.
    axis_fit = (_axis_circle_from_columns(zf, name.rsplit("/", 1)[0].lower() + "/", rows)
                if search.get("axis") else {})
    for i, row in enumerate(rows):
        try:
            fs = float(row["SlipFOS"])
            xo, yo = float(row["SlipCenterX"]), float(row["SlipCenterY"])
            radius = float(row[radius_col])
        except (KeyError, TypeError, ValueError):
            continue                     # a trial SLOPE/W could not evaluate
        if not 0.0 < fs < 100.0:
            continue                     # SLOPE/W writes a sentinel FOS for failed trials
        try:
            weight = float(row.get("SlipWeight") or 0.0)
        except ValueError:
            weight = 0.0
        code = (row.get("SlipOrigSearchMethod") or "").strip() if coded else None
        kind = by_option if (code is None or code == "6") else \
            _SEARCH_KIND.get(code, "unknown")
        if search["optimized"]:
            if search["inherited"]:
                if kind == "circle":
                    kind = "optimized"
            elif i == len(rows) - 1:
                kind = "optimized"
        if search.get("axis") and kind == "circle":
            fit = axis_fit.get((row.get("SlipNum") or "").strip())
            if fit is not None and abs(fit[2] - radius) <= 1e-3 * radius:
                xo, yo = fit[0], fit[1]
            else:
                kind = "axis"
        out.append({"fs": fs, "xo": xo, "yo": yo, "r": radius, "weight": weight,
                    "circular": kind == "circle", "kind": kind,
                    "inherited": search["inherited"]})
    return out


def _axis_circle_from_columns(zf, prefix, rows):
    """``{SlipNum: (xc, yc, R)}`` for each saved slip whose column file
    (``column_<SlipNum>.csv``) holds column bases lying on one circle.

    A least-squares circle through the columns' base corners (XLeft, YBackBotLeft)
    and (XRight, YBackBotRight), which lie on the slip surface itself, kept only when
    every corner is within 1e-4 x R of it: a surface cut short by bedrock, or
    optimized away from a circle, is not one circle and is left out.
    """
    names = {n.lower(): n for n in zf.namelist()}
    out = {}
    for row in rows:
        num = (row.get("SlipNum") or "").strip()
        name = names.get(f"{prefix}column_{num}.csv".lower())
        if not num or name is None:
            continue
        try:
            cols = list(csv.DictReader(io.StringIO(zf.read(name).decode("utf-8-sig"))))
            x = np.array([float(c[k]) for c in cols for k in ("XLeft", "XRight")])
            y = np.array([float(c[k]) for c in cols
                          for k in ("YBackBotLeft", "YBackBotRight")])
        except (KeyError, TypeError, ValueError):
            continue
        if len(x) < 3:
            continue
        a = np.c_[2.0 * x, 2.0 * y, np.ones_like(x)]
        (xc, yc, k), *_ = np.linalg.lstsq(a, x * x + y * y, rcond=None)
        r2 = k + xc * xc + yc * yc
        if r2 <= 0:
            continue
        r = math.sqrt(r2)
        if np.abs(np.hypot(x - xc, y - yc) - r).max() <= 1e-4 * r:
            out[num] = (float(xc), float(yc), float(r))
    return out


def _search_option(gsz, analysis_name):
    """How an analysis's slip surfaces were made, for classifying its saved trial rows.

    Returns ``{'found', 'option', 'optimized', 'inherited', 'axis'}``. ``axis`` is
    True when the search names a moment axis (``<Axis>``). ``option`` is the
    ``<SlipSurface><Option>`` of the search that made the surfaces: None for the default
    entry-and-exit search, which GeoStudio writes as no ``<Option>`` at all. A
    ``FromOther`` analysis solved its parent's surface, so it answers with the parent's
    search (walking up a chain of them) and ``inherited`` is True. ``optimized`` is that
    search's ``<IsOptimizationOn>``. ``found`` is False when the analysis or its
    stability entry cannot be located; the caller then treats its rows as not circular
    rather than guess.

    The XML saved in the analysis's own result folder is read first: it is the model
    SLOPE/W solved, where the root XML may have been edited since.
    """
    none = {"found": False, "option": None, "optimized": False, "inherited": False,
            "axis": False}
    zf = gsz["zip"]
    names = zf.namelist()
    folder = _result_folder(gsz, analysis_name)
    xml = None
    if folder is not None:
        prefix = folder.lower() + "/"
        xml = next((n for n in sorted(names) if n.lower().startswith(prefix)
                    and "/" not in n[len(prefix):] and n.lower().endswith(".xml")), None)
    if xml is None:
        xml = next((n for n in names if "/" not in n and n.lower().endswith(".xml")), None)
    if xml is None:
        return none
    try:
        root = ET.fromstring(zf.read(xml))
    except ET.ParseError:
        return none

    ids, parents = {}, {}
    for a in root.findall("./Analyses/Analysis"):
        aid = _text(a, "ID")
        ids[(_text(a, "Name") or "").lower()] = aid
        parents[aid] = _text(a, "ParentID")
    entries = {_text(s, "AnalysisID"): s.find("Entry")
               for s in root.findall("./StabilityItems/StabilityItem")}

    aid = ids.get(analysis_name.lower())
    inherited = False
    for _ in range(len(parents) + 1):       # bounds the walk up a FromOther chain
        entry = entries.get(aid)
        if entry is None:
            return none
        slip = entry.find("SlipSurface")
        option = _text(slip, "Option")
        if option != "FromOther":
            optimized = (_text(slip, "IsOptimizationOn", "") or "").lower() == "true"
            axis = slip is not None and slip.find("Axis") is not None
            return {"found": True, "option": option, "optimized": optimized,
                    "inherited": inherited, "axis": axis}
        inherited = True
        aid = parents.get(aid)
    return none


# GeoStudio saves an analysis's results in a folder named after the analysis. A
# character a folder name cannot hold is written as '&' and a code instead. Read off
# GeoStudio's own files: '/' is '&3' ("SLOPE/W Analysis" is saved as "SLOPE&3W
# Analysis"), ':' is '&4' ("Case 1: fixed-end anchor"), and a final '.' is '&-'
# ("... Suction St." is saved as "... Suction St&-"; Windows drops a trailing period).
# The other characters Windows refuses (\ * ? " < > |) and '&' itself must have codes
# too, but none appears in any file we have, so _result_folder lets an unknown code
# stand for any of them rather than guess which.
_FOLDER_ESCAPES = {"/": "&3", ":": "&4"}
_FOLDER_UNESCAPES = {code[1]: ch for ch, code in _FOLDER_ESCAPES.items()}
_FOLDER_RESERVED = set('\\/:*?"<>|&')


def _folder_key(text, is_folder):
    """``text`` lower-cased, in a form an analysis name and its folder name share.

    In a folder name each known code is decoded ('&3' -> '/', '&4' -> ':', a final
    '&-' -> '.') and any other '&' code becomes a placeholder. In an analysis name a
    character GeoStudio must escape but whose code is unknown becomes the same
    placeholder. A known code is therefore matched exactly, so "SLOPE:W" never takes
    the folder of "SLOPE/W"."""
    out, i, n = [], 0, len(text)
    while i < n:
        ch = text[i]
        if is_folder and ch == "&" and i + 1 < n:
            code = text[i + 1]
            if code in _FOLDER_UNESCAPES:
                out.append(_FOLDER_UNESCAPES[code])
            elif code == "-" and i + 2 == n:
                out.append(".")
            else:
                out.append("\0")
            i += 2
            continue
        unknown = (not is_folder and ch in _FOLDER_RESERVED
                   and ch not in _FOLDER_ESCAPES)
        out.append("\0" if unknown else ch)
        i += 1
    return "".join(out).lower()


def _result_folder(gsz, analysis_name):
    """The folder in the .gsz that holds an analysis's results, or None if it has none.

    Every reader of an analysis's results goes through this -- its trial surfaces, its
    saved model, its mesh and its SEEP/W field. Taking the analysis name as the folder
    name found nothing for an analysis named "SLOPE/W Analysis", whose folder is
    "SLOPE&3W Analysis" (see ``_FOLDER_ESCAPES``), so six of SLOPE/W's verification
    files imported no surface at all. The name is tried as it is, then against every
    escaped folder name (``_folder_key``); that match must be unique, or no folder is
    returned rather than another analysis's.
    """
    folders = {}
    for n in gsz["zip"].namelist():
        if "/" in n:
            f = n.split("/", 1)[0]
            folders.setdefault(f.lower(), f)
    name = analysis_name or ""
    if name.lower() in folders:
        return folders[name.lower()]
    key = _folder_key(name, is_folder=False)
    hits = [f for f in folders.values()
            if "&" in f and _folder_key(f, is_folder=True) == key]
    return hits[0] if len(hits) == 1 else None


# --------------------------------------------------------------------------------------
# SEEP/W: the pore-pressure field a parent seepage analysis computed.
#
# A stability analysis whose water option is "Parent" takes its pore pressure from a
# SEEP/W run, and SLOPE/W writes the *transferred* field into the stability analysis's
# own result folder. So there is no need to traverse to the parent, and no need to
# re-solve anything: the head field xslope needs is already sitting in the file, on the
# mesh SEEP/W solved it on.
# --------------------------------------------------------------------------------------

# PLY scalar types -> struct codes. GeoStudio writes binary_little_endian.
_PLY_TYPES = {"char": "b", "uchar": "B", "short": "h", "ushort": "H", "int": "i",
              "uint": "I", "float": "f", "double": "d",
              "int8": "b", "uint8": "B", "int16": "h", "uint16": "H",
              "int32": "i", "uint32": "I", "float32": "f", "float64": "d"}

# Node counts xslope's element library knows. Anything else is reported, not guessed at.
_MESH_ELEMENT_SIZES = {3, 4, 6, 8, 9}


def _read_ply(raw):
    """Parse a binary PLY into ``{element name: [{property: value}, ...]}``.

    Generic over the header rather than hard-coded to the layout GeoStudio happens to
    write today: offsets come from the declared property types, so a reordered or
    extended element block still reads correctly instead of silently misaligning.
    """
    marker = raw.find(b"end_header")
    if marker < 0:
        raise ValueError("mesh file is not a PLY (no end_header)")
    header = raw[:marker].decode("ascii", "replace").splitlines()
    off = marker + len("end_header")
    while raw[off:off + 1] in (b"\r", b"\n"):
        off += 1

    if not any(l.startswith("format binary_little_endian") for l in header):
        raise ValueError("only binary little-endian PLY meshes are supported; this one "
                         "declares: " + next((l for l in header if l.startswith("format")),
                                             "no format at all"))

    blocks, cur = [], None
    for line in header:
        w = line.split()
        if not w:
            continue
        if w[0] == "element":
            cur = {"name": w[1], "count": int(w[2]), "props": []}
            blocks.append(cur)
        elif w[0] == "property" and cur is not None:
            if w[1] == "list":
                cur["props"].append(("list", w[2], w[3], w[4]))
            else:
                cur["props"].append(("scalar", w[1], w[2]))

    out = {}
    for b in blocks:
        rows = []
        for _ in range(b["count"]):
            row = {}
            for p in b["props"]:
                if p[0] == "scalar":
                    _, ty, name = p
                    fmt = _PLY_TYPES[ty]
                    row[name], = struct.unpack_from("<" + fmt, raw, off)
                    off += struct.calcsize(fmt)
                else:
                    _, count_ty, item_ty, name = p
                    cf, itf = _PLY_TYPES[count_ty], _PLY_TYPES[item_ty]
                    n, = struct.unpack_from("<" + cf, raw, off)
                    off += struct.calcsize(cf)
                    row[name] = list(struct.unpack_from(f"<{n}{itf}", raw, off))
                    off += struct.calcsize(itf) * n
            rows.append(row)
        out[b["name"]] = rows
    return out


def gsz_mesh(gsz, analysis_name, polygons=None):
    """The finite element mesh an analysis was solved on, in xslope's mesh form.

    Returns ``{'nodes': (N,2), 'elements': (M,8), 'element_types': (M,),
    'element_materials': (M,)}`` — node indices 0-based and padded to 8 columns, exactly
    as ``xslope.mesh`` wants them — or None if the analysis carries no mesh.

    ``polygons`` are the imported material zones. They are what ``element_materials``
    comes from (by which zone each element's centroid falls in), since GeoStudio's mesh
    records the region an element belongs to, not xslope's material index. Without it the
    mesh is missing a key that ``fem.py`` and ``seep.py`` both require.

    A GeoStudio mesh block also holds the *line* and *point* entities that mark region
    edges and vertices. Those are not domain elements and are dropped; taking them for
    elements would leave degenerate zero-area cells that no point ever falls inside.
    """
    zf = gsz["zip"]
    folder = _result_folder(gsz, analysis_name)
    if folder is None:
        return None
    want = (folder + "/Mesh.ply").lower()
    name = next((n for n in zf.namelist() if n.lower() == want), None)
    if name is None:
        return None

    ply = _read_ply(zf.read(name))
    if "node" not in ply or "element" not in ply:
        raise ValueError(f"the mesh in '{analysis_name}' has no node/element blocks")

    nodes = np.array([[r["x"], r["y"]] for r in ply["node"]], dtype=float)

    rows, sizes, odd = [], [], set()
    for r in ply["element"]:
        ids = r.get("id") or []
        if len(ids) < 3:
            continue                       # a line or point entity, not a domain element
        if len(ids) not in _MESH_ELEMENT_SIZES:
            odd.add(len(ids))
            continue
        rows.append([i - 1 for i in ids] + [0] * (8 - len(ids)))   # PLY ids are 1-based
        sizes.append(len(ids))
    if odd:
        raise ValueError(
            f"the mesh in '{analysis_name}' contains element(s) with "
            f"{sorted(odd)} nodes, which xslope's element library does not have. "
            f"Importing them as anything else would fabricate pore pressures.")
    if not rows:
        return None

    elements = np.array(rows, dtype=int)
    element_types = np.array(sizes, dtype=int)
    ensure_ccw_elements(nodes, elements, element_types)

    # element_materials is 1-based, and every element must have one: fem.py and seep.py
    # index straight into it. An element whose centroid lands in no zone (the mesh edge
    # rounding outside a polygon) falls back to the first material rather than 0, which
    # would index the wrong soil.
    mats = np.ones(len(elements), dtype=int)
    if polygons:
        from shapely.geometry import Point
        for i, (elem, k) in enumerate(zip(elements, element_types)):
            cx = float(np.mean(nodes[elem[:k], 0]))
            cy = float(np.mean(nodes[elem[:k], 1]))
            c = Point(cx, cy)
            for p in polygons:
                if p["polygon"].contains(c):
                    mats[i] = p["mat_id"] + 1
                    break

    return {"nodes": nodes, "elements": elements, "element_types": element_types,
            "element_materials": mats}


def gsz_seep_steps(gsz, analysis_name):
    """The saved time steps of an analysis that carry a pore-pressure field.

    Returns ``[{'step': '000', 'time': 0.0}, ...]``, oldest first. A steady-state run
    has one; a transient one has a result set per saved step, and a stability analysis
    on top of it has a *different factor of safety at every one*. Empty if the file was
    saved unsolved — in which case there is no head field to import at all.
    """
    zf = gsz["zip"]
    folder = _result_folder(gsz, analysis_name)
    if folder is None:
        return []
    prefix = (folder + "/").lower()
    steps = sorted({n[len(prefix):].split("/")[0]
                    for n in zf.namelist()
                    if n.lower().startswith(prefix) and n.lower().endswith("/node.csv")})
    times = next((a["times"] for a in gsz["analyses"] if a["name"] == analysis_name), [])
    out = []
    for i, s in enumerate(steps):
        out.append({"step": s, "time": times[i] if i < len(times) else None})
    return out


def read_gsz_seep(gsz, analysis_name, step=None, polygons=None):
    """The pore-pressure field a parent SEEP/W analysis computed, on its own mesh.

    Returns ``(mesh, u, step)`` where ``mesh`` is :func:`gsz_mesh`'s dict and ``u`` is
    the pore pressure at each node — together, exactly xslope's ``seep`` inputs. Returns
    None if the analysis has no solved field.

    ``step`` selects a saved time step by its folder name ('000', '001', ...); the
    default is the *last*. Note that SLOPE/W stores the transferred field alongside the
    stability results, so this reads the stability analysis's own folder — the parent
    SEEP/W analysis never has to be re-solved, or even looked at.
    """
    steps = gsz_seep_steps(gsz, analysis_name)
    if not steps:
        return None
    chosen = step if step is not None else steps[-1]["step"]
    if not any(s["step"] == chosen for s in steps):
        raise ValueError(f"'{analysis_name}' has no time step '{chosen}' "
                         f"(it has {', '.join(s['step'] for s in steps)})")

    mesh = gsz_mesh(gsz, analysis_name, polygons)
    if mesh is None:
        return None

    zf = gsz["zip"]
    want = f"{_result_folder(gsz, analysis_name)}/{chosen}/node.csv".lower()
    name = next((n for n in zf.namelist() if n.lower() == want), None)
    if name is None:
        return None

    text = zf.read(name).decode("utf-8-sig")
    rows = list(csv.DictReader(io.StringIO(text)))
    n = len(mesh["nodes"])
    # Guard the count BEFORE indexing: a field on a different mesh must fail with a clear
    # message, not an IndexError from u[i] when a stray row's node number runs past n.
    if len(rows) != n:
        raise ValueError(
            f"the pore-pressure field for '{analysis_name}' step {chosen} has "
            f"{len(rows)} values but the mesh has {n} nodes — they do not describe the "
            f"same problem, and interpolating one onto the other would invent pressures")
    u = np.zeros(n, dtype=float)
    for row in rows:
        try:
            i = int(row["Node"]) - 1                       # 1-based, as in the mesh
            if 0 <= i < n:
                u[i] = float(row["PoreWaterPressure"])
        except (KeyError, TypeError, ValueError):
            continue
    return mesh, u, chosen


def seep_ponded_dload(ground_surface, mesh, u, gamma_water, n=200):
    """The water standing against the slope, derived from a seepage pore-pressure field.

    GeoStudio stores no ponded-water object *and* no water surface when the pore
    pressure comes from SEEP/W — yet SLOPE/W still puts the weight of the reservoir on
    the submerged face, and in a drawdown problem that load is the bigger half of the
    answer. It can only be coming from the head field itself, and it is: at any point on
    the ground surface, water stands to elevation ``y + u/gamma_w``.

    That one rule reproduces both cases exactly. Under a reservoir the surface is a
    specified-head boundary, so ``u = gamma_w (H - y)`` and the implied level comes out
    at the reservoir level H, to the millimetre, at every submerged point. On a drained
    or seepage face ``u = 0``, the implied level collapses onto the ground, and no load
    is produced. Measured against GeoStudio's own rapid-drawdown example, this puts
    water on the face at the full-reservoir step and none at any drawn-down step, which
    is what SLOPE/W's own per-slice surcharge forces do.

    Returns dload blocks in xslope's form, or [] if nothing is submerged.
    """
    x0, _, x1, _ = ground_surface.bounds
    d = max(1e-6, 1e-3 * (x1 - x0))       # sample just inside the mesh; the offset cancels
    tol = 1e-3 * (x1 - x0)                # below this, it is round-off, not a reservoir

    # Sample at the ground's OWN vertices as well as on an even grid. Without them the
    # two polylines have different breakpoints, and across a downward kink in the ground
    # the water line's straight chord rides above the corner -- inventing a puddle a few
    # centimetres deep at every convex bend, out of nothing but the sampling.
    xs = sorted(set(np.linspace(x0, x1, n)) |
                {x for x, _ in ground_surface.coords if x0 <= x <= x1})

    water, wet = [], False
    for x in xs:
        yg = _y_on(ground_surface, x)
        if yg is None:
            continue
        ys = yg - d
        val, found = interpolate_at_point(mesh["nodes"], mesh["elements"],
                                          mesh["element_types"], u, (x, ys),
                                          return_found=True)
        # The water surface implied here. Taking it from the SAMPLE elevation rather
        # than the ground is what makes the inward offset cancel out of the answer.
        level = ys + val / gamma_water if found else yg
        if level - yg <= tol:
            level = yg                     # dry face, or too shallow to be real water
        else:
            wet = True
        water.append((x, level))

    if not wet:
        return []
    return material_above_ground_dload(ground_surface, water, gamma_water)


def _xml_escape(text):
    return (str(text).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


def _export_strength(m, mat_id, polygons):
    """``(SlopeModel, [(tag, value), ...], notes)`` for one XSLOPE material.

    The reverse of :func:`_strength`. ``ratio`` is SLOPE/W's SHANSEP with a constant
    ratio, which has no minimum strength. ``cp`` is S=f(depth), which measures depth
    from the top of the material's own layer: the same strength when that top is level
    and at or below r_elev. Otherwise the export writes the strength XSLOPE has at the
    layer's highest point and the same rate, and says where the two differ.
    """
    name = m.get("name") or f"material {mat_id + 1}"
    opt = (m.get("option") or "").strip().lower()
    c, phi = float(m.get("c") or 0.0), float(m.get("phi") or 0.0)
    rate, r_elev = float(m.get("cp") or 0.0), float(m.get("r_elev") or 0.0)
    if opt == "elastic":
        return "Bedrock", [("CohesionPrime", c), ("PhiPrime", phi)], []
    if opt in ("mc", ""):
        return "MohrCoulomb", [("CohesionPrime", c), ("PhiPrime", phi)], []
    if opt == "ratio":
        notes = []
        if c > 0:
            notes.append(
                f"material '{name}' was written as SLOPE/W's SHANSEP model with ratio "
                f"{rate:g}. Its minimum strength c = {c:g} was not written, because "
                f"SHANSEP has none. Where c/p x sigma'_v is below {c:g}, SLOPE/W's "
                f"strength is lower than XSLOPE's.")
        return "SHANSEP", [("CohesionPrime", 0.0), ("PhiPrime", 0.0),
                           ("TauSigmaRatio", rate)], notes
    if opt == "cp":
        from .fileio import material_top_range
        lo, hi = material_top_range(polygons, mat_id)
        if hi is None:
            lo = hi = r_elev
        c_top = c + rate * max(0.0, r_elev - hi)
        notes = []
        tol = 1e-6 * max(1.0, abs(hi), abs(r_elev))
        if rate and hi > r_elev + tol:
            # XSLOPE holds c above r_elev; S=f(depth) has no constant part. Writing
            # the strength at r_elev minus the rate over the band keeps the strength
            # exact below r_elev and lower above it, where it can be.
            c_top = c - rate * (hi - r_elev)
            if c_top >= 0 or rate < 0:
                notes.append(
                    f"material '{name}' holds c = {c:g} above r_elev = {r_elev:g}, and "
                    f"its layer reaches y = {hi:g}. SLOPE/W's S=f(depth) has no constant "
                    f"part, so it was written as {c_top:g} at the top of the layer, "
                    f"changing {rate:g} per unit depth: the same strength below "
                    f"y = {r_elev:g}, and {'less' if rate > 0 else 'more'} than "
                    f"XSLOPE's above it.")
            else:
                extra = -c_top
                c_top = 0.0
                notes.append(
                    f"material '{name}' holds c = {c:g} above r_elev = {r_elev:g}, and "
                    f"its layer reaches y = {hi:g}. SLOPE/W's S=f(depth) has no constant "
                    f"part, so it was written as 0 at the top of the layer, rising "
                    f"{rate:g} per unit depth. Below y = {r_elev:g} that is {extra:g} "
                    f"more than XSLOPE's strength.")
        if rate and hi - lo > tol:
            notes.append(
                f"material '{name}': SLOPE/W's S=f(depth) measures depth from the top "
                f"of the layer above each slice, and this layer's top is not level "
                f"(y = {lo:g} to {hi:g}). Where it is lower, SLOPE/W's strength is "
                f"{'lower' if rate > 0 else 'higher'} than XSLOPE's, by up to "
                f"{abs(rate) * (hi - lo):g}.")
        return "SFnDepth", [("CohesionPrime", 0.0), ("PhiPrime", 0.0),
                            ("CTopOfLayer", c_top), ("CRateOfIncrease", rate)], notes
    if opt == "hb":
        why = ("SLOPE/W has a Hoek-Brown model, but the export does not write it yet. "
               "It")
    else:
        why = f"SLOPE/W has no equivalent of strength option '{opt}'. It"
    return "MohrCoulomb", [("CohesionPrime", c), ("PhiPrime", phi)], [
        f"material '{name}': {why} was written as Mohr-Coulomb with c = {c:g} and "
        f"phi = {phi:g}, which is not the same strength."]


def _has_split(m):
    """True when an XSLOPE material weighs differently below the water table."""
    g_sat = m.get("gamma_sat")
    try:
        return g_sat is not None and abs(float(g_sat) - float(m.get("gamma") or 0.0)) > 1e-12
    except (TypeError, ValueError):
        return False


def _unit_weight_out(m, fields):
    """A material's unit weight and strength fields, as SLOPE/W writes them.

    Its <UnitWeight> is the saturated weight. A different weight above the water
    table is <DryWeight>, after the strength fields, and switched on by
    <SlopeUseUnitWeightAboveWT>ConstantValue ahead of <UnitWeight>, in the order
    SLOPE/W's own Unit Weight Definition example writes them."""
    gamma = float(m.get("gamma") or 0.0)
    split = _has_split(m)
    return ((f'<SlopeUseUnitWeightAboveWT>ConstantValue</SlopeUseUnitWeightAboveWT>'
             if split else '')
            + f'<UnitWeight>{(float(m["gamma_sat"]) if split else gamma):.10g}</UnitWeight>'
            + "".join(f'<{tag}>{val:.10g}</{tag}>' for tag, val in fields)
            + (f'<DryWeight>{gamma:.10g}</DryWeight>' if split else ''))


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
    list), piles, and the power-curve and Hoek-Brown strengths. Each is reported rather
    than dropped in silence. A ``ratio`` material is written as SLOPE/W's SHANSEP and a
    ``cp`` one as S=f(depth) (see :func:`_export_strength`).

    gamma_sat, where it differs from gamma, is SLOPE/W's <UnitWeight> with gamma as its
    separate unit weight above the water table (:func:`_unit_weight_out`), and a
    phreatic piezo line switches on SLOPE/W's <IsPiezPhreatic>.

    A rapid drawdown model -- piezo lines 1 and 2, and d and psi on a Mohr-Coulomb
    material -- is written as SLOPE/W's staged rapid drawdown, the way SLOPE/W's own
    examples write one: the <IsRapidDrawdown> switch under the analysis's <Control>,
    both lines as piezometric surfaces (before drawdown first), and each d and psi as
    the <DrawdownCohesion> and <DrawdownPhi> (c_R, phi_R) that SLOPE/W converts back to
    them (:func:`_r_envelope`).

    A finite-element seepage field crosses as DATA, not as a model: its nodal pressures
    are written as SLOPE/W's spatial pore-pressure function (a set of discrete points
    GeoStudio interpolates between), with no SEEP/W analysis and no fabricated solution.
    The receiving user gets the pressures xslope solved with; re-solving the seepage
    problem means re-creating it in SEEP/W.

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

    gamma_water = require_gamma_water(slope_data, "GeoStudio export")
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
        ground_surface, _ = build_ground_surface_from_polygons(polygons)

    # A finite-element seepage field. The MODEL behind it (mesh, conductivity functions,
    # boundary conditions) has no place in a SLOPE/W-only .gsz, and fabricating a SEEP/W
    # parent with invented properties and a pretended solution would be a lie in a file
    # format. But the PRESSURES do have a place: SLOPE/W's own documented spatial
    # pore-pressure input is a set of discrete (x, y, value) points that it interpolates
    # between -- GeoStudio's <Fn3D> with OutputType PressureHead, selected on the
    # analysis by <ResultInputInfo><Option>3DFunction</Option>. So the field crosses as
    # DATA, on xslope's own seepage nodes, at full resolution.
    #
    # Substituting a piezometric line would still be wrong (it MEANS hydrostatic pressure
    # beneath it, which is exactly the assumption a seepage run exists to avoid) and is
    # not done.
    seep_pts = []
    other_u = sorted({(m.get("u") or "none") for m in materials} - {"none", "piezo"})
    if "seep" in other_u:
        mesh = slope_data.get("mesh") or {}
        nodes = mesh.get("nodes") if isinstance(mesh, dict) else None
        u_field = slope_data.get("seep_u")
        n_nodes = 0 if nodes is None else len(nodes)
        if nodes is not None and u_field is not None and len(u_field) == n_nodes > 0:
            # GeoStudio's water spatial function is a PRESSURE HEAD function in every
            # vendor file that defines one, so the pressures are divided through by the
            # unit weight of water written on the same analysis.
            seep_pts = [(float(nodes[i][0]), float(nodes[i][1]),
                         float(u_field[i]) / gamma_water) for i in range(n_nodes)]
            n_suction = sum(1 for _, _, z in seep_pts if z < 0.0)
            caveats.append(
                f"this model's pore pressure is a finite-element SEEPAGE FIELD. The "
                f"PRESSURES were written, as SLOPE/W's spatial pore-pressure function: "
                f"{len(seep_pts)} discrete points, one per node of xslope's seepage mesh, "
                f"carrying pressure head (u / {gamma_water:g}) for GeoStudio to "
                f"interpolate between. What did NOT cross is the seepage MODEL — there is "
                f"no SEEP/W analysis in this file, so the field cannot be re-solved or "
                f"its boundary conditions changed on the GeoStudio side; re-create the "
                f"seepage problem in SEEP/W if you need that. GeoStudio applies a spatial "
                f"function to every material in the analysis")
            if n_suction:
                caveats.append(
                    f"{n_suction} of those {len(seep_pts)} points are NEGATIVE (suction "
                    f"above the phreatic surface) and were written as they stand. Neither "
                    f"program credits them by default — SLOPE/W only through a phi_b, "
                    f"which is not written, and xslope clamps the effective-normal pore "
                    f"pressure at zero — so the two agree, but check it if you set phi_b "
                    f"in GeoStudio")
            if slope_data.get("seep_u2") is not None:
                caveats.append(
                    "this model carries a SECOND seepage field (rapid drawdown); only the "
                    "first was written. A .gsz analysis holds one water condition, so the "
                    "drawdown stage would have to be a second analysis, which xslope does "
                    "not write")
            if piezo_line and uses_piezo:
                caveats.append(
                    "both a piezo line and a seepage field feed pore pressure in this "
                    "model. GeoStudio picks ONE source per analysis, and the seepage field "
                    "was chosen — the piezometric surface is still written, but it is "
                    "inert unless you switch the analysis over to it")
        else:
            caveats.append(
                "THE FACTOR OF SAFETY WILL DIFFER: this model's materials take their pore "
                "pressure from a seepage field, but no solved field is loaded (the mesh "
                "or the nodal pressures are missing), so there was nothing to write — "
                "GeoStudio will open the model with NO pore pressure. Load the seepage "
                "solution and export again")
        if slope_data.get("dloads"):
            caveats.append(
                "any ponded water in this model IS written, as a surcharge, and must be "
                "KEPT: a spatial pore-pressure function gives GeoStudio no water surface, "
                "so unlike the piezometric-surface case it derives no reservoir of its "
                "own and would otherwise miss the weight of the water entirely")
    other_u = [u for u in other_u if u != "seep"]
    if other_u:
        caveats.append(
            f"pore pressure from {', '.join(repr(u) for u in other_u)} is not written "
            f"— GeoStudio will open the model with no pore pressure from that source")

    # 'elastic' HAS a SLOPE/W equivalent — the impenetrable Bedrock model, written below
    # — so it is not among the models that have no equivalent. ratio and cp have one
    # too (SHANSEP, S=f(depth)); the power curve and Hoek-Brown are written with their
    # c/phi and flagged (_export_strength).
    elastic_names = [m.get("name") or "?" for m in materials
                     if (m.get("option") or "") == "elastic"]
    if elastic_names:
        caveats.append(
            f"elastic material(s) {', '.join(repr(n) for n in elastic_names)} were "
            f"written as GeoStudio's impenetrable Bedrock — a slip surface cannot enter "
            f"them, which is xslope's elastic behavior")
    # Each material's SLOPE/W model and strength fields: ratio -> SHANSEP, cp ->
    # S=f(depth) (see _export_strength); written in the <Materials> block below.
    strength_out = []
    for k, m in enumerate(materials):
        model, fields, notes = _export_strength(m, k, polygons)
        strength_out.append((model, fields))
        caveats.extend(notes)

    # A rapid drawdown. XSLOPE's takes the water before and after drawdown as piezo
    # lines 1 and 2 and a K_c = 1 envelope (d, psi) on each Mohr-Coulomb material that
    # does not drain; SLOPE/W's staged rapid drawdown takes two piezometric surfaces
    # and an R envelope (c_R, phi_R), which it converts to the same (d, psi). An XSLOPE
    # file records no analysis type, so a model carrying both lines and the envelope is
    # written as a staged rapid drawdown, and the caveat says so.
    piezo_line2 = slope_data.get("piezo_line2") or []
    drawdown = [k for k, m in enumerate(materials)
                if (m.get("option") or "mc").strip().lower() == "mc"
                and (float(m.get("d") or 0.0) or float(m.get("psi") or 0.0))]
    rapid = (len(piezo_line) >= 2 and len(piezo_line2) >= 2 and bool(drawdown)
             and not seep_pts)
    drawdown_out = []
    if rapid:
        for k in drawdown:
            m = materials[k]
            name = m.get("name") or f"material {k + 1}"
            d, psi = float(m.get("d") or 0.0), float(m.get("psi") or 0.0)
            c, phi = float(m.get("c") or 0.0), float(m.get("phi") or 0.0)
            r = _r_envelope(d, psi, phi)
            if r is None:
                caveats.append(
                    f"material '{name}': d = {d:g}, psi = {psi:g} with phi' = {phi:g} "
                    f"has no R envelope SLOPE/W's conversion gives back, so its drawdown "
                    f"strengths were not written and GeoStudio will treat it as "
                    f"free-draining")
                continue
            c_r, phi_r = r
            strength_out[k] = (strength_out[k][0], list(strength_out[k][1])
                               + [("DrawdownCohesion", c_r), ("DrawdownPhi", phi_r)])
            drawdown_out.append(f"'{name}' d = {d:g}, psi = {psi:g} -> c_R = {c_r:.4g}, "
                                f"phi_R = {phi_r:.3g}")
            # SLOPE/W's own check on an undrained material (Slope Stability Modeling,
            # sec. 9.3): c_R above c' and phi_R below phi'.
            if not (c_r > c and phi_r < phi):
                caveats.append(
                    f"material '{name}' was written with c_R = {c_r:.4g} and phi_R = "
                    f"{phi_r:.3g}. SLOPE/W checks that a soil's c_R is above its c' "
                    f"({c:g}) and its phi_R below its phi' ({phi:g}), so GeoStudio may "
                    f"refuse to solve it. XSLOPE has no such rule")
        # SLOPE/W tells the two surfaces apart by elevation, not by order.
        p1 = list(piezo_line)
        high, _low, crossing = _before_after(p1, list(piezo_line2))
        if crossing is not None or high is not p1:
            caveats.append(
                "piezo line 2 stands above line 1"
                + (f" near x = {crossing:g}" if crossing is not None else "")
                + ". SLOPE/W takes the higher of its two piezometric surfaces as the "
                  "water before drawdown, so check the two in GeoStudio")
        caveats.append(
            "this model was written as SLOPE/W's staged rapid drawdown (Duncan, Wright "
            "and Wong 1990): piezo line 1 as the piezometric surface before drawdown, "
            "line 2 as the one after, and the d and psi of each material that does not "
            "drain as the R envelope that SLOPE/W converts back to the same d and psi"
            + (f" ({'; '.join(drawdown_out)})" if drawdown_out else "")
            + ". An XSLOPE file does not say whether a model is to be run as a rapid "
              "drawdown, so the export writes one whenever the model has both lines "
              "and a d and psi. To solve it in GeoStudio as an ordinary analysis, turn "
              "its rapid drawdown off there")
        # GeoStudio derives the lowered pool from the lower surface, as XSLOPE's
        # automatic water loads do. A stage-2 load that is anything else has no place
        # in the file: SLOPE/W applies one set of surcharges in every stage.
        ponded2 = material_above_ground_dload(ground_surface, piezo_line2, gamma_water)
        extra2 = [b for b in (slope_data.get("dloads2") or [])
                  if not all(abs(p["Normal"] - _dload_pressure_at(ponded2, p["X"]))
                             <= 1e-6 * max(1.0, abs(p["Normal"])) for p in b)]
        if extra2:
            caveats.append(
                f"{len(extra2)} distributed load(s) after drawdown (the dloads (2) "
                f"sheet) were not written: SLOPE/W's staged rapid drawdown applies the "
                f"same surcharges in every stage, so only those before drawdown were "
                f"written")
    else:
        unused = []
        if drawdown:
            unused.append(
                f"the d and psi of {', '.join(repr(materials[k].get('name')) for k in drawdown)}")
        if len(piezo_line2) >= 2:
            unused.append("piezo line 2")
        if (slope_data.get("dloads2") or []) and str(
                slope_data.get("water_loads") or "").lower() == "manual":
            unused.append("the distributed loads after drawdown (the dloads (2) sheet)")
        if unused:
            have = [n for n, ln in (("line 1", piezo_line), ("line 2", piezo_line2))
                    if len(ln) >= 2]
            if seep_pts or any((m.get("u") or "") == "seep" for m in materials):
                why = ("this model's water is a seepage field, and SLOPE/W's staged "
                       "rapid drawdown takes the water before and after drawdown as two "
                       "piezometric surfaces")
            elif len(have) < 2:
                why = ("SLOPE/W's staged rapid drawdown takes the water before and "
                       "after drawdown as piezo lines 1 and 2, and this model has "
                       + (f"only {have[0]}" if have else "neither"))
            else:
                why = (f"only a rapid drawdown reads "
                       f"{'it' if unused == ['piezo line 2'] else 'them'}, and no "
                       f"Mohr-Coulomb material in this model carries d and psi")
            caveats.append(f"{' and '.join(unused)} "
                           f"{'was' if unused == ['piezo line 2'] else 'were'} not "
                           f"written: {why}")

    # SLOPE/W's phreatic correction is one switch on the analysis, for every
    # piezometric surface in it; XSLOPE's is a Type on each piezo line.
    switches = []
    phreatic = bool(slope_data.get("piezo_phreatic"))
    if piezo_ids and phreatic:
        switches.append("IsPiezPhreatic")
        spans = _line_on_ground(piezo_line, ground_surface)
        if spans:
            caveats.append(
                "SLOPE/W's phreatic correction was switched on for the piezometric "
                "surface. SLOPE/W leaves the head uncorrected where the line runs along "
                "the ground surface ("
                + ", ".join(f"x = {a:g} to {b:g}" for a, b in spans)
                + "), where XSLOPE corrects it, so SLOPE/W's pore pressure is higher "
                  "there")
    if rapid and phreatic != bool(slope_data.get("piezo_phreatic2")):
        caveats.append(
            f"piezo line 1 is {'' if phreatic else 'not '}phreatic and line 2 is "
            f"{'' if slope_data.get('piezo_phreatic2') else 'not '}. SLOPE/W's "
            f"phreatic correction is one switch for both surfaces, so it was written "
            f"{'on' if phreatic else 'off'} for both, as line 1 has it")
    if rapid:
        switches.append("IsRapidDrawdown")

    # A unit weight below the water table. SLOPE/W's <UnitWeight> is the saturated
    # weight, and a separate one above the water table is its <DryWeight> under the
    # <SlopeUseUnitWeightAboveWT> switch: XSLOPE's gamma_sat and gamma.
    split = [f"'{m.get('name')}' {float(m['gamma']):g} above and "
             f"{float(m['gamma_sat']):g} below" for m in materials if _has_split(m)]
    if split:
        caveats.append(
            f"{'this material has' if len(split) == 1 else 'these materials have'} a "
            f"different unit weight below the water table: {'; '.join(split)}. Each was "
            f"written as SLOPE/W's unit weight (below) and its separate unit weight above "
            f"the water table, which SLOPE/W splits where the pore pressure turns "
            f"negative")

    # v16/v17 material features SLOPE/W has no place for. They must NOT corrupt the file
    # (a .gsz material has no phi_b tag, and its TensileStrength is a LEM tension limit,
    # not xslope's FEM Rankine cut), so they are dropped from the export and reported.
    if any(m.get("t_cut") is not None for m in materials):
        caveats.append(
            "a tensile-strength cutoff (t_cut) is set on one or more materials — it was "
            "NOT written: xslope's t_cut is an FEM Rankine cap, which SLOPE/W's LEM has "
            "no equivalent for")
    if any(m.get("phi_b") is not None or m.get("s_cap") is not None for m in materials):
        caveats.append(
            "matric-suction shear strength (phi_b / s_cap) is set on one or more "
            "materials — it was NOT written: SLOPE/W parameterizes unsaturated strength "
            "as a suction ceiling on the piezometric surface. It has no per-material phi_b, "
            "so re-create it in GeoStudio if you need it")

    for key, label in [("circles", "failure circles"), ("non_circ", "non-circular surface"),
                       ("pile_lines", "piles")]:
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
    elif ponded and str(slope_data.get("water_loads") or "").lower() == "auto":
        # An automatic-water model has no water block to leave out -- its reservoir is
        # the piezometric line plus Water loads = auto. The line crosses, and GeoStudio
        # derives the load from it exactly as xslope's engine does, so nothing is lost;
        # but a reader who expects to find the load somewhere should be told where it
        # went.
        caveats.append(
            "this model takes its water loads automatically (Water loads = auto on the "
            "main sheet), so the reservoir was never a distributed load to write: the "
            "piezometric surface was written and GeoStudio derives the water's weight "
            "from it, the same way xslope does. Do not add it as a surcharge")

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

    # The analysis-local point list. The piezometric surface AND the surcharges both index
    # into it, so it has to be built before either is written. A surcharge is a wedge of
    # fill of unit weight gamma between its line and the ground, so a pressure profile
    # p(x) becomes the line y_ground(x) + p(x)/gamma. Any gamma reproduces p exactly; the
    # unit weight of water is used because it is already in the file and unit-consistent.
    local_pts = list(piezo_line)                        # numbered 1..n
    # A rapid drawdown's second surface, the water after drawdown, follows the first.
    piezo2_ids = []
    if rapid:
        local_pts += list(piezo_line2)
        piezo2_ids = list(range(len(piezo_line) + 1, len(local_pts) + 1))
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

    # Reinforcement. The products live at the top level and are shared; the LINES live on
    # the analysis and index the local point list, like everything else under <Entry>.
    products, reinf_lines = [], []
    for r in (slope_data.get("reinforcement_lines") or []):
        product, notes = _reinforcement_out(r, gamma_water)
        caveats.extend(notes)
        products.append(product)
        local_pts.append((float(r["x1"]), float(r["y1"])))
        p1 = len(local_pts)
        local_pts.append((float(r["x2"]), float(r["y2"])))
        reinf_lines.append((len(products), p1, len(local_pts)))
    if reinf_lines:
        caveats.append(
            f"{len(reinf_lines)} reinforcement line(s) written. GeoStudio quotes capacity "
            f"per element with an out-of-plane spacing; xslope's are already per unit "
            f"width, so they are written with a spacing of 1")

    # Line loads. The point is a local DataPoint; the direction is a trend/plunge pair.
    load_pts = []
    for ll in (slope_data.get("line_loads") or []):
        local_pts.append((float(ll["x"]), float(ll["y"])))
        trend, plunge = _line_load_out(ll)
        load_pts.append((len(local_pts), float(ll.get("P") or 0.0), trend, plunge))
        if plunge < 0:
            caveats.append(
                f"line load {ll.get('label') or ''!r} points UPWARD "
                f"({ll.get('angle')} degrees); GeoStudio's plunge is measured below "
                f"horizontal, so it was written as {plunge:g} — check it")
    if load_pts:
        caveats.append(f"{len(load_pts)} line load(s) written")

    # Probabilistic standard deviations, the reverse of the import mapping. SLOPE/W keeps
    # them on the ANALYSIS as modifiers of the DEVIATION from the applied value, so an
    # xslope sigma is written as Normal(Mean=0,SD=sigma) against the material's ID.
    prob_mods, prob_dropped = [], []
    for i, m in enumerate(materials, start=1):
        elastic = (m.get("option") or "") == "elastic"
        undrained = (m.get("option") or "") in ("ratio", "cp")
        for field, tag in (("sigma_gamma", "Applied.UnitWeight"),
                           ("sigma_c", "Applied.EffectiveCohesion"),
                           ("sigma_phi", "Applied.EffectivePhi")):
            sigma = float(m.get(field) or 0.0)
            if sigma <= 0:
                continue
            if elastic and field != "sigma_gamma":
                # Written as Bedrock, which has no strength to vary.
                prob_dropped.append(f"'{m.get('name')}' {field} (elastic -> Bedrock)")
                continue
            if undrained and field != "sigma_gamma":
                # SHANSEP and S=f(depth) have no effective c or phi to vary.
                prob_dropped.append(f"'{m.get('name')}' {field} (no SLOPE/W "
                                    f"equivalent for option {m.get('option')})")
                continue
            prob_mods.append((i, tag, sigma))
        for field in ("sigma_cp", "sigma_d", "sigma_psi"):
            if float(m.get(field) or 0.0) > 0:
                prob_dropped.append(f"'{m.get('name')}' {field} (no SLOPE/W equivalent)")
    sample_lines = []
    if prob_mods:
        from .reliability import MC_DEFAULT_SAMPLES
        sample_lines = ['      <SampleSettings>',
                        '        <Type>Probabilistic</Type>',
                        f'        <NumMonteCarloTrials>{MC_DEFAULT_SAMPLES}</NumMonteCarloTrials>',
                        f'        <ProbabilisticModifiers Len="{len(prob_mods)}">']
        for mid, tag, sigma in prob_mods:
            sample_lines += [
                '          <Modifier>',
                f'            <Key>Materials[{mid}].{tag}</Key>',
                f'            <Entry Function="Normal(Mean=0,SD={sigma:.10g})" />',
                '          </Modifier>']
        sample_lines += ['        </ProbabilisticModifiers>', '      </SampleSettings>']
        caveats.append(
            f"{len(prob_mods)} standard deviation(s) written as a PROBABILISTIC SLOPE/W "
            f"analysis ({MC_DEFAULT_SAMPLES} Monte Carlo trials). They are written "
            f"normal, uncorrelated and untruncated — xslope states no correlation or "
            f"sampling bounds, which SLOPE/W can carry")
    if prob_dropped:
        caveats.append(
            f"standard deviation(s) not written: {', '.join(prob_dropped)}")

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
         *sample_lines,
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
         '  </Coordinates>']
    # The seepage field, as GeoStudio's spatial function. It lives at the top level
    # beside the material functions -- not on the analysis -- and the analysis points at
    # it by GGID from its WaterItem. <Elements>/<Nodes> (a cached triangulation of the
    # point cloud) are omitted: GeoStudio's own files carry a used 3DFunction with the
    # points alone, so it builds the interpolant itself.
    if seep_pts:
        L.append('  <Functions>')
        L.append('    <Func3Ds Len="1">')
        L.append('      <Fn3D>')
        L.append('        <ID>1</ID>')
        L.append('        <Name>xslope seepage field</Name>')
        L.append('        <OutputType>PressureHead</OutputType>')
        L.append('        <Model>LinearInterpolation</Model>')
        L.append(f'        <Points Len="{len(seep_pts)}">')
        for x, y, z in seep_pts:
            L.append(f'          <Points_ X="{x:.10g}" Y="{y:.10g}" Z="{z:.10g}" />')
        L.append('        </Points>')
        L.append('      </Fn3D>')
        L.append('    </Func3Ds>')
        L.append('  </Functions>')
    L += ['  <Geometries Len="1">',
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
    # own view; writing a guess left the model off-center with the axes in shot.
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
        # An xslope 'elastic' material is written as GeoStudio's Bedrock (Impenetrable)
        # model — the exact reverse of the Bedrock -> elastic import mapping. Its unit
        # weight and (harmless) c/phi ride along in the same StressStrain block, exactly
        # as a real Bedrock material carries them.
        slope_model, fields = strength_out[i - 1]
        L.append(
            f'    <Material><ID>{i}</ID><Color>RGB=({r8},{g8},{b8})</Color>'
            f'<Name>{_xml_escape(m.get("name") or f"material {i}")}</Name>'
            f'<SlopeModel>{slope_model}</SlopeModel>'
            f'<StressModel>MohrCoulomb</StressModel>'
            f'<StressStrain>'
            f'<JointEffectiveCohesion Missing="true" />'
            f'<IntactRockParam Missing="true" />'
            f'<GeologicalStrengthIndex Missing="true" />'
            f'<DisturbanceFactor Missing="true" />'
            f'<MaxConfiningStress Missing="true" />'
            f'<UseTensileStrength>false</UseTensileStrength>'
            + _unit_weight_out(m, fields)
            + f'<TensileStrength Missing="true" />'
            f'<OCRatio Missing="true" />'
            f'</StressStrain></Material>')
    L.append('  </Materials>')

    # The reinforcement PRODUCT library, shared across analyses (the lines that use it
    # live on the StabilityItem). GeoStudio stores no direction and no F-of-S dependence
    # here -- it infers both from <Type> -- so the type is the only thing carrying them.
    if products:
        L.append(f'  <Reinforcements Len="{len(products)}">')
        for i, p in enumerate(products, start=1):
            L.append(f'    <Reinforcement><ID>{i}</ID>'
                     f'<Name>{_xml_escape(p["name"])}</Name>'
                     f'<Type>{p["type"]}</Type>'
                     f'<Color>RGB=(85,0,192)</Color>'
                     f'<Spacing>{p["spacing"]:.10g}</Spacing>'
                     f'<Tensile>{p["tensile"]:.10g}</Tensile>'
                     f'<PlateCapacity>{p["plate"]:.10g}</PlateCapacity>'
                     f'<PulloutResistance>{p["pullout"]:.10g}</PulloutResistance>'
                     + (f'<InterfaceAdhesion>{p["interface"][0]:.10g}'
                        f'</InterfaceAdhesion>'
                        f'<InterfaceShearAngle>{p["interface"][1]:.10g}'
                        f'</InterfaceShearAngle>'
                        f'<SurfaceAreaFactor>2</SurfaceAreaFactor>'
                        f'<CalculatedPulloutResistance>true'
                        f'</CalculatedPulloutResistance>'
                        if p.get("interface") else '')
                     + f'</Reinforcement>')
        L.append('  </Reinforcements>')

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
    # Which source the analysis draws pore pressure FROM. One per analysis: a seepage
    # field written as a spatial function takes precedence over a piezometric surface
    # (reported above), because it is the more informative of the two.
    if seep_pts:
        water_src = ('<ResultInputInfo><Option>3DFunction</Option>'
                     '<DataGGID>3DFns-1</DataGGID></ResultInputInfo>')
    elif piezo_ids:
        water_src = '<ResultInputInfo><Option>PiezoSurface</Option></ResultInputInfo>'
    else:
        water_src = ''
    L.append('  <WaterItems Len="1">')
    L.append('    <WaterItem><AnalysisID>1</AnalysisID><Entry>'
             + water_src
             + f'<UnitWaterWeight>{gamma_water:.10g}</UnitWaterWeight>'
             + f'<WaterBulkMod>{bulk:.9g}</WaterBulkMod>'
             '</Entry></WaterItem>')
    L.append('  </WaterItems>')

    # The piezometric surface belongs to the StabilityItem, NOT the geometry, and the
    # materials that draw pore pressure from it are named in <MaterialUsesPiezs>.
    k = float(slope_data.get("k_seismic") or 0.0)

    L.append('  <StabilityItems Len="1">')
    L.append('    <StabilityItem><AnalysisID>1</AnalysisID><Entry>')
    if local_pts:
        L.append(f'      <DataPoints Len="{len(local_pts)}">')
        for i, (x, y) in enumerate(local_pts, start=1):
            L.append(f'        <DataPoint Number="{i}" X="{x:.10g}" Y="{y:.10g}" />')
        L.append('      </DataPoints>')
    # The analysis switches, under one <Control> ahead of the loads and the
    # piezometric surfaces, as SLOPE/W's own files write them: the phreatic correction
    # (<IsPiezPhreatic>, as on Greco's weak-layer slope) and the staged rapid drawdown
    # (<IsRapidDrawdown>, as on its four rapid drawdown examples).
    if switches:
        L.append('      <Control>')
        for tag in switches:
            L.append(f'        <{tag}>true</{tag}>')
        L.append('      </Control>')
    # SLOPE/W's vertical coefficient is positive downward, as xslope's is.
    kv = float(slope_data.get("k_vertical") or 0.0)
    L.append(f'      <Seismic Horizontal="{k:.10g}" Vertical="{kv:.10g}" />'
             if kv else
             f'      <Seismic Horizontal="{k:.10g}" Vertical="" />' if k else
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
    if reinf_lines:
        L.append(f'      <ReinforcementLines Len="{len(reinf_lines)}">')
        for i, (pid, a, b) in enumerate(reinf_lines, start=1):
            L.append(f'        <ReinforcementLine ID="{i}" Reinforcement="{pid}" '
                     f'Point1Id="{a}" Point2Id="{b}" />')
        L.append('      </ReinforcementLines>')
    if load_pts:
        L.append(f'      <LineLoadPoints Len="{len(load_pts)}">')
        for pid, value, trend, plunge in load_pts:
            L.append(f'        <LineLoadPoint><ID>{pid}</ID><LineLoad>'
                     f'<Value>{value:.10g}</Value>'
                     f'<Direction Trend="{trend:.10g}" Plunge="{plunge:.10g}" />'
                     f'</LineLoad></LineLoadPoint>')
        L.append('      </LineLoadPoints>')
    if piezo_ids:
        # A staged rapid drawdown has two: SLOPE/W reads the higher as the water before
        # drawdown and the lower as the water after it. Each material that takes pore
        # pressure from a line is assigned the first, as in SLOPE/W's own examples.
        surfaces = [piezo_ids] + ([piezo2_ids] if piezo2_ids else [])
        L.append(f'      <PiezometricSurfaces Len="{len(surfaces)}">')
        for sid, ids in enumerate(surfaces, start=1):
            L.append(f'        <PiezometricSurface><ID>{sid}</ID>'
                     f'<DataPoints Len="{len(ids)}">'
                     + "".join(f'<DataPoint>{i}</DataPoint>' for i in ids)
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


def import_gsz(gsz_path, template, out_path, analysis_id=None, step=None):
    """Read a .gsz and write it to an xslope .xlsx input file.

    The script-level route, mirroring :func:`xslope.cad.import_dxf`: parse the
    GeoStudio model, convert one analysis to ``slope_data``, and write it through
    a copy of an input template. Pass ``template=None`` for the current one, which
    is what an import normally wants: a GeoStudio model states its water as a
    surface rather than as a load, and only a v22-or-later template carries the
    cell that says so (the writer refuses an older one rather than dropping the
    reservoir).

    A pore-pressure FIELD from a parent SEEP/W analysis does not fit in a spreadsheet,
    so it is written beside the .xlsx as the two sidecar files xslope reads it from:
    ``<name>_mesh.json`` and ``<name>_seep.csv``. Without them the .xlsx would say every
    material takes its pore pressure from a seepage solution, and there would be none —
    a dry model that looks like a wet one.

    Returns the caveat list from :func:`gsz_to_slope_data`. Read it — a clean
    return does not mean a complete model (no failure surface is ever imported).
    """
    from .fileio import save_slope_data_to_xlsx
    from .mesh import export_mesh_to_json
    from .seep import export_seep_u

    gsz = read_gsz(gsz_path)
    slope_data, caveats = gsz_to_slope_data(gsz, analysis_id, step=step)
    save_slope_data_to_xlsx(slope_data, out_path, template=template)

    mesh, u = slope_data.get("mesh"), slope_data.get("seep_u")
    if mesh is not None and u is not None:
        base = os.path.splitext(out_path)[0]
        export_mesh_to_json(mesh, f"{base}_mesh.json")
        export_seep_u(mesh["nodes"], u, f"{base}_seep.csv",
                      require_gamma_water(slope_data, "GeoStudio seep export"))
        caveats.append(
            f"the SEEP/W pore-pressure field was written beside the input file as "
            f"{os.path.basename(base)}_mesh.json and {os.path.basename(base)}_seep.csv "
            f"— xslope reads a seepage solution from those, so keep all three together "
            f"or the model reverts to dry")
    return caveats
