"""Standing checks on six things the GeoStudio (.gsz) import and export used to drop.

  1. THE SATURATED UNIT WEIGHT ON EXPORT. A material whose gamma_sat differs from its
     gamma went to GeoStudio at gamma everywhere. SLOPE/W's <UnitWeight> is the
     saturated weight, and a separate one above the water table is its <DryWeight>
     under <SlopeUseUnitWeightAboveWT>ConstantValue: written so, in the order SLOPE/W's
     own Unit Weight Definition example writes them, and read back to the same model.
  2. THE UNIT WEIGHT ABOVE THE WATER TABLE ON IMPORT. With the switch on, <DryWeight>
     is XSLOPE's gamma and <UnitWeight> its gamma_sat; without it SLOPE/W stores a
     <DryWeight> and does not use it; a weight computed from the water content is
     reported. On SLOPE/W's own example the split reproduces its saved slice weights.
  3. THE TIE-BACK CAPACITY. An earlier GeoStudio 2025 build wrote a reinforcement's
     capacity as <TensileCapacity> with no plate capacity, which imported as zero.
     SLOPE/W's Tie-Back Wall example (verification manual 2.35) comes that way.
  4. THE PHREATIC SWITCH ON EXPORT. A phreatic piezo line went out as a plain one and
     came back plain. It is SLOPE/W's <Control><IsPiezPhreatic>.
  5. PARTIAL FACTORS. An analysis solved with a Eurocode 7 set of partial factors
     imported with only "an unrecognized 'Control'". The caveat now names the set and
     every factor in it that is not 1, and says XSLOPE does not apply them.
  6. THE SKILL. It said no import carries a failure surface; a solved .gsz carries
     SLOPE/W's critical circle. Its two copies must stay byte-identical.

SLOPE/W's own files are Seequent's copyrighted material and are not in this
repository. The checks that read them look for local copies (``GSZ_VENDOR_CORPUS``
and the 2025.2 examples zip in ``GSZ_DOWNLOADS``) and skip cleanly when they are
absent; every item also has a check on a file built here.
"""
import contextlib
import io
import os
import statistics
import sys
import tempfile
import xml.etree.ElementTree as ET
import zipfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

#: Local folders of GeoStudio files (git-ignored, not redistributable).
GSZ_VENDOR_CORPUS = os.path.expanduser("~/python_projects/vendor_files/gsz_corpus")
GSZ_DOWNLOADS = os.path.expanduser("~/python_projects/vendor_files/geostudio_downloads")
SLOPEW_EXAMPLES_ZIP = "2025.2 SLOPEW Examples.zip"
UNIT_WEIGHT_EXAMPLE = "SLOPE2D/Unit Weight Definition.gsz"

SKILL = os.path.join(_REPO, "docs", "usage", "claude", "xslope.md")
SKILL_COPY = os.path.join(_REPO, "xslope", "resources", "xslope_skill.md")

ST = "./Materials/Material"


def _close(a, b, tol=1e-9):
    return abs(float(a) - float(b)) <= tol * max(1.0, abs(float(b)))


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _vendor(name):
    p = os.path.join(GSZ_VENDOR_CORPUS, name)
    return p if os.path.isfile(p) else None


def _unit_weight_example(td):
    """SLOPE/W's Unit Weight Definition example, written to ``td``, or None."""
    return _slopew_example(td, UNIT_WEIGHT_EXAMPLE)


def _slopew_example(td, member):
    """One file of SLOPE/W's 2025.2 examples zip, written to ``td``, or None."""
    z = os.path.join(GSZ_DOWNLOADS, SLOPEW_EXAMPLES_ZIP)
    if not os.path.isfile(z):
        return None
    try:
        data = zipfile.ZipFile(z).read(member)
    except KeyError:
        return None
    p = os.path.join(td, os.path.basename(member).replace(" ", "_"))
    with open(p, "wb") as fh:
        fh.write(data)
    return p


def _critical_fs(g, analysis, sd, solver, n=30):
    """``(XSLOPE, SLOPE/W)`` factors of safety on SLOPE/W's critical circle of one
    analysis, or None when it has no circle or XSLOPE cannot solve it. A circle
    drawn through a vertex of the ground surface (an Entry and Exit search whose exit
    is a point) is moved off it by a micrometre."""
    from xslope.geostudio import read_gsz_results
    from xslope.slice import generate_slices
    rows = [r for r in read_gsz_results(g, analysis["name"]) if r["circular"]]
    if not rows:
        return None
    best = min(rows, key=lambda t: t["fs"])
    for eps in (0.0, -1e-6, 1e-6):
        circle = {"Xo": best["xo"], "Yo": best["yo"], "R": best["r"] * (1.0 + eps),
                  "Depth": None, "Y": None, "Option": "Depth",
                  "Movement": "Left to Right"}
        ok, res = _quiet(generate_slices, sd, circle=circle, num_slices=n,
                         debug=False, check_inputs=False)
        if ok:
            ok2, r = _quiet(solver, res[0])
            return (r["FS"], best["fs"]) if ok2 else None
    return None


def _built_model(gamma_sat=None, phreatic=False):
    """Two zones, a fill over a foundation, with a piezometric line."""
    from shapely.geometry import Polygon
    from xslope.geostudio import _blank_material
    from xslope.fileio import build_ground_surface_from_polygons
    fill = dict(_blank_material("fill"), option="mc", gamma=18.0, gamma_sat=gamma_sat,
                c=8.0, phi=30.0, u="piezo")
    found = dict(_blank_material("foundation"), option="mc", gamma=19.0, c=15.0,
                 phi=28.0, u="piezo")
    polygons = [
        {"mat_id": 1, "polygon": Polygon([(0, 0), (60, 0), (60, 8), (15, 8), (0, 8)])},
        {"mat_id": 0, "polygon": Polygon([(15, 8), (60, 8), (60, 20), (35, 20)])},
    ]
    ground, domain = build_ground_surface_from_polygons(polygons)
    return {"polygons": polygons, "materials": [fill, found], "gamma_water": 9.807,
            "unit_system": "si", "water_loads": "auto", "ground_surface": ground,
            "domain_polygon": domain,
            "piezo_line": [(0.0, 10.0), (20.0, 10.0), (60.0, 15.0)], "piezo_line2": [],
            "piezo_phreatic": phreatic, "piezo_phreatic2": phreatic, "k_seismic": 0.0,
            "k_vertical": 0.0, "tcrack_depth": 0.0, "tcrack_water": 0.0,
            "tcrack_surface": None, "max_depth": None, "profile_lines": [],
            "circular": True, "circles": [dict(CIRCLE)], "non_circ": [],
            "reinforcement_lines": [], "reinforce_lines": [], "pile_lines": [],
            "line_loads": [], "mesh": None, "dloads": [], "dloads2": []}


CIRCLE = {"Xo": 30.0, "Yo": 32.0, "R": 26.0, "Depth": 6.0}


def _export(sd, td, stem):
    from xslope.geostudio import export_gsz
    out = os.path.join(td, f"{stem}.gsz")
    cav = _quiet(export_gsz, sd, out, analysis_name=stem)
    return out, cav


def _root(path):
    zf = zipfile.ZipFile(path)
    name = next(n for n in zf.namelist() if "/" not in n and n.lower().endswith(".xml"))
    return ET.fromstring(zf.read(name)), name


def _mutate(src, dst, fn):
    """Copy the .gsz at ``src`` to ``dst`` with ``fn(root)`` applied to its model XML."""
    root, name = _root(src)
    fn(root)
    with zipfile.ZipFile(dst, "w") as zf:
        zf.writestr(name, ET.tostring(root, encoding="unicode"))
    return dst


def _import(path, aid=None, critical_surface=True):
    from xslope.geostudio import read_gsz, gsz_to_slope_data
    g = read_gsz(path)
    return g, _quiet(gsz_to_slope_data, g, aid, critical_surface=critical_surface)


def _mat(sd, name):
    return next((m for m in sd["materials"] if m["name"] == name), None)


def _weight(sd, circle, n=50):
    from xslope.slice import generate_slices
    ok, res = _quiet(generate_slices, sd, circle=circle, num_slices=n, debug=False,
                     check_inputs=False)
    return float(res[0]["w"].sum()) if ok else None


# --------------------------------------------------------------------------------
# 1. gamma_sat on export
# --------------------------------------------------------------------------------

def check_gamma_sat_export():
    fails = []
    with tempfile.TemporaryDirectory() as td:
        sd = _built_model(gamma_sat=20.0)
        out, _ = _export(sd, td, "split")
        root, _ = _root(out)
        ss = {m.findtext("Name"): m.find("StressStrain") for m in root.findall(ST)}
        f = ss["fill"]
        tags = [c.tag for c in f]
        if f.findtext("SlopeUseUnitWeightAboveWT") != "ConstantValue":
            fails.append("1: a material with gamma_sat was written without SLOPE/W's "
                         "<SlopeUseUnitWeightAboveWT>ConstantValue")
        if not (f.findtext("UnitWeight") and _close(f.findtext("UnitWeight"), 20.0)):
            fails.append(f"1: <UnitWeight> written as {f.findtext('UnitWeight')}, not "
                         f"gamma_sat 20")
        if not (f.findtext("DryWeight") and _close(f.findtext("DryWeight"), 18.0)):
            fails.append(f"1: <DryWeight> written as {f.findtext('DryWeight')}, not "
                         f"gamma 18")
        if all(t in tags for t in ("SlopeUseUnitWeightAboveWT", "UnitWeight", "PhiPrime",
                                   "DryWeight")):
            if not (tags.index("SlopeUseUnitWeightAboveWT") < tags.index("UnitWeight")
                    < tags.index("PhiPrime") < tags.index("DryWeight")):
                fails.append(f"1: the unit weight elements are in the order {tags}, not "
                             f"switch, UnitWeight, strengths, DryWeight")
        fo = ss["foundation"]
        if fo.find("DryWeight") is not None or fo.find("SlopeUseUnitWeightAboveWT") is not None:
            fails.append("1: a material with no gamma_sat was given a weight above the "
                         "water table")
        _g, (back, _c) = _import(out, 1)
        b = _mat(back, "fill")
        if not (_close(b["gamma"], 18.0) and b.get("gamma_sat") is not None
                and _close(b["gamma_sat"], 20.0)):
            fails.append(f"1: round trip gave gamma {b['gamma']}, gamma_sat "
                         f"{b.get('gamma_sat')}, not 18 and 20")
        if _mat(back, "foundation").get("gamma_sat") is not None:
            fails.append("1: round trip gave the foundation a gamma_sat")
        w0, w1 = _weight(sd, CIRCLE), _weight(back, CIRCLE)
        if w0 is None or w1 is None or not _close(w0, w1, 1e-9):
            fails.append(f"1: the sliding mass weighs {w0} before the round trip and "
                         f"{w1} after")

        # SLOPE/W's own order, where its example is at hand.
        ex = _unit_weight_example(td)
        if ex is not None:
            vroot, _ = _root(ex)
            v = next(m.find("StressStrain") for m in vroot.findall(ST)
                     if m.find("StressStrain").findtext("SlopeUseUnitWeightAboveWT")
                     == "ConstantValue")
            keys = ("SlopeUseUnitWeightAboveWT", "UnitWeight", "DryWeight")
            if [t for t in [c.tag for c in v] if t in keys] != \
                    [t for t in tags if t in keys]:
                fails.append("1: the unit weight elements are not in the order SLOPE/W's "
                             "Unit Weight Definition example writes them")
    return fails


# --------------------------------------------------------------------------------
# 2. The unit weight above the water table on import
# --------------------------------------------------------------------------------

def _set_fill(root, **fields):
    for m in root.findall(ST):
        if m.findtext("Name") != "fill":
            continue
        ss = m.find("StressStrain")
        for tag, text in fields.items():
            e = ss.find(tag)
            if e is None:
                e = ET.SubElement(ss, tag)
            e.text = text


def check_unit_weight_import():
    fails = []
    with tempfile.TemporaryDirectory() as td:
        base, _ = _export(_built_model(), td, "base")
        cases = {
            "switch": dict(SlopeUseUnitWeightAboveWT="ConstantValue", UnitWeight="20",
                           DryWeight="15"),
            "stored": dict(UnitWeight="20", DryWeight="15"),
            "water content": dict(SlopeUseUnitWeightAboveWT="CalculateFromVolFn",
                                  UnitWeight="20", DryWeight="20"),
        }
        got = {}
        for name, fields in cases.items():
            p = _mutate(base, os.path.join(td, f"{name}.gsz"),
                        lambda r, f=fields: _set_fill(r, **f))
            _g, got[name] = _import(p, 1)
        m = _mat(got["switch"][0], "fill")
        if not (_close(m["gamma"], 15.0) and m.get("gamma_sat") is not None
                and _close(m["gamma_sat"], 20.0)):
            fails.append(f"2: the switch on imported gamma {m['gamma']}, gamma_sat "
                         f"{m.get('gamma_sat')}, not 15 and 20")
        m, cav = _mat(got["stored"][0], "fill"), got["stored"][1]
        if not (_close(m["gamma"], 20.0) and m.get("gamma_sat") is None):
            fails.append(f"2: a stored <DryWeight> without the switch imported gamma "
                         f"{m['gamma']}, gamma_sat {m.get('gamma_sat')}; SLOPE/W does "
                         f"not use it")
        if any("above the water table" in c for c in cav):
            fails.append("2: a stored <DryWeight> that SLOPE/W does not use drew a caveat")
        m, cav = _mat(got["water content"][0], "fill"), got["water content"][1]
        if not _close(m["gamma"], 20.0) or \
                not any("volumetric water content" in c for c in cav):
            fails.append("2: a weight computed from the water content was not reported")

        # SLOPE/W's own example: Case 1 has no switch, Case 2 is 16 above and 20 below.
        # Its saved weight of every sliding mass XSLOPE can build is the check.
        ex = _unit_weight_example(td)
        if ex is not None:
            from xslope.geostudio import read_gsz, gsz_to_slope_data, read_gsz_results
            g = read_gsz(ex)
            for a in g["analyses"]:
                if not a["name"].startswith(("Case 1", "Case 2")):
                    continue
                sd, _c = _quiet(gsz_to_slope_data, g, a["id"], critical_surface=False)
                m = sd["materials"][0]
                want = (20.0, None) if a["name"].startswith("Case 1") else (16.0, 20.0)
                if not _close(m["gamma"], want[0]) or (
                        (m.get("gamma_sat") is None) != (want[1] is None)):
                    fails.append(f"2: {a['name']!r} imported gamma {m['gamma']}, "
                                 f"gamma_sat {m.get('gamma_sat')}, expected {want}")
                    continue
                errs = []
                for r in read_gsz_results(g, a["name"]):
                    if not r["circular"] or not r["weight"]:
                        continue
                    w = _weight(sd, {"Xo": r["xo"], "Yo": r["yo"], "R": r["r"],
                                     "Depth": r["yo"] - r["r"]})
                    if w is not None:
                        errs.append(100.0 * (w / r["weight"] - 1.0))
                if len(errs) < 3 or abs(statistics.median(errs)) > 0.2:
                    fails.append(f"2: {a['name']!r}: the sliding masses weigh "
                                 f"{statistics.median(errs) if errs else float('nan'):+.2f}% "
                                 f"against SLOPE/W's saved weights ({len(errs)} circles)")
    p = _vendor("Staged rapid drawdown - Pumped Storage Project Dam.gsz")
    if p is not None:
        _g, (sd, cav) = _import(p)
        if any(m.get("gamma_sat") is not None for m in sd["materials"]) or \
                any("above the water table" in c and "imported as gamma" in c for c in cav):
            fails.append("2: the pumped-storage dam's <DryWeight>, which SLOPE/W does not "
                         "use there, was applied")
    return fails


# --------------------------------------------------------------------------------
# 3. The tie-back capacity
# --------------------------------------------------------------------------------

def _add_reinforcement(root, product):
    """A reinforcement product and one line through the fill, from (40, 18) inward."""
    lib = root.find("Reinforcements")
    if lib is None:
        lib = ET.SubElement(root, "Reinforcements", Len="0")
    rf = ET.SubElement(lib, "Reinforcement")
    for tag, text in product.items():
        ET.SubElement(rf, tag).text = text
    lib.set("Len", str(len(lib)))
    entry = root.find("./StabilityItems/StabilityItem/Entry")
    dps = entry.find("DataPoints")
    n = max(int(d.get("Number")) for d in dps)
    for k, (x, y) in enumerate(((40.0, 18.0), (55.0, 12.0)), start=1):
        ET.SubElement(dps, "DataPoint", Number=str(n + k), X=str(x), Y=str(y))
    dps.set("Len", str(len(dps)))
    lines = ET.SubElement(entry, "ReinforcementLines", Len="1")
    ET.SubElement(lines, "ReinforcementLine", ID="1", Reinforcement=product["ID"],
                  Point1Id=str(n + 1), Point2Id=str(n + 2))


def check_tieback_capacity():
    fails = []
    with tempfile.TemporaryDirectory() as td:
        base, _ = _export(_built_model(), td, "base")
        old = {"Name": "older anchor", "ID": "1", "Spacing": "2", "BondLength": "10",
               "TensileCapacity": "1000", "PulloutResistance": "100"}
        p = _mutate(base, os.path.join(td, "older.gsz"),
                    lambda r: _add_reinforcement(r, old))
        _g, (sd, cav) = _import(p, 1)
        rl = (sd.get("reinforcement_lines") or [{}])[0]
        if not _close(rl.get("t_max") or 0.0, 500.0):
            fails.append(f"3: a <TensileCapacity> of 1000 at spacing 2 imported as t_max "
                         f"{rl.get('t_max')}, not 500")
        if not _close(max(rl.get("tend1") or 0.0, rl.get("tend2") or 0.0), 500.0):
            fails.append("3: the plate capacity an older file leaves out was not taken as "
                         "its tensile capacity")
        if not any("TensileCapacity" in c for c in cav):
            fails.append("3: the older capacity tag was not reported")
        new = {"Name": "anchor", "ID": "1", "Spacing": "2", "Tensile": "1000",
               "PlateCapacity": "400", "PulloutResistance": "100"}
        p = _mutate(base, os.path.join(td, "newer.gsz"),
                    lambda r: _add_reinforcement(r, new))
        _g, (sd, cav) = _import(p, 1)
        rl = (sd.get("reinforcement_lines") or [{}])[0]
        if not (_close(rl.get("t_max") or 0.0, 500.0) and
                _close(max(rl.get("tend1") or 0.0, rl.get("tend2") or 0.0), 200.0)):
            fails.append("3: <Tensile> and <PlateCapacity> no longer import as written")

    path = _vendor("Pockoski and Duncan - Tie-Back Wall.gsz")
    if path is not None:
        from xslope.geostudio import read_gsz_results
        from xslope.slice import generate_slices
        from xslope.solve import spencer
        g, (sd, _cav) = _import(path)
        # <TensileCapacity> 247343.87 over <Spacing> 4 is 61836; the anchor's bond,
        # <PulloutResistance> 4000 / 4 over <BondLength> 40, is 40000 and governs (the
        # 40,000 lb/ft of wall verification manual 2.35 quotes).
        cap = min(247343.87, 4000.0 * 40.0) / 4.0
        for rl in sd.get("reinforcement_lines") or [{}]:
            if not _close(rl.get("t_max") or 0.0, cap, 1e-6):
                fails.append(f"3: the Tie-Back Wall's anchors imported at t_max "
                             f"{rl.get('t_max')}, not {cap:.1f}")
                break
        best = min(read_gsz_results(g, g["analyses"][0]["name"]), key=lambda t: t["fs"])
        ok, res = _quiet(generate_slices, sd, circle=sd["circles"][0], num_slices=30,
                         check_inputs=False)
        ok2, r = _quiet(spencer, res[0]) if ok else (False, None)
        if not ok2 or abs(r["FS"] / best["fs"] - 1.0) > 0.01:
            fails.append(f"3: the Tie-Back Wall solves to "
                         f"{r['FS'] if ok2 else r} on SLOPE/W's circle, against SLOPE/W's "
                         f"{best['fs']:.4f}")
    return fails


# --------------------------------------------------------------------------------
# 4. The phreatic switch on export
# --------------------------------------------------------------------------------

def check_phreatic_export():
    fails = []
    with tempfile.TemporaryDirectory() as td:
        for phreatic in (True, False):
            out, _ = _export(_built_model(phreatic=phreatic), td, f"ph{int(phreatic)}")
            root, _ = _root(out)
            entry = root.find("./StabilityItems/StabilityItem/Entry")
            on = entry.findtext("./Control/IsPiezPhreatic") == "true"
            if on != phreatic:
                fails.append(f"4: a {'phreatic' if phreatic else 'plain'} line was "
                             f"written with IsPiezPhreatic {'off' if phreatic else 'on'}")
            tags = [c.tag for c in entry]
            if on and tags.index("Control") > tags.index("PiezometricSurfaces"):
                fails.append("4: <Control> is written after the piezometric surfaces")
            _g, (back, _c) = _import(out, 1)
            if bool(back.get("piezo_phreatic")) != phreatic:
                fails.append(f"4: a {'phreatic' if phreatic else 'plain'} line came back "
                             f"{'plain' if phreatic else 'phreatic'}")
        path = _vendor("Greco Weak Layer Slope.gsz")
        if path is not None:
            _g, (sd, _c) = _import(path)
            out, _ = _export(sd, td, "greco")
            _g, (back, _c) = _import(out, 1)
            if not (sd.get("piezo_phreatic") and back.get("piezo_phreatic")):
                fails.append("4: Greco's weak-layer slope, which SLOPE/W solves with the "
                             "phreatic correction, lost it on the round trip")
    return fails


# --------------------------------------------------------------------------------
# 5. Partial factors
# --------------------------------------------------------------------------------

def _add_design_factors(root):
    ds = ET.SubElement(root, "DesignFactors", Len="1")
    d = ET.SubElement(ds, "DesignFactors_")
    ET.SubElement(d, "ID").text = "1"
    ET.SubElement(d, "Name").text = "Test approach"
    fs = ET.SubElement(d, "Factors", Len="4")
    for key, fav, unfav in (("PermanentLoads", "1", "1"), ("VariableLoads", None, "1.3"),
                            ("EffectiveCohesion", "1.25", ""),
                            ("EffectiveTanPhi", "1.25", "")):
        f = ET.SubElement(fs, "Factors_")
        ET.SubElement(f, "Key").text = key
        e = ET.SubElement(f, "Entry", Unfavorable=unfav)
        if fav is not None:
            e.set("Favorable", fav)
    entry = root.find("./StabilityItems/StabilityItem/Entry")
    ctl = entry.find("Control")
    if ctl is None:
        ctl = ET.SubElement(entry, "Control")
    ET.SubElement(ctl, "DesignFactorID").text = "1"


def check_partial_factors():
    fails = []
    with tempfile.TemporaryDirectory() as td:
        base, _ = _export(_built_model(), td, "base")
        p = _mutate(base, os.path.join(td, "factored.gsz"), _add_design_factors)
        _g, (_sd, cav) = _import(p, 1)
        note = next((c for c in cav if "partial factors" in c), "")
        for want in ("'Test approach'", "c' 1.25", "tan phi' 1.25",
                     "variable loads 1.3 (unfavorable)", "applies no partial factors"):
            if want not in note:
                fails.append(f"5: the partial factor caveat does not say {want!r}")
        if "permanent loads" in note:
            fails.append("5: a factor of 1 was named")
        if any("'Control'" in c for c in cav):
            fails.append("5: the partial factors still draw the unrecognized 'Control' "
                         "caveat")
    path = _vendor("Eurocode Cutting in Clay.gsz")
    if path is not None:
        _g, (_sd, cav) = _import(path)
        note = next((c for c in cav if "partial factors" in c), "")
        if "'Eurocode 7 - DA3'" not in note or "c' 1.25" not in note:
            fails.append("5: the Eurocode cutting's DA3 factors are not named")
    return fails


# --------------------------------------------------------------------------------
# 6. The skill
# --------------------------------------------------------------------------------

def check_skill():
    fails = []
    a, b = open(SKILL, "rb").read(), open(SKILL_COPY, "rb").read()
    if a != b:
        fails.append("6: docs/usage/claude/xslope.md and xslope/resources/xslope_skill.md "
                     "differ")
    text = a.decode("utf-8")
    if "No import ever carries a failure surface" in text:
        fails.append("6: the skill still says no import carries a failure surface")
    if "SLOPE/W's critical\n  circle" not in text and "SLOPE/W's critical circle" not in text:
        fails.append("6: the skill does not say a solved .gsz brings SLOPE/W's critical "
                     "circle")
    return fails


# --------------------------------------------------------------------------------
# 7. A piezometric line that stops short of the model
# --------------------------------------------------------------------------------

def check_piezo_past_ends():
    """SLOPE/W applies no pore pressure beyond a piezometric line's ends (its Arai and
    Tagyo example saves zero at every column base left of the line's start). The
    import carries that rule as a line continued down to the model's base and along
    it, and says so."""
    from xslope.slice import generate_slices
    fails = []
    with tempfile.TemporaryDirectory() as td:
        sd0 = _built_model()
        sd0["piezo_line"] = [(20.0, 10.0), (60.0, 15.0)]
        base, _ = _export(sd0, td, "short")
        _g, (sd, cav) = _import(base, 1, critical_surface=False)
        pl = sorted(sd.get("piezo_line") or [])
        if not pl or pl[0][0] > 1e-9 or abs(pl[0][1]) > 1e-9:
            fails.append(f"7: a line starting at x = 20 imported as {pl}, not continued "
                         f"along the model's base (y = 0) to x = 0")
        if not any("does not span the model" in c and "left of x = 20" in c
                   for c in cav):
            fails.append("7: the continued line was not stated in the import notes")
        ok, res = _quiet(generate_slices, sd, circle=dict(CIRCLE), num_slices=40,
                         debug=False, check_inputs=False)
        if not ok:
            fails.append(f"7: the circle did not slice: {res}")
        else:
            df = res[0]
            left = df[df["x_c"] < 19.9]
            right = df[(df["x_c"] > 21.0) & (df["y_cb"] < 9.0)]
            if left.empty or (left["u"] != 0.0).any():
                fails.append("7: a slice left of the line's start has pore pressure")
            if right.empty or not (right["u"] > 0.0).all():
                fails.append("7: the slices under the line lost their pore pressure")

    path = _vendor("Arai and Tagyo - Pore water pressure slope.gsz")
    if path is not None:
        from xslope.solve import mprice
        g, (sd, _cav) = _import(path, critical_surface=False)
        got = _critical_fs(g, g["analyses"][0], sd,
                           lambda d: mprice(d, f_type="half_sine"))
        if got is None or abs(got[0] / got[1] - 1.0) > 0.005:
            fails.append(f"7: Arai and Tagyo's pore-water slope solves to {got} on "
                         f"SLOPE/W's critical circle (XSLOPE, SLOPE/W); expected within "
                         f"0.5%")
    return fails


# --------------------------------------------------------------------------------
# 8. Spatial Mohr-Coulomb strengths
# --------------------------------------------------------------------------------

def _set_spatial(root, points, input_param="Y", phi=None):
    """Make the foundation a Spatial Mohr-Coulomb material whose cohesion is the
    curve ``points`` [(coordinate, c), ...] against ``input_param``."""
    for m in root.findall(ST):
        if m.findtext("Name") != "foundation":
            continue
        m.find("SlopeModel").text = "MohrCoulombSpatial"
        ss = m.find("StressStrain")
        for tag in ("CohesionPrime", "PhiPrime"):
            e = ss.find(tag)
            if e is not None:
                ss.remove(e)
        ET.SubElement(ss, "SlopeCohesionModel").text = "Function"
        ET.SubElement(ss, "CohesionFnNum").text = "1"
        if phi is not None:
            ET.SubElement(ss, "PhiPrime").text = str(phi)
    fns = root.find("Functions")
    if fns is None:
        fns = ET.SubElement(root, "Functions")
    mat = ET.SubElement(fns, "Material")
    ss = ET.SubElement(mat, "StressStrain")
    cfs = ET.SubElement(ss, "CohesionFns", Len="1")
    fn = ET.SubElement(cfs, "CohesionFn")
    ET.SubElement(fn, "ID").text = "1"
    ET.SubElement(fn, "Name").text = "strength line"
    ET.SubElement(fn, "Function").text = (
        f"Spline(OutputParam=CohesionPrime,InputParam={input_param},Smooth=100,"
        f"Tension=100,SlopeBefore=0,SlopeAfter=0,LogInput=false,LogOutput=false)")
    pts = ET.SubElement(fn, "Points", Len=str(len(points)))
    for x, y in points:
        ET.SubElement(pts, "Point", X=str(x), Y=str(y))


def check_spatial_strength():
    """A Spatial Mohr-Coulomb cohesion that is a straight line rising with depth, on
    phi = 0, is XSLOPE's cp (SLOPE/W's own SHANSEP example has one). Any other
    function is replaced by its lowest value, and the notes name the material, the
    function and the value used."""
    fails = []
    with tempfile.TemporaryDirectory() as td:
        base, _ = _export(_built_model(), td, "base")
        cases = {"line": ([(0.0, 40.0), (8.0, 20.0)], "Y", None),
                 "friction": ([(0.0, 40.0), (8.0, 20.0)], "Y", 20.0),
                 "x": ([(0.0, 30.0), (30.0, 18.0), (60.0, 25.0)], "X", None)}
        got = {}
        for name, (pts, inp, phi) in cases.items():
            p = _mutate(base, os.path.join(td, f"{name}.gsz"),
                        lambda r, a=pts, b=inp, c=phi: _set_spatial(r, a, b, c))
            _g, got[name] = _import(p, 1, critical_surface=False)
        m = _mat(got["line"][0], "foundation")
        if (m["option"], m["c"], m["cp"], m["r_elev"], m["phi"]) != \
                ("cp", 20.0, 2.5, 8.0, 0.0):
            fails.append(f"8: a cohesion line 40 at y = 0 to 20 at y = 8 imported as "
                         f"{m['option']} c={m['c']} cp={m['cp']} r_elev={m['r_elev']}, "
                         f"not cp 20 / 2.5 below y = 8")
        if not any("cp option" in c for c in got["line"][1]):
            fails.append("8: the cp mapping was not stated in the import notes")
        m = _mat(got["friction"][0], "foundation")
        if (m["option"], m["c"], m["phi"]) != ("mc", 20.0, 20.0) or not any(
                "'foundation'" in c and "lowest" in c for c in got["friction"][1]):
            fails.append(f"8: a cohesion line on phi = 20 imported as {m['option']} "
                         f"c={m['c']} phi={m['phi']}, not mc with the lowest c, 20, "
                         f"and a note (XSLOPE's cp takes no friction angle)")
        m = _mat(got["x"][0], "foundation")
        if (m["option"], m["c"]) != ("mc", 18.0) or not any(
                "'foundation'" in c and "'strength line'" in c and "c = 18" in c
                for c in got["x"][1]):
            fails.append(f"8: a cohesion function of x imported as {m['option']} "
                         f"c={m['c']}, or without a note naming the material, the "
                         f"function and the value used")

        path = _slopew_example(
            td, "SLOPE2D/SHANSEP - Stress History and Normalized Soil Engineering "
                "Property.gsz")
        if path is not None:
            from xslope.solve import spencer
            from xslope.geostudio import read_gsz
            g = read_gsz(path)
            a = next(a for a in g["analyses"] if a["name"].startswith("4 "))
            _g, (sd, _cav) = _import(path, a["id"], critical_surface=False)
            got = _critical_fs(g, a, sd, spencer)
            if got is None or abs(got[0] / got[1] - 1.0) > 0.005:
                fails.append(f"8: SLOPE/W's Spatial Mohr-Coulomb example solves to "
                             f"{got} on its critical circle (XSLOPE, SLOPE/W); expected "
                             f"within 0.5%")
    return fails


# --------------------------------------------------------------------------------
# 9. Reinforcement: the untyped default, the bond length, F of S dependence
# --------------------------------------------------------------------------------

def check_anchor_import():
    """A product with no <Type> is GeoStudio's default, an anchor (axial). An anchor
    develops pullout over its <BondLength> only, so its capacity is the smaller of
    the tensile capacity and pullout x bond length. <FofSDependent>true is passive.
    A shear force, and a type XSLOPE does not map, are named in the notes."""
    fails = []
    with tempfile.TemporaryDirectory() as td:
        base, _ = _export(_built_model(), td, "base")
        products = {
            "untyped": {"Name": "bar", "ID": "1", "Spacing": "2", "BondLength": "5",
                        "Tensile": "1000", "PlateCapacity": "1000",
                        "PulloutResistance": "100"},
            "fofs": {"Name": "bar", "ID": "1", "Spacing": "2", "BondLength": "5",
                     "Tensile": "1000", "PlateCapacity": "1000",
                     "PulloutResistance": "100", "FofSDependent": "true", "Shear": "50"},
            "nail": {"Name": "bar", "ID": "1", "Type": "Nail", "Spacing": "2",
                     "BondLength": "5", "Tensile": "1000", "PlateCapacity": "1000",
                     "PulloutResistance": "100"},
            "pile": {"Name": "pile", "ID": "1", "Type": "Pile", "Shear": "50"},
        }
        got = {}
        for name, prod in products.items():
            p = _mutate(base, os.path.join(td, f"{name}.gsz"),
                        lambda r, q=prod: _add_reinforcement(r, q))
            _g, got[name] = _import(p, 1, critical_surface=False)
        rl = got["untyped"][0]["reinforcement_lines"][0]
        # rate 100 / 2 = 50 per unit length over 5 -> 250, below 1000 / 2 = 500
        if (rl["type"], rl["dir"], rl["appl"]) != ("anchor", "axial", "active") or \
                not _close(rl["t_max"], 250.0) or not _close(max(rl["lp1"], rl["lp2"]),
                                                             5.0):
            fails.append(f"9: an untyped product with a 5 m bond imported as "
                         f"{rl['type']}/{rl['dir']}/{rl['appl']} t_max={rl['t_max']} "
                         f"lp={rl['lp1']}, not an axial anchor at 250 over 5")
        cav = got["untyped"][1]
        if not any("no type" in c and "anchor" in c for c in cav) or \
                not any("bond length" in c for c in cav):
            fails.append("9: the anchor default and the bond length were not stated")
        rl = got["fofs"][0]["reinforcement_lines"][0]
        if rl["appl"] != "passive":
            fails.append("9: an F of S dependent product did not import as passive")
        if not any("shear force" in c for c in got["fofs"][1]):
            fails.append("9: a dropped shear force was not stated")
        rl = got["nail"][0]["reinforcement_lines"][0]
        if not _close(rl["t_max"], 500.0):
            fails.append(f"9: a nail's BondLength capped it at {rl['t_max']}; nails are "
                         f"bonded over their whole length")
        if not any("'pile' (Pile)" in c for c in got["pile"][1]):
            fails.append("9: a Pile product was not named as unmapped")

        from xslope.solve import janbu, spencer
        path = _vendor("Pockoski and Duncan - Reinforcement.gsz")
        if path is not None:
            from xslope.geostudio import read_gsz
            g = read_gsz(path)
            for a in g["analyses"]:
                solver = {"Janbu": janbu, "Spencer": spencer}.get(a["method"])
                if solver is None or "single" in a["name"]:
                    continue
                _g, (sd, _cav) = _import(path, a["id"], critical_surface=False)
                got_fs = _critical_fs(g, a, sd, solver)
                if got_fs is None or abs(got_fs[0] / got_fs[1] - 1.0) > 0.03:
                    fails.append(f"9: Pockoski and Duncan's tieback, {a['name']}: "
                                 f"{got_fs} on SLOPE/W's critical circle (XSLOPE, "
                                 f"SLOPE/W); expected within 3%")
        path = _slopew_example(td, "SLOPE2D/Reinforcement with Anchors.gsz")
        if path is not None:
            from xslope.geostudio import read_gsz
            g = read_gsz(path)
            for a in g["analyses"]:
                if a["name"].startswith("4 "):
                    continue                # its shear force is not imported
                _g, (sd, _cav) = _import(path, a["id"], critical_surface=False)
                got_fs = _critical_fs(g, a, sd, spencer)
                if got_fs is None or abs(got_fs[0] / got_fs[1] - 1.0) > 0.005:
                    fails.append(f"9: Reinforcement with Anchors, {a['name']}: "
                                 f"{got_fs} on SLOPE/W's critical circle (XSLOPE, "
                                 f"SLOPE/W); expected within 0.5%")
    return fails


# --------------------------------------------------------------------------------
# 10. Saved slip circles of a search with a moment axis
# --------------------------------------------------------------------------------

def check_axis_rows():
    """A search with a moment axis writes the axis as every saved row's center. The
    circle whose columns are saved gets its center back from them; the other rows
    are 'axis' rows and are not rebuilt as circles."""
    import math
    from xslope.geostudio import _axis_circle_from_columns
    fails = []
    xc, yc, r = 40.0, 25.0, 20.0
    angles = [math.radians(a) for a in range(230, 311, 10)]
    pts = [(xc + r * math.cos(t), yc + r * math.sin(t)) for t in angles]
    lines = ["XLeft,XRight,YBackBotLeft,YBackBotRight"]
    lines += [f"{a[0]},{b[0]},{a[1]},{b[1]}" for a, b in zip(pts, pts[1:])]
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr("a/001/column_7.csv", "\n".join(lines))
        zf.writestr("a/001/column_8.csv", "\n".join(lines[:3] + ["1,2,3,9", "2,3,9,3"]))
    fit = _axis_circle_from_columns(zipfile.ZipFile(buf), "a/001/",
                                    [{"SlipNum": "7"}, {"SlipNum": "8"}])
    if "7" not in fit or max(abs(fit["7"][0] - xc), abs(fit["7"][1] - yc),
                             abs(fit["7"][2] - r)) > 1e-9:
        fails.append(f"10: a saved circle's columns gave center {fit.get('7')}, not "
                     f"({xc}, {yc}, {r})")
    if "8" in fit:
        fails.append("10: columns that are not one circle were fitted as one")

    with tempfile.TemporaryDirectory() as td:
        path = _slopew_example(td, "SLOPE2D/Reinforcement with Anchors.gsz")
        if path is not None:
            from xslope.geostudio import read_gsz, read_gsz_results
            g = read_gsz(path)
            a = next(a for a in g["analyses"] if a["name"].startswith("2 "))
            rows = read_gsz_results(g, a["name"])
            circles = [t for t in rows if t["circular"]]
            if len(circles) != 1 or {t["kind"] for t in rows} != {"axis", "circle"}:
                fails.append(f"10: the anchors example's axis search gave "
                             f"{len(circles)} circles and kinds "
                             f"{sorted({t['kind'] for t in rows})}")
            elif abs(circles[0]["xo"] - 187.644) > 0.01 or \
                    abs(circles[0]["yo"] - 119.339) > 0.01:
                fails.append(f"10: its critical circle's center came back as "
                             f"({circles[0]['xo']:.3f}, {circles[0]['yo']:.3f}), not "
                             f"(187.644, 119.339)")
    return fails


CHECKS = [
    ("1. gamma_sat written on export", check_gamma_sat_export),
    ("2. unit weight above the water table", check_unit_weight_import),
    ("3. the older tie-back capacity tag", check_tieback_capacity),
    ("4. the phreatic switch written on export", check_phreatic_export),
    ("5. partial factors named", check_partial_factors),
    ("6. the skill on imported surfaces", check_skill),
    ("7. a piezometric line short of the model", check_piezo_past_ends),
    ("8. spatial Mohr-Coulomb strengths", check_spatial_strength),
    ("9. anchors: default type, bond, F of S", check_anchor_import),
    ("10. saved circles of an axis search", check_axis_rows),
]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    if not os.path.isdir(GSZ_VENDOR_CORPUS) and not os.path.isdir(GSZ_DOWNLOADS):
        print("gsz gaps: SLOPE/W's own files are not here; only the built-file checks run.")
    failures = []
    for name, fn in CHECKS:
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
    print("gsz gaps:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll gsz gap checks passed.")


if __name__ == "__main__":
    main()
