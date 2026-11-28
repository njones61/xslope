"""Standing checks on SLOPE/W's staged rapid drawdown, in and out of a .gsz.

SLOPE/W solves a rapid drawdown by the three stages of Duncan, Wright and Wong
(1990), the same method as XSLOPE's own, and it states the problem the same way:
the water before drawdown and the water after it as two piezometric surfaces, and
the strength of each soil that does not drain as an R envelope, <DrawdownCohesion>
and <DrawdownPhi> (c_R, phi_R). XSLOPE takes the two surfaces as piezo lines 1 and
2, and the envelope as d and psi after SLOPE/W's own conversion (its Eqs. 77-78).

  A. THE IMPORT, on SLOPE/W's four rapid drawdown examples. Each part of the model
     is checked against values this file reads out of the vendor's XML itself, with
     its own parser: the higher surface is line 1 and the lower is line 2, point for
     point; each undrained soil's d and psi are its c_R and phi_R converted by
     SLOPE/W's equations, written out here again rather than borrowed from the
     importer; a free-draining soil has neither; the caveat that used to say the
     drawdown was not imported is gone; and a surface that stops short of the model
     comes across as drawn, with a note on SLOPE/W's no-water rule past its end.
  B. THE ANSWER. Each example is solved as a rapid drawdown on SLOPE/W's own
     critical circle, imported from the same file, with SLOPE/W's own method, and
     the factor of safety must be the one docs/verification/geostudio.md prints.
  C. THE ROUND TRIP. Each example is imported, exported and imported again, and the
     model must come back identical, carrying the drawdown the XML states. An
     XSLOPE model built here, with no vendor file, goes round the same way.
  D. THE EXPORTED XML. The rapid drawdown is written the way SLOPE/W's own examples
     write it: the <IsRapidDrawdown> switch under the analysis's <Control>, ahead of
     the piezometric surfaces; two surfaces, the higher one first and the one the
     materials use; and <DrawdownCohesion>, <DrawdownPhi> right after <PhiPrime>, at
     the vendor's own values. Every tag path the export writes for it is one the
     vendor file has.

SLOPE/W's examples are Seequent's copyrighted files and are not in this
repository. The checks that read them (A, B, the vendor half of C, and D) look for
a local copy in ``GSZ_VENDOR_CORPUS`` and skip cleanly when it is absent.
"""
import contextlib
import io
import math
import os
import sys
import tempfile
import xml.etree.ElementTree as ET
import zipfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

#: Local folder of GeoStudio .gsz models (git-ignored, not redistributable); the
#: same folder run_tests.py reads the corpus from.
GSZ_VENDOR_CORPUS = os.path.expanduser("~/python_projects/vendor_files/gsz_corpus")

#: SLOPE/W's four staged rapid drawdown examples (its verification manual, sections
#: 2.40 to 2.43): file, the method SLOPE/W solved it with, and XSLOPE's factor of
#: safety on SLOPE/W's own critical circle as docs/verification/geostudio.md prints it.
EXAMPLES = [
    ("2.40", "Rapid drawdown with multi-stage.gsz", "spencer", 1.023),
    ("2.41", "Staged rapid drawdown - CE Example.gsz", "spencer", 1.434),
    ("2.42", "Staged rapid drawdown - Pumped Storage Project Dam.gsz", "mprice", 1.538),
    ("2.43", "Staged rapid drawdown - Pilarcitos Dam.gsz", "mprice", 1.050),
]

#: The slices the page's numbers are solved with, and how far a solve may drift
#: from the printed value (half the last printed digit, plus a little).
NUM_SLICES = 60
FS_TOL = 0.0006

#: What the import note says of a surface that stops short of the model: SLOPE/W
#: applies no pore pressure past its end, and what to do there, by where the end sits.
#: The other two examples' surfaces span their models and draw no such note.
SHORT_NOTES = {
    "2.40": ["piezometric surface before drawdown (line 1) stops short of the model",
             "SLOPE/W applies no pore pressure outside a line's x-range; XSLOPE "
             "refuses a slice there instead",
             "results were computed with no water right of x = 211",
             "right end the line is on the ground surface",
             "extend the line along the ground surface to x = 216"],
    "2.43": ["piezometric surface before drawdown (line 1) stops short of the model",
             "SLOPE/W applies no pore pressure outside a line's x-range; XSLOPE "
             "refuses a slice there instead",
             "results were computed with no water left of x = 1",
             "above the ground: a pool", "Extend the line level to x = 0"],
}

#: The caveat the import gave before it carried a staged rapid drawdown.
OLD_CAVEAT = "THE FACTOR OF SAFETY IS NOT SLOPE/W'S RAPID DRAWDOWN ANSWER"


def _close(a, b, tol=1e-9):
    return abs(float(a) - float(b)) <= tol * max(1.0, abs(float(b)))


def _vendor(name):
    path = os.path.join(GSZ_VENDOR_CORPUS, name)
    return path if os.path.isfile(path) else None


def _root(path):
    """The model XML at the root of a .gsz, parsed."""
    zf = zipfile.ZipFile(path)
    name = next(n for n in zf.namelist() if "/" not in n and n.lower().endswith(".xml"))
    return ET.fromstring(zf.read(name))


def _xml_drawdown(root):
    """What the XML states about the drawdown, read here without the importer.

    Returns ``{'switch', 'before', 'after', 'before_id', 'materials', 'uses'}``:
    the <IsRapidDrawdown> text, the two surfaces as point lists (the higher on
    average is the water before drawdown), the ID of that one, ``{name: (phi',
    c_R, phi_R)}`` for every material (c_R = phi_R = 0 for a free-draining one),
    and ``{material name: surface ID}``.
    """
    entry = root.find("./StabilityItems/StabilityItem/Entry")
    pts = {d.get("Number"): (float(d.get("X")), float(d.get("Y")))
           for d in entry.findall("./DataPoints/DataPoint") if d.get("X") is not None}
    surfaces = {}
    for s in entry.findall("./PiezometricSurfaces/PiezometricSurface"):
        surfaces[s.findtext("ID")] = [pts[d.text] for d in
                                      s.findall("./DataPoints/DataPoint")]
    mean = {k: sum(y for _, y in v) / len(v) for k, v in surfaces.items()}
    hi, lo = max(mean, key=mean.get), min(mean, key=mean.get)
    names = {m.findtext("ID"): m.findtext("Name")
             for m in root.findall("./Materials/Material")}
    used = {g.get("Entry") for g in root.findall(
        "./Contexts/Context/GeometryUsesMaterials/GeometryUsesMaterial")}
    mats = {}
    for m in root.findall("./Materials/Material"):
        if m.findtext("ID") not in used or m.findtext("SlopeModel") == "Bedrock":
            continue
        ss = m.find("StressStrain")
        mats[m.findtext("Name")] = tuple(float(ss.findtext(t) or 0.0) for t in
                                         ("PhiPrime", "DrawdownCohesion", "DrawdownPhi"))
    uses = {names[u.get("ID")]: u.get("UsesID")
            for u in entry.findall("./MaterialUsesPiezs/MaterialUsesPiez")}
    return {"switch": entry.findtext("./Control/IsRapidDrawdown"),
            "before": surfaces[hi], "after": surfaces[lo], "before_id": hi,
            "materials": mats, "uses": uses}


def _kc1(c_r, phi_r, phi):
    """SLOPE/W's Eqs. 77-78 (Slope Stability Modeling, sec. 9.3): the K_c = 1
    envelope from the R envelope drawn tangent to the Mohr circles."""
    sr = math.sin(math.radians(phi_r))
    cos_p = math.cos(math.radians(phi))
    return (c_r * math.cos(math.radians(phi_r)) * cos_p / (1.0 - sr),
            math.degrees(math.atan(sr * cos_p / (1.0 - sr))))


def _same_line(a, b):
    return len(a) == len(b) and all(_close(p, q) for u, v in zip(a, b)
                                    for p, q in zip(u, v))


def _import(path):
    from xslope.geostudio import read_gsz, gsz_to_slope_data
    with contextlib.redirect_stdout(io.StringIO()):
        return gsz_to_slope_data(read_gsz(path))


def _drawdown_problems(label, sd, want):
    """How ``sd`` departs from the drawdown the XML states (an empty list if not)."""
    # Both surfaces come across exactly as drawn, including one that stops short of
    # the model (Walter Bouldin's and Pilarcitos's surfaces before drawdown).
    out = []
    if not _same_line(sd.get("piezo_line") or [], want["before"]):
        out.append(f"{label}: piezo line 1 is {sd.get('piezo_line')}, expected the "
                   f"surface before drawdown {want['before']}")
    if not _same_line(sd.get("piezo_line2") or [], want["after"]):
        out.append(f"{label}: piezo line 2 is {sd.get('piezo_line2')}, expected the "
                   f"surface after drawdown {want['after']}")
    mats = {m["name"]: m for m in sd["materials"]}
    for name, (phi, c_r, phi_r) in want["materials"].items():
        m = mats.get(name)
        if m is None:
            out.append(f"{label}: material {name!r} is missing")
            continue
        d, psi = _kc1(c_r, phi_r, phi) if (c_r or phi_r) else (0.0, 0.0)
        if not (_close(m.get("d") or 0.0, d) and _close(m.get("psi") or 0.0, psi)):
            out.append(f"{label}: {name!r} has d = {m.get('d')}, psi = {m.get('psi')}; "
                       f"its c_R = {c_r:g}, phi_R = {phi_r:g} give d = {d:.6g}, "
                       f"psi = {psi:.6g}")
        if want["uses"] and name in want["uses"] and m.get("u") != "piezo":
            out.append(f"{label}: {name!r} takes u = {m.get('u')!r}, but SLOPE/W "
                       f"connects it to a piezometric surface")
    return out


# --------------------------------------------------------------------------------
# A. The import
# --------------------------------------------------------------------------------

def check_import():
    fails = []
    for sec, name, _meth, _fs in EXAMPLES:
        path = _vendor(name)
        if path is None:
            continue
        want = _xml_drawdown(_root(path))
        label = f"§{sec}"
        if want["switch"] != "true":
            fails.append(f"{label}: the vendor file is not a staged rapid drawdown")
            continue
        sd, cav = _import(path)
        fails += _drawdown_problems(label, sd, want)
        if any(OLD_CAVEAT in c for c in cav):
            fails.append(f"{label}: the import still says the drawdown was not imported")
        if not any("imported as XSLOPE's rapid drawdown" in c for c in cav):
            fails.append(f"{label}: no caveat says the drawdown was imported, or that "
                         f"the model is to be run as one")
        if not any("rapid drawdown FS" in c for c in cav):
            fails.append(f"{label}: the circle caveat does not say SLOPE/W's FS on it "
                         f"is its rapid drawdown answer")
        notes = [c for c in cav if "stops short of the model" in c]
        want_note = SHORT_NOTES.get(sec)
        if want_note is None and notes:
            fails.append(f"{label}: a short-line note on surfaces that span the model: "
                         f"{notes}")
        elif want_note is not None:
            missing = [w for w in want_note if not notes or w not in notes[0]]
            if len(notes) != 1 or missing:
                fails.append(f"{label}: the short-line note lacks {missing}: {notes}")
    return fails


# --------------------------------------------------------------------------------
# B. The answer, on SLOPE/W's own circle
# --------------------------------------------------------------------------------

def check_scores():
    import inspect
    from xslope.slice import generate_slices
    from xslope.advanced import rapid_drawdown
    # Slices made for a rapid drawdown also break where a surface crosses line 2.
    # (On these four circles the factor of safety is the same either way.)
    rapid_kw = ({"rapid": True} if "rapid" in inspect.signature(generate_slices).parameters
                else {})
    fails = []
    for sec, name, meth, fs_page in EXAMPLES:
        path = _vendor(name)
        if path is None:
            continue
        sd, _ = _import(path)
        if not sd.get("circles"):
            fails.append(f"§{sec}: no SLOPE/W circle was imported")
            continue
        with contextlib.redirect_stdout(io.StringIO()):
            ok, res = generate_slices(sd, circle=sd["circles"][0],
                                      num_slices=NUM_SLICES, **rapid_kw)
            if not ok:
                fails.append(f"§{sec}: SLOPE/W's circle did not slice: {res}")
                continue
            ok, r = rapid_drawdown(res[0].copy(), meth, debug_level=0)
        if not ok:
            fails.append(f"§{sec}: the rapid drawdown did not solve: {r}")
        elif abs(r["FS"] - fs_page) > FS_TOL:
            fails.append(f"§{sec}: {meth} rapid drawdown on SLOPE/W's circle gives "
                         f"{r['FS']:.4f}; the verification page prints {fs_page:.3f}")
    return fails


# --------------------------------------------------------------------------------
# C. The round trip
# --------------------------------------------------------------------------------

_MAT_KEYS = ("name", "option", "gamma", "c", "phi", "d", "psi", "u")
_MODEL_KEYS = ("gamma_water", "water_loads", "piezo_phreatic", "piezo_phreatic2",
               "k_seismic", "unit_system")


def _model_problems(label, a, b):
    """How model ``b`` differs from model ``a``: geometry, materials and water."""
    out = []
    for k in ("piezo_line", "piezo_line2"):
        if not _same_line(b.get(k) or [], a.get(k) or []):
            out.append(f"{label}: {k} {a.get(k)} came back as {b.get(k)}")
    if len(a["materials"]) != len(b["materials"]):
        out.append(f"{label}: {len(a['materials'])} materials came back as "
                   f"{len(b['materials'])}")
    for m, n in zip(a["materials"], b["materials"]):
        for k in _MAT_KEYS:
            x, y = m.get(k), n.get(k)
            same = (_close(x or 0.0, y or 0.0) if isinstance(x, (int, float))
                    and not isinstance(x, bool) else x == y)
            if not same:
                out.append(f"{label}: {m.get('name')!r} {k} {x!r} came back as {y!r}")
    if len(a["polygons"]) != len(b["polygons"]) or any(
            p["mat_id"] != q["mat_id"] or not p["polygon"].equals_exact(q["polygon"], 1e-7)
            for p, q in zip(a["polygons"], b["polygons"])):
        out.append(f"{label}: the material zones did not come back the same")
    for k in _MODEL_KEYS:
        if a.get(k) != b.get(k):
            out.append(f"{label}: {k} {a.get(k)!r} came back as {b.get(k)!r}")
    return out


def _export_back(sd, td, stem):
    from xslope.geostudio import read_gsz, gsz_to_slope_data, export_gsz
    out = os.path.join(td, f"{stem}.gsz")
    with contextlib.redirect_stdout(io.StringIO()):
        cav = export_gsz(sd, out, analysis_name=stem)
        back, _ = gsz_to_slope_data(read_gsz(out), 1)
    return out, cav, back


def _built_model():
    """An XSLOPE rapid drawdown model with no vendor file behind it: a clay core
    with a K_c = 1 envelope inside a free-draining shell, pool 20 drawn down to 6."""
    from shapely.geometry import Polygon
    from xslope.geostudio import _blank_material
    from xslope.fileio import build_ground_surface_from_polygons
    shell = dict(_blank_material("shell"), option="mc", gamma=20.0, c=0.0, phi=38.0,
                 u="piezo")
    core = dict(_blank_material("core"), option="mc", gamma=19.0, c=12.0, phi=28.0,
                d=40.0, psi=17.5, u="piezo")
    polygons = [
        {"mat_id": 0, "polygon": Polygon([(0, 0), (110, 0), (110, 5), (100, 5),
                                          (60, 5), (50, 5), (10, 5), (0, 5)])},
        {"mat_id": 0, "polygon": Polygon([(10, 5), (50, 5), (50, 25)])},
        {"mat_id": 1, "polygon": Polygon([(50, 5), (60, 5), (60, 25), (50, 25)])},
        {"mat_id": 0, "polygon": Polygon([(60, 5), (100, 5), (60, 25)])},
    ]
    ground, domain = build_ground_surface_from_polygons(polygons)
    return {"polygons": polygons, "materials": [shell, core], "gamma_water": 9.807,
            "unit_system": "si", "water_loads": "auto", "ground_surface": ground,
            "domain_polygon": domain,
            "piezo_line": [(0.0, 20.0), (55.0, 20.0), (100.0, 5.0), (110.0, 5.0)],
            "piezo_line2": [(0.0, 8.0), (40.0, 8.0), (55.0, 16.0), (100.0, 5.0),
                            (110.0, 5.0)],
            "piezo_phreatic": False, "piezo_phreatic2": False, "k_seismic": 0.0,
            "dloads": [], "dloads2": []}


def check_round_trip():
    fails = []
    with tempfile.TemporaryDirectory() as td:
        # An XSLOPE model built here: no vendor file needed.
        sd = _built_model()
        _out, _cav, back = _export_back(sd, td, "built")
        fails += _model_problems("built model", sd, back)
        for sec, name, _meth, _fs in EXAMPLES:
            path = _vendor(name)
            if path is None:
                continue
            first, _ = _import(path)
            _out, _cav, back = _export_back(first, td, f"rt{sec}")
            label = f"§{sec} round trip"
            fails += _model_problems(label, first, back)
            # Identical is not enough: a model that lost the drawdown on the way in
            # comes back identically without it. It must also carry the XML's.
            fails += _drawdown_problems(label, back, _xml_drawdown(_root(path)))
    return fails


# --------------------------------------------------------------------------------
# D. The exported XML
# --------------------------------------------------------------------------------

def _paths(elem, prefix=""):
    """Every tag path under ``elem`` (and tag@attribute), as a set."""
    here = f"{prefix}/{elem.tag}"
    out = {here} | {f"{here}@{a}" for a in elem.attrib}
    for c in elem:
        out |= _paths(c, here)
    return out


def _xml_shape_problems(label, ours, vendor):
    """How the rapid drawdown in ``ours`` departs from ``vendor``'s way of writing it."""
    out = []
    v_entry = vendor.find("./StabilityItems/StabilityItem/Entry")
    o_entry = ours.find("./StabilityItems/StabilityItem/Entry")
    # The switch, under <Control>, ahead of the piezometric surfaces.
    if o_entry.findtext("./Control/IsRapidDrawdown") != "true":
        out.append(f"{label}: no <Control><IsRapidDrawdown>true</...> on the analysis")
    v_tags, o_tags = [c.tag for c in v_entry], [c.tag for c in o_entry]
    for t in ("Control", "PiezometricSurfaces", "MaterialUsesPiezs"):
        if t in v_tags and t not in o_tags:
            out.append(f"{label}: the analysis has no <{t}>, which the vendor file has")
    common = [t for t in v_tags if t in o_tags]
    if [t for t in o_tags if t in common] != common:
        out.append(f"{label}: the analysis's {common} are in the order "
                   f"{[t for t in o_tags if t in common]}, not the vendor file's")
    # Every tag path the export writes under these elements is one the vendor writes.
    for sub in ("Control", "PiezometricSurfaces", "MaterialUsesPiezs"):
        o, v = o_entry.find(sub), v_entry.find(sub)
        if o is not None and v is not None:
            extra = sorted(_paths(o) - _paths(v))
            if extra:
                out.append(f"{label}: <{sub}> carries {extra[:3]}, which the vendor "
                           f"file does not write")
    # Two surfaces; the materials use the one before drawdown, as the vendor's do.
    want = _xml_drawdown(vendor)
    got = _xml_drawdown(ours)
    surfaces = o_entry.findall("./PiezometricSurfaces/PiezometricSurface")
    if len(surfaces) != 2 or (o_entry.find("PiezometricSurfaces") is not None and
                              o_entry.find("PiezometricSurfaces").get("Len") != "2"):
        out.append(f"{label}: {len(surfaces)} piezometric surface(s) written, not 2")
    else:
        if not _same_line(got["before"], want["before"]) or \
                not _same_line(got["after"], want["after"]):
            out.append(f"{label}: the surfaces written are not the vendor's before and "
                       f"after drawdown")
        v_uses = {n for n, u in want["uses"].items() if u == want["before_id"]}
        o_uses = {n for n, u in got["uses"].items() if u == got["before_id"]}
        if not v_uses <= o_uses:
            out.append(f"{label}: {sorted(v_uses - o_uses)} do not use the surface before "
                       f"drawdown, as they do in the vendor file")
    # <DrawdownCohesion>, <DrawdownPhi>: right after <PhiPrime>, at the vendor's values.
    v_mats = {m.findtext("Name"): m.find("StressStrain")
              for m in vendor.findall("./Materials/Material")}
    for m in ours.findall("./Materials/Material"):
        name, o_ss = m.findtext("Name"), m.find("StressStrain")
        v_ss = v_mats.get(name)
        if v_ss is None:
            continue
        for t in ("DrawdownCohesion", "DrawdownPhi"):
            vt, ot = v_ss.find(t), o_ss.find(t)
            if (vt is None) != (ot is None):
                out.append(f"{label}: {name!r} {'lacks' if ot is None else 'has'} "
                           f"<{t}>, which the vendor file {'has' if ot is None else 'does not'}")
            elif vt is not None and not _close(ot.text, vt.text):
                out.append(f"{label}: {name!r} <{t}> written as {ot.text}, the vendor's "
                           f"is {vt.text}")
        seq = ("PhiPrime", "DrawdownCohesion", "DrawdownPhi")
        if v_ss.find("DrawdownCohesion") is not None:
            vt = [c.tag for c in v_ss if c.tag in seq]
            ot = [c.tag for c in o_ss if c.tag in seq]
            o_all = [c.tag for c in o_ss]
            adjacent = all(t in o_all for t in seq) and \
                o_all.index(seq[2]) - o_all.index(seq[0]) == 2
            if ot != vt or not adjacent:
                out.append(f"{label}: {name!r} writes {ot} as {o_all}, not one after "
                           f"another in the vendor's order {vt}")
    return out


def check_exported_xml():
    fails = []
    with tempfile.TemporaryDirectory() as td:
        for sec, name, _meth, _fs in EXAMPLES:
            path = _vendor(name)
            if path is None:
                continue
            first, _ = _import(path)
            out, _cav, _back = _export_back(first, td, f"x{sec}")
            fails += _xml_shape_problems(f"§{sec} export", _root(out), _root(path))
    return fails


CHECKS = [
    ("A. import against the vendor XML", check_import),
    ("B. answer on SLOPE/W's own circle", check_scores),
    ("C. round trip", check_round_trip),
    ("D. exported XML against the vendor's", check_exported_xml),
]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    have = [n for _s, n, _m, _f in EXAMPLES if _vendor(n)]
    if not have:
        print(f"gsz rapid drawdown: SLOPE/W's examples are not in {GSZ_VENDOR_CORPUS}; "
              f"only the built model's round trip runs.")
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
    print("gsz rapid drawdown:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll gsz rapid drawdown checks passed.")


if __name__ == "__main__":
    main()
