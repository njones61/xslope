"""Build Tutorial W-2's DXF drawing and the workbook its import produces.

W-2 walks three importers. The first of them, DXF, needs a drawing — and the
drawing has to be one the reader can be handed, which rules out every vendor file
in the tutorial's other two parts. So this script draws one: a three-layer highway
embankment in feet, with a water table through the foundation and a surcharge on
the crest, written through :func:`xslope.cad.export_dxf` so the file carries
XSLOPE's own layer conventions (a layer per material, plus ``PIEZO``, ``DLOADS``
and ``SEARCH_CIRCLES``) and the import wizard's suggested targets are right on
every row.

Two files are produced:

  - ``docs/tutorials/files/w02_section.dxf`` — the drawing the tutorial opens with.
  - ``docs/tutorials/files/w02_section_imported.xlsx`` — what the reader has after
    importing it and typing in the properties a drawing cannot carry.

The second is built by running the *actual* import — ``read_dxf_layers`` +
``suggest_dxf_target`` + ``ProjectDocument.build_from_dxf_mapping``, which is the
path Studio's File → Import DXF… takes — and then filling the placeholders the
tutorial has the reader fill: the three material rows, the surcharge pressure, and
the tidier material names typed into the wizard's Material column. Nothing about
the geometry is retyped, so the shipped workbook is the drawing's own round trip
and the page's factor of safety is the one the reader's own import produces.

The source model is built here and not shipped: the tutorial's premise is that the
drawing is all the reader has.

Run:  PYTHONPATH=. python3 tools/build_w02_import_dxf.py
"""

from __future__ import annotations

import contextlib
import io
import os
import sys
import tempfile

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from shapely.geometry import Polygon                                  # noqa: E402

from xslope.cad import export_dxf, read_dxf_layers, suggest_dxf_target  # noqa: E402
from xslope.fileio import load_slope_data, save_slope_data_to_xlsx      # noqa: E402
from xslope.search import circular_search                               # noqa: E402

#: An existing polygon-geometry model in feet, read only for the global settings a
#: section does not state — unit system, unit weight of water, template layout. Its
#: geometry, materials, water and loads are all replaced below.
BASE = os.path.join(REPO_ROOT, "docs/lem/files/xslope_sloping_bottom.xlsx")

OUT_DIR = os.path.join(REPO_ROOT, "docs/tutorials/files")
DXF_OUT = os.path.join(OUT_DIR, "w02_section.dxf")
XLSX_OUT = os.path.join(OUT_DIR, "w02_section_imported.xlsx")

# --------------------------------------------------------------------------- #
# The section: a 30 ft highway embankment on 2:1 side slopes, built over 15 ft of
# silty clay on dense sand, with the ground water table 4 ft down at the left edge
# and mounded 2 ft above the base of the fill under the crest, and a 600 psf
# surcharge across a 30 ft strip of the crest.
#
# Names carry through to the DXF: export uppercases them and replaces spaces with
# underscores, so these three become layers EMBANKMENT, SILTY_CLAY and DENSE_SAND.
# --------------------------------------------------------------------------- #
MATERIALS = [
    {"name": "embankment", "gamma": 125.0, "option": "mc", "c": 250.0, "phi": 28.0,
     "u": "piezo"},
    {"name": "silty clay", "gamma": 118.0, "option": "mc", "c": 700.0, "phi": 0.0,
     "u": "piezo"},
    {"name": "dense sand", "gamma": 132.0, "option": "mc", "c": 0.0, "phi": 36.0,
     "u": "piezo"},
]

ZONES = [
    # embankment: toe at the origin, 2:1 face to a crest at elevation 30
    [(0, 0), (60, 30), (140, 30), (140, 0)],
    # silty clay: 15 ft, the full width of the section
    [(-60, 0), (140, 0), (140, -15), (-60, -15)],
    # dense sand: 20 ft below it
    [(-60, -15), (140, -15), (140, -35), (-60, -35)],
]

PIEZO = [(-60.0, -4.0), (60.0, 0.0), (140.0, 2.0)]

SURCHARGE = [{"X": 85.0, "Y": 30.0, "Normal": 600.0},
             {"X": 115.0, "Y": 30.0, "Normal": 600.0}]
SURCHARGE_DIR = "normal"

#: One circle at the base of the fill, one tangent to the top of the dense sand —
#: the two mechanisms the section allows.
CIRCLES = [{"Xo": 35.0, "Yo": 60.0, "Depth": 0.0, "R": 60.0},
           {"Xo": 45.0, "Yo": 65.0, "Depth": -15.0, "R": 80.0}]

#: What the reader types into the wizard's Material column, layer by layer. The
#: column defaults to the layer name, and the layer name is the export's uppercased
#: form; these are the names the tutorial has the reader put back.
MATERIAL_NAMES = {"EMBANKMENT": "embankment", "SILTY_CLAY": "silty clay",
                  "DENSE_SAND": "dense sand"}


def _quiet(fn, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args, **kwargs)


def build_source_model(tmpdir):
    """The model the drawing is made from, as a loaded ``slope_data``.

    Written to a workbook and read back rather than assembled in memory: the loader
    is what derives the ground surface and the domain polygon, and ``export_dxf``
    clips each search circle against that ground surface to write its arc.
    """
    sd = _quiet(load_slope_data, BASE)

    blank = dict(sd["materials"][0])
    materials = []
    for spec in MATERIALS:
        mat = dict(blank)
        mat.update({"gamma_sat": None, "cp": 0.0, "r_elev": 0.0, "d": 0.0, "psi": 0.0,
                    "ru": 0.0})
        mat.update(spec)
        materials.append(mat)
    sd["materials"] = materials

    sd["polygons"] = [{"polygon": Polygon(coords), "mat_id": i, "size": None}
                      for i, coords in enumerate(ZONES)]
    sd["profile_lines"] = []
    sd["piezo_line"] = [tuple(p) for p in PIEZO]
    sd["piezo_line2"] = []
    sd["dloads"] = [list(SURCHARGE)]
    sd["dload_dirs"] = [SURCHARGE_DIR]
    sd["dloads2"] = []
    sd["dload2_dirs"] = []
    sd["circles"] = [dict(c) for c in CIRCLES]
    sd["circular"] = True
    sd["max_depth"] = None
    sd["tcrack_depth"] = 0.0
    sd["tcrack_water"] = 0.0
    sd["water_loads"] = "auto"

    path = os.path.join(tmpdir, "w02_source.xlsx")
    save_slope_data_to_xlsx(sd, path)
    return _quiet(load_slope_data, path)


def default_mapping(layers):
    """The wizard's own defaults for a drawing XSLOPE wrote: one row per layer,
    the target ``suggest_dxf_target`` proposes, and the Material column's default
    (the layer name, with any ``PROFILE_`` prefix stripped)."""
    mapping = {}
    for name, geom in layers.items():
        target = suggest_dxf_target(name, geom)
        material = None
        if target in ("material_zone", "profile"):
            material = name[8:] if name.upper().startswith("PROFILE_") else name
            material = MATERIAL_NAMES.get(material, material)
        mapping[name] = {"target": target, "material": material}
    return mapping


def import_drawing(dxf_path):
    """Run Studio's DXF import on ``dxf_path`` and return ``(slope_data, notes)``.

    ``ProjectDocument`` is a ``QObject`` for its signals, which work without a
    running application, so this needs no Studio window and no display.
    """
    from studio.document import ProjectDocument

    layers, warnings = read_dxf_layers(dxf_path)
    doc = ProjectDocument()
    notes = doc.build_from_dxf_mapping(layers, default_mapping(layers))
    return doc.slope_data, list(notes) + list(warnings), layers


def fill_properties(sd):
    """Type in what a drawing cannot carry: the three material rows and the
    surcharge pressure. Everything else in the imported model is the drawing's."""
    by_name = {m["name"]: m for m in MATERIALS}
    for mat in sd["materials"]:
        spec = by_name.get(mat["name"])
        if spec is None:
            raise SystemExit("imported material %r is not one of the section's"
                             % mat["name"])
        mat.update(spec)

    if len(sd.get("dloads") or []) != 1:
        raise SystemExit("expected one imported load block, got %d"
                         % len(sd.get("dloads") or []))
    for point, source in zip(sd["dloads"][0], SURCHARGE):
        point["Normal"] = source["Normal"]
    sd["dload_dirs"] = [SURCHARGE_DIR]
    return sd


def build():
    os.makedirs(OUT_DIR, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        source = build_source_model(tmp)
        _quiet(export_dxf, source, DXF_OUT)
        print("-> %s" % os.path.relpath(DXF_OUT, REPO_ROOT))

        imported, notes, layers = import_drawing(DXF_OUT)

    print("   layers: %s" % ", ".join(layers))
    for note in notes:
        print("   note: %s" % note)

    # The drawing's geometry has to arrive intact — the tutorial's claim is that
    # everything but the properties comes across.
    if len(imported["polygons"]) != len(ZONES):
        raise SystemExit("imported %d zone(s), drew %d"
                         % (len(imported["polygons"]), len(ZONES)))
    if len(imported["piezo_line"]) != len(PIEZO):
        raise SystemExit("piezometric line came back with %d point(s), drew %d"
                         % (len(imported["piezo_line"]), len(PIEZO)))
    if len(imported["circles"]) != len(CIRCLES):
        raise SystemExit("imported %d circle(s), drew %d"
                         % (len(imported["circles"]), len(CIRCLES)))
    if imported.get("water_loads") != "auto":
        raise SystemExit("water loads imported %r; the surcharge does not trace the "
                         "water and nothing here is ponded" % imported.get("water_loads"))

    for circle, drawn in zip(imported["circles"], CIRCLES):
        for key in ("Xo", "Yo", "R"):
            if abs(circle[key] - drawn[key]) > 0.05:
                raise SystemExit("circle %s came back as %.4f, drew %.4f"
                                 % (key, circle[key], drawn[key]))

    fill_properties(imported)
    imported["unit_system"] = source["unit_system"]
    imported["gamma_water"] = source["gamma_water"]
    save_slope_data_to_xlsx(imported, XLSX_OUT)
    print("-> %s" % os.path.relpath(XLSX_OUT, REPO_ROOT))

    back = _quiet(load_slope_data, XLSX_OUT)
    fs_cache, converged, _path, _cache = circular_search(back, "spencer", num_slices=40)
    if not fs_cache or fs_cache[0]["FS"] >= 9999:
        raise SystemExit("the imported model's search found no valid surface")
    best = fs_cache[0]
    print("   Spencer search: FS = %.4f  (converged=%s)" % (best["FS"], converged))
    print("   critical circle: Xo = %.2f, Yo = %.2f, Depth = %.2f"
          % (best["Xo"], best["Yo"], best["Depth"]))

    # The round trip proved on the answer rather than on the vertices: the model
    # the drawing came from and the model the drawing produces have to search to
    # the same factor of safety.
    src_cache, _c, _p, _cc = circular_search(source, "spencer", num_slices=40)
    print("   source model searches to FS = %.4f (Δ %.4f)"
          % (src_cache[0]["FS"], best["FS"] - src_cache[0]["FS"]))
    if abs(best["FS"] - src_cache[0]["FS"]) > 0.002:
        raise SystemExit("the drawing's round trip moved the factor of safety")
    return DXF_OUT, XLSX_OUT, best["FS"]


if __name__ == "__main__":
    build()
