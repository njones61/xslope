"""Generate deliberately "messy" DXF variants to exercise the importer's
validation and error-handling paths.

The clean fixtures (_build_dxf_from_poly.py) are ideal CAD: one closed
LWPOLYLINE per zone, one layer per material. Real-world exports are not.
These variants reproduce the common defects the commercial packages warn
about, so the importer has something to choke on (and recover from):

  *_messy_unclosed.dxf       LWPOLYLINEs left open (close flag not set, no
                             repeated end vertex). Importer must detect and
                             auto-close, with a warning.
  *_messy_polyline.dxf       Old-style heavyweight POLYLINE entities instead
                             of LWPOLYLINE. Importer must read both types.
  *_messy_line_segments.dxf  Each zone drawn as loose, unconnected LINE
                             entities (never joined into a polyline). Importer
                             must stitch segments into rings or reject.
  *_messy_arcs.dxf           LWPOLYLINE with an arc bulge on one segment
                             (curved material interface). Importer must flatten
                             arcs to vertices.

All are derived from the sea-trench geometry. Output lands next to this file.

Run from the repo root:  python3 poly_test/_build_dxf_messy.py
"""
import os

import ezdxf

from _build_dxf_from_poly import (
    _read_polygons,
    _read_materials,
    _layer_name,
    _dedupe_closing_vertex,
)
import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCE = os.path.join(HERE, "xslope_sea_trench_poly.xlsx")
_ACI = [40, 30, 8, 5, 3, 1, 6, 2]


def _load():
    wb = openpyxl.load_workbook(SOURCE, data_only=False)
    polys = _read_polygons(wb)
    mat_names = _read_materials(wb)
    layer_for, doc_layers = {}, {}
    for mid in sorted({p["mat_id"] for p in polys}):
        layer_for[mid] = _layer_name(mat_names.get(mid), mid)
    return polys, mat_names, layer_for


def _new_doc(layer_for):
    doc = ezdxf.new("R2010")
    for i, (mid, lname) in enumerate(sorted(layer_for.items())):
        if lname not in doc.layers:
            doc.layers.add(name=lname, color=_ACI[i % len(_ACI)])
    return doc


def build_unclosed(polys, layer_for, path):
    """LWPOLYLINEs with close=False and no repeated closing vertex."""
    doc = _new_doc(layer_for)
    msp = doc.modelspace()
    for p in polys:
        coords = _dedupe_closing_vertex(p["coords"])
        msp.add_lwpolyline(coords, close=False,
                           dxfattribs={"layer": layer_for[p["mat_id"]]})
    doc.saveas(path)


def build_polyline(polys, layer_for, path):
    """Old-style heavyweight POLYLINE entities (add_polyline2d), closed."""
    doc = _new_doc(layer_for)
    msp = doc.modelspace()
    for p in polys:
        coords = _dedupe_closing_vertex(p["coords"])
        pl = msp.add_polyline2d(coords, dxfattribs={"layer": layer_for[p["mat_id"]]})
        pl.close(True)
    doc.saveas(path)


def build_line_segments(polys, layer_for, path):
    """Each polygon edge as a separate, unconnected LINE entity."""
    doc = _new_doc(layer_for)
    msp = doc.modelspace()
    for p in polys:
        coords = _dedupe_closing_vertex(p["coords"])
        layer = layer_for[p["mat_id"]]
        n = len(coords)
        for i in range(n):
            a = coords[i]
            b = coords[(i + 1) % n]  # wrap to close the ring edge-by-edge
            msp.add_line(a, b, dxfattribs={"layer": layer})
    doc.saveas(path)


def build_arcs(polys, layer_for, path):
    """Closed LWPOLYLINE with an arc bulge on one segment per polygon.

    Bulge = tan(theta/4); 0.25 ~= a gentle 56-deg arc. Applied to the first
    segment of each ring to make a curved material interface.
    """
    doc = _new_doc(layer_for)
    msp = doc.modelspace()
    for p in polys:
        coords = _dedupe_closing_vertex(p["coords"])
        # format "xyseb": x, y, start_width, end_width, bulge
        pts = []
        for i, (x, y) in enumerate(coords):
            bulge = 0.25 if i == 0 else 0.0
            pts.append((x, y, 0.0, 0.0, bulge))
        msp.add_lwpolyline(pts, format="xyseb", close=True,
                           dxfattribs={"layer": layer_for[p["mat_id"]]})
    doc.saveas(path)


def _summarize(path):
    doc = ezdxf.readfile(path)
    msp = doc.modelspace()
    counts = {}
    for e in msp:
        counts[e.dxftype()] = counts.get(e.dxftype(), 0) + 1
    bulged = 0
    for e in msp.query("LWPOLYLINE"):
        if any(abs(pt[4]) > 1e-9 for pt in e.get_points("xyseb")):
            bulged += 1
    closed_lw = sum(1 for e in msp.query("LWPOLYLINE") if e.closed)
    extra = f", {bulged} with arc bulges" if bulged else ""
    extra += f", {closed_lw} closed" if msp.query("LWPOLYLINE") else ""
    ents = ", ".join(f"{v}x {k}" for k, v in sorted(counts.items()))
    return f"{ents}{extra}"


if __name__ == "__main__":
    polys, mat_names, layer_for = _load()
    stem = os.path.join(HERE, "xslope_sea_trench_messy_")
    jobs = [
        ("unclosed",      build_unclosed),
        ("polyline",      build_polyline),
        ("line_segments", build_line_segments),
        ("arcs",          build_arcs),
    ]
    for tag, fn in jobs:
        path = f"{stem}{tag}.dxf"
        fn(polys, layer_for, path)
        print(f"{os.path.basename(path)}")
        print(f"   {_summarize(path)}")
