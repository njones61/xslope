"""RS2's own meshes for RS2 #68 (Part III), "Stability of seismically loaded slopes"
(Loukidis, Bandini & Salgado 2003), written as the corpus files' mesh sidecars.

RS2's finite element k_c on this problem comes off a coarse mesh (2,616 to 4,193
six-node triangles), and at that density the answer depends on how the elements are
laid out, not only on their size: XSLOPE's own mesher at RS2's element count gives a
different k_c from RS2's mesh on Case 1. The FEM k_c rows therefore solve on RS2's
mesh, so that the mesh cannot be the difference and the comparison tests the solver.

This module reads each case's vendor model (the ``.fea`` inside
``slope stability #068_0N.fez``) and writes ``rs2_68{a,b,c}_mesh.json`` beside the
workbook that ``build_rs2.py`` builds from hard-coded geometry. ``load_slope_data``
picks the sidecar up as the model's mesh, and the rows' tags say ``mesh=file``.

Before writing, each mesh is checked against the workbook it will sit beside, not
against itself:

* every midside node lies at the midpoint of its edge (the RS2 -> XSLOPE tri6 slot
  reorder, ``rs2_transient_seep.validate_tri6_midsides``);
* each element's material, mapped by NAME, is the material of the workbook polygon
  that contains the element's centroid;
* the mesh covers exactly the workbook's region (equal area, every node inside it).

The sidecars are committed corpus artifacts; this generator is dev-side tooling run
against the local vendor archive, which is never committed. With the archive absent it
prints a note and leaves the committed sidecars in place.

Run:

    python benchmarks/rocscience/rs2_68_vendor_mesh.py [VENDOR_DIR]
"""

import json
import os
import re
import sys
import tempfile
import zipfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

from rs2_transient_seep import parse_computed_fea, validate_tri6_midsides  # noqa: E402

#: The folder holding the RS2 native verification models; override with
#: $XSLOPE_RS2_SLOPE_VENDOR_DIR or the first argument.
_DEFAULT_VENDOR_DIR = os.environ.get(
    "XSLOPE_RS2_SLOPE_VENDOR_DIR",
    os.path.expanduser("~/python_projects/vendor_files/rocscience_downloads/work/"
                       "RS2_Slope Stability Verification RS2 and Slide2 Import/"
                       "Slope Stability Verification"))

OUT = os.path.join(ROOT, "docs", "verification", "files", "rocscience")

_CASES = [
    ("slope stability #068_01.fez", "rs2_68a"),   # Case 1, r_u = 0.5
    ("slope stability #068_02.fez", "rs2_68b"),   # Case 2, dry
    ("slope stability #068_03.fez", "rs2_68c"),   # Case 3, three layers
]


def _read_fea(fez_path, tmp):
    """Extract the single ``.fea`` member of a ``.fez`` into ``tmp``; return its path."""
    with zipfile.ZipFile(fez_path) as z:
        names = [n for n in z.namelist() if n.lower().endswith(".fea")]
        if len(names) != 1:
            raise ValueError(f"{os.path.basename(fez_path)}: expected one .fea, "
                             f"found {names}")
        out = os.path.join(tmp, "model.fea")
        with z.open(names[0]) as src, open(out, "wb") as dst:
            dst.write(src.read())
    return out


def _element_material_names(fea_path):
    """The ``material:`` name of each element row, in file order."""
    lines = open(fea_path, encoding="latin-1").read().splitlines()
    ei = next(i for i, l in enumerate(lines) if l.startswith("elements:"))
    row = re.compile(r"\s*\d+\s+type:\s*\S+\s+nodes:\s*\[[^\]]*\]\s+material:\s*(\S+)")
    names = []
    for l in lines[ei + 1:]:
        m = row.match(l)
        if not m:
            break
        names.append(m.group(1))
    return names


def _check_against_workbook(mesh, sd, base):
    """Material by centroid and region coverage, against the built workbook."""
    from shapely.geometry import Point, Polygon
    from shapely.ops import unary_union
    from xslope.mesh import get_material_polygons
    # {coords, mat_id (0-based)} per material polygon; element_materials is 1-based
    shapes = [(Polygon(p["coords"]), int(p["mat_id"]) + 1) for p in get_material_polygons(sd)]
    nodes, els = mesh["nodes"], mesh["elements"]
    corners = els[:, :3]
    cent = nodes[corners].mean(axis=1)
    bad = 0
    for (x, y), mat in zip(cent, mesh["element_materials"]):
        hit = [mid for poly, mid in shapes if poly.contains(Point(x, y))]
        if hit != [mat]:
            bad += 1
    if bad:
        raise ValueError(f"{base}: {bad} elements carry a material other than the "
                         f"workbook polygon at their centroid")
    region = unary_union([poly for poly, _ in shapes])
    p = nodes[corners]
    area = 0.5 * np.abs((p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1])
                        - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1])).sum()
    if abs(area - region.area) > 1e-6 * region.area:
        raise ValueError(f"{base}: mesh area {area:.3f} != model area {region.area:.3f}")
    # The builder transcribes vertices to 1 mm (Case 3's face corner (74.693, 41.077) is
    # RS2's (74.69273, 41.07709)), so face nodes can sit a fraction of a millimeter out.
    grown = region.buffer(1e-5 * max(region.bounds[2] - region.bounds[0], 1.0))
    outside = sum(not grown.contains(Point(x, y)) for x, y in nodes)
    if outside:
        raise ValueError(f"{base}: {outside} mesh nodes lie outside the model")
    return area


def generate(vendor_dir=_DEFAULT_VENDOR_DIR, out_dir=OUT, verbose=True):
    """Write the three mesh sidecars; skip (return []) where the archive is absent."""
    from xslope.fileio import load_slope_data
    written = []
    for fez_name, base in _CASES:
        fez = os.path.join(vendor_dir, fez_name)
        if not os.path.exists(fez):
            if verbose:
                print(f"  [skip] {fez_name} not found; keeping committed {base}_mesh.json")
            continue
        with tempfile.TemporaryDirectory() as tmp:
            fea = _read_fea(fez, tmp)
            mesh = parse_computed_fea(fea)
            rs2_names = _element_material_names(fea)
        sd = load_slope_data(os.path.join(out_dir, f"{base}.xlsx"))
        names = [m["name"] for m in sd["materials"]]
        unknown = sorted(set(rs2_names) - set(names))
        if unknown or len(rs2_names) != len(mesh["elements"]):
            raise ValueError(f"{base}: RS2 materials {unknown} not in the workbook "
                             f"{names}, or the element rows do not line up")
        mesh["element_materials"] = np.array([names.index(n) + 1 for n in rs2_names])
        worst_mid = validate_tri6_midsides(mesh)
        area = _check_against_workbook(mesh, sd, base)
        with open(os.path.join(out_dir, f"{base}_mesh.json"), "w") as f:
            json.dump({"nodes": mesh["nodes"].tolist(),
                       "elements": mesh["elements"].tolist(),
                       "element_types": mesh["element_types"].tolist(),
                       "element_materials": mesh["element_materials"].tolist()}, f)
        written.append(base)
        if verbose:
            print(f"  [{base}] nodes={len(mesh['nodes'])} elements={len(mesh['elements'])} "
                  f"area={area:.1f} midside_dev={worst_mid:.1e}")
    return written


if __name__ == "__main__":
    vdir = sys.argv[1] if len(sys.argv) > 1 else _DEFAULT_VENDOR_DIR
    print(f"RS2 #68 vendor meshes (vendor dir: {vdir})")
    generate(vdir)
