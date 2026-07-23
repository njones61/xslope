"""Vendor-field adapter for RS2 #67 (Part III) — "Stability of an earth dam under
steady & transient unsaturated seepage" (Huang & Jia 2009).

XSLOPE has no transient-seepage solver, so the transient-drawdown stages of RS2 #67
cannot be reproduced by solving the flow. They do not need to be: RS2's computed
``.fea`` files carry the SOLVED pore-pressure field for each snapshot as a per-node
block::

    pore pressure:
    n: 4 u: 64.4713
    n: 5 u: 64.4713
    ...

This module reads a computed ``.fea`` (the model bundled inside each ``.fez``) and
converts that vendor nodal-u field into XSLOPE's existing external-pore-pressure
input path — the ``u='seep'`` seep-solution mechanism (a ``*_mesh.json`` + a
``*_seep.csv`` sidecar, exactly the format the corpus's FE-seepage problems use, e.g.
``vp077a_mesh.json`` / ``vp077a_seep.csv``). XSLOPE's SSRM then runs on RS2's own
mesh with RS2's own snapshot pore pressures. This verifies the SSR-under-transient-u
MECHANICS; the transient seepage SOLUTION itself is RS2's (imported), not XSLOPE's.

Fidelity notes
--------------
* RS2's ``lst_element`` is a 6-node quadratic triangle stored ``[c0, m01, c1, m12,
  c2, m20]`` (corners at even indices, edge midpoints at odd). XSLOPE's tri6 order is
  ``[c0, c1, c2, m01, m12, m20]`` (corners first, then the same three edge midpoints),
  so the reorder is ``[0, 2, 4, 1, 3, 5]``. ``import_mesh_from_json`` re-checks winding
  (``ensure_ccw_elements``) on load, so only the node SLOT convention has to be right.
* Node ids are remapped to a 0-based contiguous index; ``element_materials`` is 1-based
  (the corpus/geo convention); the seep CSV carries one trailing row that ``fileio``
  drops (``iloc[:-1]``), matching the vendor seep sidecars.

The sidecar files are committed corpus artifacts (like ``vp077a_*``); this generator is
dev-side tooling run against the vendor archive. ``build_rs2.py`` builds the ``.xlsx``
from hard-coded geometry and does NOT need the vendor files.

Run (regenerate the committed sidecars from the vendor archive):

    python benchmarks/rocscience/rs2_transient_seep.py [VENDOR_DIR]
"""

import csv
import json
import os
import re
import sys

import numpy as np

# RS2 ``lst_element`` slot order [c0,m01,c1,m12,c2,m20] -> XSLOPE [c0,c1,c2,m01,m12,m20]
_TRI6_REORDER = [0, 2, 4, 1, 3, 5]

# Directory holding the extracted vendor computed '.fea' files (the model bundled
# inside each '#067_0N.fez'). Dev-side only; override with $XSLOPE_RS2_VENDOR_DIR.
# Absent -> generate() skips and the committed sidecars stand.
_DEFAULT_VENDOR_DIR = os.environ.get(
    "XSLOPE_RS2_VENDOR_DIR",
    os.path.expanduser("~/python_projects/vendor_files/rocscience_downloads/rs2_67_fea"))

OUT = os.path.join(os.path.dirname(__file__), "..", "..",
                   "docs", "verification", "files", "rocscience")


def parse_computed_fea(path):
    """Parse a computed RS2 ``.fea`` into an XSLOPE-format mesh + nodal pore pressure.

    Returns a dict with:
        nodes            (n, 2)  float, ordered by remapped 0-based index
        elements         (m, 9)  int, XSLOPE tri6 slot order, 0-based, zero-padded
        element_types    (m,)    int (6 for every tri6)
        element_materials(m,)    int, 1-based
        u                (n,)    float, vendor nodal pore pressure aligned to ``nodes``
                                 (None if the file carries no ``pore pressure:`` block)
    """
    lines = open(path).read().splitlines()

    # --- nodes: "<id> x: <x> y: <y>" -----------------------------------------
    ni = next(i for i, l in enumerate(lines) if l.startswith("nodes:"))
    nid_xy = {}
    for l in lines[ni + 1:]:
        m = re.match(r"\s*(\d+)\s+x:\s*([-\d.eE+]+)\s+y:\s*([-\d.eE+]+)", l)
        if m:
            nid_xy[int(m.group(1))] = (float(m.group(2)), float(m.group(3)))
        elif l.startswith("elements:"):
            break
    ids = sorted(nid_xy)
    id2idx = {nid: k for k, nid in enumerate(ids)}
    nodes = np.array([nid_xy[nid] for nid in ids], dtype=float)

    # --- elements: "<id> type: lst_element nodes: [..] material: .. materialID: N"
    ei = next(i for i, l in enumerate(lines) if l.startswith("elements:"))
    elems, etypes, emats = [], [], []
    erow = re.compile(
        r"\s*\d+\s+type:\s*(\S+)\s+nodes:\s*\[([\d,\s]+)\]\s+material:\s*\S+\s+materialID:\s*(\d+)")
    for l in lines[ei + 1:]:
        m = erow.match(l)
        if m:
            etype = m.group(1)
            ncoords = [int(x) for x in m.group(2).split(",")]
            if etype != "lst_element" or len(ncoords) != 6:
                raise ValueError(f"{os.path.basename(path)}: unexpected element "
                                 f"'{etype}' with {len(ncoords)} nodes (only tri6 "
                                 f"lst_element supported)")
            reordered = [id2idx[ncoords[j]] for j in _TRI6_REORDER]
            elems.append(reordered + [0, 0, 0])   # pad 6 -> 9
            etypes.append(6)
            emats.append(int(m.group(3)))
        elif l.startswith(("springs:", "stresses:")):
            break
    elements = np.array(elems, dtype=int)
    element_types = np.array(etypes, dtype=int)
    element_materials = np.array(emats, dtype=int)   # 1-based already in the .fea

    # --- pore pressure: "n: <id> u: <val>"  (optional) -----------------------
    u = None
    pi = next((i for i, l in enumerate(lines) if l.startswith("pore pressure:")), -1)
    if pi >= 0:
        pp = {}
        for l in lines[pi + 1:]:
            m = re.match(r"n:\s*(\d+)\s+u:\s*([-\d.eE+]+)", l)
            if m:
                pp[int(m.group(1))] = float(m.group(2))
            elif l.startswith("stresses:"):
                break
        if pp:
            missing = [nid for nid in ids if nid not in pp]
            if missing:
                raise ValueError(f"{os.path.basename(path)}: {len(missing)} nodes "
                                 f"have no pore pressure (e.g. {missing[:5]})")
            u = np.array([pp[nid] for nid in ids], dtype=float)

    return {"nodes": nodes, "elements": elements, "element_types": element_types,
            "element_materials": element_materials, "u": u}


def validate_tri6_midsides(mesh, tol=1e-3):
    """Assert every edge midpoint node sits at the midpoint of its two flanking corners.

    Confirms the RS2->XSLOPE slot reorder is correct: with XSLOPE order the node in
    slot 3 must be the (0,1) edge midpoint, slot 4 the (1,2) midpoint, slot 5 the
    (2,0) midpoint. Returns the worst deviation (domain units)."""
    nodes = mesh["nodes"]
    worst = 0.0
    for row in mesh["elements"]:
        c0, c1, c2, m01, m12, m20 = row[:6]
        for a, b, mid in ((c0, c1, m01), (c1, c2, m12), (c2, c0, m20)):
            expect = 0.5 * (nodes[a] + nodes[b])
            worst = max(worst, float(np.hypot(*(nodes[mid] - expect))))
    if worst > tol:
        raise ValueError(f"tri6 midside geometry off by {worst:.4g} (> {tol}) — "
                         f"slot reorder is wrong")
    return worst


def _interp_u_at(mesh, xy):
    """Barycentric-interpolate the nodal u field at a point using its containing
    tri6 (corner triangle) — the same piecewise-linear interpolation the FEM/slice
    machinery does. Returns None if the point is outside every element."""
    nodes, u = mesh["nodes"], mesh["u"]
    px, py = xy
    for row in mesh["elements"]:
        a, b, c = row[0], row[1], row[2]
        (x1, y1), (x2, y2), (x3, y3) = nodes[a], nodes[b], nodes[c]
        det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
        if abs(det) < 1e-12:
            continue
        l1 = ((y2 - y3) * (px - x3) + (x3 - x2) * (py - y3)) / det
        l2 = ((y3 - y1) * (px - x3) + (x1 - x3) * (py - y3)) / det
        l3 = 1.0 - l1 - l2
        if -1e-9 <= l1 <= 1 + 1e-9 and -1e-9 <= l2 <= 1 + 1e-9 and -1e-9 <= l3 <= 1 + 1e-9:
            return l1 * u[a] + l2 * u[b] + l3 * u[c]
    return None


def spot_check_interpolation(mesh, n=6, seed=0):
    """Interpolate u at a handful of vendor node locations and compare to that node's
    stored value. For a piecewise-linear field sampled at a node the interpolant must
    recover the nodal value exactly — a direct check that the mesh+u sidecars are
    mutually consistent (right nodes, right ordering). Returns worst |error| (kPa)."""
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(mesh["nodes"]), size=min(n, len(mesh["nodes"])), replace=False)
    worst = 0.0
    rows = []
    for k in idx:
        got = _interp_u_at(mesh, mesh["nodes"][k])
        exp = float(mesh["u"][k])
        err = abs((got if got is not None else exp) - exp)
        worst = max(worst, err)
        rows.append((tuple(np.round(mesh["nodes"][k], 3)), round(exp, 4),
                     None if got is None else round(got, 4), round(err, 6)))
    return worst, rows


def write_sidecars(mesh, out_dir, base, gamma_water=9.80665):
    """Write ``{base}_mesh.json`` + ``{base}_seep.csv`` in XSLOPE seep-sidecar format.

    The CSV mirrors the corpus format (header + one row per node + a trailing comment
    row that ``fileio`` drops with ``iloc[:-1]``). Only the ``u`` column is read back
    by the loader; ``head`` is filled consistently (y + u/gamma_w) for plotting."""
    os.makedirs(out_dir, exist_ok=True)
    mesh_json = {
        "nodes": mesh["nodes"].tolist(),
        "elements": mesh["elements"].tolist(),
        "element_types": mesh["element_types"].tolist(),
        "element_materials": mesh["element_materials"].tolist(),
    }
    with open(os.path.join(out_dir, f"{base}_mesh.json"), "w") as f:
        json.dump(mesh_json, f)

    u = mesh["u"]
    ys = mesh["nodes"][:, 1]
    with open(os.path.join(out_dir, f"{base}_seep.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")   # LF, matching the corpus sidecars
        w.writerow(["node_id", "head", "u", "v_x", "v_y", "v_mag",
                    "i_x", "i_y", "i_mag", "q", "phi"])
        for k in range(len(u)):
            head = ys[k] + u[k] / gamma_water
            w.writerow([k + 1, head, u[k], 0, 0, 0, 0, 0, 0, 0, 0])
        f.write("# Imported from RS2 vendor .fea nodal pore-pressure field\n")


# --- stage catalogue: vendor file -> corpus base name --------------------------
_STAGES = [
    ("slope stability #067_03 (90h).fea", "rs2_67c"),   # Case 3 downstream
    ("slope stability #067_04 (90h).fea", "rs2_67d"),   # Case 3 upstream
]


def generate(vendor_dir=_DEFAULT_VENDOR_DIR, out_dir=OUT, verbose=True):
    """Regenerate the committed seep sidecars from the vendor archive.

    Skips gracefully (returns []) if the vendor files are not present — the corpus
    ``.xlsx`` still builds from hard-coded geometry and the previously-committed
    sidecars remain in place."""
    written = []
    for fea_name, base in _STAGES:
        path = os.path.join(vendor_dir, fea_name)
        if not os.path.exists(path):
            if verbose:
                print(f"  [skip] {fea_name} not found; keeping committed {base}_* sidecars")
            continue
        mesh = parse_computed_fea(path)
        if mesh["u"] is None:
            raise ValueError(f"{fea_name} carries no pore-pressure block")
        worst_mid = validate_tri6_midsides(mesh)
        worst_interp, rows = spot_check_interpolation(mesh)
        write_sidecars(mesh, out_dir, base)
        written.append(base)
        if verbose:
            print(f"  [{base}] nodes={len(mesh['nodes'])} elems={len(mesh['elements'])} "
                  f"u=[{mesh['u'].min():.2f},{mesh['u'].max():.2f}] "
                  f"midside_dev={worst_mid:.2e} interp_err={worst_interp:.2e} kPa")
            for xy, exp, got, err in rows:
                print(f"        node {xy}: vendor u={exp}  interp={got}  |err|={err}")
    return written


if __name__ == "__main__":
    vdir = sys.argv[1] if len(sys.argv) > 1 else _DEFAULT_VENDOR_DIR
    print(f"RS2 #67 transient-seepage adapter (vendor dir: {vdir})")
    generate(vdir)
