"""SEEP-2 benchmark: xslope vs the actual SEEP2D code on identical mesh + BCs.

Exports the Johnson Reservoir problem (docs/seep/files/xslope_johnson_res.xlsx)
to a SEEP2D .s2d input file — same tri3 mesh topology, same boundary
conditions, same material parameters — runs the original USACE/WES SEEP2D
Fortran program on it, and compares nodal heads and total discharge against
xslope's seepage solution on the very same arrays.

SEEP2D binary: compiled from ../xslope_private/ref_docs/ref_docs_seep/seep2d_fortran/src/seep2d.f
with gfortran (-std=legacy), after three mechanical patches for gfortran:
  - comment out `use DFPORT` / `use DFLIB` (DEC/Intel portability modules;
    getarg/iargc are gfortran intrinsics)
  - logical comparisons `.eq./.ne. .true./.false.` -> `.eqv./.neqv.`
Build:  cd <builddir> && SDKROOT=$(xcrun --show-sdk-path) \
        gfortran -o seep2d seep2d.f -std=legacy -O1

The mesh is RCM-renumbered before export (SEEP2D's banded solver has
MXBNDW=700; gmsh node ordering is not banded).

Run from the repo root:
  PYTHONPATH=. python3 benchmarks/run_seep2d_compare.py [path/to/seep2d-binary]
"""
import io
import os
import re
import subprocess
import sys
import contextlib

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

from xslope.fileio import load_slope_data
from xslope.mesh import get_material_polygons, build_mesh_from_polygons
from xslope.seep import build_seep_data, run_seepage_analysis

MODEL = "docs/seep/files/xslope_johnson_res.xlsx"
SEEP2D_BIN = sys.argv[1] if len(sys.argv) > 1 else "/tmp/seep2d_build/seep2d"
WORKDIR = "/tmp/seep2d_johnson"


def build_tri3_seep_data():
    """Johnson Reservoir on a tri3 mesh (SEEP2D = linear triangles)."""
    sd = load_slope_data(MODEL)
    polys = get_material_polygons(sd)
    xs = [x for x, _ in sd["ground_surface"].coords]
    target = (max(xs) - min(xs)) / 120.0   # same rule as the samples harness
    mesh = build_mesh_from_polygons(polys, target, "tri3")
    seep = build_seep_data(mesh, sd, seep_bc=1)
    return sd, seep


def rcm_renumber(seep):
    """Permute nodes with reverse Cuthill-McKee to minimize bandwidth."""
    nodes = np.asarray(seep["nodes"])
    elements = np.asarray(seep["elements"])
    etypes = np.asarray(seep["element_types"])
    n = len(nodes)
    ii, jj = [], []
    for e, et in zip(elements, etypes):
        en = e[:et]
        for a in en:
            for b in en:
                ii.append(a); jj.append(b)
    A = coo_matrix((np.ones(len(ii)), (ii, jj)), shape=(n, n)).tocsr()
    perm = reverse_cuthill_mckee(A, symmetric_mode=True)   # new_order[i] = old index
    inv = np.empty(n, dtype=int)
    inv[perm] = np.arange(n)                               # old -> new
    out = dict(seep)
    out["nodes"] = nodes[perm]
    out["bc_type"] = np.asarray(seep["bc_type"])[perm]
    out["bc_values"] = np.asarray(seep["bc_values"])[perm]
    new_el = elements.copy()
    for k, (e, et) in enumerate(zip(elements, etypes)):
        new_el[k, :et] = inv[e[:et]]
        if et == 3 and new_el.shape[1] > 3:
            new_el[k, 3:] = 0
    out["elements"] = new_el
    # bandwidth check
    bw = 0
    for e, et in zip(new_el, etypes):
        en = e[:et]
        bw = max(bw, int(en.max() - en.min()))
    return out, bw


def write_s2d(seep, gamma_w, path):
    nodes = np.asarray(seep["nodes"])
    elements = np.asarray(seep["elements"])
    etypes = np.asarray(seep["element_types"])
    emats = np.asarray(seep["element_materials"])       # 1-based
    bc_type = np.asarray(seep["bc_type"])
    bc_val = np.asarray(seep["bc_values"])
    k1 = np.asarray(seep["k1_by_mat"]); k2 = np.asarray(seep["k2_by_mat"])
    ang = np.asarray(seep["angle_by_mat"])
    kr0 = np.asarray(seep["kr0_by_mat"]); h0 = np.asarray(seep["h0_by_mat"])
    nmat = len(k1)
    with open(path, "w") as f:
        f.write("xslope Johnson Reservoir export for SEEP2D cross-check\n")
        # format(4i5, 1x, a4, f10.0, 4x, a1, f10.0, i5)
        f.write(f"{len(nodes):5d}{len(elements):5d}{nmat:5d}{0:5d} PLNE"
                f"{0.0:10.1f}    F{gamma_w:10.2f}{1:5d}\n")
        for m in range(nmat):                            # format(i5, 5f15.0)
            f.write(f"{m+1:5d}{k1[m]:15.6f}{k2[m]:15.6f}{ang[m]:15.6f}"
                    f"{kr0[m]:15.6f}{h0[m]:15.6f}\n")
        for i, (x, y) in enumerate(nodes):               # format(i5, i2, i3, 3f15.0)
            nbc = int(bc_type[i])
            fx = float(bc_val[i]) if nbc == 1 else (float(y) if nbc == 2 else 0.0)
            f.write(f"{i+1:5d}{0:2d}{nbc:3d}{float(x):15.6f}{float(y):15.6f}"
                    f"{fx:15.6f}\n")
        for k, (e, et) in enumerate(zip(elements, etypes)):  # format(6i5)
            n1, n2, n3 = (int(v) + 1 for v in e[:3])
            n4 = n3 if et == 3 else int(e[3]) + 1
            f.write(f"{k+1:5d}{n1:5d}{n2:5d}{n3:5d}{n4:5d}{int(emats[k]):5d}\n")


def run_seep2d(s2d_path, workdir):
    sup = os.path.join(workdir, "model.sup")
    out = os.path.join(workdir, "model.out")
    if os.path.exists(out):
        os.remove(out)
    with open(sup, "w") as f:
        f.write(f"SEEPSUP\nSEEP {os.path.basename(s2d_path)}\n"
                f"ODAT {os.path.basename(out)}\n")
    # the program PAUSEs at the end; a newline on stdin terminates it cleanly
    subprocess.run([SEEP2D_BIN, os.path.basename(sup)], cwd=workdir,
                   input="\n", text=True, capture_output=True, timeout=600)
    return out


def parse_seep2d_out(path, n_nodes):
    """Extract nodal heads and boundary flows from the ODAT file."""
    heads = np.full(n_nodes, np.nan)
    flows = np.zeros(n_nodes)
    efmt = re.compile(r"^[-+]?\d*\.?\d+E[-+]\d+$", re.I)
    in_table = False
    for line in open(path, encoding="latin-1"):
        if "Nodal Flows and Heads" in line:
            in_table = True
            continue
        if in_table and "Element Flowrates" in line:
            break
        if not in_table:
            continue
        toks = line.split()
        if len(toks) >= 2 and toks[0].isdigit() and efmt.match(toks[1]):
            nid = int(toks[0]) - 1
            if 0 <= nid < n_nodes:
                heads[nid] = float(toks[1])
                # flow value, if present, is the next E-format token after '%'
                if "%" in toks:
                    rest = toks[toks.index("%") + 1:]
                    for t in rest:
                        if efmt.match(t):
                            flows[nid] = float(t)
                            break
    return heads, flows


def main():
    os.makedirs(WORKDIR, exist_ok=True)
    with contextlib.redirect_stdout(io.StringIO()):
        sd, seep = build_tri3_seep_data()
    seep_r, bw = rcm_renumber(seep)
    print(f"mesh: {len(seep_r['nodes'])} nodes, {len(seep_r['elements'])} "
          f"tri3 elements, RCM half-bandwidth {bw}")
    assert bw < 700, "exceeds SEEP2D MXBNDW"

    # xslope solution on the SAME (renumbered) arrays
    with contextlib.redirect_stdout(io.StringIO()):
        sol = run_seepage_analysis(dict(seep_r), tol=1e-8)
    q_xs = sol["flowrate"]
    h_xs = sol["head"]

    # SEEP2D on the exported file
    s2d = os.path.join(WORKDIR, "model.s2d")
    write_s2d(seep_r, sd["gamma_water"], s2d)
    out = run_seep2d(s2d, WORKDIR)
    h_s2d, fl = parse_seep2d_out(out, len(seep_r["nodes"]))
    q_s2d = fl[fl > 0].sum()

    ok = np.isfinite(h_s2d)
    dh = h_xs[ok] - h_s2d[ok]
    print(f"\nSEEP-2  Johnson Reservoir, xslope vs SEEP2D (identical mesh/BCs)")
    print(f"  nodes compared:        {ok.sum()}/{len(h_s2d)}")
    print(f"  total discharge q:     xslope {q_xs:.4f}   SEEP2D {q_s2d:.4f}"
          f"   diff {100*(q_xs-q_s2d)/q_s2d:+.2f}%")
    print(f"  nodal head error:      max |dh| = {np.abs(dh).max():.4f}, "
          f"RMS = {np.sqrt((dh**2).mean()):.4f}  (head range "
          f"{h_xs.min():.1f}-{h_xs.max():.1f})")


if __name__ == "__main__":
    main()
