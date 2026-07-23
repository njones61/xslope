"""SEEP-2 benchmark: xslope vs the actual SEEP2D code on identical mesh + BCs.

Exports an xslope problem to a SEEP2D .s2d input file — same tri3 mesh topology,
same boundary conditions, same material parameters, same unsaturated law — runs
the original USACE/WES SEEP2D Fortran program on it, and compares nodal heads,
total discharge and the free-surface RELEASE POINT against xslope's seepage
solution on the very same arrays.

Both unsaturated models are exported (SEEP2D's `iuntyp`: 1 = linear front, with
uspar = kr0, h0; 2 = van Genuchten, with uspar = alpha, n). The two codes use the
identical Mualem-van Genuchten form — Se = [1+(a|psi|)^n]^-m, m = 1-1/n,
kr = sqrt(Se)[1-(1-Se^(1/m))^m]^2 — so alpha and n carry straight across. xslope
additionally floors kr at kr_min (SEEP2D does not); that floor has been checked
not to matter here. iuntyp is GLOBAL in SEEP2D, so a model that mixes unsaturated
laws across materials cannot be exported and is refused rather than fudged.

SEEP2D binary: compiled from ../xslope_private/ref_docs/ref_docs_seep/seep2d_fortran/src/seep2d.f
(with seep.inc alongside it) using gfortran (-std=legacy), after two mechanical
patches:
  - comment out `use DFPORT` / `use DFLIB` (DEC/Intel portability modules;
    getarg/iargc are gfortran intrinsics)
  - logical comparisons `.eq./.ne. .true./.false.` -> `.eqv./.neqv.`
Build:  cd <builddir> && SDKROOT=$(xcrun --show-sdk-path) \
        gfortran -o seep2d seep2d.f -std=legacy -O1

The mesh is RCM-renumbered before export (SEEP2D's banded solver has
MXBNDW=700; gmsh node ordering is not banded).

Run from the repo root:
  PYTHONPATH=. python3 benchmarks/run_seep2d_compare.py            # Johnson Reservoir
  PYTHONPATH=. python3 benchmarks/run_seep2d_compare.py --model docs/verification/files/rocscience_gw/gw004.xlsx \
      --target-size 0.147 --max-iter 400
  PYTHONPATH=. python3 benchmarks/run_seep2d_compare.py --gw       # the whole GW free-surface panel
"""
import argparse
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
SEEP2D_BIN = "/tmp/seep2d_build/seep2d"
WORKDIR = "/tmp/seep2d_johnson"

# The groundwater-corpus free-surface panel (task #30): every GW problem that has
# an exit face, so the release point is actually exercised. (stem, target_size, max_iter)
GW_PANEL = [("gw004", 0.147, 400), ("gw006a", 1.0, 2000), ("gw009a", None, 400),
            ("gw010", 0.25, 1500), ("gw012", 1.0, 1500), ("gw013", 1.0, 1500)]

IUNTYP = {"lf": 1, "vg": 2}


def unsat_code(slope_data):
    """SEEP2D's global iuntyp for this model. Raises if the materials disagree —
    SEEP2D carries one unsaturated law for the whole problem, so a mixed model
    cannot be exported faithfully."""
    models = {str(m.get("unsat", "lf") or "lf") for m in slope_data["materials"]}
    if len(models) > 1:
        raise ValueError(f"mixed unsaturated models {models}: SEEP2D's iuntyp is global")
    (model,) = models
    if model not in IUNTYP:
        raise ValueError(f"unsaturated model {model!r} has no SEEP2D equivalent")
    return IUNTYP[model]


def build_tri3_seep_data(model_path=MODEL, target_size=None):
    """Load a problem onto a tri3 mesh (SEEP2D = linear triangles)."""
    sd = load_slope_data(model_path)
    polys = get_material_polygons(sd)
    if target_size is None:
        xs = [x for x, _ in sd["ground_surface"].coords]
        target_size = (max(xs) - min(xs)) / 120.0   # same rule as the samples harness
    mesh = build_mesh_from_polygons(polys, target_size, "tri3")
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
    if seep.get("flux_nodal") is not None:
        # Must be permuted with everything else, or the nodal flows land on the wrong nodes.
        out["flux_nodal"] = np.asarray(seep["flux_nodal"], dtype=float)[perm]
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


def write_s2d(seep, gamma_w, path, iuntyp=1, title="xslope export for SEEP2D cross-check",
              flux_edges=None):
    """flux_edges: optional list of (i, j, q) — a boundary EDGE (0-based node pair) carrying
    a uniform normal Darcy velocity q, positive into the domain.

    SEEP2D expresses a Neumann BC as "flowrate cards" read after the element cards, with
    their count in the header's 4th field (nflcd). It does its OWN consistent-load assembly
    (seep2d.f:2023):

        sij   = |x_j - x_i|                       ! edge length
        fx(i) += 0.5 * sij * flrt                 ! q*L/2 at each end
        fx(j) += 0.5 * sij * flrt
        nbc(i) = -1 ;  nbc(j) = -1                ! and sets the flags itself

    which is the same formula XSLOPE uses. So handing SEEP2D the RAW (edge, q) — rather than
    our own assembled nodal vector — makes the comparison a genuine test of one assembly
    against the other. It also picks up SEEP2D's internal `flrt / sck` scaling (sck = max k),
    which a hand-written nbc=-1 node card does NOT get: writing fx directly on the node card
    silently produces a head field ~1/sck too small.
    """
    nodes = np.asarray(seep["nodes"])
    elements = np.asarray(seep["elements"])
    etypes = np.asarray(seep["element_types"])
    emats = np.asarray(seep["element_materials"])       # 1-based
    bc_type = np.asarray(seep["bc_type"])
    bc_val = np.asarray(seep["bc_values"])
    k1 = np.asarray(seep["k1_by_mat"]); k2 = np.asarray(seep["k2_by_mat"])
    ang = np.asarray(seep["angle_by_mat"])
    # SEEP2D's two per-material unsaturated slots, uspar(1..2): the linear front
    # reads them as (kr0, h0), van Genuchten as (alpha, n).
    if iuntyp == 2:
        us1 = np.asarray(seep["vg_a_by_mat"]); us2 = np.asarray(seep["vg_n_by_mat"])
    else:
        us1 = np.asarray(seep["kr0_by_mat"]); us2 = np.asarray(seep["h0_by_mat"])
    nmat = len(k1)
    fedges = list(flux_edges or [])
    with open(path, "w") as f:
        f.write(f"{title}\n")
        # format(4i5, 1x, a4, f10.0, 4x, a1, f10.0, i5) — 4th int is nflcd (flow-card count),
        # trailing field is iuntyp
        f.write(f"{len(nodes):5d}{len(elements):5d}{nmat:5d}{len(fedges):5d} PLNE"
                f"{0.0:10.1f}    F{gamma_w:10.2f}{iuntyp:5d}\n")
        for m in range(nmat):                            # format(i5, 5f15.0)
            f.write(f"{m+1:5d}{k1[m]:15.6f}{k2[m]:15.6f}{ang[m]:15.6f}"
                    f"{us1[m]:15.6f}{us2[m]:15.6f}\n")
        # Node cards carry ONLY the Dirichlet conditions (nbc 1 = fixed head, 2 = exit face).
        # Flux nodes stay nbc = 0 here: the flow cards below set nbc = -1 themselves.
        # Do NOT hand-write nbc = -1 with a nodal load on the node card — that path skips
        # SEEP2D's `flrt / sck` scaling and yields a head field ~1/max(k) too small.
        for i, (x, y) in enumerate(nodes):               # format(i5, i2, i3, 3f15.0)
            nbc = int(bc_type[i])
            fx = float(bc_val[i]) if nbc == 1 else (float(y) if nbc == 2 else 0.0)
            # Not %15.6f: heads and flows can be far below 1e-6 and six decimals would
            # silently truncate them to zero. Fortran's F descriptor accepts E-notation.
            f.write(f"{i+1:5d}{0:2d}{nbc:3d}{float(x):15.6f}{float(y):15.6f}"
                    f"{fx:15.7E}\n")
        for k, (e, et) in enumerate(zip(elements, etypes)):  # format(6i5)
            n1, n2, n3 = (int(v) + 1 for v in e[:3])
            n4 = n3 if et == 3 else int(e[3]) + 1
            f.write(f"{k+1:5d}{n1:5d}{n2:5d}{n3:5d}{n4:5d}{int(emats[k]):5d}\n")
        for (i, j, q) in fedges:
            # format(2i5, f10.0) — the flowrate field is EXACTLY 10 columns. %10.3E is the
            # widest E-form that still fits with a leading minus sign ("-4.400E-06"); a wider
            # one overruns into the node numbers and SEEP2D silently reads garbage.
            f.write(f"{int(i)+1:5d}{int(j)+1:5d}{float(q):10.3E}\n")


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


def release_point(nodes, head, exit_nodes, tol=1e-6):
    """Elevation of the top of the saturated seepage face — where the phreatic
    surface daylights on the exit face. This is the quantity task #30 is about.

    An ACTIVE seepage-face node is one the solver has toggled to p = 0 (it is
    exfiltrating); above the release point the face is unsaturated and p < 0. So
    the release point is the highest exit node with p >= 0 to within tolerance —
    NOT p strictly positive, which would exclude the p == 0 nodes that define it."""
    if len(exit_nodes) == 0:
        return float("nan")
    p = head[exit_nodes] - nodes[exit_nodes, 1]
    wet = exit_nodes[p >= -tol]
    return float(nodes[wet, 1].max()) if len(wet) else float("nan")


def compare(model_path, target_size, max_iter, workdir, label):
    os.makedirs(workdir, exist_ok=True)
    with contextlib.redirect_stdout(io.StringIO()):
        sd, seep = build_tri3_seep_data(model_path, target_size)
    iuntyp = unsat_code(sd)
    seep_r, bw = rcm_renumber(seep)
    if bw >= 700:
        print(f"{label}: SKIP — RCM half-bandwidth {bw} exceeds SEEP2D MXBNDW")
        return

    with contextlib.redirect_stdout(io.StringIO()):
        sol = run_seepage_analysis(dict(seep_r), tol=1e-8, max_iter=max_iter)
    h_xs, q_xs = np.asarray(sol["head"]), sol["flowrate"]

    s2d = os.path.join(workdir, "model.s2d")
    write_s2d(seep_r, sd["gamma_water"], s2d, iuntyp=iuntyp, title=f"xslope {label}")
    h_s2d, fl = parse_seep2d_out(run_seep2d(s2d, workdir), len(seep_r["nodes"]))
    q_s2d = fl[fl > 0].sum()

    nodes = np.asarray(seep_r["nodes"])
    exit_nodes = np.where(np.asarray(seep_r["bc_type"]) == 2)[0]
    ok = np.isfinite(h_s2d)
    dh = h_xs[ok] - h_s2d[ok]
    r_xs, r_s2d = (release_point(nodes, h, exit_nodes) for h in (h_xs, h_s2d))
    law = {1: "linear front", 2: "van Genuchten"}[iuntyp]

    print(f"\n=== {label}  [{law}]  {len(nodes)} nodes, "
          f"{len(exit_nodes)} exit-face nodes, bandwidth {bw} ===")
    dq = 100 * (q_xs - q_s2d) / q_s2d if q_s2d else float("nan")
    print(f"  discharge q   : xslope {q_xs:.5e}   SEEP2D {q_s2d:.5e}   diff {dq:+.2f}%")
    print(f"  nodal head    : max|dh| {np.abs(dh).max():.4f}   RMS "
          f"{np.sqrt((dh**2).mean()):.4f}   (head range {h_xs.min():.2f}-{h_xs.max():.2f})")
    print(f"  release point : xslope y={r_xs:.3f}   SEEP2D y={r_s2d:.3f}   "
          f"diff {r_xs - r_s2d:+.3f}")
    sys.stdout.flush()


def main():
    global SEEP2D_BIN
    ap = argparse.ArgumentParser()
    ap.add_argument("bin", nargs="?", default=SEEP2D_BIN)
    ap.add_argument("--model", default=MODEL)
    ap.add_argument("--target-size", type=float, default=None)
    ap.add_argument("--max-iter", type=int, default=400)
    ap.add_argument("--gw", action="store_true",
                    help="run the whole groundwater free-surface panel (task #30)")
    a = ap.parse_args()
    SEEP2D_BIN = a.bin

    if a.gw:
        for stem, ts, mi in GW_PANEL:
            try:
                compare(f"docs/verification/files/rocscience_gw/{stem}.xlsx", ts, mi,
                        f"/tmp/seep2d_{stem}", stem)
            except Exception as e:                                   # noqa: BLE001
                print(f"\n=== {stem} ===\n  FAILED: {e}")
        return

    compare(a.model, a.target_size, a.max_iter, WORKDIR,
            os.path.basename(a.model))


if __name__ == "__main__":
    main()
