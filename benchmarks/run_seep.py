"""Run the FE-seepage benchmark cases and compare to the reference.

SEEP-1: Kozeny homogeneous dam with horizontal toe drain. Compares xslope's
total discharge to the closed form q = k * y0 (y0 = sqrt(d^2+h^2) - d).

Run from the repo root:  PYTHONPATH=. python3 benchmarks/run_seep.py
"""
import math
import numpy as np
from xslope.fileio import load_slope_data
from xslope.mesh import get_material_polygons, build_mesh_from_polygons
from xslope.seep import build_seep_data, run_seepage_analysis
from benchmarks.build_seep import (KOZENY, CONFINED, T_SP, H1_SP, H2_SP,
                                   sheetpile_q_exact)


def run_sheetpile(s_over_T, target_size=1.0, element_type="tri6", tol=1e-10):
    """Partially penetrating sheetpile. Returns (q, n_nodes, head error below tip).

    Second exact check: by antisymmetry the head on the wall plane below the tip
    is exactly (H1+H2)/2.
    """
    sfx = str(int(round(100 * s_over_T)))
    sd = load_slope_data(f"docs/seep/files/xslope_sheetpile_s{sfx}.xlsx")
    polygons = get_material_polygons(sd)
    mesh = build_mesh_from_polygons(polygons, target_size=target_size,
                                    element_type=element_type)
    seep = build_seep_data(mesh, sd, seep_bc=1)
    sol = run_seepage_analysis(seep, tol=tol)
    nodes = np.array(mesh["nodes"])
    # nodes on the wall plane below the tip (x ~ 0, y < T - s)
    tip_y = T_SP - s_over_T * T_SP
    m = (np.abs(nodes[:, 0]) < 1e-6) & (nodes[:, 1] < tip_y - 1e-6)
    h_mid = 0.5 * (H1_SP + H2_SP)
    err = float(np.max(np.abs(sol["head"][m] - h_mid))) if m.any() else float("nan")
    return sol["flowrate"], len(nodes), err


def run_confined(target_size=1.0, element_type="tri6", tol=1e-10):
    """Confined radial flow: returns (q_xslope, n_nodes, max head error vs exact)."""
    sd = load_slope_data("docs/seep/files/xslope_confined_radial.xlsx")
    polygons = get_material_polygons(sd)
    mesh = build_mesh_from_polygons(polygons, target_size=target_size,
                                    element_type=element_type)
    seep = build_seep_data(mesh, sd, seep_bc=1)
    sol = run_seepage_analysis(seep, tol=tol)
    nodes = np.array(mesh["nodes"])
    r = np.hypot(nodes[:, 0], nodes[:, 1])
    h_exact = (CONFINED["h1"] + (CONFINED["h2"] - CONFINED["h1"])
               * np.log(r / CONFINED["r1"]) / math.log(CONFINED["r2"] / CONFINED["r1"]))
    head_err = float(np.max(np.abs(sol["head"] - h_exact)))
    return sol["flowrate"], len(nodes), head_err


def run_kozeny(target_size=2.0, element_type="tri6", tol=1e-8):
    sd = load_slope_data("docs/seep/files/xslope_kozeny_dam.xlsx")
    polygons = get_material_polygons(sd)
    mesh = build_mesh_from_polygons(polygons, target_size=target_size,
                                    element_type=element_type)
    seep = build_seep_data(mesh, sd, seep_bc=1)
    sol = run_seepage_analysis(seep, tol=tol)
    q = sol["flowrate"]
    return q, len(mesh["nodes"])


if __name__ == "__main__":
    import io
    import contextlib

    # SEEP-1 (analytical anchor): confined radial flow -- exact, saturated.
    qc_ref = CONFINED["q"]
    conf = []
    with contextlib.redirect_stdout(io.StringIO()):
        for et, ts in [("tri6", 2.0), ("tri6", 1.0), ("tri3", 0.5)]:
            try:
                q, n, he = run_confined(target_size=ts, element_type=et)
                conf.append((et, ts, n, q, 100.0 * (q - qc_ref) / qc_ref, he))
            except Exception as e:
                conf.append((et, ts, None, None, None, str(e)))
    print(f"\nSEEP-1  Confined radial flow  (analytical anchor, saturated)")
    print(f"  exact  q = k*alpha*(h1-h2)/ln(r2/r1) = {qc_ref:.4f}")
    print(f"  {'mesh':<16}{'nodes':>8}{'xslope q':>12}{'diff %':>9}{'max|dh|':>10}")
    for et, ts, n, q, diff, he in conf:
        if q is None:
            print(f"  {et+' ts='+str(ts):<16} ERROR: {he}")
        else:
            print(f"  {et+' ts='+str(ts):<16}{n:>8}{q:>12.4f}{diff:>9.2f}{he:>10.4f}")

    # SEEP-1c: partially penetrating sheetpile (Pavlovsky exact form factor).
    sp = []
    with contextlib.redirect_stdout(io.StringIO()):
        for sT in [0.5, 0.75]:
            try:
                q, n, he = run_sheetpile(sT, target_size=1.0)
                qe = sheetpile_q_exact(sT)
                sp.append((sT, n, q, qe, 100.0 * (q - qe) / qe, he))
            except Exception as e:
                sp.append((sT, None, None, None, None, str(e)))
    print(f"\nSEEP-1c  Sheetpile cutoff (Pavlovsky exact, confined)")
    print(f"  {'s/T':<8}{'nodes':>8}{'xslope q':>12}{'exact q':>10}{'diff %':>9}"
          f"{'|dh| below tip':>16}")
    for sT, n, q, qe, diff, he in sp:
        if q is None:
            print(f"  {sT:<8} ERROR: {he}")
        else:
            print(f"  {sT:<8}{n:>8}{q:>12.4f}{qe:>10.4f}{diff:>9.2f}{he:>16.4f}")

    # Kozeny dam (free-surface sample, retained for scrutiny -- see samples.md).
    qk_ref = KOZENY["q"]
    koz = []
    with contextlib.redirect_stdout(io.StringIO()):
        for et, ts in [("tri6", 1.5), ("tri6", 0.8)]:
            try:
                q, n = run_kozeny(target_size=ts, element_type=et)
                koz.append((et, ts, n, q, 100.0 * (q - qk_ref) / qk_ref))
            except Exception as e:
                koz.append((et, ts, None, None, str(e)))
    print(f"\nKozeny dam  (free-surface sample; discharge sensitive to drain-start"
          f" tangency)")
    print(f"  basic-parabola  q = k*y0 = {qk_ref:.4f}")
    print(f"  {'mesh':<16}{'nodes':>8}{'xslope q':>12}{'diff %':>9}")
    for et, ts, n, q, diff in koz:
        if q is None:
            print(f"  {et+' ts='+str(ts):<16} ERROR: {diff}")
        else:
            print(f"  {et+' ts='+str(ts):<16}{n:>8}{q:>12.4f}{diff:>9.1f}")
