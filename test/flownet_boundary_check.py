"""The flow lines follow the solver's velocity, hold one value along every no-flow
boundary, and span the total flow.

A no-flow boundary is a flow line, so the stream function phi takes one value along
it, and between the two bounding flow lines phi changes by exactly the total flow Q.
Phi's boundary values are built by walking the boundary and adding up the flow
across each edge, read from the adjacent element's velocity. Water crosses the
boundary only at nodes with a specified head, on an active exit face, or with a
specified flux, and the finite element solution spreads each such node's flow over
the edges beside it.

Two rules have failed here. Counting every edge took in the leak across no-flow
edges, where the solution meets the no-flow condition only on average: on
sheetpile_s50 at 4 m, phi along the pile faces spread over 18% of Q, the flow lines
crossed the pile above its tip, and phi's range was 0.74 Q. The fix for that
(2026-10-27, shipped in 1.0.6 to 1.0.9) counted only edges with BOTH ends carrying
flow. Where the phreatic surface reaches a dam's toe with no seepage face, most of
the flow leaves across the dry-slope edge beside the toe node, so that rule dropped
83% of the outflow on the book's earth dam, and the rescale to Q then sent all the
flow out through the tailwater ground. This check then tested only A to C, which
that rule satisfies by construction. D tests the flow lines against the velocity.
An edge now counts when at least one end carries flow, and an edge with one such end
keeps its flow only when it runs the same way as that node's own.

  A. NO-FLOW EDGES HOLD PHI. On sheetpile_s50 (tri3 and tri6, 4 m) and earth_dam1
     (tri3, its committed mesh), no boundary edge with neither end carrying flow
     changes phi by more than 1e-9 Q.
  B. THE RANGE IS THE FLOW. On the same three solves, max(phi) - min(phi) = Q to 1e-9.
  C. A FLUX EDGE STILL CARRIES ITS FLOW. A rainfall flux on the earth dam's crest
     enters through the flux edges: phi changes along them, and the range is still Q.
  D. THE FLOW LINES FOLLOW THE VELOCITY. In each element, grad(phi) is the Darcy
     velocity turned 90 degrees. Summed over the mesh, area-weighted, the mismatch
     |grad(phi) - Jv| is under 10% of |v| on earth_dam_down_seep (the book's earth
     dam, which exits at the toe; its committed mesh), johnson_res (its committed
     mesh), sheetpile_s50 (tri6, 4 m), and earth_dam1. Under the both-ends rule the
     book's dam read 35% and johnson_res 18%; counting every edge, 18% and 5%.
"""
import contextlib
import io
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

FILES = os.path.join(_REPO, "docs", "seep", "files")


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        return fn(*a, **k)


def _solve(name, element_type=None, size=None, flux=None, folder=FILES):
    from xslope.fileio import load_slope_data
    from xslope.mesh import build_mesh_from_polygons, get_material_polygons, import_mesh_from_json
    from xslope.seep import build_seep_data, run_seepage_analysis
    stem = os.path.join(folder, name)
    sd = _quiet(load_slope_data, stem + ".xlsx")
    if size is None:
        mesh = import_mesh_from_json(stem + "_mesh.json")
    else:
        mesh = _quiet(build_mesh_from_polygons, get_material_polygons(sd), size, element_type)
    if flux is not None:
        bcs = dict(sd["seepage_bc"])
        bcs["specified_fluxes"] = list(bcs.get("specified_fluxes") or []) + [flux]
        sd = dict(sd, seepage_bc=bcs)
    seep = _quiet(build_seep_data, mesh, sd)
    # the docs' solve settings (tools/make_seep_sidecars.py): the earth dam closes in 525 sweeps
    return mesh, seep, _quiet(run_seepage_analysis, seep, tol=1e-4, max_iter=1000)


def _carrying(seep, sol):
    # an exit-face node is active when the solve pinned it at its elevation (p = 0)
    bc = np.asarray(seep["bc_type"])
    p = np.asarray(sol["head"]) - np.asarray(seep["nodes"])[:, 1]
    dirichlet = (bc == 1) | ((bc == 2) & (np.abs(p) < 1e-9))
    flux = np.asarray(sol.get("flux_nodal", np.zeros(len(bc)))) != 0.0
    return dirichlet, flux


def _edges(mesh):
    from xslope.seep import _boundary_edge_map
    edges, _ = _boundary_edge_map(mesh["elements"], mesh["element_types"])
    return list(edges)


def _mismatch(seep, sol):
    """Area-weighted |grad(phi) - J v| / |J v| over the mesh, with J v = (-vy, vx).

    Each element's velocity and phi gradient come from its own shape functions: the
    mesh is exploded so no node is shared, and the library's nodal averaging then
    returns each element's own values at its nodes."""
    from xslope.seep import _elem_props, compute_gradient, compute_velocity
    nodes = np.asarray(seep["nodes"], float)
    E = np.asarray(seep["elements"])
    et = np.asarray(seep.get("element_types", np.full(len(E), 3))).astype(int)
    used = np.arange(E.shape[1])[None, :] < et[:, None]
    newE = np.zeros_like(E)
    newE[used] = np.arange(used.sum())
    idx = E[used]
    X = nodes[idx]
    _, _, k1, k2, ang, kr0, h0, vg_a, vg_n, vg_l, model = _elem_props(seep)
    if not sol.get("unconfined", False):
        kr0 = h0 = vg_a = vg_n = vg_l = model = None
    V = compute_velocity(X, newE, np.asarray(sol["head"], float)[idx], k1, k2, ang, kr0, h0, et,
                         vg_a=vg_a, vg_n=vg_n, vg_l=vg_l, model=model)
    G = -compute_gradient(X, newE, np.asarray(sol["phi"], float)[idx], et)
    nc = np.where(np.isin(et, (4, 8, 9)), 4, 3)
    corner = (np.arange(E.shape[1])[None, :] < nc[:, None]).astype(float)
    ve = (V[newE] * corner[..., None]).sum(1) / nc[:, None]
    ge = (G[newE] * corner[..., None]).sum(1) / nc[:, None]
    ring = E[:, :4].copy()
    ring[nc == 3, 3] = ring[nc == 3, 0]
    px, py = nodes[ring, 0], nodes[ring, 1]
    area = 0.5 * np.abs(np.sum(px * np.roll(py, -1, 1) - np.roll(px, -1, 1) * py, 1))
    Jv = np.column_stack([-ve[:, 1], ve[:, 0]])
    s = 1.0 if np.sum(area * np.sum(ge * Jv, 1)) >= 0 else -1.0
    r = np.linalg.norm(ge - s * Jv, axis=1)
    return float(np.sum(area * r) / np.sum(area * np.linalg.norm(Jv, axis=1)))


def run():
    failures = []
    for label, args in [("sheetpile_s50 tri3 4 m", ("xslope_sheetpile_s50", "tri3", 4.0)),
                        ("sheetpile_s50 tri6 4 m", ("xslope_sheetpile_s50", "tri6", 4.0)),
                        ("earth_dam1 committed mesh", ("xslope_earth_dam1",))]:
        mesh, seep, sol = _solve(*args)
        phi, Q = np.asarray(sol["phi"]), float(sol["flowrate"])
        dirichlet, flux = _carrying(seep, sol)
        carries = dirichlet | flux
        steps = [abs(phi[a] - phi[b]) for a, b in _edges(mesh) if not (carries[a] or carries[b])]
        worst = max(steps) / Q if steps else 0.0
        span = (phi.max() - phi.min()) / Q
        print(f"  {label:<28} no-flow step {worst:.1e} Q, range {span:.10f} Q")
        if worst > 1e-9:
            failures.append(f"A. {label}: phi changes by {worst:.3g} Q along a no-flow edge")
        if abs(span - 1.0) > 1e-9:
            failures.append(f"B. {label}: phi spans {span:.6f} Q")

    # C: rain on the crest (x 51 to 59 at elevation 22)
    flux = {"coords": [(51.0, 22.0), (59.0, 22.0)], "flux": 0.5}
    mesh, seep, sol = _solve("xslope_earth_dam1", flux=flux)
    phi, Q = np.asarray(sol["phi"]), float(sol["flowrate"])
    _, fl = _carrying(seep, sol)
    crest = [abs(phi[a] - phi[b]) for a, b in _edges(mesh) if fl[a] and fl[b]]
    span = (phi.max() - phi.min()) / Q
    print(f"  {'earth_dam1 + crest rain':<28} {len(crest)} flux edges, "
          f"phi change along them {sum(crest) / Q:.3f} Q, range {span:.10f} Q")
    if not crest or sum(crest) <= 1e-6 * Q:
        failures.append("C. the crest flux edges carry no flow")
    if abs(span - 1.0) > 1e-9:
        failures.append(f"C. with crest rain, phi spans {span:.6f} Q")

    # D: the flow lines against the velocity
    lem = os.path.join(_REPO, "docs", "lem", "files")
    for label, args in [("earth_dam_down_seep committed mesh", ("xslope_earth_dam_down_seep", None, None, None, lem)),
                        ("johnson_res committed mesh", ("xslope_johnson_res",)),
                        ("sheetpile_s50 tri6 4 m", ("xslope_sheetpile_s50", "tri6", 4.0)),
                        ("earth_dam1 committed mesh", ("xslope_earth_dam1",))]:
        mesh, seep, sol = _solve(*args)
        if not sol.get("converged", True):
            failures.append(f"D. {label}: the seepage solve did not converge")
            continue
        m = _mismatch(seep, sol)
        print(f"  {label:<36} flow lines vs velocity: mismatch {m:.3f}")
        if m > 0.10:
            failures.append(f"D. {label}: the flow lines miss the velocity by {m:.1%}")
    return failures


if __name__ == "__main__":
    f = run()
    print("\n".join(["FAILED:"] + f) if f else "\nAll flow-net boundary checks passed.")
    sys.exit(1 if f else 0)
