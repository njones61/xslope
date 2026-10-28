"""The stream function is constant along every no-flow boundary, and spans the total flow.

A no-flow boundary is a flow line, so the stream function phi takes one value along
it, and between the two bounding flow lines phi changes by exactly the total flow Q.
Phi's boundary values are built by walking the boundary and adding up the flow
across each edge. That flow was read from the adjacent element's velocity on every
edge, no-flow edges included, where the finite element solution leaks: it meets a
no-flow condition only on average. Near a singularity the leak was large. On
sheetpile_s50 meshed at 4 m, phi along the pile faces spread over 18% of Q, the flow
lines crossed the pile above its tip, and phi's range was 0.74 Q. Only edges between
two fixed-head or active exit-face nodes, or two specified-flux nodes, now carry flow.

  A. NO-FLOW EDGES HOLD PHI. On sheetpile_s50 (tri3 and tri6, 4 m) and earth_dam1
     (tri3, its committed mesh), no boundary edge outside the flow-carrying set
     changes phi by more than 1e-9 Q.
  B. THE RANGE IS THE FLOW. On the same three solves, max(phi) - min(phi) = Q to 1e-9.
  C. A FLUX EDGE STILL CARRIES ITS FLOW. A rainfall flux on the earth dam's crest
     enters through the flux edges: phi changes along them, and the range is still Q.
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


def _solve(name, element_type=None, size=None, flux=None):
    from xslope.fileio import load_slope_data
    from xslope.mesh import build_mesh_from_polygons, get_material_polygons, import_mesh_from_json
    from xslope.seep import build_seep_data, run_seepage_analysis
    stem = os.path.join(FILES, name)
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
    return mesh, seep, _quiet(run_seepage_analysis, seep)


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


def run():
    failures = []
    for label, args in [("sheetpile_s50 tri3 4 m", ("xslope_sheetpile_s50", "tri3", 4.0)),
                        ("sheetpile_s50 tri6 4 m", ("xslope_sheetpile_s50", "tri6", 4.0)),
                        ("earth_dam1 committed mesh", ("xslope_earth_dam1",))]:
        mesh, seep, sol = _solve(*args)
        phi, Q = np.asarray(sol["phi"]), float(sol["flowrate"])
        dirichlet, flux = _carrying(seep, sol)
        steps = [abs(phi[a] - phi[b]) for a, b in _edges(mesh)
                 if not ((dirichlet[a] and dirichlet[b]) or (flux[a] and flux[b]))]
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
    return failures


if __name__ == "__main__":
    f = run()
    print("\n".join(["FAILED:"] + f) if f else "\nAll flow-net boundary checks passed.")
    sys.exit(1 if f else 0)
