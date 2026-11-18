"""The vertical seismic coefficient k_v, in every LEM method and in the FEM loads.

k_v is signed, positive downward, and adds k_v·W to each slice's weight through
its center of gravity. Three checks:

* **Every LEM method, on a circle.** On a dry slope with no cohesion and no
  external loads, every force in every method's equilibrium scales with the body
  forces, so a slope under (k, k_v) has exactly the factor of safety of the same
  slope under k/(1 + k_v) with no vertical force. The planar wedge on the
  verification page checks the force-equilibrium methods against a closed form;
  this checks the moment methods (Ordinary, Bishop, Janbu) and the rest on a
  circle, where there is none.
* **The cohesion term does not scale.** With cohesion the identity must fail, or
  the vertical force is being applied to the strength as well as the weight.
* **The FEM load vector.** The vertical loads are (1 + k_v) times the weight, and
  the horizontal ones are unchanged.

Runs in a few seconds.
"""
import os
import sys
import warnings

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, ROOT)

from xslope import solve  # noqa: E402
from xslope.fileio import load_slope_data  # noqa: E402
from xslope.slice import generate_slices  # noqa: E402

MODEL = os.path.join(ROOT, "docs", "lem", "files", "xslope_acads_simple.xlsx")
METHODS = ("oms", "bishop", "janbu", "corps", "lowe", "spencer", "mprice")
K, KV = 0.15, -0.075


def _dry(cohesion):
    sd = load_slope_data(MODEL)
    for m in sd["materials"]:
        m["c"] = cohesion
        m["phi"] = max(float(m.get("phi") or 0.0), 30.0)
        m["u"] = "none"
    sd["piezo_line"] = []
    sd["dloads"] = []
    sd["tcrack_depth"] = 0.0
    sd["tcrack_water"] = 0.0
    return sd


def _fs(sd, k, kv, method):
    sd = dict(sd, k_seismic=k, k_vertical=kv)
    ok, out = generate_slices(sd, circle=sd["circles"][0], num_slices=40)
    if not ok:
        raise RuntimeError(out)
    df, _ = out
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        good, res = getattr(solve, method)(df)
    return float(res["FS"]) if good else None


def check_every_method_scales_on_a_circle():
    fails = []
    sd = _dry(0.0)
    for method in METHODS:
        a = _fs(sd, K, KV, method)
        b = _fs(sd, K / (1.0 + KV), 0.0, method)
        if a is None or b is None:
            fails.append(f"{method} found no solution on the dry frictional circle")
        elif abs(a - b) > 1e-6 * max(1.0, abs(b)):
            fails.append(f"{method}: F(k={K}, k_v={KV}) = {a:.6f}, but "
                         f"F(k={K / (1 + KV):.4f}, k_v=0) = {b:.6f}")
    return fails


def check_cohesion_does_not_scale():
    sd = _dry(10.0)
    a = _fs(sd, K, KV, "bishop")
    b = _fs(sd, K / (1.0 + KV), 0.0, "bishop")
    if a is None or b is None:
        return ["Bishop found no solution on the c-phi circle"]
    if abs(a - b) < 1e-4:
        return [f"with cohesion, F(k, k_v) = {a:.6f} still equals F(k/(1 + k_v), 0) "
                f"= {b:.6f}: the vertical force is scaling the strength too"]
    return []


def check_fem_loads():
    from xslope.fem import build_fem_data
    from xslope.mesh import build_mesh_from_polygons, get_material_polygons

    sd = load_slope_data(MODEL)
    mesh = build_mesh_from_polygons(get_material_polygons(sd), target_size=2.0,
                                    element_type="tri6", debug=False)
    loads = {}
    for kv in (0.0, KV):
        fd = build_fem_data(dict(sd, k_seismic=K, k_vertical=kv), mesh)
        from xslope.fem import build_gravity_loads
        loads[kv] = build_gravity_loads(fd["nodes"], fd["elements"],
                                        fd["element_types"], fd["element_materials"],
                                        fd["gamma_by_mat"], fd["k_seismic"], fem_data=fd)
    n = len(mesh["nodes"])
    y = 2 * np.arange(n) + 1
    x = 2 * np.arange(n)
    fails = []
    if not np.allclose(loads[KV][y], (1.0 + KV) * loads[0.0][y]):
        fails.append("the FEM vertical loads are not (1 + k_v) times the weight")
    if not np.allclose(loads[KV][x], loads[0.0][x]):
        fails.append("k_v changed the FEM horizontal loads")
    return fails


def run():
    """Failures as a list, for run_tests.py."""
    fails = []
    for check in (check_every_method_scales_on_a_circle, check_cohesion_does_not_scale,
                  check_fem_loads):
        try:
            fails += check()
        except Exception as exc:          # a crash is a failure, not a pass
            fails.append(f"{check.__name__} raised {exc!r}")
    return fails


if __name__ == "__main__":
    failures = run()
    if failures:
        print("FAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("Vertical seismic coefficient checks passed.")
