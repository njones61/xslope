# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Guard: the dry-buoyant still-water oracle for xslope's ponded-water dload.

A FULLY submerged slope under STILL water must read IDENTICALLY to the same slope
analyzed DRY with buoyant unit weight gamma' = gamma_sat - gamma_w (no water, no
distributed load, no piezo). This is the classical exact effective-stress
equivalence, and it is the oracle that validates the perpendicular-distributed-load
water convention that every vendor importer now converts ponded water INTO:

    (a) submerged: gamma_sat, a hydrostatic perpendicular dload on the submerged
        face, and piezometric pore pressure u on the base;
    (c) dry buoyant: gamma' = gamma_sat - gamma_w, u = 'none', no dload, no piezo.

If (a) did not reduce to (c), the reservoir loads the Slide2 and RS2 importers now
synthesize would be putting the wrong statics on the slope. The diagnostic that
motivated the import fix measured the residual at < 0.006 FS; this guard freezes it.

Asserts |FS_submerged - FS_buoyant| < 0.01 for Bishop AND Spencer, across three
strength combinations, on BOTH a right-descending and a left-descending slope.
Deterministic, no network; builds its slopes on docs/lem/files/xslope_acads_simple.xlsx.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/submerged_oracle_guard.py
Exits non-zero on any failure.
"""
import os
import sys
import warnings

warnings.filterwarnings("ignore")

from shapely.geometry import Polygon

from xslope.fileio import load_slope_data, build_ground_surface_from_polygons
from xslope.slice import generate_slices
from xslope.solve import bishop, spencer

GW = 9.81
TOL = 0.01
_BASE_XLSX = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "docs", "lem", "files", "xslope_acads_simple.xlsx")

# Two fully submerged slopes (reservoir el. 40, crest 30), one descending each way.
_CASES = [
    ("CASE R descend-right",
     [(0, -10), (120, -10), (120, 10), (60, 10), (20, 30), (0, 30)],
     [(0, 30), (20, 30), (60, 10), (120, 10)],
     {"Xo": 50.0, "Yo": 42.0, "Depth": 38.0, "R": 38.0}),
    ("CASE L descend-left",
     [(0, -10), (120, -10), (120, 30), (100, 30), (60, 10), (0, 10)],
     [(0, 10), (60, 10), (100, 30), (120, 30)],
     {"Xo": 70.0, "Yo": 42.0, "Depth": 38.0, "R": 38.0}),
]
_STRENGTHS = [(5.0, 30.0), (0.0, 35.0), (15.0, 20.0)]
_RES = 40.0


def _base():
    return load_slope_data(_BASE_XLSX)


def _make_material(base, name, c, phi, gamma, u):
    m = dict(base["materials"][0])
    m.update(name=name, c=c, phi=phi, gamma=gamma, gamma_sat=gamma,
             option="mc", u=u, ru=0.0)
    return m


def _slope_data(base, materials, polygons_coords, piezo_line, dloads, circle):
    sd = dict(base)
    sd["materials"] = materials
    polys = [{"polygon": Polygon(pts), "mat_id": mid} for mid, pts in polygons_coords]
    sd["polygons"] = polys
    gs, dom = build_ground_surface_from_polygons(polys)
    sd["ground_surface"] = gs
    sd["domain_polygon"] = dom
    sd["piezo_line"] = piezo_line
    sd["piezo_line2"] = []
    sd["piezo_phreatic"] = False
    sd["dloads"] = dloads
    sd["dloads2"] = []
    # The oracle STATES its own water load: the whole point is to price a known
    # hydrostatic block against the buoyant-unit-weight equivalent, so the block
    # above is the quantity under test. The base file is a v22 model in automatic
    # mode, and inheriting that mode would have the engine derive the same
    # reservoir from the piezometric line and apply it on top -- the submerged
    # cases read FS ~ 161 instead of ~ 2.05, which is the double count this
    # guard would otherwise be measuring. Set directly, not derived: a model
    # built in memory declares its own mode.
    sd["water_loads"] = "manual"
    sd["circular"] = True
    sd["circles"] = [circle]
    sd["non_circ"] = []
    sd["reinforce_lines"] = []
    sd["reinforcement_lines"] = []
    sd["pile_lines"] = []
    sd["line_loads"] = []
    sd["k_seismic"] = 0.0
    sd["tcrack_depth"] = 0.0
    sd["tcrack_water"] = 0.0
    sd["mesh"] = None
    sd.pop("_water_table_profile", None)
    return sd


def _hydro_dload(surface_pts, res_level):
    """Hydrostatic normal intensity gamma_w*(res - y) at each surface vertex,
    clamped to the submerged span — the perpendicular-dload water convention."""
    block = [{"X": float(x), "Y": float(y), "Normal": float(GW * max(0.0, res_level - y))}
             for (x, y) in surface_pts]
    return [block]


def _fs(sd, num_slices=40):
    ok, res = generate_slices(sd, circle=sd["circles"][0], num_slices=num_slices)
    if not ok:
        return {"bishop": f"slices FAIL: {res}", "spencer": None}
    slice_df, _ = res
    out = {}
    for name, fn in (("bishop", bishop), ("spencer", spencer)):
        try:
            s, r = fn(slice_df.copy())
            out[name] = r["FS"] if s else f"FAIL:{r}"
        except Exception as e:                    # pragma: no cover - defensive
            out[name] = f"EXC:{e}"
    return out


def check():
    """Run the oracle; return a list of failure strings (empty on success)."""
    failures = []
    base = _base()
    for title, soil, surf, circ in _CASES:
        for c, phi in _STRENGTHS:
            mat_a = _make_material(base, "soil", c, phi, gamma=20.0, u="piezo")
            piezo = [(-50, _RES), (400, _RES)]
            sd_a = _slope_data(base, [mat_a], [(0, soil)], piezo,
                               _hydro_dload(surf, _RES), circ)
            fa = _fs(sd_a)

            mat_c = _make_material(base, "soilb", c, phi, gamma=20.0 - GW, u="none")
            sd_c = _slope_data(base, [mat_c], [(0, soil)], [], [], circ)
            fc = _fs(sd_c)

            for method in ("bishop", "spencer"):
                a, cc = fa[method], fc[method]
                if not isinstance(a, float) or not isinstance(cc, float):
                    failures.append(f"{title} c={c} phi={phi} {method}: "
                                    f"submerged={a!r} buoyant={cc!r}")
                    continue
                if abs(a - cc) >= TOL:
                    failures.append(f"{title} c={c} phi={phi} {method}: "
                                    f"|{a:.4f}-{cc:.4f}| = {abs(a - cc):.4f} >= {TOL}")
    return failures


def main():
    if not os.path.exists(_BASE_XLSX):
        print(f"SKIP: base template not found ({_BASE_XLSX})")
        return 0
    failures = check()
    if failures:
        print("FAILED:")
        for f in failures:
            print("  -", f)
        return 1
    print("OK: a fully submerged still-water slope reads identically to the same slope "
          "run dry with gamma' = gamma_sat - gamma_w (|dFS| < 0.01, Bishop + Spencer, "
          "both facings, 3 strengths) — the perpendicular-dload water convention is exact.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
