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

"""Guard: the opt-in matric-suction apparent-cohesion option (Fredlund extended
Mohr-Coulomb) in generate_slices / the LEM solvers.

A simple homogeneous slope with a piezometric line drawn WELL BELOW the ground
surface, so the upper part of the circular failure surface runs above the water
table and carries matric suction. The guard freezes four invariants of the
``suction_phi_b`` / ``suction_cap`` run options:

  1. OFF is bit-identical: ``suction_phi_b=None`` and ``phi_b=0`` both leave FS
     EXACTLY equal to the clamped baseline (c_suction = s*tan(0) = 0.0, and adding
     0.0 is exact in IEEE-754) — the default-off invariance the corpus locks rely on.
  2. phi_b > 0 STRICTLY raises FS: apparent cohesion s*tan(phi_b) on the slices
     above the water table can only add resisting strength.
  3. A finite ``suction_cap`` gives an FS strictly between baseline and uncapped:
     capping the suction s caps the apparent cohesion, so 0 < capped-gain < full-gain.
  4. Hand-check: on a slice whose base is above the piezo line, the stored
     ``c_suction`` equals gamma_w*(y_piezo - y_base) clamped to >= 0, times
     tan(phi_b) — the exact (u_a - u_w) tan(phi_b) term with u_a = 0.

Deterministic, no network; builds its slope on docs/lem/files/xslope_acads_simple.xlsx.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/suction_guard.py
Exits non-zero on any failure.
"""
import math
import os
import sys
import warnings

warnings.filterwarnings("ignore")

from shapely.geometry import Polygon

from xslope.fileio import load_slope_data, build_ground_surface_from_polygons
from xslope.slice import generate_slices
from xslope.solve import oms, bishop, janbu, spencer

GW = 9.81
PHI_B = 20.0                 # suction-strength angle for the ON runs
SUCTION_CAP = 15.0           # stress units; below the peak suction so the cap bites
MAT_NAME = "suction soil"
_BASE_XLSX = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "docs", "lem", "files", "xslope_acads_simple.xlsx")

# Homogeneous slope descending left-to-right: crest bench at y=30 (x 0..20),
# a 28-ish deg face down to a toe bench at y=10, extended to x=120. The circle
# (center high, base near y=4) cuts up through the crest region, so the upper
# slices' bases sit above the piezo line (y=6) and carry suction.
_SOIL_POLY = [(0, -10), (120, -10), (120, 10), (60, 10), (20, 30), (0, 30)]
_SURFACE = [(0, 30), (20, 30), (60, 10), (120, 10)]
_PIEZO = [(-50.0, 6.0), (400.0, 6.0)]
_CIRCLE = {"Xo": 50.0, "Yo": 42.0, "Depth": 38.0, "R": 38.0}
_NUM_SLICES = 40


def _base():
    return load_slope_data(_BASE_XLSX)


def _slope_data(base):
    m = dict(base["materials"][0])
    m.update(name=MAT_NAME, c=5.0, phi=30.0, gamma=20.0, gamma_sat=20.0,
             option="mc", u="piezo", ru=0.0)
    polys = [{"polygon": Polygon(_SOIL_POLY), "mat_id": 0}]
    gs, dom = build_ground_surface_from_polygons(polys)
    sd = dict(base)
    sd["materials"] = [m]
    sd["polygons"] = polys
    sd["ground_surface"] = gs
    sd["domain_polygon"] = dom
    sd["piezo_line"] = _PIEZO
    sd["piezo_line2"] = []
    sd["piezo_phreatic"] = False
    sd["dloads"] = []
    sd["dloads2"] = []
    sd["circular"] = True
    sd["circles"] = [_CIRCLE]
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


def _slices(sd, **kw):
    ok, res = generate_slices(sd, circle=sd["circles"][0], num_slices=_NUM_SLICES,
                              debug=False, **kw)
    if not ok:
        raise RuntimeError(f"generate_slices failed: {res}")
    return res[0]


def _fs(slice_df):
    out = {}
    for name, fn in (("oms", oms), ("bishop", bishop), ("janbu", janbu), ("spencer", spencer)):
        s, r = fn(slice_df.copy())
        out[name] = r["FS"] if s else None
    return out


def check():
    """Run the guard; return a list of failure strings (empty on success)."""
    failures = []
    sd = _slope_data(_base())

    df_off = _slices(sd)                                            # default (None)
    df_none = _slices(sd, suction_phi_b=None)                       # explicit None
    df_zero = _slices(sd, suction_phi_b={MAT_NAME: 0.0})            # phi_b = 0
    df_on = _slices(sd, suction_phi_b={MAT_NAME: PHI_B})            # phi_b = 20
    df_cap = _slices(sd, suction_phi_b={MAT_NAME: PHI_B}, suction_cap=SUCTION_CAP)

    # The problem itself must exercise suction, else the guard proves nothing.
    n_susc = int((df_on["c_suction"] > 0).sum())
    if n_susc == 0:
        failures.append("setup: no slice carries matric suction (piezo line above "
                        "the whole failure surface) — guard would be vacuous")
        return failures

    fs_off = _fs(df_off)
    fs_none = _fs(df_none)
    fs_zero = _fs(df_zero)
    fs_on = _fs(df_on)
    fs_cap = _fs(df_cap)

    for method in ("oms", "bishop", "janbu", "spencer"):
        b = fs_off[method]
        if b is None:
            failures.append(f"{method}: baseline solve failed")
            continue

        # (1) OFF is bit-identical (None and phi_b=0 both).
        if fs_none[method] != b:
            failures.append(f"{method}: suction_phi_b=None changed FS "
                            f"({fs_none[method]!r} != baseline {b!r})")
        if fs_zero[method] != b:
            failures.append(f"{method}: phi_b=0 changed FS "
                            f"({fs_zero[method]!r} != baseline {b!r})")

        # (2) phi_b > 0 strictly raises FS.
        on = fs_on[method]
        if not (on > b):
            failures.append(f"{method}: phi_b={PHI_B} did not raise FS "
                            f"(on={on:.6f} <= baseline={b:.6f})")

        # (3) capped is strictly between baseline and uncapped.
        cap = fs_cap[method]
        if not (b < cap < on):
            failures.append(f"{method}: capped FS not strictly between baseline and "
                            f"uncapped (baseline={b:.6f}, cap={cap:.6f}, on={on:.6f})")

    # (4) Hand-check c_suction on the suction slice with the largest suction.
    row = df_on.loc[df_on["c_suction"].idxmax()]
    y_base = float(row["y_cb"])
    y_piezo = float(row["piezo_y"])
    s_expected = max(0.0, GW * (y_piezo - y_base) * -1.0)  # gamma_w*(y_base - y_piezo)
    c_expected = s_expected * math.tan(math.radians(PHI_B))
    c_got = float(row["c_suction"])
    if abs(c_got - c_expected) > 1e-9:
        failures.append(f"hand-check: c_suction={c_got:.6f} != s*tan(phi_b)="
                        f"{c_expected:.6f} (s={s_expected:.4f}, y_base={y_base:.3f}, "
                        f"y_piezo={y_piezo:.3f})")

    # (4b) The same slice under the cap must equal min(s, cap)*tan(phi_b).
    row_c = df_cap.iloc[int(row["slice #"]) - 1]
    c_cap_expected = min(s_expected, SUCTION_CAP) * math.tan(math.radians(PHI_B))
    if abs(float(row_c["c_suction"]) - c_cap_expected) > 1e-9:
        failures.append(f"hand-check(cap): c_suction={float(row_c['c_suction']):.6f} "
                        f"!= min(s,cap)*tan(phi_b)={c_cap_expected:.6f}")

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
    print("OK: matric-suction apparent cohesion is off-by-default bit-identical "
          "(None and phi_b=0), strictly raises FS for phi_b>0, respects suction_cap "
          "(baseline < capped < uncapped), and c_suction = max(0, gamma_w*(y_base - "
          "y_piezo))*tan(phi_b) exactly — OMS, Bishop, Janbu, Spencer.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
