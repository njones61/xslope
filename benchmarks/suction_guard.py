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

import numpy as np

warnings.filterwarnings("ignore")

from shapely.geometry import Polygon

from xslope.fileio import load_slope_data, build_ground_surface_from_polygons
from xslope.mesh import interpolate_at_point
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


# === SEEP branch: matric suction delivered from a seepage u-field ===============
# The piezo guard above proves the strength side for a hand-drawn piezo line. This
# second guard proves the SAME apparent cohesion is delivered when the base pore
# pressure comes from a SEEPAGE solution (material u='seep'), which reads the mesh
# via interpolate_at_point. That read must hand the SIGNED field to the suction
# option (negative u = matric suction above the water table); if it clamps at 0
# first, the seep u-source can never raise suction and c_suction stays 0 for every
# slice. Regression target for exactly that latent clamp.
#
# The fixture puts a HYDROSTATIC field u(x, y) = gamma_w*(y_wt - y) onto a synthetic
# triangular mesh. On linear triangles that field interpolates EXACTLY, so the seep
# u-source must reproduce the piezo u-source (piezo line at y_wt) slice-for-slice:
# same c_suction, same FS. gamma_sat is left unset so the gamma/gamma_sat weight
# split (a different mesh consumer) stays out of the comparison.
Y_WT_SEEP = 6.0              # water table elevation for the synthetic seep field
                            # (matched to the piezo guard's _PIEZO y so the two
                            # u-sources are directly comparable)


def _synthetic_hydrostatic_mesh(y_wt, gamma_w, x0=-5.0, x1=125.0, y0=-5.0, y1=46.0,
                                 nx=54, ny=22):
    """A structured triangular mesh over [x0,x1] x [y0,y1] carrying the exact
    hydrostatic field u = gamma_w*(y_wt - y) at every node (negative above y_wt).
    Returns (mesh_dict, seep_u) in the format generate_slices expects."""
    xs = np.linspace(x0, x1, nx)
    ys = np.linspace(y0, y1, ny)
    XX, YY = np.meshgrid(xs, ys)                      # (ny, nx)
    nodes = np.column_stack([XX.ravel(), YY.ravel()])  # row-major: idx = j*nx + i
    elements = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            a = j * nx + i
            b = j * nx + (i + 1)
            c = (j + 1) * nx + i
            d = (j + 1) * nx + (i + 1)
            elements.append([a, b, d])                 # lower triangle
            elements.append([a, d, c])                 # upper triangle
    n_elem = len(elements)
    elem_arr = np.zeros((n_elem, 9), dtype=int)
    elem_arr[:, :3] = np.asarray(elements, dtype=int)
    element_types = np.full(n_elem, 3, dtype=int)
    mesh = {"nodes": nodes, "elements": elem_arr, "element_types": element_types}
    seep_u = gamma_w * (y_wt - nodes[:, 1])            # signed: < 0 above y_wt
    return mesh, seep_u


def _seep_slope_data(base):
    """The piezo guard's slope, but the u-source is the synthetic seepage mesh."""
    m = dict(base["materials"][0])
    m.update(name=MAT_NAME, c=5.0, phi=30.0, gamma=20.0, gamma_sat=None,
             option="mc", u="seep", ru=0.0)
    polys = [{"polygon": Polygon(_SOIL_POLY), "mat_id": 0}]
    gs, dom = build_ground_surface_from_polygons(polys)
    mesh, seep_u = _synthetic_hydrostatic_mesh(Y_WT_SEEP, GW)
    sd = dict(base)
    sd["materials"] = [m]
    sd["polygons"] = polys
    sd["ground_surface"] = gs
    sd["domain_polygon"] = dom
    sd["piezo_line"] = []
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
    sd["mesh"] = mesh
    sd["seep_u"] = seep_u
    sd.pop("_water_table_profile", None)
    return sd


def check_seep():
    """Guard the seep-u-source -> matric-suction delivery. Returns failures list."""
    failures = []
    base = _base()
    sd = _seep_slope_data(base)

    df_off = _slices(sd)                                     # default (None)
    df_on = _slices(sd, suction_phi_b={MAT_NAME: PHI_B})     # phi_b = 20

    # The seep u-source must actually raise suction on the slices above y_wt --
    # this is exactly what the pre-fix clamp defeated (c_suction == 0 everywhere).
    n_susc = int((df_on["c_suction"] > 0).sum())
    if n_susc == 0:
        failures.append("seep: no slice carries matric suction from the seep u-source "
                        "(the interpolated field was clamped at 0 before reaching the "
                        "suction option) -- the seep->suction delivery is broken")
        return failures

    # OFF is bit-identical for the seep source too.
    if not (df_off["c_suction"] == 0.0).all():
        failures.append("seep: suction_phi_b=None left a nonzero c_suction")
    fs_off = _fs(df_off)
    fs_on = _fs(df_on)
    for method in ("oms", "bishop", "janbu", "spencer"):
        if fs_off[method] is None or fs_on[method] is None:
            failures.append(f"seep/{method}: solve failed")
            continue
        if not (fs_on[method] > fs_off[method]):
            failures.append(f"seep/{method}: phi_b did not raise FS from the seep "
                            f"u-source (on={fs_on[method]:.6f} <= off={fs_off[method]:.6f})")

    # Hand-check on the largest-suction slice: the field is exactly hydrostatic, so
    # the SIGNED interpolation equals gamma_w*(y_wt - y_base), and c_suction must be
    # max(0, gamma_w*(y_base - y_wt)) * tan(phi_b). Also confirm the CLAMPED read
    # (the old default) would have returned 0 there -- i.e. the sign is load-bearing.
    row = df_on.loc[df_on["c_suction"].idxmax()]
    x_c, y_base = float(row["x_c"]), float(row["y_cb"])
    mesh, seep_u = _synthetic_hydrostatic_mesh(Y_WT_SEEP, GW)
    u_signed = interpolate_at_point(mesh["nodes"], mesh["elements"],
                                    mesh["element_types"], seep_u, (x_c, y_base),
                                    signed=True)
    u_clamped = interpolate_at_point(mesh["nodes"], mesh["elements"],
                                     mesh["element_types"], seep_u, (x_c, y_base))
    s_expected = max(0.0, GW * (y_base - Y_WT_SEEP))
    c_expected = s_expected * math.tan(math.radians(PHI_B))
    c_got = float(row["c_suction"])
    if abs(u_signed - GW * (Y_WT_SEEP - y_base)) > 1e-6:
        failures.append(f"seep: signed interpolation {u_signed:.6f} != hydrostatic "
                        f"{GW*(Y_WT_SEEP - y_base):.6f} at the check slice")
    if u_clamped != 0.0:
        failures.append(f"seep: clamped read at the suction slice was {u_clamped:.6f}, "
                        f"expected 0 -- the sign carries the suction, so signed=True "
                        f"is what makes c_suction nonzero")
    if abs(c_got - c_expected) > 1e-6:
        failures.append(f"seep hand-check: c_suction={c_got:.6f} != "
                        f"max(0, gamma_w*(y_base - y_wt))*tan(phi_b)={c_expected:.6f} "
                        f"(y_base={y_base:.3f}, y_wt={Y_WT_SEEP})")

    # Full per-slice check against the exact hydrostatic field. The single field
    # feeds BOTH terms and both must be right on every slice:
    #   effective-normal u = max(0, gamma_w*(y_wt - y_base))  (clamped, saturated only)
    #   c_suction         = max(0, gamma_w*(y_base - y_wt)) * tan(phi_b)  (signed)
    tanb = math.tan(math.radians(PHI_B))
    yb = df_on["y_cb"].to_numpy()
    u_exp = np.maximum(0.0, GW * (Y_WT_SEEP - yb))
    cs_exp = np.maximum(0.0, GW * (yb - Y_WT_SEEP)) * tanb
    du_all = float(np.max(np.abs(df_on["u"].to_numpy() - u_exp)))
    dcs_all = float(np.max(np.abs(df_on["c_suction"].to_numpy() - cs_exp)))
    if du_all > 1e-6:
        failures.append(f"seep per-slice: effective-normal u off by {du_all:.2e} from "
                        f"max(0, gamma_w*(y_wt - y_base)) -- the caller's own clamp is wrong")
    if dcs_all > 1e-6:
        failures.append(f"seep per-slice: c_suction off by {dcs_all:.2e} from "
                        f"max(0, gamma_w*(y_base - y_wt))*tan(phi_b) on some slice")

    return failures


# === v17 auto-wiring branch: file-carried phi_b/s_cap == explicit kwargs =========
# load_slope_data (v17) carries phi_b/s_cap on each material dict. generate_slices
# auto-wires them into the suction_phi_b / suction_cap options when the caller passes
# no explicit kwarg (t_cut override semantics: an explicit kwarg wins). This guard
# freezes that equivalence -- the whole point of the template column is that opening
# a file with phi_b=20 must give EXACTLY the answer the corpus locks got by passing
# suction_phi_b={name:20} as a tag -- plus the default-off and override paths, and the
# per-material suction_cap DICT path (the template is per-material; the old kwarg was a
# single scalar, so the dict path is new and must match the scalar for one material).
def check_autowire():
    """Return failures for the v17 file-carried auto-wiring equivalence."""
    failures = []
    base = _base()
    sd = _slope_data(base)   # material carries phi_b=None, s_cap=None (loader default)

    # (A) Default path (no phi_b in the file) is off-by-default bit-identical: no
    # kwargs derives an empty dict -> None, so every c_suction is exactly 0.0.
    df_default = _slices(sd)
    if not (df_default["c_suction"] == 0.0).all():
        failures.append("autowire: a material with no phi_b produced nonzero c_suction "
                        "(default path is not bit-identical to pre-v17)")

    # (B) File-carried == kwarg-driven, bit-for-bit. Put phi_b/s_cap on the material
    # dict (the load_slope_data path) and pass NO kwargs; compare to leaving them off
    # the dict and passing the identical values as explicit kwargs.
    sd_file = _slope_data(base)
    sd_file["materials"][0]["phi_b"] = PHI_B
    sd_file["materials"][0]["s_cap"] = SUCTION_CAP
    df_file = _slices(sd_file)                                    # auto-wired from file
    df_kwarg = _slices(sd, suction_phi_b={MAT_NAME: PHI_B}, suction_cap=SUCTION_CAP)
    dcs = float((df_file["c_suction"] - df_kwarg["c_suction"]).abs().max())
    if dcs != 0.0:
        failures.append(f"autowire: file-carried c_suction != kwarg-driven "
                        f"(max |delta| = {dcs:.3e}, must be exactly 0)")
    if _fs(df_file) != _fs(df_kwarg):
        failures.append(f"autowire: file-carried FS != kwarg-driven FS "
                        f"({_fs(df_file)} vs {_fs(df_kwarg)})")

    # (C) An explicit kwarg OVERRIDES the file (t_cut semantics): a file phi_b=PHI_B
    # with an explicit empty suction_phi_b forces suction off.
    df_override = _slices(sd_file, suction_phi_b={})
    if not (df_override["c_suction"] == 0.0).all():
        failures.append("autowire: an explicit empty suction_phi_b did not override "
                        "the file's phi_b (kwarg must win)")

    # (D) The per-material suction_cap DICT path equals the scalar for a single
    # material -- the minimal extension that lets the per-material template column work.
    df_cap_scalar = _slices(sd, suction_phi_b={MAT_NAME: PHI_B}, suction_cap=SUCTION_CAP)
    df_cap_dict = _slices(sd, suction_phi_b={MAT_NAME: PHI_B},
                          suction_cap={MAT_NAME: SUCTION_CAP})
    dcap = float((df_cap_scalar["c_suction"] - df_cap_dict["c_suction"]).abs().max())
    if dcap != 0.0:
        failures.append(f"autowire: per-material suction_cap dict != scalar cap "
                        f"(max |delta| = {dcap:.3e}, must be exactly 0)")
    # A material ABSENT from the cap dict is uncapped (matches suction_cap=None).
    df_cap_none = _slices(sd, suction_phi_b={MAT_NAME: PHI_B}, suction_cap={})
    df_uncapped = _slices(sd, suction_phi_b={MAT_NAME: PHI_B})
    if float((df_cap_none["c_suction"] - df_uncapped["c_suction"]).abs().max()) != 0.0:
        failures.append("autowire: a material absent from the suction_cap dict was not "
                        "left uncapped")

    return failures


def main():
    if not os.path.exists(_BASE_XLSX):
        print(f"SKIP: base template not found ({_BASE_XLSX})")
        return 0
    failures = check()
    failures += check_seep()
    failures += check_autowire()
    if failures:
        print("FAILED:")
        for f in failures:
            print("  -", f)
        return 1
    print("OK: matric-suction apparent cohesion is off-by-default bit-identical "
          "(None and phi_b=0), strictly raises FS for phi_b>0, respects suction_cap "
          "(baseline < capped < uncapped), and c_suction = max(0, gamma_w*(y_base - "
          "y_piezo))*tan(phi_b) exactly — OMS, Bishop, Janbu, Spencer.\n"
          "OK: the SEEP u-source delivers the same suction — a hydrostatic seepage "
          "field raises c_suction on the slices above the water table (signed read; "
          "the clamped read returns 0 there), matches max(0, gamma_w*(y_base - y_wt))"
          "*tan(phi_b) on every slice while the effective-normal u keeps its own "
          "clamp, and is off-by-default bit-identical.")
    print("OK: v17 file-carried phi_b/s_cap auto-wire to EXACTLY the explicit-kwarg "
          "answer (c_suction and FS bit-for-bit); the default path (no phi_b) stays "
          "off-by-default bit-identical; an explicit kwarg overrides the file; and the "
          "per-material suction_cap dict matches the scalar cap.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
