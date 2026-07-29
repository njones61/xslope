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

"""Guard: the two ways pore pressure used to reach a solver as a silent zero.

Both defects shared a signature — the model asks for pore pressure, the code
cannot produce it, and u = 0 is used instead. Zero pore pressure below a water
table over-predicts the factor of safety, so neither may ever be silent.

  A. An unrecognized ``u`` option on the mat sheet. A typo ('peizo') fell through
     every branch to u = 0. ``load_slope_data`` now rejects anything outside
     none / piezo / seep / ru, naming the material and the bad string, and the
     slicer's own else-branch backs it up if the two lists ever drift.

  B. A point that samples a piezometric line from OUTSIDE the line's x-extent.
     The interpolation returns NaN there (a piezometric line assigns pressure
     where it exists and nothing beyond it — the Slide2/RS2 convention), and NaN
     was read as zero. Both sampling layers now raise instead: the LEM slicer per
     slice base (naming the surface, the slice, its x, and the extent), and the
     FEM at the mesh node and Gauss point.

  There is deliberately NO requirement that a piezometric line span the domain.
  A line that stops short is a legitimate model (an upstream pool with no
  downstream pond, say) and must keep solving as long as nothing samples it
  beyond its end — checked here as the short-but-covering case, which must give
  the SAME factor of safety as a line carried on to the domain edge.

Deterministic, no network; builds its slope on docs/lem/files/xslope_acads_simple.xlsx.
The FEM half needs gmsh and is skipped cleanly without it.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/piezo_extent_guard.py
Exits non-zero on any failure.
"""
import os
import sys
import tempfile
import warnings

import numpy as np

warnings.filterwarnings("ignore")

from shapely.geometry import Polygon

from xslope.fileio import (load_slope_data, save_slope_data_to_xlsx,
                           build_ground_surface_from_polygons)
from xslope.slice import generate_slices
from xslope.solve import bishop, spencer

MAT_NAME = "piezo soil"
BAD_U = "peizo"                    # the transposition that used to zero u silently
_BASE_XLSX = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "docs", "lem", "files", "xslope_acads_simple.xlsx")

# Same homogeneous slope the suction guard uses: crest bench at y=30 (x 0..20),
# a face down to a toe bench at y=10 out to x=120. The circle's base bottoms out
# near y=4, so the piezo line at y=12 puts real pore pressure on the downslope
# half of the arc — the half a too-short line cuts off.
_SOIL_POLY = [(0, -10), (120, -10), (120, 10), (60, 10), (20, 30), (0, 30)]
_CIRCLE = {"Xo": 50.0, "Yo": 42.0, "Depth": 38.0, "R": 38.0}
_PIEZO_Y = 12.0                    # over the lower half of the arc -> real u on the slices a short line would cut off
_NUM_SLICES = 40
_FEM_TARGET = 8.0                  # coarse: this guard never solves, it only builds


def _base():
    return load_slope_data(_BASE_XLSX)


def _slope_data(base, piezo_x, piezo2_x=None):
    """The slope with its piezometric line spanning ``piezo_x`` = (x_lo, x_hi)."""
    m = dict(base["materials"][0])
    m.update(name=MAT_NAME, c=15.0, phi=30.0, gamma=20.0, gamma_sat=20.0,
             option="mc", u="piezo", ru=0.0, E=30000.0, nu=0.33)
    polys = [{"polygon": Polygon(_SOIL_POLY), "mat_id": 0}]
    gs, dom = build_ground_surface_from_polygons(polys)
    sd = dict(base)
    sd["materials"] = [m]
    sd["polygons"] = polys
    sd["ground_surface"] = gs
    sd["domain_polygon"] = dom
    sd["piezo_line"] = [(piezo_x[0], _PIEZO_Y), (piezo_x[1], _PIEZO_Y)]
    sd["piezo_line2"] = ([] if piezo2_x is None else
                         [(piezo2_x[0], _PIEZO_Y), (piezo2_x[1], _PIEZO_Y)])
    sd["piezo_phreatic"] = False
    sd["piezo_phreatic2"] = False
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


def _slices(sd):
    ok, res = generate_slices(sd, circle=sd["circles"][0], num_slices=_NUM_SLICES,
                              debug=False)
    if not ok:
        raise RuntimeError(f"generate_slices failed: {res}")
    return res[0]


def _fs(slice_df):
    out = {}
    for name, fn in (("bishop", bishop), ("spencer", spencer)):
        s, r = fn(slice_df.copy())
        out[name] = r["FS"] if s else None
    return out


def _raises(fn, *needles):
    """Call fn; return (ok, message). ok when it raised a ValueError whose text
    contains every needle."""
    try:
        fn()
    except ValueError as exc:
        text = str(exc)
        missing = [n for n in needles if n not in text]
        if missing:
            return False, f"raised but the message omits {missing}: {text[:220]}"
        return True, text
    except Exception as exc:                      # noqa: BLE001 - any other error is a failure
        return False, f"raised {type(exc).__name__} instead of ValueError: {exc}"
    return False, "did not raise"


# --------------------------------------------------------------------------- #
# A. the mat-sheet u option
# --------------------------------------------------------------------------- #

def check_u_option():
    """An unrecognized ``u`` string must be refused by name, at load."""
    failures = []
    base = _base()

    # Through the real Excel path: write the typo into the mat sheet and reload.
    sd = _slope_data(base, (-50.0, 400.0))
    sd["materials"][0]["u"] = BAD_U
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "bad_u.xlsx")
        save_slope_data_to_xlsx(sd, path)
        ok, msg = _raises(lambda: load_slope_data(path), MAT_NAME, BAD_U)
        if not ok:
            failures.append(f"u option: load_slope_data on u='{BAD_U}' {msg}")

        # And each allowed value must still load. 'PIEZO ' also proves the
        # normalization (strip + lower) happens before the check, so a stray
        # capital or trailing space is not a hard error.
        for good in ("none", "piezo", "seep", "ru", "PIEZO "):
            sd["materials"][0]["u"] = good
            gpath = os.path.join(tmp, f"good_{good.strip().lower()}.xlsx")
            save_slope_data_to_xlsx(sd, gpath)
            try:
                got = load_slope_data(gpath)["materials"][0].get("u")
            except Exception as exc:              # noqa: BLE001
                failures.append(f"u option: u='{good}' must load, but raised "
                                f"{type(exc).__name__}: {str(exc)[:160]}")
                continue
            if got != good.strip().lower():
                failures.append(f"u option: u='{good}' loaded as '{got}'")

    # The slicer's own backstop, in case the two allowed-value lists ever drift.
    sd2 = _slope_data(base, (-50.0, 400.0))
    sd2["materials"][0]["u"] = BAD_U
    ok, msg = _raises(lambda: _slices(sd2), MAT_NAME, BAD_U)
    if not ok:
        failures.append(f"u option: generate_slices on u='{BAD_U}' {msg}")
    return failures


# --------------------------------------------------------------------------- #
# B1. the LEM slicer
# --------------------------------------------------------------------------- #

def check_lem_extent():
    """A slice base that samples the line from outside it must raise; one inside
    a SHORT line must solve, and give the full-line answer."""
    failures = []
    base = _base()

    # Full-width line: the reference answer.
    df_full = _slices(_slope_data(base, (-50.0, 400.0)))
    fs_full = _fs(df_full)

    # Where the failure surface actually lives, so the two short lines below can be
    # placed on either side of it deliberately rather than by luck.
    x_lo = float(df_full["x_l"].min())
    x_hi = float(df_full["x_r"].max())
    x_stop = 0.5 * (x_lo + x_hi)

    # The half of the arc a too-short line would cut off must carry REAL pore
    # pressure, else the guard proves nothing: silently zeroing it has to be a
    # change the guard is actually preventing.
    cut = df_full[df_full["x_c"] > x_stop]
    if len(cut) == 0 or not (cut["u"] > 0).all():
        failures.append("setup: the slices a short line would cut off carry no pore "
                        "pressure — the extent guard would be vacuous")
        return failures

    # (B1a) SHORT but covering: the line stops well short of the domain (x=120) yet
    # still spans the whole failure surface. No domain-span requirement — this must
    # solve, and to the same factor of safety as the full-width line.
    df_short = _slices(_slope_data(base, (x_lo - 2.0, x_hi + 2.0)))
    fs_short = _fs(df_short)
    du = float(np.max(np.abs(df_short["u"].to_numpy() - df_full["u"].to_numpy())))
    if du != 0.0:
        failures.append(f"lem: a line that stops short of the domain but covers the "
                        f"failure surface changed the base pore pressure "
                        f"(max|du| = {du:.3e})")
    for k in fs_full:
        if fs_full[k] is None or fs_short[k] is None:
            failures.append(f"lem: {k} failed to solve on the short-but-covering line")
        elif abs(fs_full[k] - fs_short[k]) > 1e-9:
            failures.append(f"lem: a line that stops short of the domain but covers the "
                            f"failure surface changed {k} FS "
                            f"{fs_full[k]:.6f} -> {fs_short[k]:.6f}; there is no "
                            f"domain-span requirement")

    # (B1b) TOO short: the line ends inside the failure surface, so the downslope
    # slices sample it from outside. That used to be u = 0 and a higher FS.
    ok, msg = _raises(lambda: _slices(_slope_data(base, (x_lo - 2.0, x_stop))),
                      MAT_NAME, "x-extent", f"{x_stop:.6g}")
    if not ok:
        failures.append(f"lem: a slice base beyond the line's right end {msg}")

    # ... and on the other side.
    ok, msg = _raises(lambda: _slices(_slope_data(base, (x_stop, x_hi + 2.0))),
                      MAT_NAME, "x-extent", f"{x_stop:.6g}")
    if not ok:
        failures.append(f"lem: a slice base before the line's left end {msg}")

    # (B1c) The SECOND line (rapid drawdown) is sampled by the same slices and gets
    # the same treatment — named separately so the message points at the right line.
    ok, msg = _raises(
        lambda: _slices(_slope_data(base, (x_lo - 2.0, x_hi + 2.0),
                                    piezo2_x=(x_lo - 2.0, x_stop))),
        "second piezometric line", f"{x_stop:.6g}")
    if not ok:
        failures.append(f"lem: a slice base beyond the SECOND line's right end {msg}")

    # An absent second line is not "outside" anything — it must stay silent.
    try:
        _slices(_slope_data(base, (x_lo - 2.0, x_hi + 2.0)))
    except ValueError as exc:
        failures.append(f"lem: no second piezo line must not trip the guard: {exc}")
    return failures


# --------------------------------------------------------------------------- #
# B2. the FEM
# --------------------------------------------------------------------------- #

def _mesh(sd):
    from xslope.mesh import build_mesh_from_polygons
    polys = [{"coords": list(p["polygon"].exterior.coords)[:-1], "mat_id": p["mat_id"]}
             for p in sd["polygons"]]
    return build_mesh_from_polygons(polys, target_size=_FEM_TARGET, element_type="tri6")


def gmsh_available():
    try:
        import gmsh                                # noqa: F401
    except Exception:                              # noqa: BLE001
        return False
    return True


def check_fem_extent():
    """A mesh node — and a Gauss point — outside the line must raise."""
    from xslope.fem import build_fem_data, _prepare_fem_model
    failures = []
    base = _base()

    sd_full = _slope_data(base, (-50.0, 400.0))
    mesh = _mesh(sd_full)
    x_min = float(np.min(mesh["nodes"][:, 0]))
    x_max = float(np.max(mesh["nodes"][:, 0]))

    # (B2a) In-extent: builds, and really does carry pore pressure.
    fem_full = build_fem_data(sd_full, mesh)
    if not np.any(fem_full["u"] > 0):
        failures.append("setup: the full-width line put no pore pressure on the mesh — "
                        "the FEM extent guard would be vacuous")
        return failures

    # A line that stops short of the mesh but still covers it is fine (no
    # domain-span requirement), and gives the identical nodal field.
    fem_snug = build_fem_data(_slope_data(base, (x_min, x_max)), mesh)
    du = float(np.max(np.abs(fem_snug["u"] - fem_full["u"])))
    if du != 0.0:
        failures.append(f"fem: a line trimmed to the mesh's own x-extent changed the "
                        f"nodal pore pressure (max|du| = {du:.3e})")

    # (B2b) Out of extent at the node layer.
    x_stop = 0.5 * (x_min + x_max)
    ok, msg = _raises(lambda: build_fem_data(_slope_data(base, (x_min, x_stop)), mesh),
                      "mesh node", "x-extent", f"{x_stop:.6g}")
    if not ok:
        failures.append(f"fem: a mesh node beyond the line's right end {msg}")
    ok, msg = _raises(lambda: build_fem_data(_slope_data(base, (x_stop, x_max)), mesh),
                      "mesh node", "x-extent", f"{x_stop:.6g}")
    if not ok:
        failures.append(f"fem: a mesh node before the line's left end {msg}")

    # (B2c) The Gauss-point layer is guarded independently of the nodal one, so a
    # fem_data assembled elsewhere cannot slip a short line past it.
    short = dict(fem_full)
    short["piezo_line_coords"] = [(x_min, _PIEZO_Y), (x_stop, _PIEZO_Y)]
    ok, msg = _raises(lambda: _prepare_fem_model(short, debug_level=0),
                      "Gauss point", "x-extent", f"{x_stop:.6g}")
    if not ok:
        failures.append(f"fem: a Gauss point beyond the line's right end {msg}")

    # And the signed field the suction option reads takes the same guard.
    ok, msg = _raises(
        lambda: _prepare_fem_model(short, suction_phi_b={MAT_NAME: 20.0}, debug_level=0),
        "Gauss point", "x-extent", f"{x_stop:.6g}")
    if not ok:
        failures.append(f"fem: a Gauss point of the signed (suction) field beyond the "
                        f"line's right end {msg}")
    return failures


def check():
    """Run every guard; return a list of failure strings (empty on success)."""
    failures = []
    failures += check_u_option()
    failures += check_lem_extent()
    if gmsh_available():
        failures += check_fem_extent()
    return failures


if __name__ == "__main__":
    if not os.path.exists(_BASE_XLSX):
        print(f"SKIP: {_BASE_XLSX} not found")
        sys.exit(0)
    fails = check()
    for f in fails:
        print("FAIL:", f)
    print("piezo/u guard:", "PASS" if not fails else f"{len(fails)} FAILURE(S)")
    sys.exit(1 if fails else 0)
