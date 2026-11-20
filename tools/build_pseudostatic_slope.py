"""Build Tutorial LEM-14's file pair (pseudostatic analysis, one stage and two).

The slope
---------
A 15 ft slope at 2H:1V with level ground 45 ft deep behind the crest and in front
of the toe, in one soil to a base 15 ft below the toe, facing right (the crest on
the left), so that a positive seismic coefficient pushes the finite element model
down the slope as well as the limit equilibrium one.

The soil is soil 1 of Duncan, Wright, and Brandon (2014), Table 10.3: the
Pilarcitos Dam fill, gamma = 135 pcf. The table gives it two strength sets from the
same consolidated-undrained tests, and the tutorial enters each in turn:

  - the R envelope, c_R = 60 psf and phi_R = 23 degrees, entered as c and phi, for
    the single-stage analysis;
  - c' = 0 and phi' = 45 degrees with the K_c = 1 envelope's d = 64 psf and
    psi = 24.4 degrees, for the two-stage analysis.

The same table is the source of the verification models in
benchmarks/published/build_dwb_t10_4.py (the book's own slope); this slope is a
different one.

The water table
---------------
The water case has a water table from Dupuit's solution for a slope receiving
recharge, with the bottom of the model impermeable: h^2 = h1^2 - (h1^2 - h2^2) s/L +
(W/k) s (L - s), measured up from the base, with s the distance from the back
edge. It stands 6 ft below the crest at the back edge (h1 = 24 ft) and comes out at
the toe (h2 = 15 ft), and W/k = (h1^2 - h2^2)/L^2 makes the back edge a groundwater
divide, where the water table is level. It steepens all the way to the toe and lies
at the ground surface in front of the slope. It is entered as points every 5 ft,
rounded to 0.01 ft.

The two files
-------------
  xslope_pseudostatic_start.xlsx   dry, the R envelope, no seismic coefficient: the
                                   file the tutorial opens.
  xslope_pseudostatic.xlsx         the state the tutorial ends in: the water table,
                                   the R envelope, k = 0.15.

Both carry E and nu for the strength reduction run at the end of the page, and the
finite element run it is made with (tri6 at TARGET_SIZE, the bracket F_MIN to
F_MAX, no K0 initial stress), so Build Mesh and Run FEM open on the tutorial's own
settings. The elastic pair is FEM-1's, E = 2,088,500 psf (100 MPa) and nu = 0.3.
Neither file carries a mesh sidecar.

Run:  PYTHONPATH=. python3 tools/build_pseudostatic_slope.py
"""

from __future__ import annotations

import math
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from xslope.fileio import (default_template_path, load_slope_data,  # noqa: E402
                           save_slope_data_to_xlsx)

TUTORIAL_FILES = os.path.join(REPO_ROOT, "docs", "tutorials", "files")
SEED = os.path.join(REPO_ROOT, "docs", "lem", "files", "xslope_acads_simple.xlsx")
START_OUT = "xslope_pseudostatic_start.xlsx"
DONE_OUT = "xslope_pseudostatic.xlsx"

H = 15.0                     # slope height, ft
RUN = 2.0                    # 2H:1V
BACK = 45.0                  # level ground behind the crest and in front of the toe, ft
X_CREST, X_TOE = 0.0, RUN * H
X_LEFT, X_RIGHT = X_CREST - BACK, X_TOE + BACK
Y_BASE = -H                  # the bottom of the model, one slope height below the toe
K = 0.15                     # the seismic coefficient of the tutorial's runs

GAMMA = 135.0                # pcf, Table 10.3
C_R, PHI_R = 60.0, 23.0      # the R envelope, Table 10.3 soil 1

#: FEM-1's elastic pair (Griffiths and Lane's 100 MPa, in psf).
E, NU = 2088500.0, 0.3
ELEMENT_TYPE = "tri6"
TARGET_SIZE = 1.5
F_MIN, F_MAX = 0.5, 1.5
#: No K0 initial stress: at k = 0.15 the slope does not stand at full strength
#: (F < 1), so no in-situ state can be equilibrated to carry into the trials.
K0 = None

#: The water table: 6 ft below the crest at the back edge, out at the toe.
WT_BACK = H - 6.0
WT_STEP = 5.0


def water_table():
    """The Dupuit water table with recharge, as (x, y) points every WT_STEP ft
    from the back edge to the toe, then along the ground to the right edge."""
    h1, h2 = WT_BACK - Y_BASE, 0.0 - Y_BASE
    L = X_TOE - X_LEFT
    wk = (h1 ** 2 - h2 ** 2) / L ** 2          # the back edge is a divide

    def y(x):
        s = x - X_LEFT
        return math.sqrt(h1 ** 2 - (h1 ** 2 - h2 ** 2) * s / L + wk * s * (L - s)) + Y_BASE

    n = int(round(L / WT_STEP))
    pts = [(X_LEFT + i * WT_STEP, round(y(X_LEFT + i * WT_STEP), 2)) for i in range(n)]
    return pts + [(X_TOE, 0.0), (X_RIGHT, 0.0)]


def _material(u):
    """The fill with the R envelope entered as c and phi, d and psi blank."""
    return {
        "name": "fill", "gamma": GAMMA, "gamma_sat": None, "option": "mc",
        "c": C_R, "phi": PHI_R, "cp": 0.0, "r_elev": 0.0, "d": None, "psi": None,
        "t_cut": 0.0, "phi_b": None, "s_cap": None, "Ss": None, "Sy": None,
        "pow_a": 0.0, "pow_b": 0.0, "pow_c": 0.0, "pow_d": 0.0,
        "u": u, "ru": 0.0,
        "sigma_gamma": 0.0, "sigma_c": 0.0, "sigma_phi": 0.0,
        "sigma_cp": 0.0, "sigma_d": 0.0, "sigma_psi": 0.0,
        "k1": 0.0, "k2": 0.0, "alpha": 0.0, "unsat": "lf",
        "kr0": 0.0, "h0": 0.0, "vg_a": 0.0, "vg_n": 0.0,
        "E": E, "nu": NU,
        "hb_sci": 0.0, "hb_gsi": 0.0, "hb_mi": 0.0, "hb_d": 0.0,
    }


def _model(wet):
    sd = load_slope_data(SEED)
    sd.pop("mesh", None)
    sd["unit_system"] = "imperial"
    sd["gamma_water"] = 62.4
    sd["tcrack_depth"] = 0.0
    sd["tcrack_water"] = 0.0
    sd["k_seismic"] = K if wet else 0.0
    sd["k_vertical"] = 0.0
    sd["materials"] = [_material("piezo" if wet else "none")]
    sd["profile_lines"] = [{"coords": [(X_LEFT, H), (X_CREST, H), (X_TOE, 0.0),
                                       (X_RIGHT, 0.0)],
                            "mat_id": 0, "size": None}]
    sd["polygons"] = []
    sd["max_depth"] = Y_BASE
    sd["piezo_line"] = water_table() if wet else []
    sd["piezo_line2"] = []
    sd["water_loads"] = "auto"
    for key in ("dloads", "dload_dirs", "dloads2", "dload2_dirs", "reinforcement_lines",
                "reinforce_lines", "pile_lines", "line_loads", "non_circ"):
        sd[key] = []
    sd["circular"] = True
    sd["search_window"] = {}
    # One starting circle through the toe: centered over it, tangent to the toe's
    # elevation, and coming out on the crest just behind the top of the face.
    sd["circles"] = [{"Xo": X_TOE, "Yo": 40.0, "Depth": 0.0, "R": 40.0}]
    sd["element_type"] = ELEMENT_TYPE
    sd["target_size"] = TARGET_SIZE
    sd["ssrm_f_min"] = F_MIN
    sd["ssrm_f_max"] = F_MAX
    sd["k0"] = K0
    return sd


def _write(sd, filename):
    out = os.path.join(TUTORIAL_FILES, filename)
    save_slope_data_to_xlsx(sd, out, template=default_template_path())
    for sidecar in (f"{os.path.splitext(out)[0]}_mesh.json",
                    f"{os.path.splitext(out)[0]}_seep.csv"):
        if os.path.exists(sidecar):
            raise SystemExit(f"tutorial downloads are sidecar-free, but "
                             f"{os.path.basename(sidecar)} sits beside {filename}")
    print(f"wrote {os.path.relpath(out, REPO_ROOT)}")
    return out


def build():
    return [_write(_model(wet=False), START_OUT), _write(_model(wet=True), DONE_OUT)]


if __name__ == "__main__":
    build()
