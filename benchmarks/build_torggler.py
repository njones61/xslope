"""Build the Torggler (2016) benchmark input files.

Torggler, N. (2016), "Numerical Studies of Embedded Beam Row in Safety Analysis
in PLAXIS 2D", MSc thesis, Graz University of Technology.
https://diglib.tugraz.at/download.php?id=5891c94c5ba8d&location=browse

Two models are transcribed, each in an unsupported and a plate-supported
variant.  Both are dry (the thesis states "without water loads" for §3 and "No
water loads are applied again" for §4), plane strain, Mohr-Coulomb, SI.

  xslope_torggler_3a_*   thesis §3, homogeneous slope, 7.5 m plate
  xslope_torggler_3b_*   thesis §4, slope with a 1 m weak layer, 15 m plate

GEOMETRY (thesis Fig. 9 and Fig. 63)
------------------------------------
Both models are a 10 m high slope at 30 deg on 20 m of the same soil.  The
figures dimension the domain as a chain of horizontal runs summing to the total
width, and a chain of vertical runs summing to 30.0 m:

  Fig. 9  (§3):  10.0 + 5.0 + 5.0 + 17.32 + 4.68 + 5.0 + 10.0 = 57.0
  Fig. 63 (§4):  10.0 + 10.0 + 17.32 + 17.68 + 10.0           = 65.0
  both:          10.0 + 5.0 + 5.0 + 10.0                      = 30.0

so in both models the toe sits at (20, 20), the 17.32 m slope run reaches the
crest at (37.32, 30), and the base of the mesh is y = 0.  §4 is §3 "expanded on
the right side" (thesis §4.1) — same slope, plateau carried out to x = 65.

PLATE (thesis Table 4, Fig. 17)
-------------------------------
The plate is elastic with EA = 2.0e6 kN/m and EI = 4.0e4 kNm^2/m.  Thesis
Table 5 gives the section those stiffnesses come from — E = 2.0e8 kPa,
A = 0.01 m^2, I = 2.0e-4 m^4 at a spacing of 1.0 m — and thesis Eq. (21)
confirms the pair by deriving D_eq = sqrt(12 I / A) = 0.49 m from them.

XSLOPE's FEM smears a pile row over its spacing as EA = E*A/S, EI = E*I/S
(xslope/fem.py), so entering Torggler's own E, A, I and S reproduces his plate
stiffnesses exactly and with no derived quantity in between:

  EA = 2.0e8 * 0.01 / 1.0 = 2.0e6 kN/m     (Table 4: 2.0E6)
  EI = 2.0e8 * 2.0e-4 / 1.0 = 4.0e4 kNm^2/m (Table 4: 4.0E4)

S = 1.0 m is a continuous wall, which is what a 2D plate element is, so no
out-of-plane arching question arises.  D is left blank: it is only used to
derive A and I when those are absent, and here they are given.  The plate has
no rotational restraint at its head, so Fixity = free.

The plate station is dimensioned in Fig. 17 as 8.66 m from the crest and 8.66 m
from the toe — mid-slope, x = 20 + 8.66 = 28.66, where the slope surface is at
elevation 20 + 8.66*tan(30) = 25.0.  §4 does not print the station; see
docs/verification/ssrm.md for the figure read that places it at the same
station.

WEAK LAYER (thesis §4.1, Table 18)
-----------------------------------
"A weak layer with a thickness of one meter is added.  The layer is based on
the failure line from SLIDE in Appendix C, which is offset by 0.5 m in both
directions (coordinates in Table 18)."  Table 18 lists 32 points of that
failure line in slope-foot coordinates (its caption: X = 0 m, Y = 0 m at the
slope foot), transcribed verbatim below and shifted to the model frame by the
toe position (20, 20).  The band is the polyline buffered 0.5 m each side; the
centreline is extended past both daylight points before buffering so the band
is cut off by the ground surface rather than by a rounded polyline cap.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/build_torggler.py
"""
import math
import os

from shapely.geometry import LineString, Polygon
from shapely.ops import unary_union

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'docs', 'fem', 'files')
# A blank template carries no failure surface and so cannot be loaded; the
# corpus builders seed from an existing metric model and overwrite every field.
SEED = os.path.join(ROOT, 'docs', 'lem', 'files', 'xslope_acads_simple.xlsx')

# --- domain (metres) ---------------------------------------------------------
X_TOE, Y_TOE = 20.0, 20.0
SLOPE_RUN = 17.32                     # Fig. 9 / Fig. 63 dimension chain
SLOPE_H = 10.0
X_CREST, Y_CREST = X_TOE + SLOPE_RUN, Y_TOE + SLOPE_H     # (37.32, 30.0)
X_RIGHT_3A, X_RIGHT_3B = 57.0, 65.0
Y_BASE = 0.0

# --- plate -------------------------------------------------------------------
X_PLATE = X_TOE + SLOPE_RUN / 2.0                          # 28.66
Y_PLATE_HEAD = Y_TOE + SLOPE_H / 2.0                       # 25.0
E_PLATE, A_PLATE, I_PLATE, S_PLATE = 2.0e8, 0.01, 2.0e-4, 1.0

# --- Table 18: failure line, slope-foot coordinates --------------------------
TABLE_18 = [
    (-1.84, 0.00), (-1.18, -0.59), (-0.80, -0.92), (0.46, -1.79),
    (1.39, -2.21), (2.55, -2.54), (3.93, -2.73), (4.76, -2.75),
    (5.60, -2.77), (6.61, -2.70), (7.63, -2.63), (8.66, -2.50),
    (9.77, -2.03), (10.87, -1.57), (11.88, -1.07), (12.89, -0.57),
    (13.77, -0.05), (14.65, 0.46), (15.42, 0.99), (16.18, 1.51),
    (17.30, 2.40), (17.89, 2.94), (18.48, 3.47), (19.08, 4.01),
    (19.67, 4.54), (20.51, 5.33), (21.11, 5.94), (21.70, 6.54),
    (22.75, 7.65), (23.64, 8.62), (24.17, 9.16), (25.03, 10.00),
]
BAND_HALF = 0.5


def _ground(x_right):
    return [(0.0, Y_TOE), (X_TOE, Y_TOE), (X_CREST, Y_CREST), (x_right, Y_CREST)]


def _domain(x_right):
    return Polygon(_ground(x_right) + [(x_right, Y_BASE), (0.0, Y_BASE)])


def _material(name, gamma, c, phi, E, nu):
    """A dry Mohr-Coulomb material row; psi = 0 throughout (thesis Tables 1, 10)."""
    return {
        'name': name, 'gamma': gamma, 'gamma_sat': None, 'option': 'mc',
        'c': c, 'phi': phi, 'cp': 0.0, 'r_elev': 0.0, 'd': 0.0, 'psi': 0.0,
        't_cut': None, 'phi_b': None, 's_cap': None, 'Ss': None, 'Sy': None,
        'pow_a': 0.0, 'pow_b': 0.0, 'pow_c': 0.0, 'pow_d': 0.0,
        'u': 'none', 'ru': 0.0,
        'sigma_gamma': 0.0, 'sigma_c': 0.0, 'sigma_phi': 0.0,
        'sigma_cp': 0.0, 'sigma_d': 0.0, 'sigma_psi': 0.0,
        'k1': 0.0, 'k2': 0.0, 'alpha': 0.0, 'unsat': 'lf',
        'kr0': 0.0, 'h0': 0.0, 'vg_a': 0.0, 'vg_n': 0.0,
        'E': E, 'nu': nu,
        'hb_sci': 0.0, 'hb_gsi': 0.0, 'hb_mi': 0.0, 'hb_d': 0.0,
    }


def _plate(y_tip, label):
    """Torggler's plate as an XSLOPE pile row: E, A, I and S entered as his own
    Table 5 section so that E*A/S and E*I/S are his Table 4 EA and EI exactly."""
    return {
        'x1': X_PLATE, 'y1': Y_PLATE_HEAD, 'x2': X_PLATE, 'y2': y_tip,
        'H': None, 'theta_p': 0.0,
        'D_pile': None, 'S': S_PLATE,
        'E': E_PLATE, 'I': I_PLATE, 'area': A_PLATE,
        'V_cap': None, 'M_cap': None, 'fixity': 'free', 'appl': 'active',
        'label': label,
    }


def _base(x_right, target_size):
    """A loaded template stripped to the fields these models set."""
    sd = load_slope_data(SEED)
    sd['gamma_water'] = 9.81
    sd['unit_system'] = 'si'
    sd['tcrack_depth'] = 0.0
    sd['tcrack_water'] = 0.0
    sd['k_seismic'] = 0.0
    sd['max_depth'] = Y_BASE
    sd['piezo_line'] = []
    sd['piezo_line2'] = []
    sd['dloads'] = []
    sd['dloads2'] = []
    sd['dload_dirs'] = []
    sd['dload2_dirs'] = []
    sd['reinforce_lines'] = []
    sd['reinforcement_lines'] = []
    sd['pile_lines'] = []
    sd['line_loads'] = []
    sd['non_circ'] = []
    sd['circular'] = True
    sd['element_type'] = 'tri6'
    sd['target_size'] = target_size
    sd['ssrm_f_min'] = None
    sd['ssrm_f_max'] = None
    sd['k0'] = None
    sd['tension_srf'] = None
    sd['lem_method'] = None
    sd['num_slices'] = None
    sd['water_loads'] = 'manual'
    sd['seepage_bc'] = {'specified_heads': [], 'specified_fluxes': [], 'exit_face': []}
    sd['seepage_bc2'] = {'specified_heads': [], 'specified_fluxes': [], 'exit_face': []}
    sd['has_seepage_bc2'] = False
    # Search seed: centre above mid-slope at the toe elevation + 2H, tangent
    # elevations spanning toe-level to below it (see docs/usage modelling rules).
    sd['circles'] = [
        {'Xo': X_PLATE, 'Yo': Y_TOE + 2 * SLOPE_H, 'Depth': Y_TOE, 'R': Y_TOE + 2 * SLOPE_H - Y_TOE},
        {'Xo': X_PLATE, 'Yo': Y_TOE + 2 * SLOPE_H, 'Depth': Y_TOE - 3.0,
         'R': Y_TOE + 2 * SLOPE_H - (Y_TOE - 3.0)},
    ]
    return sd


# =============================================================================
# 3a — homogeneous slope (thesis §3)
# =============================================================================
def _sd_3a(with_plate, target_size):
    sd = _base(X_RIGHT_3A, target_size)
    # Table 1: gamma = gamma_sat = 16, E = 2000, nu = 0.4, c = 10, phi = 15, psi = 0
    sd['materials'] = [_material('Klei', 16.0, 10.0, 15.0, 2000.0, 0.4)]
    sd['profile_lines'] = [{'mat_id': 0, 'coords': _ground(X_RIGHT_3A), 'size': None}]
    if with_plate:
        sd['pile_lines'] = [_plate(Y_PLATE_HEAD - 7.5, 'plate 7.5 m')]   # Table 4
    return sd


def build_torggler_3a_nopile():
    """Thesis §3.1, unsupported homogeneous slope.  PLAXIS phi/c reduction
    SumMsf = 1.111 (Table 2 at mesh factor 0.01, quoted in Table 3)."""
    sd = _sd_3a(False, 1.0)
    dst = os.path.join(OUT, 'xslope_torggler_3a_nopile.xlsx')
    save_slope_data_to_xlsx(sd, dst)
    return dst


def build_torggler_3a_plate():
    """Thesis §3.2.1, same slope with the 7.5 m plate and no interfaces —
    the variant that matches XSLOPE's shared-node beam.  SumMsf = 1.175."""
    sd = _sd_3a(True, 1.0)
    dst = os.path.join(OUT, 'xslope_torggler_3a_plate.xlsx')
    save_slope_data_to_xlsx(sd, dst)
    return dst


# =============================================================================
# 3b — slope with a weak layer (thesis §4)
# =============================================================================
def _weak_band():
    """The 1 m weak band: Table 18 shifted to the model frame and buffered
    0.5 m each side.  The centreline is extended 2 m beyond each end along its
    end-segment direction so the buffer's flat caps fall outside the domain and
    the band is trimmed by the ground surface, not by its own cap."""
    pts = [(x + X_TOE, y + Y_TOE) for x, y in TABLE_18]

    def _extend(p_in, p_end, length):
        dx, dy = p_end[0] - p_in[0], p_end[1] - p_in[1]
        n = math.hypot(dx, dy)
        return (p_end[0] + dx / n * length, p_end[1] + dy / n * length)

    line = LineString([_extend(pts[1], pts[0], 2.0)] + pts
                      + [_extend(pts[-2], pts[-1], 2.0)])
    return line.buffer(BAND_HALF, cap_style=2, join_style=2)


def _sd_3b(with_plate, target_size):
    sd = _base(X_RIGHT_3B, target_size)
    # Table 10: body E = 15000, nu = 0.3, c = 10, phi = 25;
    #           weak layer E = 5000, nu = 0.3, c = 0.01, phi = 20; gamma = 16 both.
    sd['materials'] = [
        _material('Soil body', 16.0, 10.0, 25.0, 15000.0, 0.3),
        _material('Weak layer', 16.0, 0.01, 20.0, 5000.0, 0.3),
    ]
    domain = _domain(X_RIGHT_3B)
    band = _weak_band().intersection(domain)
    body = domain.difference(band)
    parts = list(getattr(body, 'geoms', [body]))
    parts.sort(key=lambda p: (-p.bounds[3], p.bounds[0]))     # upper wedge first
    polys = [{'polygon': _round(p), 'mat_id': 0, 'size': None} for p in parts]
    polys.append({'polygon': _round(band), 'mat_id': 1, 'size': 0.5})
    sd['polygons'] = polys
    sd['profile_lines'] = []
    # Slide's own failure line is the band's centreline, so it seeds the
    # non-circular search directly.  Both ends daylight on the ground surface
    # (the flat in front of the toe, and the plateau), so both are Free.
    sd['non_circ'] = [{'X': round(x + X_TOE, 6), 'Y': round(y + Y_TOE, 6),
                       'Movement': 'Free'} for x, y in TABLE_18]
    if with_plate:
        sd['pile_lines'] = [_plate(Y_PLATE_HEAD - 15.0, 'plate 15 m')]   # §4.2
    return sd


def _round(poly, dp=6):
    """Round a polygon's vertices so the written file is byte-stable.  The band
    and the body pieces come from one shapely difference, so their shared
    vertices are identical before rounding and stay identical after it."""
    return Polygon([(round(x, dp), round(y, dp)) for x, y in poly.exterior.coords])


def build_torggler_3b_nopile():
    """Thesis §4.1, unsupported slope with the weak layer.  PLAXIS phi/c
    reduction SumMsf = 1.045 (Table 11 at mesh factor 0.03, Table 12)."""
    sd = _sd_3b(False, 1.0)
    dst = os.path.join(OUT, 'xslope_torggler_3b_nopile.xlsx')
    save_slope_data_to_xlsx(sd, dst)
    return dst


def build_torggler_3b_plate():
    """Thesis §4.2, same slope with the 15 m plate.  SumMsf = 1.725, reported
    for the plate with interfaces and the plate without interfaces alike."""
    sd = _sd_3b(True, 1.0)
    dst = os.path.join(OUT, 'xslope_torggler_3b_plate.xlsx')
    save_slope_data_to_xlsx(sd, dst)
    return dst


BUILDERS = [
    build_torggler_3a_nopile,
    build_torggler_3a_plate,
    build_torggler_3b_nopile,
    build_torggler_3b_plate,
]


def main():
    os.makedirs(OUT, exist_ok=True)
    for fn in BUILDERS:
        print('built', fn())


if __name__ == '__main__':
    main()
