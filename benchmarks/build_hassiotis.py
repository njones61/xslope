"""Build the Hassiotis et al. (1997) pile-stabilized slope input files.

Hassiotis, S., Chameau, J. L., and Gunaratne, M. (1997), "Design method for
stabilization of slopes with piles", Journal of Geotechnical and
Geoenvironmental Engineering, 123(4), 314-323.
https://doi.org/10.1061/(ASCE)1090-0241(1997)123:4(314)

A homogeneous 30-degree slope, 13.7 m high, in a soil with c = 23.94 kPa,
phi = 10 degrees, gamma = 19.63 kN/m3, dry.  One row of 1.0 m piles at 2.5 m
centres is placed 13.7 m and (in a second case) 23.1 m horizontally from the
toe.  The published clear-to-centre spacing ratio D2/D1 = 0.6 is reproduced
exactly by D = 1.0 with S = 2.5 (clear opening 1.5 m, 1.5/2.5 = 0.6).

Geometry.  The paper fixes the slope (height and angle) and the pile stations
as horizontal distances from the toe; it does not dimension the mesh box.  The
domain here carries 30 m of flat ground in front of the toe, 20 m of the same
soil below it, and 36 m of crest platform, all far enough from the slope that
they do not enter the critical mechanism.

Pile force.  H is left blank on the piles sheet, so XSLOPE computes the
Ito & Matsui (1975) lateral force for every trial surface from D and S and the
soil above the surface at the pile, and divides the per-pile force by the
2.5 m centre spacing to get the per-unit-width value the slice equations use.
V_cap and M_cap are left blank: the paper specifies no structural capacities,
and the published factors of safety are the full soil force, so a cap would
change the quantity being compared.  Appl = active, i.e. the force is applied
as published rather than divided by the computed factor of safety.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/build_hassiotis.py
"""
import math
import os

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'docs', 'lem', 'files')
SEED = os.path.join(OUT, 'xslope_acads_simple.xlsx')

H_SLOPE = 13.7
BETA = 30.0
GAMMA, COH, PHI = 19.63, 23.94, 10.0
D_PILE, S_PILE, L_PILE = 1.0, 2.5, 17.0

X_TOE, Y_TOE = 30.0, 20.0
RUN = H_SLOPE / math.tan(math.radians(BETA))          # 23.7291 m
X_CREST, Y_CREST = X_TOE + RUN, Y_TOE + H_SLOPE
X_RIGHT, Y_BASE = 90.0, 0.0

GROUND = [(0.0, Y_TOE), (X_TOE, Y_TOE),
          (round(X_CREST, 4), Y_CREST), (X_RIGHT, Y_CREST)]


def _pile_at(distance_from_toe):
    """A vertical pile row whose head sits on the slope face `distance_from_toe`
    metres downslope-horizontally from the toe.  theta_p = 0: the Ito & Matsui
    reaction acts horizontally, which is the direction the plastic-deformation
    theory derives it in."""
    x = X_TOE + distance_from_toe
    y_head = Y_TOE + distance_from_toe * math.tan(math.radians(BETA))
    return {
        'x1': round(x, 4), 'y1': round(y_head, 4),
        'x2': round(x, 4), 'y2': round(y_head - L_PILE, 4),
        'H': None, 'theta_p': 0.0,
        'D_pile': D_PILE, 'S': S_PILE,
        'E': None, 'I': None, 'area': None,
        'V_cap': None, 'M_cap': None, 'fixity': 'free', 'appl': 'active',
        'label': 'pile row %.1f m from toe' % distance_from_toe,
    }


def _slope_data():
    sd = load_slope_data(SEED)
    m = sd['materials'][0]
    m.update(name='Soil', gamma=GAMMA, gamma_sat=None, option='mc',
             c=COH, phi=PHI, u='none', ru=0.0, psi=0.0)
    sd['materials'] = [m]
    sd['gamma_water'] = 9.81
    sd['unit_system'] = 'si'
    sd['tcrack_depth'] = 0.0
    sd['tcrack_water'] = 0.0
    sd['k_seismic'] = 0.0
    sd['profile_lines'] = [{'mat_id': 0, 'coords': GROUND, 'size': None}]
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
    sd['water_loads'] = 'manual'
    sd['element_type'] = None
    sd['target_size'] = None
    # Seed circles: centre above mid-slope at the toe elevation + 2H, one
    # tangent to the toe elevation and one 6 m below it.
    x_mid = round(X_TOE + RUN / 2.0, 4)
    y_c = Y_TOE + 2 * H_SLOPE
    sd['circles'] = [
        {'Xo': x_mid, 'Yo': y_c, 'Depth': Y_TOE, 'R': y_c - Y_TOE},
        {'Xo': x_mid, 'Yo': y_c, 'Depth': Y_TOE - 6.0, 'R': y_c - (Y_TOE - 6.0)},
    ]
    return sd


def build_hassiotis():
    """The unreinforced slope.  Published: 1.08 (Hassiotis, friction circle),
    1.12 (Hull & Poulos, Bishop), 1.11 (Ausilio et al., friction circle)."""
    sd = _slope_data()
    dst = os.path.join(OUT, 'xslope_hassiotis.xlsx')
    save_slope_data_to_xlsx(sd, dst)
    return dst


def _windowed(sd, pile):
    """Restrict the search to surfaces the pile row actually stabilizes.

    An unrestricted search on this slope finds a deep surface that passes
    below the pile tip and collects no pile force at all — the published
    comparisons are for the mechanism through the row, and the source
    comparisons restrict the search the same way (a range at the crest and a
    range at the toe).  Stated once here for both rows: the surface daylights
    within a few metres of the toe, enters behind the crest, and its lowest
    point stays above the pile tip.
    """
    sd['search_window'] = {
        'entry_x_min': 54.0, 'entry_x_max': 75.0,
        'exit_x_min': 25.0, 'exit_x_max': 32.0,
        'max_tangent_depth': round(pile['y2'] + 0.1, 2),
    }
    return sd


def build_hassiotis_p1():
    """Pile row 13.7 m from the toe.  Published: 1.82 (Hassiotis),
    1.45 (Hull & Poulos)."""
    sd = _slope_data()
    pile = _pile_at(13.7)
    sd['pile_lines'] = [pile]
    _windowed(sd, pile)
    dst = os.path.join(OUT, 'xslope_hassiotis_p1.xlsx')
    save_slope_data_to_xlsx(sd, dst)
    return dst


def build_hassiotis_p2():
    """Pile row 23.1 m from the toe.  Published: 1.64 (Hassiotis)."""
    sd = _slope_data()
    pile = _pile_at(23.1)
    sd['pile_lines'] = [pile]
    _windowed(sd, pile)
    dst = os.path.join(OUT, 'xslope_hassiotis_p2.xlsx')
    save_slope_data_to_xlsx(sd, dst)
    return dst


BUILDERS = [build_hassiotis, build_hassiotis_p1, build_hassiotis_p2]


def main():
    for fn in BUILDERS:
        print('built', fn())


if __name__ == '__main__':
    main()
