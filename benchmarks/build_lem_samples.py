"""Build the limit-equilibrium sample workbooks under docs/lem/files that are
variants of another shipped sample model.

  xslope_gsat_rapid   docs/lem/samples.md, problem 16 (the gamma_sat rapid-
                      drawdown case): the earth dam of problem 8 with a saturated
                      unit weight for each zone and a post-drawdown
                      piezometric line

The dam is not transcribed here. It is read from problem 8's own workbook
(``xslope_earth_dam_up.xlsx``) and only the departures that make it problem 16
are written on top, so the two models cannot drift apart: a change to the dam
in problem 8 shows up here as a rebuild difference
(``benchmarks/verify_rebuild.py --group lem_samples``) rather than going
unnoticed. Every departure below is the value the shipped workbook carried when
this builder was written (2026-10-01), except the base of the model: the shipped
file declared Max depth = 0, where problem 8 declares El. 182 (the bottom of the
sand). That one value is now inherited from problem 8 like the rest of the dam.

Each model is written through :func:`xslope.fileio.save_slope_data_to_xlsx`,
which copies the current template and never writes a formula cell.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/build_lem_samples.py
"""
import copy
import os

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'docs', 'lem', 'files')
# The source models are read from the shipped tree, never from OUT, so a rebuild
# into a scratch directory still starts from the committed problem 8.
SRC = os.path.join(ROOT, 'docs', 'lem', 'files')

#: Problem 8's workbook: the earth dam every variant here starts from.
EARTH_DAM = 'xslope_earth_dam_up.xlsx'


def _earth_dam():
    """Problem 8's model as the loader returns it, with the geometry the loader
    derives (and the writer recomputes) removed."""
    sd = copy.deepcopy(load_slope_data(os.path.join(SRC, EARTH_DAM)))
    for k in ('ground_surface', 'domain_polygon', 'tcrack_surface', 'polygons',
              'filename', 'template_version'):
        sd.pop(k, None)
    return sd


def _write(sd, name):
    dst = os.path.join(OUT, name)
    save_slope_data_to_xlsx(sd, dst)
    return dst


# =============================================================================
# docs/lem/samples.md, problem 16: gamma_sat under rapid drawdown
# =============================================================================
#: Saturated unit weight (pcf) per zone, in problem 8's material order.
_GSAT_RAPID_GAMMA_SAT = {'Shell': 130.0, 'Core': 127.0, 'Clay': 128.0,
                         'Sand': 132.0}

#: The post-drawdown piezometric line: pool at El. 250, the phreatic surface
#: falling through the dam to the downstream ground at El. 227. The
#: pre-drawdown line (pool at El. 302) is problem 8's own piezometric line.
_GSAT_RAPID_PIEZO2 = [(-150.0, 250.0), (251.5, 250.0), (315.0, 242.0),
                      (347.0, 233.0), (590.0, 227.0), (740.0, 227.0)]


def build_gsat_rapid():
    """Problem 8's dam, reservoir at El. 302 before drawdown and El. 250 after.

    Departures from problem 8, each as the shipped file carries it:
      * gamma_sat on every zone (``_GSAT_RAPID_GAMMA_SAT``);
      * a second piezometric line, the post-drawdown water
        (``_GSAT_RAPID_PIEZO2``);
      * E and nu entered as 0 and t_cut left blank on every zone. This file is
        solved by limit equilibrium only, which reads none of the three.
    """
    sd = _earth_dam()
    for m in sd['materials']:
        m['gamma_sat'] = _GSAT_RAPID_GAMMA_SAT[m['name']]
        m['E'] = 0.0
        m['nu'] = 0.0
        m['t_cut'] = None
    sd['piezo_line2'] = list(_GSAT_RAPID_PIEZO2)
    return _write(sd, 'xslope_gsat_rapid.xlsx')


BUILDERS = [
    build_gsat_rapid,
]


def main():
    os.makedirs(OUT, exist_ok=True)
    for fn in BUILDERS:
        print('built', fn())


if __name__ == '__main__':
    main()
