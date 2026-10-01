"""Build the earth-dam sample workbooks under docs/lem/files from one definition
of the dam.

  xslope_earth_dam_up       docs/lem/samples.md, problem 8 (upstream slope,
                            steady piezometric line)
  xslope_earth_dam_down     docs/lem/samples.md, problem 8 (downstream slope)
  xslope_gsat_seep          docs/lem/samples.md, problem 16 (gamma_sat with the
                            water table from a steady seepage solution)
  xslope_gsat_rapid         docs/lem/samples.md, problem 16 (gamma_sat under
                            rapid drawdown, staged piezometric lines)
  xslope_earth_dam_rapid    the rapid-drawdown model with two seepage stages
                            (seep-bc and seep-bc (2)); a fixture of the seepage
                            and report checks and of the Studio screenshots

The dam is Duncan, Wright & Brandon's earth dam (Shear Strength and Slope
Stability, 2nd ed., p. 121): a shell and a clay core on a clay layer over sand,
crest at El. 317, foundation surface at El. 227, clay to El. 197 and sand to
El. 182, the base of the model. All five files are built from ``_dam()`` below,
so the section, the strengths and the base cannot differ between them.

The rapid-drawdown strength parameters d and psi (the Duncan, Wright & Wong
total-stress envelope of a low-permeability material) are written only into the
two rapid-drawdown runs, with the values of the course problem the model comes
from (CE 544, unit 2, rapid drawdown homework): core d = 300 psf, psi = 20 deg;
clay d = 100 psf, psi = 18 deg; shell and sand blank, as free-draining
materials. A blank d or psi loads as 0, and a slice whose material has d = 0 and
psi = 0 is treated as drained by the rapid-drawdown solver, which is how the
shell and sand have always been analyzed. The three steady runs leave d and psi
blank on every material: nothing in a steady analysis reads them.

The seepage meshes and solutions beside the two ``u = seep`` workbooks
(``{base}_mesh.json``, ``{base}_seep.csv``, ``{base}_seep2.csv``) are not written
here; tools/make_seep_sidecars.py solves them. Nothing written here changes what
they depend on (section, permeabilities, boundary conditions).

Every value below is the value the shipped workbooks carried on 2026-10-01
except: the base of xslope_gsat_rapid (Max depth 0 -> El. 182, problem 8's
base); the clay's psi in xslope_gsat_rapid (19 -> 18, the course problem's
value); and d and psi on the three steady runs (removed). Each model is
written through :func:`xslope.fileio.save_slope_data_to_xlsx`, which copies the
current template and never writes a formula cell.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/build_lem_samples.py
"""
import os

from xslope.fileio import save_slope_data_to_xlsx

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'docs', 'lem', 'files')

_EMPTY_BC = {'specified_heads': [], 'specified_fluxes': [], 'exit_face': []}

# =============================================================================
# The dam
# =============================================================================
#: Profile lines in material order: shell, core, clay, sand.
_PROFILE = [
    ([(0.0, 227.0), (270.0, 317.0), (320.0, 317.0), (590.0, 227.0)], 0),
    ([(240.0, 227.0), (280.0, 307.0), (310.0, 307.0), (350.0, 227.0)], 1),
    ([(-150.0, 227.0), (740.0, 227.0)], 2),
    ([(-150.0, 197.0), (740.0, 197.0)], 3),
]
_MAX_DEPTH = 182.0

#: name, gamma (pcf), c' (psf), phi' (deg), gamma_sat (pcf), k (ft/day)
_MATERIALS = [
    ('Shell', 125.0, 0.0, 34.0, 130.0, 864.0),
    ('Core', 122.0, 100.0, 26.0, 127.0, 0.0864),
    ('Clay', 123.0, 0.0, 24.0, 128.0, 0.864),
    ('Sand', 127.0, 0.0, 32.0, 132.0, 86.4),
]

#: Rapid-drawdown total-stress envelope (d psf, psi deg) of the low-permeability
#: materials; every other material is free-draining (d and psi blank).
_RAPID_D_PSI = {'Core': (300.0, 20.0), 'Clay': (100.0, 18.0)}

#: Reservoir at El. 302 (before drawdown) and El. 250 (after).
_PIEZO_FULL = [(-150.0, 302.0), (277.5, 302.0), (315.0, 275.0),
               (343.5, 240.0), (590.0, 227.0), (740.0, 227.0)]
_PIEZO_DRAWN = [(-150.0, 250.0), (251.5, 250.0), (315.0, 242.0),
                (347.0, 233.0), (590.0, 227.0), (740.0, 227.0)]

#: Seepage boundaries: the reservoir face as a specified head, the downstream
#: face and ground as a potential exit face.
_EXIT_FACE = [(320.0, 317.0), (590.0, 227.0), (740.0, 227.0)]
_SEEP_BC_FULL = {
    'specified_heads': [{'head': 302.0, 'kind': 'head',
                         'coords': [(-150.0, 227.0), (0.0, 227.0), (225.0, 302.0)]}],
    'specified_fluxes': [], 'exit_face': list(_EXIT_FACE)}
_SEEP_BC_DRAWN = {
    'specified_heads': [{'head': 250.0, 'kind': 'head',
                         'coords': [(-150.0, 227.0), (0.0, 227.0), (69.0, 250.0)]}],
    'specified_fluxes': [], 'exit_face': list(_EXIT_FACE)}

#: Starting circles (Depth option): one tangent to the top of the sand, one to
#: the base. Upstream and downstream differ only in the centre's x.
_YO = 500.0
_DEPTHS = (197.0, 182.0)
_XO_UP, _XO_DOWN = 100.0, 490.0


def _material(name, gamma, c, phi, gamma_sat, k, *, u, rapid, seep_props,
              E=700000.0, nu=0.3, t_cut=0.0):
    d, psi = _RAPID_D_PSI.get(name, (None, None)) if rapid else (None, None)
    return {
        'name': name, 'gamma': gamma, 'gamma_sat': gamma_sat, 'option': 'mc',
        'c': c, 'phi': phi, 'cp': 0.0, 'r_elev': 0.0, 'd': d, 'psi': psi,
        't_cut': t_cut, 'phi_b': None, 's_cap': None, 'Ss': None, 'Sy': None,
        'pow_a': 0.0, 'pow_b': 0.0, 'pow_c': 0.0, 'pow_d': 0.0,
        'u': u, 'ru': 0.0,
        'sigma_gamma': 0.0, 'sigma_c': 0.0, 'sigma_phi': 0.0,
        'sigma_cp': 0.0, 'sigma_d': 0.0, 'sigma_psi': 0.0,
        'k1': k if seep_props else 0.0, 'k2': k if seep_props else 0.0,
        'alpha': 0.0, 'unsat': 'lf',
        'kr0': 0.0001 if seep_props else 0.0, 'h0': -1.0 if seep_props else 0.0,
        'vg_a': 0.0, 'vg_n': 0.0, 'vg_l': 0.5,
        'E': E, 'nu': nu,
        'hb_sci': 0.0, 'hb_gsi': 0.0, 'hb_mi': 0.0, 'hb_d': 0.0,
    }


def _dam(*, u, rapid, gsat, seep_props, xo=_XO_UP, **mat_extra):
    """The dam with the water and analysis options of one run.

    u           the pore-pressure option of every material ('piezo' or 'seep')
    rapid       write the rapid-drawdown d and psi (only a rapid run reads them)
    gsat        give each zone its saturated unit weight
    seep_props  write the permeabilities and the linear-front unsaturated law
                (the runs that solve their own seepage)
    """
    return {
        'unit_system': 'imperial',
        'gamma_water': 62.4,
        'tcrack_depth': 0.0,
        'tcrack_water': 0.0,
        'k_seismic': 0.0,
        'water_loads': 'auto',
        'materials': [_material(n, g, c, f, gs if gsat else None, k, u=u,
                                rapid=rapid, seep_props=seep_props, **mat_extra)
                      for n, g, c, f, gs, k in _MATERIALS],
        'profile_lines': [{'coords': list(c), 'mat_id': m, 'size': None}
                          for c, m in _PROFILE],
        'max_depth': _MAX_DEPTH,
        'piezo_line': [],
        'piezo_line2': [],
        'piezo_phreatic': False,
        'piezo_phreatic2': False,
        'circles': [{'Xo': xo, 'Yo': _YO, 'Depth': dep, 'R': _YO - dep}
                    for dep in _DEPTHS],
        'non_circ': [],
        'seepage_bc': {'specified_heads': [dict(h, coords=list(h['coords']))
                                           for h in _SEEP_BC_FULL['specified_heads']],
                       'specified_fluxes': [],
                       'exit_face': list(_EXIT_FACE)},
        'seepage_bc2': dict(_EMPTY_BC),
        'has_seepage_bc2': False,
    }


def _write(sd, name):
    dst = os.path.join(OUT, name)
    save_slope_data_to_xlsx(sd, dst)
    return dst


# =============================================================================
# docs/lem/samples.md, problem 8: steady piezometric line, both slopes
# =============================================================================
def build_earth_dam_up():
    sd = _dam(u='piezo', rapid=False, gsat=False, seep_props=False)
    sd['piezo_line'] = list(_PIEZO_FULL)
    return _write(sd, 'xslope_earth_dam_up.xlsx')


def build_earth_dam_down():
    sd = _dam(u='piezo', rapid=False, gsat=False, seep_props=False, xo=_XO_DOWN)
    sd['piezo_line'] = list(_PIEZO_FULL)
    return _write(sd, 'xslope_earth_dam_down.xlsx')


# =============================================================================
# docs/lem/samples.md, problem 16: gamma_sat with a seepage water table
# =============================================================================
def build_gsat_seep():
    sd = _dam(u='seep', rapid=False, gsat=True, seep_props=True)
    sd['time_unit'] = 'day'
    return _write(sd, 'xslope_gsat_seep.xlsx')


# =============================================================================
# docs/lem/samples.md, problem 16: gamma_sat under rapid drawdown
# =============================================================================
def build_gsat_rapid():
    """Staged piezometric lines (pool at El. 302, then El. 250). E and nu are
    entered as 0 and t_cut is blank, as the shipped file has them: this file is
    solved by limit equilibrium only, which reads none of the three."""
    sd = _dam(u='piezo', rapid=True, gsat=True, seep_props=False,
              E=0.0, nu=0.0, t_cut=None)
    sd['piezo_line'] = list(_PIEZO_FULL)
    sd['piezo_line2'] = list(_PIEZO_DRAWN)
    return _write(sd, 'xslope_gsat_rapid.xlsx')


# =============================================================================
# Rapid drawdown from two seepage stages
# =============================================================================
def build_earth_dam_rapid():
    sd = _dam(u='seep', rapid=True, gsat=False, seep_props=True)
    sd['seepage_bc2'] = {
        'specified_heads': [dict(h, coords=list(h['coords']))
                            for h in _SEEP_BC_DRAWN['specified_heads']],
        'specified_fluxes': [], 'exit_face': list(_EXIT_FACE)}
    sd['has_seepage_bc2'] = True
    sd['time_unit'] = 'day'
    return _write(sd, 'xslope_earth_dam_rapid.xlsx')


BUILDERS = [
    build_earth_dam_up,
    build_earth_dam_down,
    build_gsat_seep,
    build_gsat_rapid,
    build_earth_dam_rapid,
]


def main():
    os.makedirs(OUT, exist_ok=True)
    for fn in BUILDERS:
        print('built', fn())


if __name__ == '__main__':
    main()
