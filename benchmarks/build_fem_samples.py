"""Build the finite element sample workbooks under docs/fem/files that no other
builder writes.

  xslope_noncircular_fem          docs/fem/samples.md, problem 3 (the thin weak
                                  clay layer; fem_ssrm row)
  xslope_simple_mult_layers_fem   docs/fem/samples.md, problem 4 (two-layer c-phi
                                  slope; fem_reliability row REL-FEM)
  xslope_piles_fem                the two-pile slope; a fixture of report_check,
                                  nr_ssrm_check, fem_1d_details_check, the
                                  round-trip set and the FEM pile preflight rules
  xslope_piles_fem_nopile         the same slope with the two pile rows cleared,
                                  the before of tools/make_fem_docs_sidecars.py's
                                  before and after
  xslope_griffiths1_load          the Griffiths & Lane Example 1 slope (imperial)
                                  with an 800 psf strip load on the crest; the
                                  loaded-slope fixture of report_check and the
                                  v21 round-trip

Every one of them is solved by strength reduction (or the reliability analysis
built on it) and by nothing else: no page runs a limit-equilibrium analysis on
any of them. So none carries a failure surface. The loader does not require one,
and preflight asks for one only of a run that reads it (docs/usage/preflight.md).

The meshes beside these workbooks (``{base}_mesh.json``) and the solved-run
sidecars are not written here; they are produced by the runs that use them.

Every value below is the value the shipped workbook carried when this builder
was written (2026-10-01): it transcribes the files, it does not re-derive them.
Each model is written through :func:`xslope.fileio.save_slope_data_to_xlsx`,
which copies the current template and never writes a formula cell.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/build_fem_samples.py
"""
import os

from xslope.fileio import save_slope_data_to_xlsx

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'docs', 'fem', 'files')

_NO_SEEP_BC = {'specified_heads': [], 'specified_fluxes': [], 'exit_face': []}


def _material(name, gamma, c, phi, E, nu, t_cut=0.0, u='none',
              sigma_gamma=0.0, sigma_c=0.0, sigma_phi=0.0):
    """One Mohr-Coulomb row of the mat sheet, every other column at the value an
    unfilled cell loads as. ``t_cut`` is None where the shipped cell is blank and
    0.0 where it holds a zero; both load as no tensile cap."""
    return {
        'name': name, 'gamma': gamma, 'gamma_sat': None, 'option': 'mc',
        'c': c, 'phi': phi, 'cp': 0.0, 'r_elev': 0.0, 'd': 0.0, 'psi': 0.0,
        't_cut': t_cut, 'phi_b': None, 's_cap': None, 'Ss': None, 'Sy': None,
        'pow_a': 0.0, 'pow_b': 0.0, 'pow_c': 0.0, 'pow_d': 0.0,
        'u': u, 'ru': 0.0,
        'sigma_gamma': sigma_gamma, 'sigma_c': sigma_c, 'sigma_phi': sigma_phi,
        'sigma_cp': 0.0, 'sigma_d': 0.0, 'sigma_psi': 0.0,
        'k1': 0.0, 'k2': 0.0, 'alpha': 0.0, 'unsat': 'lf',
        'kr0': 0.0, 'h0': 0.0, 'vg_a': 0.0, 'vg_n': 0.0, 'vg_l': 0.5,
        'E': E, 'nu': nu,
        'hb_sci': 0.0, 'hb_gsi': 0.0, 'hb_mi': 0.0, 'hb_d': 0.0,
    }


def _model(materials, profile_lines, max_depth, **extra):
    """An imperial, dry-by-default model with no failure surface. The water-load
    mode is automatic, as every one of these files carries it."""
    sd = {
        'unit_system': 'imperial',
        'gamma_water': 62.4,
        'tcrack_depth': 0.0,
        'tcrack_water': 0.0,
        'k_seismic': 0.0,
        'water_loads': 'auto',
        'materials': materials,
        'profile_lines': [{'coords': list(c), 'mat_id': m, 'size': None}
                          for c, m in profile_lines],
        'max_depth': max_depth,
        'piezo_phreatic': False,
        'piezo_phreatic2': False,
        'circles': [],
        'non_circ': [],
        'seepage_bc': dict(_NO_SEEP_BC),
        'seepage_bc2': dict(_NO_SEEP_BC),
        'has_seepage_bc2': False,
    }
    sd.update(extra)
    return sd


def _write(sd, name):
    dst = os.path.join(OUT, name)
    save_slope_data_to_xlsx(sd, dst)
    return dst


# =============================================================================
# docs/fem/samples.md, problem 3: thin weak clay layer
# =============================================================================
def build_noncircular_fem():
    """The LEM non-circular sample's section (UTEXASED manual problem) with E and
    nu added: sand fill over sand, a 2 ft soft clay layer, dense sand to y = -10.
    The water is a horizontal piezometric line at y = -2 read by the sand only."""
    sd = _model(
        materials=[
            _material('Sand Fill', 120.0, 0.0, 37.0, 1.0e6, 0.30),
            _material('Sand', 123.0, 0.0, 33.0, 7.0e5, 0.30, u='piezo'),
            _material('Soft Clay', 118.0, 200.0, 0.0, 6.0e4, 0.40),
            _material('Dense Sand', 131.0, 0.0, 37.0, 1.5e6, 0.28),
        ],
        profile_lines=[
            ([(0.0, 0.0), (30.0, 10.0), (50.0, 10.0)], 0),
            ([(-20.0, 0.0), (50.0, 0.0)], 1),
            ([(-20.0, -4.0), (50.0, -4.0)], 2),
            ([(-20.0, -6.0), (50.0, -6.0)], 3),
        ],
        max_depth=-10.0,
        piezo_line=[(-20.0, -2.0), (50.0, -2.0)],
    )
    return _write(sd, 'xslope_noncircular_fem.xlsx')


# =============================================================================
# docs/fem/samples.md, problem 4: two-layer c-phi slope (reliability)
# =============================================================================
def build_simple_mult_layers_fem():
    """Embankment over a foundation layer, with the standard deviations the
    Taylor-series reliability analysis perturbs (c, phi, gamma on both)."""
    sd = _model(
        materials=[
            _material('embankment', 130.0, 70.0, 20.0, 5.0e5, 0.35, t_cut=None,
                      sigma_gamma=6.5, sigma_c=18.0, sigma_phi=2.0),
            _material('foundation', 135.0, 140.0, 20.0, 5.0e5, 0.35, t_cut=None,
                      sigma_gamma=6.75, sigma_c=35.0, sigma_phi=2.0),
        ],
        profile_lines=[
            ([(0.0, 0.0), (40.0, 20.0), (90.0, 20.0)], 0),
            ([(-30.0, 0.0), (90.0, 0.0)], 1),
        ],
        max_depth=-10.0,
    )
    return _write(sd, 'xslope_simple_mult_layers_fem.xlsx')


# =============================================================================
# The two-pile slope and its unstabilized twin
# =============================================================================
_PILES_PROFILE = [([(-30.0, 0.0), (0.0, 0.0), (20.0, 20.0), (80.0, 20.0)], 0)]


def _pile(x, y_head, y_tip):
    """A vertical 2 ft pile row at 6 ft centres, free at both ends. E is
    518,400,000 psf; the section (area, I) is derived from D by the engine."""
    return {
        'x1': x, 'y1': y_head, 'x2': x, 'y2': y_tip,
        'H': None, 'theta_p': 0.0, 'D_pile': 2.0, 'S': 6.0,
        'E': 518400000.0, 'I': None, 'area': None,
        'V_cap': 46000.0, 'M_cap': 60000.0,
        'head_fixity': 'free', 'tip_fixity': 'free', 'appl': 'active',
        'label': 'pile',
    }


def _piles_model(with_piles):
    extra = {}
    if with_piles:
        # Each head on the slope face (y = x on the 1:1 face), each tip at the
        # base of the domain.
        extra['pile_lines'] = [_pile(5.0, 5.0, -10.0), _pile(10.0, 10.0, -10.0)]
    return _model(
        materials=[_material('soil', 120.0, 200.0, 20.0, 2.0e6, 0.30)],
        profile_lines=_PILES_PROFILE,
        max_depth=-10.0,
        **extra,
    )


def build_piles_fem():
    return _write(_piles_model(True), 'xslope_piles_fem.xlsx')


def build_piles_fem_nopile():
    return _write(_piles_model(False), 'xslope_piles_fem_nopile.xlsx')


# =============================================================================
# Griffiths & Lane Example 1 slope, imperial, with a crest strip load
# =============================================================================
def build_griffiths1_load():
    """A 50 ft, 2:1 homogeneous slope with c/(gamma H) = 0.05 and phi = 20 deg
    (the Griffiths & Lane Example 1 ratios), firm base at the toe, and an
    800 psf normal strip load on the crest from x = 30 to 45."""
    sd = _model(
        materials=[_material('soil', 125.0, 312.5, 20.0, 7.0e5, 0.30, t_cut=None)],
        profile_lines=[([(0.0, 50.0), (60.0, 50.0), (160.0, 0.0)], 0)],
        max_depth=0.0,
        dloads=[[{'X': 30.0, 'Y': 50.0, 'Normal': 800.0},
                 {'X': 45.0, 'Y': 50.0, 'Normal': 800.0}]],
        dload_dirs=['normal'],
    )
    return _write(sd, 'xslope_griffiths1_load.xlsx')


BUILDERS = [
    build_noncircular_fem,
    build_simple_mult_layers_fem,
    build_piles_fem,
    build_piles_fem_nopile,
    build_griffiths1_load,
]


def main():
    os.makedirs(OUT, exist_ok=True)
    for fn in BUILDERS:
        print('built', fn())


if __name__ == '__main__':
    main()
