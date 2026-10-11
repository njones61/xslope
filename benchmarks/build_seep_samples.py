"""Build the steady seepage sample workbooks under docs/seep/files that no other
builder writes.

  xslope_sea_trench   docs/seep/samples.md, the sea trench (two sheetpile walls
                      through silt into silty clay, the sea on both sides, the
                      pumped trench between them); type=seep row
  xslope_earth_dam2   docs/seep/samples.md, the earth dam with a clay core, a
                      chimney filter and a blanket drain; type=seep row

Both models are seepage only: no failure surface, no strength, no unit weight.
The loader does not require a surface, and preflight asks for one only of a run
that reads it (docs/usage/preflight.md).

Both carry their conductivities in ft/day and declare it: the Time selector
(main sheet, Time) is "day". The length unit comes from the Units selector
(Imperial), so the conductivity unit the program reports is ft/day and the flow
rate ft³/day per ft.

Every other value below is the value the shipped workbook carried when this
builder was written (2026-10-01): it transcribes the files, it does not re-derive
them. The time unit is the one departure; the files shipped with the Time
selector blank.

The companions beside xslope_earth_dam2.xlsx (``_mesh.json``, ``_seep.csv``) are
not written here. The mesh is the one its shipped field was solved on, and the
field is re-solved on that mesh by ``tools/make_seep_sidecars.py``. The sea trench
has no companion.

Each model is written through :func:`xslope.fileio.save_slope_data_to_xlsx`,
which copies the current template and never writes a formula cell.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/build_seep_samples.py
"""
import os

from xslope.fileio import save_slope_data_to_xlsx

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'docs', 'seep', 'files')

#: The time base of every conductivity below.
TIME_UNIT = 'day'

_NO_SEEP_BC = {'specified_heads': [], 'specified_fluxes': [], 'exit_face': []}


def _material(name, k):
    """One seepage-only row of the mat sheet: an isotropic, saturated-only
    conductivity ``k`` (ft/day) with the linear-front unsaturated defaults the
    files carry, and every strength and weight column at the value an unfilled
    cell loads as."""
    return {
        'name': name, 'gamma': 0.0, 'gamma_sat': None, 'option': '',
        'c': 0.0, 'phi': 0.0, 'cp': 0.0, 'r_elev': 0.0, 'd': 0.0, 'psi': 0.0,
        't_cut': None, 'phi_b': None, 's_cap': None, 'Ss': None, 'Sy': None,
        'pow_a': 0.0, 'pow_b': 0.0, 'pow_c': 0.0, 'pow_d': 0.0,
        'u': 'none', 'ru': 0.0,
        'sigma_gamma': 0.0, 'sigma_c': 0.0, 'sigma_phi': 0.0,
        'sigma_cp': 0.0, 'sigma_d': 0.0, 'sigma_psi': 0.0,
        'k1': k, 'k2': k, 'alpha': 0.0, 'unsat': 'lf',
        'kr0': 0.001, 'h0': -1.0, 'vg_a': 0.0, 'vg_n': 0.0, 'vg_l': 0.5,
        'E': 0.0, 'nu': 0.0,
        'hb_sci': 0.0, 'hb_gsi': 0.0, 'hb_mi': 0.0, 'hb_d': 0.0,
    }


def _model(materials, profile_lines, seepage_bc, water_loads):
    """An imperial seepage model in ft/day with no failure surface and no base
    depth (the profile lines close the domain)."""
    return {
        'unit_system': 'imperial',
        'time_unit': TIME_UNIT,
        'gamma_water': 62.4,
        'tcrack_depth': 0.0,
        'tcrack_water': 0.0,
        'k_seismic': 0.0,
        'water_loads': water_loads,
        'materials': materials,
        'profile_lines': [{'coords': list(c), 'mat_id': m, 'size': None}
                          for c, m in profile_lines],
        'max_depth': 0.0,
        'piezo_line': [],
        'piezo_line2': [],
        'piezo_phreatic': False,
        'piezo_phreatic2': False,
        'circles': [],
        'non_circ': [],
        'seepage_bc': seepage_bc,
        'seepage_bc2': dict(_NO_SEEP_BC),
        'has_seepage_bc2': False,
    }


def _head(h, coords):
    return {'head': h, 'coords': list(coords), 'kind': 'head'}


def _write(sd, name):
    dst = os.path.join(OUT, name)
    save_slope_data_to_xlsx(sd, dst)
    return dst


# =============================================================================
# docs/seep/samples.md: the sea trench
# =============================================================================
def build_sea_trench():
    """Silt banks (k = 0.5 ft/day) to El. 26 over silty clay (k = 0.1 ft/day);
    each sheetpile is a 0.4 ft wide notch from El. 26 to El. 22 at x = 43 and
    x = 67. The sea stands at El. 56 on both banks; the trench floor between the
    walls is pumped to El. 26.

    The water-load mode is manual, as the shipped file carries it, with no
    distributed loads: no limit-equilibrium run reads this model."""
    sd = _model(
        materials=[_material('Silt', 0.5), _material('Silty Clay', 0.1)],
        profile_lines=[
            ([(0.0, 40.0), (43.0, 40.0)], 0),
            ([(67.0, 40.0), (110.0, 40.0)], 0),
            ([(0.0, 26.0), (42.8, 26.0), (43.0, 22.0), (43.2, 26.0),
              (66.8, 26.0), (67.0, 22.0), (67.2, 26.0), (110.0, 26.0)], 1),
        ],
        seepage_bc={
            'specified_heads': [_head(56.0, [(0.0, 40.0), (43.0, 40.0)]),
                                _head(26.0, [(43.2, 26.0), (66.8, 26.0)]),
                                _head(56.0, [(67.0, 40.0), (110.0, 40.0)])],
            'specified_fluxes': [],
            'exit_face': [],
        },
        water_loads='manual',
    )
    return _write(sd, 'xslope_sea_trench.xlsx')


# =============================================================================
# docs/seep/samples.md: the earth dam with a filter
# =============================================================================
def build_earth_dam2():
    """A 72 ft dam on an impervious base: shell, a clay core, and a chimney
    filter that turns into a blanket drain under the downstream shell. The
    conductivities are 1e-3, 1e-4 and 1e-5 cm/s converted to ft/day
    (1e-3 cm/s = 0.864 m/day = 2.8346456693 ft/day), entered to ten decimals as
    the shipped file carries them. Pool at El. 60 on the upstream face; the
    downstream face from the crest to the toe is an exit face.

    The head boundary ends at x = 138.3, not at the face's vertex
    x = 138.3333333333 where the pool meets it; transcribed as shipped."""
    sd = _model(
        materials=[_material('shell', 2.8346456693),
                   _material('filter', 0.2834645669),
                   _material('core', 0.0283464567)],
        profile_lines=[
            ([(0.0, 0.0), (138.3333333333, 60.0), (166.0, 72.0), (194.0, 72.0),
              (346.1666666667, 6.0)], 0),
            ([(186.0, 60.0), (192.0, 60.0), (211.9, 6.0), (346.1666666667, 6.0),
              (360.0, 0.0)], 1),
            ([(152.0, 0.0), (174.0, 60.0), (186.0, 60.0), (208.0, 0.0)], 2),
        ],
        seepage_bc={
            'specified_heads': [_head(60.0, [(0.0, 0.0), (138.3, 60.0)])],
            'specified_fluxes': [],
            'exit_face': [(194.0, 72.0), (360.0, 0.0)],
        },
        water_loads='auto',
    )
    return _write(sd, 'xslope_earth_dam2.xlsx')


BUILDERS = [
    build_sea_trench,
    build_earth_dam2,
]


def main():
    os.makedirs(OUT, exist_ok=True)
    for fn in BUILDERS:
        print('built', fn())


if __name__ == '__main__':
    main()
